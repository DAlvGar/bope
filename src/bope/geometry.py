"""Geometry-driven bond-order perception.

In-house perception used when no CCD template applies: planar rings whose
bonds fall in aromatic bond-length envelopes are scored by a Hückel (4n+2)
electron-count judge over per-atom pi assignments; the winning mask sets
aromatic bonds and pyrrole-type N-H's, and the remaining bonds get orders
from length thresholds with chemistry fixups (amidines, carbonyls,
exocyclic N/S).  Over-valent non-aromatic atoms are demoted before
sanitization.

Validated to reproduce the exact molecular formulas of the 16
eval-dataset crystal ligands, including N-rich fused heteroaromatics
where OpenBabel's ``PerceiveBondOrders`` corrupts ring systems
(pentavalent carbons).
"""

from __future__ import annotations

import math
from itertools import combinations
from typing import Any

from bope._deps import Chem, rdGeometry
from bope.helpers import _planarity_rms, _sym_pair
from bope.tables import (
    _AROMATIC_ENVELOPE,
    _AROMATIC_RESCUE_RMS_MAX,
    _RESCUE_RING_BOND_MAX,
    _AROMATIC_SLACK_BOND,
    _AROMATIC_SLACK_DROP,
    _AROMATIC_SLACK_RING,
    _AROMATIC_SLACK_SHORT,
    _BOND_ORDER_TABLE,
    _CRYSTAL_CARBONYL,
    _HUCKEL,
    _MAX_VALENCE,
    _TRIPLE_BOND_TABLE,
)


def perceive_bond_orders_geometric(
    elements: list[str],
    coords: list[tuple[float, float, float]],
    graph: set[tuple[int, int]],
) -> tuple[Any | None, str | None]:
    """Perceive bond orders from geometry alone (no CCD template).

    Args:
        elements: upper-case element symbols (index i).
        coords: 3-D coordinates as ``(x, y, z)`` tuples.
        graph: heavy-atom connectivity as a set of ``(i, j)`` pairs.

    Returns:
        A 2-tuple ``(mol, err)``; *mol* is an RDKit Mol with 3-D
        coordinates and perceived bond orders (or ``None``), *err* is a
        failure reason when *mol* is ``None``.
    """
    def rdk_el(el: str) -> str:
        return el.capitalize() if len(el) == 2 else el

    n = len(elements)
    rw = Chem.RWMol()  # type: ignore[attr-defined]
    for el in elements:
        rw.AddAtom(Chem.Atom(rdk_el(el)))  # type: ignore[attr-defined]
    conf = Chem.Conformer(n)  # type: ignore[attr-defined]
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, rdGeometry.Point3D(float(x), float(y), float(z)))
    rw.AddConformer(conf)

    blen: dict[tuple[int, int], float] = {}
    for i, j in graph:
        p1, p2 = coords[i], coords[j]
        blen[_sym_pair(i, j)] = math.sqrt(
            (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2
        )

    # candidate aromatic rings: planar + all ring bonds in the envelope
    rw0 = Chem.RWMol(rw)  # type: ignore[attr-defined]
    for i, j in graph:
        rw0.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
    rings_info = Chem.GetSymmSSSR(rw0.GetMol())  # type: ignore[attr-defined]

    graph_nbrs_all: dict[int, list[int]] = {}
    for a, b in graph:
        graph_nbrs_all.setdefault(a, []).append(b)
        graph_nbrs_all.setdefault(b, []).append(a)

    cand_rings: list[list[int]] = []
    ring_excess: dict[tuple[int, ...], float] = {}
    # atoms shared with other rings: a bond from a ring atom to one of these
    # is a fusion bond or an N-aryl bond of the ring system, never an exo
    # double (the indole 2,3-bond, the purine imidazole C-N, an N-aryl
    # pyrrole) - those are single bonds that stay inside the fused aromatic
    # system or are policed by the aniline rule and the valence demotion
    # pass at assemble time.
    other_ring_atoms: set[int] = set()
    for r2 in rings_info:
        other_ring_atoms.update(r2)
    for ring in rings_info:
        rl = list(ring)
        if len(rl) < 5:
            continue  # 4-membered (and smaller) rings are never aromatic:
                      # the square faces of a saturated cage like cubane are
                      # perfectly planar and their fused 6-atom unions score
                      # a Huckel 6-pi rescue, which would mark the cage
                      # aromatic and strip its hydrogens (C8H8 -> C8).
                      # Anti-aromatic or non-aromatic, 4-rings (cubane,
                      # cyclobutadiene, biphenylene's central ring, the
                      # beta-lactam azetidinone) never take part in aromatic
                      # perception.
        rms = _planarity_rms([coords[i] for i in rl])
        if rms > 0.12:
            continue
        # a 5-membered aromatic ring cannot carry an exocyclic double bond:
        # in any valid kekule pattern the lone-pair atom (N-H / O / S) takes
        # two singles and every other atom takes exactly one ring double, so
        # an atom with a genuine exo double (valence already spent) makes
        # the ring unkekulizable.  The 8TO5 azadiene ring
        # (C2'=N-C38-C32=C2''-C2' with the exo C38=C9b) is such a
        # non-aromatic diene whose N-H form scores 6 pi and would be marked
        # aromatic - then SanitizeMol fails.  Six-membered rings are
        # unaffected (uracil / styrene / coumarin carbonyl and vinyl
        # substituents take ring singles and kekulize).
        # The length test alone cannot tell a real exo double from a short
        # single (the PKK benzofuran Ar-C(=O) bond refines to 1.36 A): a
        # bond is a genuine blocker only when the exo atom stays within
        # valence with it as a double - otherwise the valence demotion pass
        # fixes it at assemble time and the ring stays aromatic.
        if len(rl) == 5:
            rset = set(rl)
            blocked = False
            for a in rl:
                for j in graph_nbrs_all.get(a, ()):
                    if j in rset or j in other_ring_atoms:
                        continue
                    dlen = blen[_sym_pair(a, j)]
                    pair = _sym_pair(elements[a], elements[j])
                    if pair in _TRIPLE_BOND_TABLE and dlen <= _TRIPLE_BOND_TABLE[pair]:
                        continue  # a nitrile / alkyne substituent is a
                                  # single-atom triple, not a double
                    dmax = _BOND_ORDER_TABLE.get(pair, (None, None))[1]
                    if dmax is None or dlen > dmax:
                        continue
                    vmax = _MAX_VALENCE.get(elements[j])
                    if vmax is not None:
                        val = 2.0  # this bond as a double
                        for k in graph_nbrs_all.get(j, ()):
                            if k == a:
                                continue
                            dk = blen[_sym_pair(j, k)]
                            pk = _sym_pair(elements[j], elements[k])
                            if pk in _TRIPLE_BOND_TABLE and dk <= _TRIPLE_BOND_TABLE[pk]:
                                val += 3.0
                            else:
                                dkmax = _BOND_ORDER_TABLE.get(pk, (None, None))[1]
                                val += 2.0 if dkmax is not None and dk <= dkmax else 1.0
                        if val > vmax:
                            continue  # demotable - leave it to the demotion pass
                    blocked = True
                    break
                if blocked:
                    break
            if blocked:
                continue
        ok = True
        excess = 0.0
        for k in range(len(rl)):
            i, j = rl[k], rl[(k + 1) % len(rl)]
            d = blen[_sym_pair(i, j)]
            lo, hi = _AROMATIC_ENVELOPE.get(
                _sym_pair(elements[i], elements[j]), (1.27, 1.50)
            )
            if d < lo - _AROMATIC_SLACK_SHORT:
                ok = False  # short bonds stay strict: a double-embedded
                break       # ring is a diene, not an aromatic candidate
            if d > hi:
                over = d - hi
                if over > _AROMATIC_SLACK_BOND:
                    ok = False
                    break
                excess += over
        if not ok or excess > _AROMATIC_SLACK_RING:
            continue
        # a ring whose lengths need slack is a low-resolution aromatic
        # candidate: demand better planarity than the generous 0.12, or a
        # saturated ether ring (1T3R/017 THF: C-O/C-C at 1.40-1.51, rms
        # 0.115) is let through - 07L's coumarin at rms 0.053 qualifies,
        # its two C-O at 1.444/1.462 and one C-C at 1.514 riding the slack
        # (each at most +0.022, within the per-bond cap).
        if excess > 0.0 and rms > 0.08:
            continue
        ring_excess[tuple(rl)] = excess
        cand_rings.append(rl)

    cand_atoms = set()
    for rl in cand_rings:
        cand_atoms.update(rl)

    # ring sigma from candidate-ring EDGES only: an exo substituent that is
    # itself a candidate ring (triazolyl-phenyl) must not inflate an atom's
    # ring sigma - the N-C(phenyl) bond is NOT part of the triazole
    ring_edges: set[tuple[int, int]] = set()
    for rl in cand_rings:
        for k in range(len(rl)):
            i, j = rl[k], rl[(k + 1) % len(rl)]
            ring_edges.add(_sym_pair(i, j))
    ring_sigma: dict[int, int] = {a: 0 for a in cand_atoms}
    for i, j in graph:
        if _sym_pair(i, j) in ring_edges:
            ring_sigma[i] += 1
            ring_sigma[j] += 1

    def exo_sigma(a: int) -> int:
        return sum(1 for i, j in graph if a in (i, j)) - ring_sigma[a]

    # total heavy-atom degree of every atom, and the exo neighbours of each
    # candidate atom (any graph neighbour not joined via a candidate ring
    # edge).  Used by the pi count for ring carbons bearing a
    # directly-attached terminal O (carbonyl): their p orbital is in the C=O
    # pi bond.
    deg: dict[int, int] = {}
    for a, b in graph:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    exo_nbrs: dict[int, list[int]] = {a: [] for a in cand_atoms}
    ring_nbrs: dict[int, list[int]] = {a: [] for a in cand_atoms}
    graph_nbrs: dict[int, list[int]] = {a: [] for a in range(n)}
    for a, b in graph:
        graph_nbrs[a].append(b)
        graph_nbrs[b].append(a)
        if _sym_pair(a, b) in ring_edges:
            ring_nbrs[a].append(b)
            ring_nbrs[b].append(a)
            continue
        if a in exo_nbrs:
            exo_nbrs[a].append(b)
        if b in exo_nbrs:
            exo_nbrs[b].append(a)

    def carbonyl_c(c: int, dmax: float | None = None) -> bool:
        """True for a C with a directly-attached terminal O at C=O length.
        The default bound is the strict 1.30 (length-rule double cutoff):
        its p orbital sits in the C=O pi bond, so such a ring carbon
        contributes 0 pi.  The N-pyridone/amide discriminator passes the
        generous crystal-carbonyl bound (1.40) - a caffeine C=O refines to
        1.34-1.36 A and its N-methyl is amide-type regardless."""
        if elements[c] != "C":
            return False
        if dmax is None:
            dmax = _BOND_ORDER_TABLE[("C", "O")][1]
        for n in graph_nbrs[c]:
            if elements[n] == "O" and deg[n] == 1 and blen[_sym_pair(c, n)] <= dmax:
                return True
        return False

    def pi_of(a: int, h_choice: bool) -> int | None:
        """pi electrons contributed by candidate atom *a*."""
        el = elements[a]
        if el == "C":
            # a ring carbon with a directly-attached terminal O at C=O
            # distance (<= 1.30, the length-rule double cutoff) has its p
            # orbital in the C=O pi bond and contributes 0 pi to the ring.
            # An O attached through another atom (aryl ketones), an O-H /
            # O-R substituent, or a long C-O (phenol/enol) leaves 1 pi.
            # Without this, uracil-type rings score 8 pi with pyrrole N-H
            # (not Huckel) and 6 pi with the wrong pyridine-type N-H - the
            # correct 6-pi count needs the 2 carbonyl carbons at 0 and
            # both N's at 2.
            return 0 if carbonyl_c(a) else 1
        if el == "N":
            rs, exo = ring_sigma[a], exo_sigma(a)
            if rs == 2 and exo == 1:
                # 5-ring: N-methyl/aryl-pyrrole, lone pair in the ring
                # plane-free p orbital: 2 pi.  6-ring with a plain N-alkyl:
                # pyridinium, the alkyl consumes the lone pair: 1 pi.  An
                # N-alkyl 6-ring N bonded to (or exo to) a carbonyl is
                # pyridone/amide-type: lone pair delocalised, 2 pi.
                # Without the 6-ring rule the nicotinamide rings of
                # NAD/NAP/NDP score 7 pi, fail Huckel and come out
                # saturated (+6 H in the formula).
                if max((len(rl) for rl in cand_rings if a in rl), default=0) <= 5:
                    return 2
                if any(carbonyl_c(j, _CRYSTAL_CARBONYL) for j in ring_nbrs[a]):
                    return 2
                if any(carbonyl_c(e, _CRYSTAL_CARBONYL) for e in exo_nbrs[a]):
                    return 2
                # An N-alkyl 6-ring N is pyridinium only in a pyridine-like
                # ring (exactly one N).  A ring holding a second N is either
                # a saturated lactam/amine (1D1: the sibling N is amide-type
                # with a heavy exo substituent, so the alkyl N has no p
                # orbital - None) or a fused N-heteroaromatic where the alkyl
                # N is a neutral pyrrole-type N with its lone pair in the
                # ring p orbital (flavin N10 in FMN/RS3/FAD: all sibling N's
                # are pyridine-type with no exo - 2 pi; the fused system
                # supplies the pi the per-ring count lacks).  The amide-type
                # sibling is the discriminator: pyridine-type N's carry no
                # heavy exo substituent.
                for rl in cand_rings:
                    if a in rl and len(rl) > 5:
                        others = [b for b in rl if elements[b] == "N" and b != a]
                        if not others:
                            return 1  # exactly one N: plain pyridinium
                        if any(exo_sigma(b) > 0 for b in others):
                            return None  # amide-type sibling: saturated amine
                        return 2  # pyridine-type siblings: pyrrole-like N
                return 1
            if rs == 3 and exo == 0:
                # fusion N shared by two rings (triazolo-triazine bridgehead,
                # etc.): 3 ring sigma, no H - RDKit gives it Two electrons
                # (countAtomElec: dv 3, degree 3, nlp 2 -> 2).  The old 1-pi
                # assignment made RNL/QUP/9KI triazoles score 5-7, failing
                # Huckel and gaining +1 H in the formula; 2 pi makes the
                # triazole itself Huckel (6) and unblocks kekulization.
                return 2
            if rs == 2 and exo == 0:
                return 2 if h_choice else 1  # pyrrole(1H) vs pyridine(0H)
            return None
        if el in ("O", "S"):
            return 2 if ring_sigma[a] == 2 else None
        return None

    def attempt(cand_rings: list[list[int]]) -> tuple[Any | None, str | None]:
        # flexible N's: 2 ring sigma, no exo (pyrrole-vs-pyridine tautomer bit)
        flex_n = [
            a for a in sorted(cand_atoms)
            if elements[a] == "N" and ring_sigma[a] == 2 and exo_sigma(a) == 0
        ]
        flex_idx = {a: i for i, a in enumerate(flex_n)}

        # fused systems: rings sharing an edge (>= 2 atoms) belong to one
        # conjugated system, aromaticity decided by the exact RDKit
        # subset-union rule in arom_for(): a ring is aromatic iff it is in a
        # fused subset whose 1-2-ring atom union is Huckel.  Isoalloxazine's
        # pyrazine ring is 7 pi per-ring and its uracil ring 5 pi, yet the
        # {uracil,pyrazine} union is 10 (and the tricycle 14) - both 4n+2;
        # the triazolo-triazines of RNL/QUP/9KI pass the same way (6-ring
        # 7 pi, 5-ring 6 pi, union 10).
        ring_sys = list(range(len(cand_rings)))
        for i in range(len(cand_rings)):
            for j in range(i):
                if len(set(cand_rings[i]) & set(cand_rings[j])) >= 2:
                    si, sj = ring_sys[i], ring_sys[j]
                    for k in range(len(cand_rings)):
                        if ring_sys[k] == sj:
                            ring_sys[k] = si
        sys_rings: dict[int, list[int]] = {}
        for i, s in enumerate(ring_sys):
            sys_rings.setdefault(s, []).append(i)

        def pi_val(a: int, mask: int) -> int | None:
            return pi_of(a, (mask >> flex_idx[a]) & 1 if a in flex_idx else 0)

        def huckel_subset(mask: int, rings: list[int]) -> bool:
            """Exact RDKit fused-ring rule (applyHuckelToFused): the union of
            the subset's ring atoms - each atom present in 1-2 rings of the
            subset counted once, atoms in 3+ rings excluded (they have no
            p-orbital contribution, the acepentalene fix) - must be a Huckel
            4n+2 pi set.  The subset must also be fused (connected), checked
            by the caller."""
            cnt: dict[int, int] = {}
            for ri in rings:
                for a in cand_rings[ri]:
                    cnt[a] = cnt.get(a, 0) + 1
            unon = [a for a, c in cnt.items() if c in (1, 2)]
            if len(unon) < 3:
                return False
            tot = 0
            for a in unon:
                p = pi_val(a, mask)
                if p is None:
                    return False  # an sp3-type N poisons the subset, not a
                                  # 0-pi atom (1D1's lactam alkyl N)
                tot += p
            return tot in _HUCKEL

        def subset_fused(rings: list[int]) -> bool:
            """checkFused: the subset's rings must form one connected
            component via shared edges (>= 2 shared atoms)."""
            n = len(rings)
            if n <= 1:
                return True
            parent = list(range(n))

            def find(x: int) -> int:
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            for i in range(n):
                for j in range(i):
                    if len(set(cand_rings[rings[i]]) & set(cand_rings[rings[j]])) >= 2:
                        ri, rj = find(i), find(j)
                        if ri != rj:
                            parent[rj] = ri
            return len({find(i) for i in range(n)}) == 1

        def arom_for(mask: int) -> tuple[set[int], set[int]]:
            """(aromatic ring indices, per-ring aromatic ring indices) per
            the exact RDKit fused-subset rule: a ring is aromatic iff it
            belongs to a fused subset (sizes 1..min(n,6)) whose 1-2-ring
            atom union is Huckel.  Trying subsets in size order and marking
            every ring of a passing subset reproduces the C++ loop."""
            arom: set[int] = set()
            per_ring: set[int] = set()
            def rescue_eligible(ri: int) -> bool:
                """A ring may be rescued through a fused subset only if it
                shows no reduction and contains no in-ring ketone-type C=O.
                A reduced ring (an sp3 atom: the FADH2 pyrazine of 7VKD is
                puckered with 1.47-1.48 A edges) must stay out of the
                rescue, which would rewrite reduced hydrogens back into the
                formula.  Two independent reduction tests, OR'd because
                neither works alone: at fine resolution the oxidized RS3
                pyrazine is flat but refines to the same 1.48-1.49 A edges,
                while at 2.5-3.0 A the oxidized 8QIN flavin rings are bent
                by pure coordinate noise (rms 0.096) but keep short
                1.36-1.40 A edges.  A ring passes if it is essentially
                planar (rms <= _AROMATIC_RESCUE_RMS_MAX) or free of long
                edges (<= _RESCUE_RING_BOND_MAX, C=O-adjacent excepted: a
                normal aromatic uracil carries ~1.47 A C-C(=O) singles).
                The ketone gate keeps a cyclohexadienone ring (in-ring
                C-C(=O)-C, e.g. the 537/1UKI ketone ring) out of the
                rescue: its C=O carbon's p orbital sits in the exo double
                and the ring is a non-aromatic dienone, and a ring carbon
                carrying an exo double cannot take the alternating kekule
                pattern - the mol is then unkekulizable.  Amide-type C=O
                (a ring neighbour N - the flavin uracil, the RNL
                pyrimidinone) stays rescuable: those rings are genuinely
                heteroaromatic."""
                if (_planarity_rms([coords[i] for i in cand_rings[ri]])
                        > _AROMATIC_RESCUE_RMS_MAX):
                    rl = cand_rings[ri]
                    for k in range(len(rl)):
                        i, j = rl[k], rl[(k + 1) % len(rl)]
                        if blen[_sym_pair(i, j)] <= _RESCUE_RING_BOND_MAX:
                            continue
                        if (elements[i] == "C" and carbonyl_c(i)) or (
                                elements[j] == "C" and carbonyl_c(j)):
                            continue
                        return False
                rset = set(cand_rings[ri])
                for a in cand_rings[ri]:
                    if elements[a] != "C" or not carbonyl_c(a):
                        continue
                    rnbrs = [b for b in graph_nbrs[a] if b in rset]
                    if len(rnbrs) == 2 and all(elements[b] == "C" for b in rnbrs):
                        return False
                return True

            for rings in sys_rings.values():
                n = len(rings)
                for size in range(1, min(n, 6) + 1):
                    for comb in combinations(rings, size):
                        if not subset_fused(comb):
                            continue
                        if size > 1 and not all(rescue_eligible(ri) for ri in comb):
                            continue
                        if not huckel_subset(mask, comb):
                            continue
                        for ri in comb:
                            arom.add(ri)
                        if size == 1:
                            per_ring.add(comb[0])
            return arom, per_ring

        def score(mask: int) -> tuple[int, int, int]:
            """(aromatic rings, per-ring aromatic rings, amide-H count).
            The total counts rings aromatic per-ring OR via any passing fused
            subset (isoalloxazine's pyrazine is 7 pi per-ring and its uracil
            5, yet the {uracil,pyrazine} union is 10 and the tricycle 14 -
            both 4n+2); the per-ring count breaks ties toward the mask that
            needs no system rescue (purine's N7-H tautomer over the
            pyrimidine-N1-H tautomer, which is aromatic only via the
            system); the amide-H count then breaks flavin-type ties: the
            pyrrole-H prefers the N with the most carbonyl ring neighbours
            (flavin N3 between the two C=O's, matching the ccd, over N1)."""
            arom, per_ring = arom_for(mask)
            ok = len(arom)
            per_ring_n = len(per_ring)
            amide_h = 0
            for a in flex_n:
                if (mask >> flex_idx[a]) & 1:
                    amide_h += sum(
                        1 for j in ring_nbrs[a]
                        if carbonyl_c(j, _CRYSTAL_CARBONYL)
                    )
            return ok, per_ring_n, amide_h

        best_mask, best_score = None, (-1, -1, -1)
        for mask in range(1 << len(flex_n)):
            s = score(mask)
            if s > best_score:
                best_score, best_mask = s, mask

        arom, _ = arom_for(best_mask)
        arom_rings = [cand_rings[ri] for ri in sorted(arom)]
        arom_atoms = set()
        for rl in arom_rings:
            arom_atoms.update(rl)
        arom_bonds = set()
        for rl in arom_rings:
            for k in range(len(rl)):
                i, j = rl[k], rl[(k + 1) % len(rl)]
                arom_bonds.add(_sym_pair(i, j))
        def assemble(exo_force_all: bool) -> tuple[Any | None, str | None]:
            """Build the mol with bond orders.  exo_force_all=False: only degree-1
            N/S exo to an aromatic atom (NH2 / amino-pyridine) is forced single;
            substituted N/S (C=N imines, N-aryl triazoles) stay on the length rule
            so a triazine C=N survives.  exo_force_all=True: every N/S exo to an
            aromatic atom is single (aniline fallback - fixes kekulization when the
            length rule puts a double on a benzene-ring carbon)."""
            m = Chem.RWMol(rw)  # type: ignore[attr-defined]
            for i, j in graph:
                if _sym_pair(i, j) in arom_bonds:
                    m.AddBond(i, j, Chem.BondType.AROMATIC)  # type: ignore[attr-defined]
                else:
                    d = blen[_sym_pair(i, j)]
                    pair = _sym_pair(elements[i], elements[j])
                    if pair in _TRIPLE_BOND_TABLE and d <= _TRIPLE_BOND_TABLE[pair]:
                        # nitrile / alkyne: unmistakably short, checked first so an
                        # N/S exo triple on an aromatic ring (benzonitrile) survives
                        m.AddBond(i, j, Chem.BondType.TRIPLE)  # type: ignore[attr-defined]
                    elif (i in arom_atoms) != (j in arom_atoms):
                        # exactly one aromatic endpoint.  A non-aromatic atom
                        # bonded to two or more aromatic atoms is a macrocycle
                        # bridge (the porphyrin methine carbons of HEM/HEC/
                        # ZNH, each sitting between two pyrrole rings): its
                        # bonds to them are single, never double - the length
                        # rule reads the 1.37-1.40 A bridges as C=C doubles,
                        # and a bridge double saturates the pyrrole ring
                        # carbon (valence 5 with its ring bonds), which then
                        # cannot take its ring double and the kekulization
                        # dies.  A genuine exo double substituent (vinyl,
                        # carbonyl, imine) has exactly one aromatic neighbour.
                        narm = i if i not in arom_atoms else j
                        arom_nbrs = sum(
                            1 for a, b in graph
                            if (a == narm and b in arom_atoms)
                            or (b == narm and a in arom_atoms)
                        )
                        if arom_nbrs >= 2:
                            m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                        elif elements[i] in ("N", "S") or elements[j] in ("N", "S"):
                            if exo_force_all:
                                m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                            else:
                                ns = i if elements[i] in ("N", "S") else j
                                deg_ns = sum(1 for a, b in graph if a == ns or b == ns)
                                if deg_ns == 1:
                                    # aniline / amino-pyridine NH2: single, always
                                    m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                                elif pair in _BOND_ORDER_TABLE:
                                    _, dmax = _BOND_ORDER_TABLE[pair]
                                    m.AddBond(
                                        i, j,
                                        Chem.BondType.DOUBLE if (dmax is not None and d <= dmax)
                                        else Chem.BondType.SINGLE,
                                    )  # type: ignore[attr-defined]
                                else:
                                    m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                        elif pair in _BOND_ORDER_TABLE:
                            # exo C-C / C-O / C-S substituent (methyl, vinyl,
                            # carbonyl) of a single aromatic neighbour: on the
                            # length rule like the non-aromatic pairs, so a
                            # short conjugated C=C or C=O stays double while
                            # the 1.43-1.55 A ring-substituent singles stay
                            # single.
                            _, dmax = _BOND_ORDER_TABLE[pair]
                            m.AddBond(
                                i, j,
                                Chem.BondType.DOUBLE if (dmax is not None and d <= dmax)
                                else Chem.BondType.SINGLE,
                            )  # type: ignore[attr-defined]
                        else:
                            m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                    elif i in arom_atoms:
                        # both endpoints aromatic but the bond is not a ring
                        # edge: an inter-ring linkage (N-aryl, S-aryl,
                        # biphenyl, a fused-ring junction).  These are
                        # single, never on the length rule - a double
                        # between two aromatic atoms has no valid kekule,
                        # and N-aryl C-N refines to 1.33-1.43 A in crystals
                        # (the N-aryl pyrrole of QUP/5AEP measures 1.337),
                        # inside the C=N cutoff.  Ring-internal doubles live
                        # in the AROMATIC branch above.
                        m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                    elif pair in _BOND_ORDER_TABLE:
                        _, dmax = _BOND_ORDER_TABLE[pair]
                        m.AddBond(
                            i, j,
                            Chem.BondType.DOUBLE if (dmax is not None and d <= dmax)
                            else Chem.BondType.SINGLE,
                        )  # type: ignore[attr-defined]
                    else:
                        m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]

            # --- fixup pass: corrections the pure length rule gets wrong ---------
            # (1) carbonyl rescue: non-aromatic C with a terminal O single bond and
            #     at least one N neighbor -> C=O.  Crystal carbonyls often refine
            #     to 1.34-1.36 A (caffeine C2/C6 in 3RFM); aromatic ring C's are
            #     excluded, and internal O's (ether/alcohol C-O) stay single.
            # (2) amidine/imine: non-aromatic C with exactly 2 N single neighbors
            #     and no O/S neighbor: the shorter C-N becomes double (delocalized
            #     amidines measure 1.31-1.33 for BOTH bonds; neutral benzamidine
            #     needs one).
            for i in range(n):
                if elements[i] != "C" or i in arom_atoms:
                    continue
                nbrs = []
                for a, b in graph:
                    if a == i:
                        nbrs.append(b)
                    elif b == i:
                        nbrs.append(a)
                n_nbrs = [j for j in nbrs if elements[j] == "N"]
                o_nbrs = [j for j in nbrs if elements[j] == "O"]
                if o_nbrs and n_nbrs:
                    for o in o_nbrs:
                        deg_o = sum(1 for a, b in graph if a == o or b == o)
                        if deg_o == 1 and blen[_sym_pair(i, o)] <= 1.40:
                            bo = m.GetBondBetweenAtoms(i, o)  # type: ignore[attr-defined]
                            if bo.GetBondType() == Chem.BondType.SINGLE:  # type: ignore[attr-defined]
                                bo.SetBondType(Chem.BondType.DOUBLE)  # type: ignore[attr-defined]
                elif len(n_nbrs) == 2 and not any(
                    elements[j] in ("O", "S") for j in nbrs
                ):
                    d1, d2 = blen[_sym_pair(i, n_nbrs[0])], blen[_sym_pair(i, n_nbrs[1])]
                    # 1.40 (not 1.36): delocalized amidines / imines measure
                    # 1.31-1.36 (ETKDG +0.03 bias, plus noise), while a genuine
                    # C-N single pair sits at 1.45+ - the 0.05 gap keeps the
                    # raised threshold safe.
                    if max(d1, d2) <= 1.40:
                        tgt = n_nbrs[0] if d1 <= d2 else n_nbrs[1]
                        bo = m.GetBondBetweenAtoms(i, tgt)  # type: ignore[attr-defined]
                        if bo.GetBondType() == Chem.BondType.SINGLE:  # type: ignore[attr-defined]
                            bo.SetBondType(Chem.BondType.DOUBLE)  # type: ignore[attr-defined]

            # aromatic atoms: flags + pyrrole-N explicit H
            for a in arom_atoms:
                at = m.GetAtomWithIdx(a)  # type: ignore[attr-defined]
                at.SetIsAromatic(True)  # type: ignore[attr-defined]
                at.SetHybridization(Chem.HybridizationType.SP2)  # type: ignore[attr-defined]
                if (
                    elements[a] == "N" and ring_sigma[a] == 2 and exo_sigma(a) == 0
                    and ((best_mask >> flex_idx[a]) & 1)
                ):
                    at.SetNoImplicit(True)  # type: ignore[attr-defined]
                    at.SetNumExplicitHs(1)  # type: ignore[attr-defined]

            # pyridinium: a 6-ring aromatic N with a heavy exo substituent and no
            # carbonyl adjacency contributes 1 pi and cannot take a ring double
            # (its lone pair is consumed by the exo bond): the neutral form is
            # unkekulizable, the +1 charge is forced by valence.  This is exactly
            # the nicotinamide N of NAD/NAP/NDP (and any N-alkyl pyridinium).
            # N-alkyl 5-ring pyrroles (2 pi) and pyridone/amide N's (2 pi) stay
            # neutral - the same tests as the pi-count's pyridinium branch.
            for a in arom_atoms:
                if (
                    elements[a] == "N"
                    and ring_sigma[a] == 2
                    and exo_sigma(a) == 1
                    and max((len(rl) for rl in cand_rings if a in rl), default=0) > 5
                    and all(
                        sum(1 for b in rl if elements[b] == "N") == 1
                        for rl in cand_rings if a in rl and len(rl) > 5
                    )
                    and not any(carbonyl_c(j, _CRYSTAL_CARBONYL) for j in ring_nbrs[a])
                    and not any(carbonyl_c(e, _CRYSTAL_CARBONYL) for e in exo_nbrs[a])
                ):
                    m.GetAtomWithIdx(a).SetFormalCharge(1)  # type: ignore[attr-defined]

            mol = m.GetMol()

            # valence demotion pass on non-aromatic atoms (before sanitize)
            changed = True
            while changed:
                changed = False
                for a in mol.GetAtoms():
                    vmax = _MAX_VALENCE.get(a.GetSymbol())
                    if vmax is None or a.GetIsAromatic() or a.GetFormalCharge() != 0:
                        continue
                    val = sum(b.GetBondTypeAsDouble() for b in a.GetBonds())
                    if val == vmax and a.GetSymbol() == "C":
                        # exactly-valent sp3-looking carbon holding a C=N
                        # double: the length rule reads a noisy amine C-N as
                        # a short imine (the STU sugar 4'-N-methylamino C-N
                        # at 1.416 A measures 1.355 under 0.03-A bond-RMS
                        # noise - inside the raised 1.37 imine cutoff that
                        # the 8TO5 azadiene needs).  A genuine imine carbon
                        # is sp2 with at most one carbon single (the 8TO5
                        # ring C2'=N, PLP's pyridoxal C=N against an
                        # aromatic ring); a carbon with two saturated
                        # single-bonded carbon neighbours cannot hold a
                        # double, so demote it.  Real C=N's survive:
                        # oximes bond to O, amidines/guanidines to N's,
                        # N-H ketimines have no heavy single on the N, and
                        # conjugated imines fail the all-single-neighbour
                        # test.
                        db = [
                            b for b in a.GetBonds()
                            if b.GetBondType() == Chem.BondType.DOUBLE
                            and not b.GetIsAromatic()
                        ]
                        if len(db) == 1:
                            nbr = db[0].GetOtherAtom(a)
                            if nbr.GetSymbol() == "N" and not nbr.GetIsAromatic():
                                # exclude the double itself by endpoints - RDKit
                                # wraps the same C++ bond in a fresh Python
                                # object per GetBonds() call, so `is not` fails
                                # and the double would be counted as a single.
                                dbl_ends = (
                                    db[0].GetBeginAtomIdx(), db[0].GetEndAtomIdx()
                                )
                                cn_singles = [
                                    b for b in a.GetBonds()
                                    if (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) != dbl_ends
                                    and (b.GetEndAtomIdx(), b.GetBeginAtomIdx()) != dbl_ends
                                ]
                                if (
                                    len(cn_singles) >= 2
                                    and all(
                                        b.GetBondType() == Chem.BondType.SINGLE  # type: ignore[attr-defined]
                                        and b.GetOtherAtom(a).GetSymbol() == "C"
                                        and not b.GetOtherAtom(a).GetIsAromatic()
                                        and all(
                                            bb.GetBondType() == Chem.BondType.SINGLE  # type: ignore[attr-defined]
                                            for bb in b.GetOtherAtom(a).GetBonds()
                                        )
                                        for b in cn_singles
                                    )
                                    and any(
                                        b.GetBondType() == Chem.BondType.SINGLE  # type: ignore[attr-defined]
                                        and b.GetOtherAtom(nbr).GetSymbol() == "C"
                                        for b in nbr.GetBonds()
                                    )
                                ):
                                    mol.GetBondBetweenAtoms(  # type: ignore[attr-defined]
                                        db[0].GetBeginAtomIdx(), db[0].GetEndAtomIdx()
                                    ).SetBondType(Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                                    changed = True
                                    break
                    if val > vmax:
                        db = [
                            b for b in a.GetBonds()
                            if b.GetBondType() == Chem.BondType.DOUBLE
                            and not b.GetIsAromatic()
                        ]
                        if not db:
                            # An N with exactly four single bonds and no double to
                            # demote is a quaternary ammonium: the geometry alone
                            # cannot know about the formal charge, but the +1 is
                            # forced by valence (the only neutral alternative
                            # fails sanitization).  Charged groups are otherwise
                            # out of scope (the API's ``charge`` argument is
                            # unused).  An N over-valent without four singles
                            # (e.g. a spurious nitrile triple from a bad
                            # geometry) is not resolvable - report the failure.
                            if a.GetSymbol() == "N" and (
                                sum(1 for b in a.GetBonds()) == 4
                                and all(
                                    b.GetBondType() == Chem.BondType.SINGLE  # type: ignore[attr-defined]
                                    for b in a.GetBonds()
                                )
                            ):
                                a.SetFormalCharge(1)
                                changed = True
                                break
                            return None, f"over-valent {a.GetSymbol()} {a.GetIdx()} no double to demote"
                        longest = max(
                            db, key=lambda b: blen[
                                _sym_pair(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                            ]
                        )
                        mol.GetBondBetweenAtoms(  # type: ignore[attr-defined]
                            longest.GetBeginAtomIdx(), longest.GetEndAtomIdx()
                        ).SetBondType(Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                        changed = True
                        break
            try:
                Chem.SanitizeMol(mol)  # type: ignore[attr-defined]
            except Exception as exc:  # noqa: BLE001
                return None, f"santize: {exc}"
            return mol, None

        # pass 1: keep length-rule doubles on substituted N/S (imines, triazine
        # C=N).  If RDKit cannot kekulize (a length double on a benzene-ring
        # carbon, aniline-type), pass 2 forces every N/S exo single.
        mol, err = assemble(exo_force_all=False)
        if mol is None:
            mol, err = assemble(exo_force_all=True)
        if mol is None:
            return None, err
        return mol, None

    # pass 1/2 failed on the full candidate set: a ring admitted only
    # through noise-level slack (1D1's pyridinone ring at +0.003, fused
    # so its forced double makes the neighbour ring unkekulizable) can be
    # the odd one out - drop slack rings at or below the drop line and
    # retry once.  True low-resolution aromatics (07L +0.041, NDP +0.012)
    # sit well above it.
    mol, err = attempt(cand_rings)
    if mol is None:
        drop = [rl for rl in cand_rings
                if 0.0 < ring_excess.get(tuple(rl), 0.0) <= _AROMATIC_SLACK_DROP]
        if drop:
            mol, err = attempt([rl for rl in cand_rings if rl not in drop])
    return mol, err
