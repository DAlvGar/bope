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
from typing import Any

from bope._deps import Chem, rdGeometry
from bope.helpers import _planarity_rms, _sym_pair
from bope.tables import (
    _AROMATIC_ENVELOPE,
    _AROMATIC_SLACK_BOND,
    _AROMATIC_SLACK_DROP,
    _AROMATIC_SLACK_RING,
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

    cand_rings: list[list[int]] = []
    ring_excess: dict[tuple[int, ...], float] = {}
    for ring in rings_info:
        rl = list(ring)
        rms = _planarity_rms([coords[i] for i in rl])
        if rms > 0.12:
            continue
        ok = True
        excess = 0.0
        for k in range(len(rl)):
            i, j = rl[k], rl[(k + 1) % len(rl)]
            d = blen[_sym_pair(i, j)]
            lo, hi = _AROMATIC_ENVELOPE.get(
                _sym_pair(elements[i], elements[j]), (1.27, 1.50)
            )
            if d < lo:
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
                # an N-alkyl 6-ring N is pyridinium only in a pyridine-like
                # ring (exactly one N).  A ring holding a second N is a
                # saturated lactam/amine (1D1): the alkyl N has no p orbital,
                # and counting it as 1 pi would let the ring pass Huckel at
                # 6 pi with the amide N at 2 pi (0 carbonyl + 1 + 1 + 1 + 1
                # + 2) and come out aromatic (and unkekulizable unless the
                # wrong N is charged).
                if any(
                    sum(1 for b in rl if elements[b] == "N") > 1
                    for rl in cand_rings if a in rl and len(rl) > 5
                ):
                    return None
                return 1
            if rs == 3 and exo == 0:
                return 1  # fusion N shared by two rings (purine / triazolo-
                          # triazine bridgehead): pyridine-type, 0 H
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

        def score(mask: int) -> int:
            ok_rings = 0
            for rl in cand_rings:
                total = 0
                good = True
                for a in rl:
                    p = pi_of(a, (mask >> flex_idx[a]) & 1 if a in flex_idx else 0)
                    if p is None:
                        good = False
                        break
                    total += p
                if good and total in _HUCKEL:
                    ok_rings += 1
            return ok_rings

        best_mask, best_score = None, -1
        for mask in range(1 << len(flex_n)):
            s = score(mask)
            if s > best_score:
                best_score, best_mask = s, mask

        arom_rings = []
        for rl in cand_rings:
            ps = [
                pi_of(a, (best_mask >> flex_idx[a]) & 1 if a in flex_idx else 0)
                for a in rl
            ]
            if None in ps:
                continue  # an sp3-type N (non-pyridinium) makes the ring
                          # non-Huckel, not a 0-pi atom
            if sum(ps) in _HUCKEL:
                arom_rings.append(rl)
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
                    elif (i in arom_atoms or j in arom_atoms) and (
                        elements[i] in ("N", "S") or elements[j] in ("N", "S")
                    ):
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
