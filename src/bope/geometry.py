"""Geometry-driven bond-order perception.

In-house perception used when no CCD template applies: planar rings whose
bonds fall in aromatic bond-length envelopes are scored by a Hückel (4n+2)
electron-count judge over per-atom pi assignments; the winning mask sets
aromatic bonds and pyrrole-type N-H's, and the remaining bonds get orders
from length thresholds with chemistry fixups (carbonyls, amidines,
nitro, phosphate P=O, sulfonamides).  Over-valent non-aromatic atoms
are demoted before sanitization.

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
from bope.fixups import FixupEngine
from bope.helpers import _planarity_rms, _sym_pair
from bope.tables import (
    _AROMATIC_ENVELOPE,
    _AROMATIC_HARD_MAX,
    _AROMATIC_RESCUE_RMS_MAX,
    _AROMATIC_SLACK_BOND,
    _AROMATIC_SLACK_DROP,
    _AROMATIC_SLACK_RING,
    _AROMATIC_SLACK_SHORT,
    _BOND_ORDER_TABLE,
    _CRYSTAL_CARBONYL,
    _HUCKEL,
    _MAX_VALENCE,
    _RESCUE_RING_BOND_MAX,
    _TRIPLE_BOND_TABLE,
)


def _find(parent: list[int], x: int) -> int:
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


class _GeometricPerceiver:
    """Geometric bond-order perception for one molecule.

    Private state container for :func:`perceive_bond_orders_geometric`: the
    former closure stack of the monolith function as a small object, so the
    perception pipeline (ring candidacy, the Huckel judge, the fixup passes)
    can be read and tested in isolation.  Instances are single-use.
    """

    def __init__(
        self,
        elements: list[str],
        coords: list[tuple[float, float, float]],
        graph: set[tuple[int, int]],
    ) -> None:
        self.elements = elements
        self.coords = coords
        self.graph = graph
        self.n = len(elements)

    def perceive(self) -> tuple[Any | None, str | None]:
        """Perceive bond orders for the molecule.

        Runs the pipeline: candidate-ring filter over the GetSymmSSSR
        rings, then the Huckel-judge attempt with its three-pass
        """
        elements, coords, graph = self.elements, self.coords, self.graph
        n = len(elements)
        rw = Chem.RWMol()  # type: ignore[attr-defined]
        for el in elements:
            rw.AddAtom(Chem.Atom(self._rdk_el(el)))  # type: ignore[attr-defined]
        conf = Chem.Conformer(n)  # type: ignore[attr-defined]
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, rdGeometry.Point3D(float(x), float(y), float(z)))
        rw.AddConformer(conf)
        self.rw = rw

        self.blen: dict[tuple[int, int], float] = {}
        for i, j in graph:
            p1, p2 = coords[i], coords[j]
            self.blen[_sym_pair(i, j)] = math.sqrt(
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

        self.cand_rings: list[list[int]] = []
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
        # edges shared by two rings: a fusion edge (indole 2,3, the
        # benzofuran C4-C5).  The double-edge sanity gate below must not
        # count them - every aromatic 5-ring fused to a benzene shares a
        # short aromatic edge, and V0X's dihydrobenzofuran ring leans on
        # exactly that one.
        edge_ring_count: dict[tuple[int, int], int] = {}
        for r2 in rings_info:
            rl2 = list(r2)
            for k in range(len(rl2)):
                e = _sym_pair(rl2[k], rl2[(k + 1) % len(rl2)])
                edge_ring_count[e] = edge_ring_count.get(e, 0) + 1
        # a SINGLE-ring molecule whose ring is broken in both directions is
        # refused so the caller can fall back to OpenBabel: an edge beyond the
        # envelope + per-bond slack AND an edge below the envelope lo - slack
        # together mean the ring does not hold together in these coordinates
        # (the ETKDG thiazole embed refines S-C to 1.90 A and C-N to 1.16 A -
        # a fragment, not a 5-membered aromatic).  Real data never shows the
        # combination: low resolution stretches rings, but never also
        # compresses an edge below ~1.25 A (planar rings across the crystal
        # datasets carry over-long edges at 2.5-3.0 A resolution - porphyrins,
        # nucleotides - with no under-short partner, and the length rule keeps
        # those rings out of the aromatic set without poisoning the molecule).
        # Rejecting the broken ring from candidacy alone is not enough: the
        # length rule would rebuild it as a valid-looking ring in the output
        # mol (stamping the 1.16 A C-N as a ring-internal triple), silently
        # beating the OpenBabel fallback the corpus reserves for exactly these
        # S-heterocycle embeds.
        broken_ring = False
        if len(rings_info) == 1:
            rl0 = list(next(iter(rings_info)))
            if len(rl0) >= 5 and _planarity_rms([coords[i] for i in rl0]) <= 0.12:
                over_hit = under_hit = False
                for k in range(len(rl0)):
                    i, j = rl0[k], rl0[(k + 1) % len(rl0)]
                    d = self.blen[_sym_pair(i, j)]
                    lo, hi = _AROMATIC_ENVELOPE.get(
                        _sym_pair(elements[i], elements[j]), (1.27, 1.50)
                    )
                    if d > hi + _AROMATIC_SLACK_BOND:
                        over_hit = True
                    if d < lo - _AROMATIC_SLACK_SHORT:
                        under_hit = True
                if over_hit and under_hit:
                    broken_ring = True
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
                        dlen = self.blen[_sym_pair(a, j)]
                        pair = _sym_pair(elements[a], elements[j])
                        if pair in _TRIPLE_BOND_TABLE and dlen <= _TRIPLE_BOND_TABLE[pair]:
                            continue  # a nitrile / alkyne substituent is a
                                      # single-atom triple, not a double
                        dmax = _BOND_ORDER_TABLE.get(pair, (None, None))[1]
                        if dmax is None or dlen > dmax:
                            continue
                        # an exocyclic C=N at 1.34-1.37 A is aniline-type, not an
                        # imine: the N's lone pair conjugates into the ring and
                        # the bond is a single (3QTX X43's 2-aminothiazole N at
                        # 1.353 A).  A genuine exo imine refines to 1.28-1.32 A,
                        # so an N exo needs the tight bound; the ring-internal
                        # C=N of a diene (the 8TO5 azadiene at 1.360) keeps the
                        # generous 1.37 because this gate never sees it.
                        if elements[j] == "N" and dlen > 1.33:
                            continue
                        vmax = _MAX_VALENCE.get(elements[j])
                        # S and P exo bonds never block: the valence test is
                        # vacuous for them (vmax 6/5 absorbs any double), and
                        # a thioether C-S refines to 1.68-1.76 A inside the
                        # C=S cutoff - 7FOY W5C's 1.687 A C-S blocked a
                        # textbook planar triazole.  Leave S/P exo bonds to
                        # the length rule and the demotion pass.
                        if vmax is None or elements[j] in ("S", "P"):
                            continue
                        val = 2.0  # this bond as a double
                        for k in graph_nbrs_all.get(j, ()):
                            if k == a:
                                continue
                            dk = self.blen[_sym_pair(j, k)]
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
                d = self.blen[_sym_pair(i, j)]
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
            # a ring edge longer than any real aromatic bond of its pair is not
            # aromatic: the envelope's slack exists for low-res benzenes (C-C
            # elongated to ~1.52), but must not admit an sp3 ring riding one
            # "good" edge (5KYA 6Y4's pyrrolidine passes every gate with C-N
            # 1.471/1.488 and 6 pi - the longest aromatic C-N in the tuning set
            # is 1XQS AMP's 1.450).  C-C and C-O are deliberately uncapped
            # (07L's coumarin lactone ring C-O refines to 1.462).
            for k in range(len(rl)):
                i, j = rl[k], rl[(k + 1) % len(rl)]
                hard = _AROMATIC_HARD_MAX.get(
                    _sym_pair(elements[i], elements[j])
                )
                if hard is not None and self.blen[_sym_pair(i, j)] >= hard:
                    ok = False
                    break
            if not ok:
                continue
            # a ring whose lengths need slack is a low-resolution aromatic
            # candidate: demand better planarity than the generous 0.12, or a
            # saturated ether ring (1T3R/017 THF: C-O/C-C at 1.40-1.51, rms
            # 0.115) is let through - 07L's coumarin at rms 0.053 qualifies,
            # its two C-O at 1.444/1.462 and one C-C at 1.514 riding the slack
            # (each at most +0.022, within the per-bond cap).
            if excess > 0.0 and rms > 0.08:
                continue
            # a 5-ring whose pi count leans on an O/S lone pair (furan /
            # thiophene-type 2 pi) must carry at least one genuine
            # double-length edge not incident to it: V0X's
            # 2,3-dihydrobenzofuran 5-ring scores O 2pi + 4C = 6pi and
            # passes Huckel, but every C-C edge is a 1.47-1.52 A single -
            # there is no conjugation for the count to reward, and marking
            # the ring aromatic strips real hydrogens (BT5's saturated
            # THF-fused ring is the same signature).  Pure-carbon and
            # pyridine-type rings pass unchecked: low-resolution benzenes
            # refine to uniform 1.45-1.52 A edges.  6-rings are exempt (a
            # saturated 6-ring needs 4n+2 only with two heteroatoms, and
            # the coumarin pyranone at 1.416 A C=C must survive).
            if len(rl) == 5 and any(elements[b] in ("O", "S") for b in rl):
                short_ok = False
                for k in range(len(rl)):
                    i, j = rl[k], rl[(k + 1) % len(rl)]
                    pair = _sym_pair(elements[i], elements[j])
                    if elements[i] == "O" or elements[j] == "O":
                        # an aryl C-O is short even in a saturated ring
                        # (phenol / aryl ether refines to 1.31-1.44), so O
                        # edges prove nothing about conjugation
                        continue
                    if elements[i] == "S" or elements[j] == "S":
                        # a thiophene-type C-S is genuinely short (1.58-1.76
                        # across the dataset; 4L9Q 9TP's thiophene at 2.6 A
                        # refines its C-C edges to 1.46-1.55 - too long for
                        # the C-C test below - but its C-S stays at
                        # 1.64-1.66), while a thioether ring C-S refines to
                        # 1.78+ (8J7D BTI's tetrahydrothiophene 1.81, 3UDI
                        # PNM's 1.84 - all of them non-planar anyway).  A
                        # short C-S edge is real conjugation evidence.
                        dmax = _BOND_ORDER_TABLE.get(pair, (None, None))[1]
                        if (dmax is not None
                                and self.blen[_sym_pair(i, j)] <= dmax + _AROMATIC_SLACK_SHORT):
                            short_ok = True
                            break
                        continue
                    # fusion edge: shared with another ring, it proves nothing
                    # about this ring's own conjugation (V0X's
                    # dihydrobenzofuran borrows the benzene's 1.356 A C4-C5)
                    if edge_ring_count.get(_sym_pair(i, j), 0) > 1:
                        continue
                    dmax = _BOND_ORDER_TABLE.get(pair, (None, None))[1]
                    if (dmax is not None
                            and self.blen[_sym_pair(i, j)] <= dmax + _AROMATIC_SLACK_BOND):
                        short_ok = True
                        break
                if not short_ok:
                    continue
            ring_excess[tuple(rl)] = excess
            self.cand_rings.append(rl)

        if broken_ring:
            return None, ("broken ring: the molecule's only ring has an edge "
                          "beyond the aromatic envelope plus slack and an edge "
                          "below it - the ring does not hold together in these "
                          "coordinates, refusing the molecule so the caller "
                          "can fall back to OpenBabel")

        cand_atoms = set()
        for rl in self.cand_rings:
            cand_atoms.update(rl)
        self.cand_atoms = cand_atoms

        # ring sigma from candidate-ring EDGES only: an exo substituent that is
        # itself a candidate ring (triazolyl-phenyl) must not inflate an atom's
        # ring sigma - the N-C(phenyl) bond is NOT part of the triazole
        ring_edges: set[tuple[int, int]] = set()
        for rl in self.cand_rings:
            for k in range(len(rl)):
                i, j = rl[k], rl[(k + 1) % len(rl)]
                ring_edges.add(_sym_pair(i, j))
        self.ring_sigma: dict[int, int] = {a: 0 for a in cand_atoms}
        for i, j in graph:
            if _sym_pair(i, j) in ring_edges:
                self.ring_sigma[i] += 1
                self.ring_sigma[j] += 1

        # total heavy-atom degree of every atom, and the exo neighbours of each
        # candidate atom (any graph neighbour not joined via a candidate ring
        # edge).  Used by the pi count for ring carbons bearing a
        # directly-attached terminal O (carbonyl): their p orbital is in the C=O
        # pi bond.
        self.deg: dict[int, int] = {}
        for a, b in graph:
            self.deg[a] = self.deg.get(a, 0) + 1
            self.deg[b] = self.deg.get(b, 0) + 1
        self.exo_nbrs: dict[int, list[int]] = {a: [] for a in cand_atoms}
        self.ring_nbrs: dict[int, list[int]] = {a: [] for a in cand_atoms}
        self.graph_nbrs: dict[int, list[int]] = {a: [] for a in range(n)}
        for a, b in graph:
            self.graph_nbrs[a].append(b)
            self.graph_nbrs[b].append(a)
            if _sym_pair(a, b) in ring_edges:
                self.ring_nbrs[a].append(b)
                self.ring_nbrs[b].append(a)
                continue
            if a in self.exo_nbrs:
                self.exo_nbrs[a].append(b)
            if b in self.exo_nbrs:
                self.exo_nbrs[b].append(a)

        # functional-group fixup engine (the carbonyl / amidine / nitro /
        # phosphate / sulfonamide rulebook, see bope.fixups): constructed
        # once from the static graph data, applied per assembly pass with
        # that pass's aromatic set.
        self._fixup_engine = FixupEngine(self.elements, self.graph, self.deg, self.blen)

        # pass 1/2 failed on the full candidate set: a ring admitted only
        # through noise-level slack (1D1's pyridinone ring at +0.003, fused
        # so its forced double makes the neighbour ring unkekulizable) can be
        # the odd one out - drop slack rings at or below the drop line and
        # retry once.  True low-resolution aromatics (07L +0.041, NDP +0.012)
        # sit well above it.
        mol, err = self._attempt(self.cand_rings)
        if mol is None:
            drop = [rl for rl in self.cand_rings
                    if 0.0 < ring_excess.get(tuple(rl), 0.0) <= _AROMATIC_SLACK_DROP]
            if drop:
                mol, err = self._attempt([rl for rl in self.cand_rings if rl not in drop])
        return mol, err

    def _rdk_el(self, el: str) -> str:
        return el.capitalize() if len(el) == 2 else el

    def _exo_sigma(self, a: int) -> int:
        return sum(1 for i, j in self.graph if a in (i, j)) - self.ring_sigma[a]

    def _carbonyl_c(self, c: int, dmax: float | None = None) -> bool:
        """True for a C with a directly-attached terminal O at C=O length.
        The default bound is the strict 1.30 (length-rule double cutoff):
        its p orbital sits in the C=O pi bond, so such a ring carbon
        contributes 0 pi.  The N-pyridone/amide discriminator passes the
        generous crystal-carbonyl bound (1.40) - a caffeine C=O refines to
        1.34-1.36 A and its N-methyl is amide-type regardless."""
        if self.elements[c] != "C":
            return False
        if dmax is None:
            dmax = _BOND_ORDER_TABLE[("C", "O")][1]
        for n in self.graph_nbrs[c]:
            if self.elements[n] == "O" and self.deg[n] == 1 and self.blen[_sym_pair(c, n)] <= dmax:
                return True
        return False

    def _pi_of(self, a: int, h_choice: bool) -> int | None:
        """pi electrons contributed by candidate atom *a*."""
        el = self.elements[a]
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
            return 0 if self._carbonyl_c(a) else 1
        if el == "N":
            rs, exo = self.ring_sigma[a], self._exo_sigma(a)
            if rs == 2 and exo == 1:
                # 5-ring: N-methyl/aryl-pyrrole, lone pair in the ring
                # plane-free p orbital: 2 pi.  6-ring with a plain N-alkyl:
                # pyridinium, the alkyl consumes the lone pair: 1 pi.  An
                # N-alkyl 6-ring N bonded to (or exo to) a carbonyl is
                # pyridone/amide-type: lone pair delocalised, 2 pi.
                # Without the 6-ring rule the nicotinamide rings of
                # NAD/NAP/NDP score 7 pi, fail Huckel and come out
                # saturated (+6 H in the formula).
                # NOTE: reads self.cand_rings - the FULL candidate list, never the
                # attempt's filtered subset: the pi count is defined over all rings.
                if max((len(rl) for rl in self.cand_rings if a in rl), default=0) <= 5:
                    # 5-ring with an exo substituent: the N is sp2 with
                    # three sigma bonds, so its lone pair sits in the
                    # ring p orbital: 2 pi (pyrrole-type: N-methylpyrrole,
                    # N-alkyl imidazole / triazole, the N-glycosides of
                    # nucleosides).  A neutral 5-ring N contributes 1 pi
                    # only when the exo bond is itself a double - an
                    # imine N=X, detected by the length rule.  The old
                    # discriminator (any other heteroatom in the ring,
                    # no N-N bond) mislabeled 7-methylguanine N7-CH3 and
                    # triazole N1-R: at 1 pi the N must double-bond and
                    # either strangles the ring (no kekule: 7FOY W5C,
                    # 5MUY MGT) or turns a 7-pi neutral ring aromatic
                    # (an N-alkyl thiazole whose neutral form has no
                    # conjugation to reward).
                    for e in self.exo_nbrs[a]:
                        # an exo carbon that is itself ring-bound is an
                        # N-aryl (5KYA 6Y4's pyrazole N-phenyl at 1.327 A
                        # sits inside the imine window but is a plain
                        # single bond): the aryl ring's own candidate
                        # status decides it.  A genuine exo imine carbon
                        # is never in a candidate ring.
                        if self.elements[e] == "C" and any(
                            e in rl for rl in self.cand_rings
                        ):
                            continue
                        pair = _sym_pair(self.elements[a], self.elements[e])
                        dmax = _BOND_ORDER_TABLE.get(pair, (None, None))[1]
                        if dmax is None:
                            continue
                        # an exocyclic C=N at 1.34-1.37 A is an N-aryl or
                        # N-alkyl aniline-type single (5AEP QUP's pyrrole
                        # N-aryl at 1.36-1.37 A), not an imine: a genuine
                        # exo imine refines to 1.28-1.32 A, so a C exo
                        # needs the tight bound (non-C exo N=O / N=S keep
                        # the table cutoff).
                        if self.elements[e] == "C" and self.blen[_sym_pair(a, e)] > 1.33:
                            continue
                        if self.blen[_sym_pair(a, e)] <= dmax:
                            return 1  # exo imine: the p pair is spent
                    return 2
                if any(self._carbonyl_c(j, _CRYSTAL_CARBONYL) for j in self.ring_nbrs[a]):
                    return 2
                if any(self._carbonyl_c(e, _CRYSTAL_CARBONYL) for e in self.exo_nbrs[a]):
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
                for rl in self.cand_rings:
                    if a in rl and len(rl) > 5:
                        others = [b for b in rl if self.elements[b] == "N" and b != a]
                        if not others:
                            return 1  # exactly one N: plain pyridinium
                        if any(self._exo_sigma(b) > 0 for b in others):
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
            return 2 if self.ring_sigma[a] == 2 else None
        return None

    def _pi_val(self, a: int, mask: int) -> int | None:
        return self._pi_of(a, (mask >> self._flex_idx[a]) & 1 if a in self._flex_idx else 0)

    def _huckel_subset(self, mask: int, rings: list[int]) -> bool:
        """Exact RDKit fused-ring rule (applyHuckelToFused): the union of
        the subset's ring atoms - each atom present in 1-2 rings of the
        subset counted once, atoms in 3+ rings excluded (they have no
        p-orbital contribution, the acepentalene fix) - must be a Huckel
        4n+2 pi set.  The subset must also be fused (connected), checked
        by the caller."""
        cnt: dict[int, int] = {}
        for ri in rings:
            for a in self._cur_cand_rings[ri]:
                cnt[a] = cnt.get(a, 0) + 1
        unon = [a for a, c in cnt.items() if c in (1, 2)]
        if len(unon) < 3:
            return False
        tot = 0
        for a in unon:
            p = self._pi_val(a, mask)
            if p is None:
                return False  # an sp3-type N poisons the subset, not a
                              # 0-pi atom (1D1's lactam alkyl N)
            tot += p
        return tot in _HUCKEL

    def _subset_fused(self, rings: list[int]) -> bool:
        """checkFused: the subset's rings must form one connected
        component via shared edges (>= 2 shared atoms)."""
        n = len(rings)
        if n <= 1:
            return True
        parent = list(range(n))

        for i in range(n):
            for j in range(i):
                if len(set(self._cur_cand_rings[rings[i]]) & set(self._cur_cand_rings[rings[j]])) >= 2:
                    ri, rj = _find(parent, i), _find(parent, j)
                    if ri != rj:
                        parent[rj] = ri
        return len({_find(parent, i) for i in range(n)}) == 1

    def _arom_for(self, mask: int) -> tuple[set[int], set[int]]:
        """(aromatic ring indices, per-ring aromatic ring indices) per
        the exact RDKit fused-subset rule: a ring is aromatic iff it
        belongs to a fused subset (sizes 1..min(n,6)) whose 1-2-ring
        atom union is Huckel.  Trying subsets in size order and marking
        every ring of a passing subset reproduces the C++ loop."""
        arom: set[int] = set()
        per_ring: set[int] = set()
        for rings in self._sys_rings.values():
            n = len(rings)
            for size in range(1, min(n, 6) + 1):
                for comb in combinations(rings, size):
                    if not self._subset_fused(comb):
                        continue
                    if size > 1 and not all(self._rescue_eligible(ri) for ri in comb):
                        continue
                    if not self._huckel_subset(mask, comb):
                        continue
                    for ri in comb:
                        arom.add(ri)
                    if size == 1:
                        per_ring.add(comb[0])
        return arom, per_ring

    def _rescue_eligible(self, ri: int) -> bool:
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
            if (_planarity_rms([self.coords[i] for i in self._cur_cand_rings[ri]])
                    > _AROMATIC_RESCUE_RMS_MAX):
                rl = self._cur_cand_rings[ri]
                for k in range(len(rl)):
                    i, j = rl[k], rl[(k + 1) % len(rl)]
                    if self.blen[_sym_pair(i, j)] <= _RESCUE_RING_BOND_MAX:
                        continue
                    if (self.elements[i] == "C" and self._carbonyl_c(i)) or (
                            self.elements[j] == "C" and self._carbonyl_c(j)):
                        continue
                    return False
            rset = set(self._cur_cand_rings[ri])
            for a in self._cur_cand_rings[ri]:
                if self.elements[a] != "C" or not self._carbonyl_c(a):
                    continue
                rnbrs = [b for b in self.graph_nbrs[a] if b in rset]
                if len(rnbrs) == 2 and all(self.elements[b] == "C" for b in rnbrs):
                    return False
            return True

    def _score(self, mask: int) -> tuple[int, int, int]:
        """(aromatic rings, per-ring aromatic rings, amide-H count).
        The total counts rings aromatic per-ring OR via any passing
        fused subset (isoalloxazine's pyrazine is 7 pi per-ring and
        its uracil 5, yet the {uracil,pyrazine} union is 10 and the
        tricycle 14 - both 4n+2).  A mask whose rings are ALL
        rescue-only (no ring stands alone on its own Huckel count)
        scores 0: the rescue is then an accident of the mask
        enumeration - 7-methylguanine's pyridine tautomer puts the
        6-ring at 5 pi and the {5-ring,6-ring} union lands on 10,
        marking the 7-pi imidazole aromatic and making the mol
        unkekulizable, where the N1-H tautomer stands alone at 6 pi
        and the imidazole stays non-aromatic, exactly as the ccd
        draws it.  The flavin's rescue is never lost: its uracil /
        pyrazine pair is rescued alongside the per-ring benzene, so
        the mask with the N3-H tautomer (uracil 5 pi, pyrazine 7 pi,
        3 rings) beats the both-H mask (uracil 6 pi per-ring but the
        pyrazine dropped, 2 rings) on the total itself.  The
        per-ring count then breaks ties toward the mask that needs
        no system rescue (purine's N7-H tautomer over the
        pyrimidine-N1-H tautomer), and the amide-H count breaks
        flavin-type ties: the pyrrole-H prefers the N with the most
        carbonyl ring neighbours (flavin N3 between the two C=O's,
        matching the ccd, over N1)."""
        arom, per_ring = self._arom_for(mask)
        ok = len(arom) if per_ring else 0
        per_ring_n = len(per_ring)
        amide_h = 0
        for a in self._flex_n:
            if (mask >> self._flex_idx[a]) & 1:
                amide_h += sum(
                    1 for j in self.ring_nbrs[a]
                    if self._carbonyl_c(j, _CRYSTAL_CARBONYL)
                )
        return ok, per_ring_n, amide_h

    def _attempt(self, cand_rings: list[list[int]]) -> tuple[Any | None, str | None]:
        elements, cand_atoms = self.elements, self.cand_atoms
        self._cur_cand_rings = cand_rings  # per-attempt view: the slack-drop retry passes a filtered subset
        # flexible N's: 2 ring sigma, no exo (pyrrole-vs-pyridine tautomer bit)
        flex_n = [
            a for a in sorted(cand_atoms)
            if elements[a] == "N" and self.ring_sigma[a] == 2 and self._exo_sigma(a) == 0
        ]
        flex_idx = {a: i for i, a in enumerate(flex_n)}
        self._flex_n = flex_n
        self._flex_idx = flex_idx

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
        self._sys_rings = sys_rings

        self._best_mask, best_score = None, (-1, -1, -1)
        for mask in range(1 << len(flex_n)):
            s = self._score(mask)
            if s > best_score:
                best_score, self._best_mask = s, mask

        arom, _ = self._arom_for(self._best_mask)
        arom_rings = [cand_rings[ri] for ri in sorted(arom)]
        arom_atoms = set()
        for rl in arom_rings:
            arom_atoms.update(rl)
        self._arom_atoms = arom_atoms
        arom_bonds = set()
        for rl in arom_rings:
            for k in range(len(rl)):
                i, j = rl[k], rl[(k + 1) % len(rl)]
                arom_bonds.add(_sym_pair(i, j))
        self._arom_bonds = arom_bonds

        # pass 1: keep length-rule doubles on substituted N/S (imines, triazine
        # C=N).  If RDKit cannot kekulize (a length double on a benzene-ring
        # carbon, aniline-type), pass 2 forces every N/S exo single.  Pass 3
        # drops the aromatic set entirely (pure length-rule molecule, demotion
        # pass still applied) so perception never returns None on valid input:
        # a kekulization failure is a wrong aromatic call, not a reason to
        # report no molecule at all.
        mol, err = self._assemble(exo_force_all=False)
        if mol is None:
            mol, err = self._assemble(exo_force_all=True)
        if mol is None:
            mol, err = self._assemble(
                exo_force_all=False,
                arom_atoms_arg=set(),
                arom_bonds_arg=set(),
            )
        if mol is None:
            return None, err
        return mol, None

    def _assemble(
        self,
        exo_force_all: bool,
        arom_atoms_arg: set[int] | None = None,
        arom_bonds_arg: set[tuple[int, int]] | None = None,
    ) -> tuple[Any | None, str | None]:
        """Build the mol with bond orders.  exo_force_all=False: only degree-1
        N/S exo to an aromatic atom (NH2 / amino-pyridine) is forced single;
        substituted N/S (C=N imines, N-aryl triazoles) stay on the length rule
        so a triazine C=N survives.  exo_force_all=True: every exo bond of an
        aromatic atom is single, N/S and C-C alike (the aniline fallback,
        generalised - fixes kekulization when the length rule puts a double
        on the direct exo bond of a ring atom: 7RS8 7EI's phenol ring carbons
        hold 1.348 / 1.353 A aryl-alkenyl singles inside the 1.38 C=C cutoff,
        and such a double leaves the ring atom no pi for its ring double).
        Passing empty aromatic sets builds a pure length-rule molecule (the
        final never-return-None fallback)."""
        a_atoms = self._arom_atoms if arom_atoms_arg is None else arom_atoms_arg
        a_bonds = self._arom_bonds if arom_bonds_arg is None else arom_bonds_arg
        m = Chem.RWMol(self.rw)  # type: ignore[attr-defined]
        for i, j in self.graph:
            if _sym_pair(i, j) in a_bonds:
                m.AddBond(i, j, Chem.BondType.AROMATIC)  # type: ignore[attr-defined]
            else:
                d = self.blen[_sym_pair(i, j)]
                pair = _sym_pair(self.elements[i], self.elements[j])
                if pair in _TRIPLE_BOND_TABLE and d <= _TRIPLE_BOND_TABLE[pair]:
                    # nitrile / alkyne: unmistakably short, checked first so an
                    # N/S exo triple on an aromatic ring (benzonitrile) survives
                    m.AddBond(i, j, Chem.BondType.TRIPLE)  # type: ignore[attr-defined]
                elif (i in a_atoms) != (j in a_atoms):
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
                    narm = i if i not in a_atoms else j
                    arom_nbrs = sum(
                        1 for a, b in self.graph
                        if (a == narm and b in a_atoms)
                        or (b == narm and a in a_atoms)
                    )
                    if arom_nbrs >= 2:
                        m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                    elif exo_force_all:
                        # pass 2: no double on the direct exo bond of an
                        # aromatic atom - the ring atom's pi belongs to the
                        # ring (7RS8 7EI, 5KYA 6Y4's alkene-bearing rings).
                        m.AddBond(i, j, Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                    elif self.elements[i] in ("N", "S") or self.elements[j] in ("N", "S"):
                        ns = i if self.elements[i] in ("N", "S") else j
                        deg_ns = sum(1 for a, b in self.graph if a == ns or b == ns)
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
                elif i in a_atoms:
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
        # The functional-group rulebook (carbonyl rescue + ester-O
        # protection, amidine, nitro, phosphate, sulfonamide) lives in
        # bope.fixups as data; the engine runs it in order on the molecule
        # as it stands.
        self._fixup_engine.apply(m, a_atoms)
        self._remove_phantom_edges(m, a_atoms)
        self._mark_aromatic_pyrrole_h(m, a_atoms)
        self._set_pyridinium(m, a_atoms)
        mol = m.GetMol()

        self._demote_valences(mol)
        try:
            Chem.SanitizeMol(mol)  # type: ignore[attr-defined]
        except Exception as exc:  # noqa: BLE001
            return None, f"santize: {exc}"
        return mol, None

    def _remove_phantom_edges(self, m: Chem.RWMol, a_atoms: set[int]) -> None:  # type: ignore[attr-defined]
        # phantom-edge removal: an atom over-valent with only single
        #     bonds (plus nothing reducible - no double, no triple, not
        #     the quaternary-N case the demotion pass charges) carries a
        #     spurious graph edge.  The bond graph's 0.40-A tolerance
        #     admits 1.69-1.78 A C-O edges (4L9Q 9TP O38-C39, 5GHV X5G
        #     O12-C10) that would make O trivalent and crash sanitize;
        #     the phantom is always the longest bond on the atom, so
        #     drop it.  Aromatic atoms are excluded - their exo singles
        #     are load-bearing (the pyridinium charge and pyrrole-N-H
        #     logic below runs on them), and their over-valence is
        #     handled by the aromatic rules.
        while True:
            dropped = False
            for i in range(self.n):
                if i in a_atoms:
                    continue
                vmax = _MAX_VALENCE.get(self.elements[i])
                if vmax is None:
                    continue
                atom = m.GetAtomWithIdx(i)  # type: ignore[attr-defined]
                if atom.GetFormalCharge() != 0:  # type: ignore[attr-defined]
                    continue
                val = 0.0
                singles = []
                has_multi = False
                for b in atom.GetBonds():  # type: ignore[attr-defined]
                    bt = b.GetBondType()
                    if bt == Chem.BondType.SINGLE:
                        val += 1.0
                        singles.append(b)
                    elif bt == Chem.BondType.AROMATIC:
                        val += 1.5
                    else:
                        has_multi = True
                        val += b.GetBondTypeAsDouble()
                if val <= vmax or has_multi or not singles:
                    continue
                # quaternary ammonium: leave to the demotion pass,
                # which knows the charge is the resolution
                if self.elements[i] == "N" and len(singles) == 4 and val == 4:
                    continue
                longest = max(
                    singles,
                    key=lambda b: self.blen[  # type: ignore[arg-type]
                        _sym_pair(b.GetBeginAtomIdx(),  # type: ignore[attr-defined]
                                  b.GetEndAtomIdx())  # type: ignore[attr-defined]
                    ],
                )
                m.RemoveBond(longest.GetBeginAtomIdx(),  # type: ignore[attr-defined]
                             longest.GetEndAtomIdx())  # type: ignore[attr-defined]
                dropped = True
                break
            if not dropped:
                break

    def _mark_aromatic_pyrrole_h(self, m: Chem.RWMol, a_atoms: set[int]) -> None:  # type: ignore[attr-defined]
        # aromatic atoms: flags + pyrrole-N explicit H
        for a in a_atoms:
            at = m.GetAtomWithIdx(a)  # type: ignore[attr-defined]
            at.SetIsAromatic(True)  # type: ignore[attr-defined]
            at.SetHybridization(Chem.HybridizationType.SP2)  # type: ignore[attr-defined]
            if (
                self.elements[a] == "N" and self.ring_sigma[a] == 2 and self._exo_sigma(a) == 0
                and ((self._best_mask >> self._flex_idx[a]) & 1)
            ):
                at.SetNoImplicit(True)  # type: ignore[attr-defined]
                at.SetNumExplicitHs(1)  # type: ignore[attr-defined]

    def _set_pyridinium(self, m: Chem.RWMol, a_atoms: set[int]) -> None:  # type: ignore[attr-defined]
        # pyridinium: a 6-ring aromatic N with a heavy exo substituent and no
        # carbonyl adjacency contributes 1 pi and cannot take a ring double
        # (its lone pair is consumed by the exo bond): the neutral form is
        # unkekulizable, the +1 charge is forced by valence.  This is exactly
        # the nicotinamide N of NAD/NAP/NDP (and any N-alkyl pyridinium).
        # N-alkyl 5-ring pyrroles (2 pi) and pyridone/amide N's (2 pi) stay
        # neutral - the same tests as the pi-count's pyridinium branch.
        for a in a_atoms:
            if (
                self.elements[a] == "N"
                and self.ring_sigma[a] == 2
                and self._exo_sigma(a) == 1
                and max((len(rl) for rl in self._cur_cand_rings if a in rl), default=0) > 5
                and all(
                    sum(1 for b in rl if self.elements[b] == "N") == 1
                    for rl in self._cur_cand_rings if a in rl and len(rl) > 5
                )
                and not any(self._carbonyl_c(j, _CRYSTAL_CARBONYL) for j in self.ring_nbrs[a])
                and not any(self._carbonyl_c(e, _CRYSTAL_CARBONYL) for e in self.exo_nbrs[a])
            ):
                m.GetAtomWithIdx(a).SetFormalCharge(1)  # type: ignore[attr-defined]

    def _demote_valences(self, mol: Chem.Mol) -> None:  # type: ignore[attr-defined]
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
                        if nbr.GetSymbol() == "N" and not nbr.GetIsAromatic():
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
                        elif (nbr.GetSymbol() == "C"
                              and not nbr.GetIsAromatic()
                              # the C=C analog: a noisy C-C single refined
                              # short (4L9Q 9TP's lactone CH2-C2 at 1.286 A
                              # vs the 1.52 A true single) lands inside the
                              # 1.38 double cutoff and gives an
                              # exactly-valent sp3-looking carbon a spurious
                              # C=C.  A real alkene carbon never carries two
                              # saturated carbon singles AND a double
                              # partner with a non-C single (that partner is
                              # an enol-ether/ester carbon); tetrasubstituted
                              # alkenes have all-carbon partners and survive.
                              and any(
                                  b.GetBondType() == Chem.BondType.SINGLE  # type: ignore[attr-defined]
                                  and b.GetOtherAtom(nbr).GetSymbol() != "C"
                                  and not b.GetOtherAtom(nbr).GetIsAromatic()
                                  for b in nbr.GetBonds()
                              )
                              and len(cn_singles) >= 2
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
                        # no double to demote: a spurious nitrile /
                        # alkyne triple can be reduced stepwise (the
                        # 7FOZ WD0 / 2QES ADE class: an N or C with a
                        # phantom triple).  With nothing reducible the
                        # atom stays over-valent and sanitize decides;
                        # attempt()'s slack-ring and length-rule
                        # fallbacks give the molecule every chance
                        # before a final failure is reported.
                        tb = [
                            b for b in a.GetBonds()
                            if b.GetBondType() == Chem.BondType.TRIPLE
                        ]
                        if tb:
                            tlong = max(
                                tb, key=lambda b: self.blen[
                                    _sym_pair(b.GetBeginAtomIdx(),
                                              b.GetEndAtomIdx())
                                ]
                            )
                            mol.GetBondBetweenAtoms(  # type: ignore[attr-defined]
                                tlong.GetBeginAtomIdx(),
                                tlong.GetEndAtomIdx()
                            ).SetBondType(Chem.BondType.DOUBLE)  # type: ignore[attr-defined]
                            changed = True
                            break
                        continue
                    # prefer demoting a C=N double over a C=C double
                    # when the C=N is at amide-plausible length: amide
                    # / aniline C-N singles refine to 1.31-1.37 A
                    # inside the C=N cutoff, while a C-C single never
                    # refines below ~1.45 - the C=C is the stronger
                    # double signal (8JZ7 FI7: C=N 1.327 vs C=C
                    # 1.361; keeping the C=N and demoting the C=C
                    # put a double at 1.327 A that cannot be).
                    cn = [
                        b for b in db
                        if {self.elements[b.GetBeginAtomIdx()],
                            self.elements[b.GetEndAtomIdx()]} == {"C", "N"}
                        and self.blen[_sym_pair(b.GetBeginAtomIdx(),
                                           b.GetEndAtomIdx())] >= 1.30
                    ]
                    if cn and len(db) > 1:
                        longest = max(
                            cn, key=lambda b: self.blen[
                                _sym_pair(b.GetBeginAtomIdx(),
                                          b.GetEndAtomIdx())
                            ]
                        )
                    else:
                        longest = max(
                            db, key=lambda b: self.blen[
                                _sym_pair(b.GetBeginAtomIdx(),
                                          b.GetEndAtomIdx())
                            ]
                        )
                    mol.GetBondBetweenAtoms(  # type: ignore[attr-defined]
                        longest.GetBeginAtomIdx(), longest.GetEndAtomIdx()
                    ).SetBondType(Chem.BondType.SINGLE)  # type: ignore[attr-defined]
                    changed = True
                    break


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

    return _GeometricPerceiver(elements, coords, graph).perceive()

