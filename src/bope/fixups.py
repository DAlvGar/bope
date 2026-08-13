"""Functional-group fixup rules for the geometry perception.

The length-rule bond orders assembled by the geometry perception are
corrected by a small rulebook of functional-group fixups: carbonyl
rescue + ester-O protection, amidine/imine, nitro, phosphate P=O and
sulfonamide.  Each rule is declarative data - a trigger (center
element, neighbor composition, measured bond lengths) and an action
(bond-order and charge assignments) - evaluated by the generic
:class:`FixupEngine`.

The rules run in order and each acts on the molecule as it stands, so a
later rule sees an earlier rule's edits - exactly the semantics of the
hand-written passes they replaced.  Adding a functional group is one
data entry in :data:`FIXUP_RULES` plus a test; no new code path.

The engine is constructed once per perception from the static graph
data and applied per assembly pass with that pass's perceived
aromatic-atom set, so a rulebook run reflects the assembly pass that
produced the molecule.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from bope._deps import Chem
from bope.helpers import _sym_pair


@dataclass(frozen=True)
class NbrGroup:
    """A named neighbor set of a rule's center atom.

    Members are selected from the static graph (elements, degrees,
    perceived aromaticity); bond lengths are looked up per pair at
    evaluation time.  ``min`` / ``max`` / ``exact`` bound the member
    count; ``exact`` is mutually exclusive with the other two.
    """

    name: str
    element: str = "*"              # element symbol, or "*" for any element
    aromatic: bool | None = None    # None=any, True=member of a_atoms, False=not
    terminal: bool | None = None    # None=any degree, True=degree 1, False=degree > 1
    min: int = 0
    max: int | None = None          # None = unbounded
    exact: int | None = None


@dataclass(frozen=True)
class FixupRule:
    """One functional-group correction.

    The engine iterates the atoms matching ``center`` and ``aromatic``,
    builds the named neighbor groups from ``groups``, evaluates the gate
    predicates (all must hold for the rule to fire), then applies the
    actions in order: ``make_single`` bonds, ``make_double`` bonds,
    ``charges``.

    Gate predicates:

    - group count bounds - a group's ``min`` / ``max`` / ``exact`` bound
      is enforced through the ``require`` / ``require_or`` clauses that
      reference it (a clause passes only at the bounded count); a group
      outside both clauses is bounded unconditionally.
    - ``require`` - named groups that must be non-empty (at their
      bounded count).
    - ``require_or`` - clauses of group names; at least one group per
      clause must be non-empty (at its bounded count).
    - ``exclude_nbrs`` - elements that must not appear among ANY
      neighbor of the center.
    - ``max_len`` - every member of the named group must sit within the
      cutoff (gate semantics: the rule is skipped otherwise).
    - ``no_double_to`` - the rule is skipped if any center bond to
      these groups is already DOUBLE (a molecule-state guard).

    Actions:

    - ``make_double`` - ``"shortest:<group>"`` or ``"all:<group>"``.
      ``action_len`` filters the members before the selection (filter
      semantics: members beyond the cutoff are left untouched).
    - ``make_single`` - ``"group:<name>"`` (all members),
      ``"<group>_others"`` (members minus the double target), or
      ``"non_single_non_arom"`` (the center's non-single, non-aromatic
      bonds).
    - ``charges`` - ``("center", n)`` or ``("<group>_others", n)``.
    - ``only_if_single`` - ``make_double`` touches only
      currently-SINGLE bonds.
    """

    name: str
    center: str
    note: str = ""  # crystal citations + reasoning for the thresholds
    aromatic: bool = False
    min_degree: int = 0
    groups: tuple[NbrGroup, ...] = ()
    require: tuple[str, ...] = ()
    require_or: tuple[tuple[str, ...], ...] = ()
    exclude_nbrs: tuple[str, ...] = ()
    max_len: dict[str, float] = field(default_factory=dict)
    action_len: dict[str, float] = field(default_factory=dict)
    no_double_to: tuple[str, ...] = ()
    make_double: str = ""
    make_single: tuple[str, ...] = ()
    only_if_single: bool = False
    charges: tuple[tuple[str, int], ...] = ()


# ---------------------------------------------------------------------------
# The rulebook.  Order is load-bearing: rules run in order, each on the
# molecule as the previous rule left it (identical to the pass sequence
# they replace).  The notes carry the crystal citations that justify the
# thresholds - they are the scientific record, keep them.
# ---------------------------------------------------------------------------
FIXUP_RULES: tuple[FixupRule, ...] = (
    FixupRule(
        name="carbonyl_rescue",
        center="C",
        groups=(
            NbrGroup("any_o", "O", min=1),
            NbrGroup("n", "N", min=1),
            NbrGroup("arom", aromatic=True, min=2),
            NbrGroup("term_o", "O", terminal=True),
        ),
        require=("any_o",),
        require_or=(("n", "arom"),),
        action_len={"term_o": 1.40},
        make_double="all:term_o",
        only_if_single=True,
        note="non-aromatic C with a terminal O single and at least one N "
             "neighbor, or two aromatic neighbors (a diaryl / aryl-heteroaryl "
             "ketone bridge - 3QTX X43's C(=O) at 1.358 A between the thiazole "
             "and the phenyl) -> C=O.  Crystal carbonyls often refine to "
             "1.34-1.36 A (caffeine C2/C6 in 3RFM); aromatic ring C's are "
             "excluded, and internal O's (ether/alcohol C-O) stay single.",
    ),
    FixupRule(
        name="ester_o_single",
        center="C",
        groups=(NbrGroup("bridged_o", "O", terminal=False, min=1),),
        require=("bridged_o",),
        make_single=("group:bridged_o",),
        note="a C-O double can never land on a bridged O (an O with a second "
             "heavy neighbour - ester, lactone, acyl phosphate): the O's "
             "valence 2 is spent on its two sigma bonds.  When a crystal "
             "refines the ester C-O shorter than the C=O (BT5/1WQW measures "
             "1.247 vs 1.272), the length rule doubles the wrong bond, the "
             "bridged O goes over-valent, and the demotion pass destroys the "
             "carbonyl.  Force bridged C-O single; the rescue above and the "
             "length rule keep the terminal O double.",
    ),
    FixupRule(
        name="amidine",
        center="C",
        groups=(NbrGroup("n", "N", exact=2),),
        require=("n",),
        exclude_nbrs=("O", "S"),
        max_len={"n": 1.40},
        make_double="shortest:n",
        only_if_single=True,
        note="amidine/imine: non-aromatic C with exactly 2 N neighbors and no "
             "O/S neighbor: the shorter C-N becomes double.  Delocalized "
             "amidines measure 1.31-1.33 for BOTH bonds; neutral benzamidine "
             "needs one.  1.40 (not 1.36): delocalized amidines / imines "
             "measure 1.31-1.36 (ETKDG +0.03 bias, plus noise), while a "
             "genuine C-N single pair sits at 1.45+ - the 0.05 gap keeps the "
             "raised threshold safe.",
    ),
    FixupRule(
        name="nitro",
        center="N",
        min_degree=3,
        groups=(NbrGroup("term_o", "O", terminal=True, exact=2),),
        require=("term_o",),
        max_len={"term_o": 1.45},
        make_double="shortest:term_o",
        make_single=("term_o_others", "non_single_non_arom"),
        charges=(("center", 1), ("term_o_others", -1)),
        note="nitro: N with two terminal O's at N-O <= 1.45 A (short nitro "
             "refines to 1.18-1.22, delocalised nitro to 1.36-1.43, 3B67 B67 "
             "rides 1.420/1.420 and 6SUH LVE 1.364/1.428) is the "
             "charge-separated [N+](=O)[O-] of the CCD record.  N vmax 3 "
             "cannot hold two N=O, so the length rule alone reads a single "
             "plus an O-H (6SUH LVE, 4EK8 16K at 1.214) and the demotion pass "
             "destroys even textbook-length nitro; the +1 charge exempts the N "
             "from the demotion pass and reproduces the CCD form exactly (the "
             "benchmark neutralises before comparing).  The only neutral "
             "mislabel the gate guards against is N(OH)2 (two N-O "
             "hydroxylamine singles at 1.40-1.47), which does not occur in "
             "crystal ligands.",
    ),
    FixupRule(
        name="phosphate",
        center="P",
        groups=(NbrGroup("term_o", "O", terminal=True, min=2),),
        require=("term_o",),
        no_double_to=("term_o",),
        make_double="shortest:term_o",
        note="phosphate P=O: P with 2+ terminal O's and no P=O double yet.  "
             "Crystal phosphate P-O refines to 1.55-1.70 A (1TPB PGH 1.701, "
             "2VF5 GLP 1.591, 2I22 I22 1.604), above the 1.55 P=O cutoff, so "
             "the length rule leaves every P-O single and P falls back to "
             "P-H; forcing the shortest P-O to double puts P at exactly "
             "valence 5, no charge needed.",
    ),
    FixupRule(
        name="sulfonamide",
        center="S",
        groups=(
            NbrGroup("term_o", "O", terminal=True, min=2),
            NbrGroup("n", "N", min=1),
        ),
        require=("term_o", "n"),
        make_double="all:term_o",
        make_single=("group:n",),
        note="sulfonamide S-N: an S with two terminal O's and an N neighbor is "
             "the sulfonyl of a sulfonamide R-S(=O)(=O)-NH2: both terminal O's "
             "take the double (crystal S=O refines to 1.44-1.64 A, above the "
             "1.55 cutoff - 3QTX X43's O at 1.639) and the S-N stays single "
             "(sulfonamide S-N measures 1.56-1.63 A, inside the (N,S) double "
             "cutoff 1.70, so the length rule misorders the pair).  A genuine "
             "S=N double occurs only in sulfoximines / sulfonimidamides, "
             "which carry a single S=O - never two terminal O's.",
    ),
)


class FixupEngine:
    """Evaluate :class:`FixupRule` instances against a partially
    assembled molecule.

    Constructed once per perception from the static graph data
    (``elements``, ``graph`` edges, per-atom ``deg``, per-pair
    ``blen``); :meth:`apply` is called per assembly pass with that
    pass's perceived aromatic-atom set.
    """

    def __init__(
        self,
        elements: list[str],
        graph: set[tuple[int, int]],
        deg: dict[int, int],
        blen: dict[tuple[int, int], float],
    ) -> None:
        self.elements = elements
        self.deg = deg
        self.blen = blen
        self.adj: list[list[int]] = [[] for _ in elements]
        for a, b in graph:
            self.adj[a].append(b)
            self.adj[b].append(a)

    def apply(
        self,
        m: Chem.RWMol,  # type: ignore[attr-defined]
        a_atoms: set[int],
        rules: tuple[FixupRule, ...] = FIXUP_RULES,
    ) -> list[str]:
        """Run ``rules`` on ``m``; return the names of the rules fired.

        ``m`` must be an RWMol with all graph edges added as bonds (the
        pre-fixup state of the assembly pass).  ``a_atoms`` is the
        pass's perceived aromatic-atom set.  The returned names let
        tests and diagnostics report exactly which corrections fired.
        """
        fired: list[str] = []
        for rule in rules:
            for i in range(len(self.elements)):
                if self.elements[i] != rule.center or (i in a_atoms) != rule.aromatic:
                    continue
                if len(self.adj[i]) < rule.min_degree:
                    continue
                groups = self._groups(i, a_atoms, rule)
                if not self._gate(i, rule, groups, m):
                    continue
                self._apply(i, rule, groups, m)
                fired.append(rule.name)
        return fired

    def _groups(
        self, i: int, a_atoms: set[int], rule: FixupRule
    ) -> dict[str, list[int]]:
        groups: dict[str, list[int]] = {}
        for g in rule.groups:
            members: list[int] = []
            for j in self.adj[i]:
                if g.element != "*" and self.elements[j] != g.element:
                    continue
                if g.aromatic is not None and ((j in a_atoms) != g.aromatic):
                    continue
                if g.terminal is True and self.deg[j] != 1:
                    continue
                if g.terminal is False and self.deg[j] <= 1:
                    continue
                members.append(j)
            groups[g.name] = members
        return groups

    @staticmethod
    def _count_ok(rule: FixupRule, groups: dict[str, list[int]], name: str) -> bool:
        """Non-empty AND the group's exact/min/max bounds hold."""
        g = next(g for g in rule.groups if g.name == name)
        cnt = len(groups[name])
        if cnt == 0:
            return False
        if g.exact is not None:
            return cnt == g.exact
        return cnt >= g.min and (g.max is None or cnt <= g.max)

    def _gate(
        self,
        i: int,
        rule: FixupRule,
        groups: dict[str, list[int]],
        m: Chem.RWMol,  # type: ignore[attr-defined]
    ) -> bool:
        # group count bounds are gate predicates, enforced through the
        # require / require_or clauses that reference the group (the
        # amidine's exactly-two-N gate, the nitro's exactly-two-O, the
        # phosphate's two-terminal-O floor, the carbonyl's two-aromatic
        # fallback trigger ...).  A group outside both clauses is still
        # bounded unconditionally.
        referenced = set(rule.require)
        for clause in rule.require_or:
            referenced.update(clause)
        for g in rule.groups:
            if g.name in referenced:
                continue
            cnt = len(groups[g.name])
            if g.exact is not None and cnt != g.exact:
                return False
            if cnt < g.min or (g.max is not None and cnt > g.max):
                return False
        for name in rule.require:
            if not self._count_ok(rule, groups, name):
                return False
        for clause in rule.require_or:
            if not any(self._count_ok(rule, groups, name) for name in clause):
                return False
        if any(self.elements[j] in rule.exclude_nbrs for j in self.adj[i]):
            return False
        for name, dmax in rule.max_len.items():
            if not all(self.blen[_sym_pair(i, j)] <= dmax for j in groups[name]):
                return False
        for name in rule.no_double_to:
            for j in groups[name]:
                if m.GetBondBetweenAtoms(i, j).GetBondType() == Chem.BondType.DOUBLE:  # type: ignore[attr-defined]
                    return False
        return True

    def _apply(
        self,
        i: int,
        rule: FixupRule,
        groups: dict[str, list[int]],
        m: Chem.RWMol,  # type: ignore[attr-defined]
    ) -> None:
        # --- resolve the double-target selection first: single actions
        # may reference "the group minus the double target" ---------------
        dbl_group: str | None = None
        dbl_targets: list[int] = []
        if rule.make_double:
            gname = rule.make_double.split(":", 1)[1]
            dbl_group = gname
            dmax = rule.action_len.get(gname)
            eligible = [
                j for j in groups[gname]
                if dmax is None or self.blen[_sym_pair(i, j)] <= dmax
            ]
            if rule.make_double.startswith("shortest:"):
                if eligible:
                    dbl_targets = [
                        min(eligible, key=lambda j: self.blen[_sym_pair(i, j)])
                    ]
            else:
                dbl_targets = list(eligible)

        # --- singles ------------------------------------------------------
        for spec in rule.make_single:
            if spec == "non_single_non_arom":
                for b in m.GetAtomWithIdx(i).GetBonds():  # type: ignore[attr-defined]
                    if b.GetBondType() not in (  # type: ignore[attr-defined]
                        Chem.BondType.SINGLE, Chem.BondType.AROMATIC  # type: ignore[attr-defined]
                    ):
                        b.SetBondType(Chem.BondType.SINGLE)  # type: ignore[attr-defined]
            elif spec.endswith("_others"):
                if dbl_group is None:
                    raise ValueError(
                        f"make_single {spec!r} needs make_double in rule {rule.name!r}"
                    )
                for j in groups[spec[: -len("_others")]]:
                    if j not in dbl_targets:
                        m.GetBondBetweenAtoms(i, j).SetBondType(  # type: ignore[attr-defined]
                            Chem.BondType.SINGLE  # type: ignore[attr-defined]
                        )
            elif spec.startswith("group:"):
                for j in groups[spec.split(":", 1)[1]]:
                    m.GetBondBetweenAtoms(i, j).SetBondType(  # type: ignore[attr-defined]
                        Chem.BondType.SINGLE  # type: ignore[attr-defined]
                    )
            else:
                raise ValueError(f"unknown make_single spec {spec!r} in rule {rule.name!r}")

        # --- doubles ------------------------------------------------------
        for j in dbl_targets:
            bo = m.GetBondBetweenAtoms(i, j)  # type: ignore[attr-defined]
            if not rule.only_if_single or bo.GetBondType() == Chem.BondType.SINGLE:  # type: ignore[attr-defined]
                bo.SetBondType(Chem.BondType.DOUBLE)  # type: ignore[attr-defined]

        # --- charges ------------------------------------------------------
        for spec, charge in rule.charges:
            if spec == "center":
                m.GetAtomWithIdx(i).SetFormalCharge(charge)  # type: ignore[attr-defined]
            elif spec.endswith("_others"):
                if dbl_group is None:
                    raise ValueError(
                        f"charge spec {spec!r} needs make_double in rule {rule.name!r}"
                    )
                for j in groups[spec[: -len("_others")]]:
                    if j not in dbl_targets:
                        m.GetAtomWithIdx(j).SetFormalCharge(charge)  # type: ignore[attr-defined]
            else:
                raise ValueError(f"unknown charge spec {spec!r} in rule {rule.name!r}")
