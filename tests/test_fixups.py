"""Fixup-rulebook tests: schema validation and per-rule behavior.

The geometry perception's functional-group corrections (carbonyl rescue
+ ester-O protection, amidine, nitro, phosphate P=O, sulfonamide) are
declarative :class:`FixupRule` instances evaluated by
:class:`FixupEngine`.  Every rule gets a focused unit test on a
hand-built molecule - trigger gates, action selection (shortest vs all),
charge placement - plus a schema-validation test that catches malformed
rulebook entries.  The full-pipeline parity that proves the rulebook
reproduces the hand-written passes is the benchmark run.
"""

from __future__ import annotations

import json
import os
import sys

# Ensure the src directory is on the path when running without install
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import pytest

import bope as perception
from bope._deps import Chem
from bope.fixups import FIXUP_RULES, FixupEngine
from bope.helpers import _sym_pair

pytestmark = pytest.mark.skipif(
    not perception._RDKIT_AVAILABLE, reason="rdkit not installed"
)


# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------
def _make_mol(elements, edges):
    """RWMol with the graph edges added as SINGLE bonds - the pre-fixup
    assembly state for a molecule whose length rule left everything
    single."""
    m = Chem.RWMol()
    for el in elements:
        m.AddAtom(Chem.Atom(el))
    for a, b in edges:
        m.AddBond(a, b, Chem.BondType.SINGLE)
    return m


def _make_engine(elements, edges, blen):
    deg = {a: 0 for a in range(len(elements))}
    for a, b in edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    return FixupEngine(elements, set(edges), deg, blen)


def _blen(edges, lengths):
    """blen map from (edge, length) pairs."""
    return {_sym_pair(a, b): d for (a, b), d in zip(edges, lengths)}


def _bonds(m, i):
    """{neighbor: BondType} of atom i's bonds."""
    atom = m.GetAtomWithIdx(i)
    return {b.GetOtherAtom(atom).GetIdx(): b.GetBondType() for b in atom.GetBonds()}


def _rule(name):
    for r in FIXUP_RULES:
        if r.name == name:
            return r
    raise AssertionError(f"no rule named {name!r}")


# ---------------------------------------------------------------------------
# Rulebook schema
# ---------------------------------------------------------------------------
def test_rulebook_schema():
    """Every rule references only its own groups; exact is exclusive
    with min/max; make_double / make_single / charges are well-formed."""
    for r in FIXUP_RULES:
        names = {g.name for g in r.groups}
        assert r.name, "rules must be named"
        assert r.center, f"{r.name}: no center element"
        for g in r.groups:
            assert g.name, f"{r.name}: unnamed group"
            assert not (g.exact is not None and (g.min or g.max is not None)), (
                f"{r.name}: group {g.name!r} mixes exact with min/max"
            )
        for name in (*r.require, *r.no_double_to):
            assert name in names, f"{r.name}: unknown group {name!r} in gate"
        for clause in r.require_or:
            assert all(name in names for name in clause), (
                f"{r.name}: unknown group in require_or {clause!r}"
            )
        for name, dmax in (*r.max_len.items(), *r.action_len.items()):
            assert name in names, f"{r.name}: unknown group {name!r} in length cap"
            assert dmax > 0, f"{r.name}: non-positive length cap for {name!r}"
        if r.make_double:
            kind, sep, gname = r.make_double.partition(":")
            assert kind in ("shortest", "all") and sep, (
                f"{r.name}: bad make_double {r.make_double!r}"
            )
            assert gname in names, f"{r.name}: unknown group {gname!r} in make_double"
        for spec in r.make_single:
            if spec == "non_single_non_arom":
                continue
            if spec.startswith("group:"):
                assert spec.split(":", 1)[1] in names, f"{r.name}: bad make_single {spec!r}"
            elif spec.endswith("_others"):
                assert r.make_double, f"{r.name}: {spec!r} needs make_double"
                assert spec[: -len("_others")] in names, f"{r.name}: bad make_single {spec!r}"
            else:
                raise AssertionError(f"{r.name}: bad make_single {spec!r}")
        for spec, _charge in r.charges:
            if spec == "center":
                continue
            if spec.endswith("_others"):
                assert r.make_double, f"{r.name}: {spec!r} needs make_double"
                assert spec[: -len("_others")] in names, f"{r.name}: bad charge spec {spec!r}"
            else:
                raise AssertionError(f"{r.name}: bad charge spec {spec!r}")


def test_plain_alkane_fires_nothing():
    """The whole rulebook leaves a plain alkane untouched."""
    elements = ["C", "C", "C"]
    edges = [(0, 1), (1, 2)]
    blen = _blen(edges, [1.53, 1.53])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(m, set()) == []


# ---------------------------------------------------------------------------
# carbonyl_rescue
# ---------------------------------------------------------------------------
def test_carbonyl_rescue_with_n_neighbor():
    # C(=O)-NH2-like: C-O at 1.36 (crystal-short, caffeine C2/C6 in
    # 3RFM), N neighbor -> C=O
    elements = ["C", "N", "O"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.36, 1.36])
    m = _make_mol(elements, edges)
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("carbonyl_rescue"),)
    )
    assert fired == ["carbonyl_rescue"]
    assert _bonds(m, 0) == {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE}


def test_carbonyl_rescue_skips_long_o():
    # C-O at 1.45 > 1.40: rule fires (O + N present) but leaves it single
    elements = ["C", "N", "O"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.36, 1.45])
    m = _make_mol(elements, edges)
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("carbonyl_rescue"),)
    )
    assert fired == ["carbonyl_rescue"]
    assert _bonds(m, 0)[2] == Chem.BondType.SINGLE


def test_carbonyl_rescue_keeps_existing_double():
    # the length rule already doubled a short C=O: rescue is a no-op
    elements = ["C", "N", "O"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.36, 1.30])
    m = _make_mol(elements, edges)
    m.GetBondBetweenAtoms(0, 2).SetBondType(Chem.BondType.DOUBLE)
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("carbonyl_rescue"),)
    )
    assert fired == ["carbonyl_rescue"]
    assert _bonds(m, 0)[2] == Chem.BondType.DOUBLE


def test_carbonyl_rescue_needs_o_and_trigger():
    # no O neighbor: nothing to rescue
    elements = ["C", "N"]
    edges = [(0, 1)]
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, {}).apply(
        m, set(), (_rule("carbonyl_rescue"),)
    ) == []
    # O but no N and fewer than two aromatic neighbors: no fire
    elements = ["C", "O", "C"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.36, 1.40])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("carbonyl_rescue"),)
    ) == []


def test_carbonyl_rescue_aromatic_bridge_trigger():
    # diaryl / aryl-heteroaryl ketone bridge (3QTX X43's C(=O) at
    # 1.358): two aromatic neighbors fire the rescue
    elements = ["C", "O", "C", "C"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.358, 1.40, 1.40])
    m = _make_mol(elements, edges)
    fired = _make_engine(elements, edges, blen).apply(
        m, {2, 3}, (_rule("carbonyl_rescue"),)
    )
    assert fired == ["carbonyl_rescue"]
    assert _bonds(m, 0)[1] == Chem.BondType.DOUBLE


def test_carbonyl_rescue_skips_aromatic_center():
    # aromatic ring C's are never rescue targets
    elements = ["C", "O", "N"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.36, 1.36])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, {0}, (_rule("carbonyl_rescue"),)
    ) == []


# ---------------------------------------------------------------------------
# ester_o_single
# ---------------------------------------------------------------------------
def test_ester_o_single_forces_bridged_o_single():
    # ester C-O-C: the length rule doubled the wrong bond (crystal ester
    # C-O 1.247 < C=O 1.272, BT5/1WQW); the bridged O must go back to
    # single while the terminal O stays double
    elements = ["C", "O", "O", "C"]
    edges = [(0, 1), (0, 2), (1, 3)]  # O1 bridged (O-C), O2 terminal
    blen = _blen(edges, [1.247, 1.272, 1.36])
    m = _make_mol(elements, edges)
    m.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.DOUBLE)
    m.GetBondBetweenAtoms(0, 2).SetBondType(Chem.BondType.DOUBLE)
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("ester_o_single"),)
    )
    # fires per center atom: C0 (the ester carbon) and C3 (whose O1 is
    # also bridged - a no-op there)
    assert "ester_o_single" in fired
    assert _bonds(m, 0) == {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE}


# ---------------------------------------------------------------------------
# amidine
# ---------------------------------------------------------------------------
def test_amidine_doubles_shorter_n():
    # benzamidine-like: C with two N's at 1.33 / 1.35 -> shorter C=N
    elements = ["C", "N", "N", "C"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.33, 1.35, 1.50])
    m = _make_mol(elements, edges)
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("amidine"),)
    )
    assert fired == ["amidine"]
    assert _bonds(m, 0)[1] == Chem.BondType.DOUBLE
    assert _bonds(m, 0)[2] == Chem.BondType.SINGLE


def test_amidine_needs_exactly_two_n():
    # a single short C-N (a plain amine, e.g. the drug side-chain of
    # carazolol) is NOT an amidine: the exact=2 gate must refuse it
    elements = ["C", "N", "C"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.35, 1.53])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("amidine"),)
    ) == []
    # and three N's (a guanidine-like center) is out of scope too
    elements = ["C", "N", "N", "N"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.33, 1.33, 1.33])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("amidine"),)
    ) == []


def test_amidine_skips_long_pair_and_hetero_nbr():
    # true C-N single pair (1.50 / 1.51): no double
    elements = ["C", "N", "N"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.50, 1.51])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("amidine"),)
    ) == []
    # an O neighbor (amide C) excludes the rule
    elements = ["C", "N", "N", "O"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.33, 1.35, 1.36])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("amidine"),)
    ) == []


# ---------------------------------------------------------------------------
# nitro
# ---------------------------------------------------------------------------
def test_nitro_charge_separated_form():
    # N with two terminal O's (1.20 / 1.42) + a carbon: [N+](=O)[O-]
    elements = ["N", "C", "O", "O"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.30, 1.20, 1.42])
    m = _make_mol(elements, edges)
    m.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.DOUBLE)  # spurious nitro-aryl double
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("nitro"),)
    )
    assert fired == ["nitro"]
    bonds = _bonds(m, 0)
    assert bonds[1] == Chem.BondType.SINGLE  # nitro-aryl C-N cleared
    assert bonds[2] == Chem.BondType.DOUBLE  # shortest O
    assert bonds[3] == Chem.BondType.SINGLE
    assert m.GetAtomWithIdx(0).GetFormalCharge() == 1
    assert m.GetAtomWithIdx(2).GetFormalCharge() == 0
    assert m.GetAtomWithIdx(3).GetFormalCharge() == -1


def test_nitro_needs_exactly_two_term_o():
    # one terminal O (an N-oxide / hydroxylamine N-O) is not nitro
    elements = ["N", "C", "O"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.30, 1.20])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("nitro"),)
    ) == []


def test_nitro_skips_long_o_and_low_degree():
    # one O at 1.50 > 1.45: the gate refuses the whole rule
    elements = ["N", "C", "O", "O"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.30, 1.20, 1.50])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("nitro"),)
    ) == []
    # degree 2 (no carbon neighbor): below min_degree 3
    elements = ["N", "O", "O"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.20, 1.42])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("nitro"),)
    ) == []


# ---------------------------------------------------------------------------
# phosphate
# ---------------------------------------------------------------------------
def test_phosphate_doubles_shortest_term_o():
    # crystal P-O refines to 1.55-1.70 (1TPB PGH 1.701); the length rule
    # leaves everything single, the fixup doubles the shortest
    elements = ["P", "O", "O", "C"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.60, 1.70, 1.80])
    m = _make_mol(elements, edges)
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("phosphate"),)
    )
    assert fired == ["phosphate"]
    assert _bonds(m, 0)[1] == Chem.BondType.DOUBLE
    assert _bonds(m, 0)[2] == Chem.BondType.SINGLE


def test_phosphate_needs_two_term_o():
    # a single terminal O (a phosphinate P-O) is not a phosphate P=O
    elements = ["P", "O", "C"]
    edges = [(0, 1), (0, 2)]
    blen = _blen(edges, [1.60, 1.80])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("phosphate"),)
    ) == []


def test_phosphate_skips_when_double_exists():
    # a P=O already present (the length rule caught a short one): no fire
    elements = ["P", "O", "O", "C"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.60, 1.70, 1.80])
    m = _make_mol(elements, edges)
    m.GetBondBetweenAtoms(0, 2).SetBondType(Chem.BondType.DOUBLE)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("phosphate"),)
    ) == []


# ---------------------------------------------------------------------------
# sulfonamide
# ---------------------------------------------------------------------------
def test_sulfonamide_both_o_double_n_single():
    # R-S(=O)(=O)-NH2: both terminal O's double, S-N forced single even
    # when the length rule misordered it (S-N measures 1.56-1.63, inside
    # the (N,S) double cutoff 1.70)
    elements = ["S", "O", "O", "N", "C"]
    edges = [(0, 1), (0, 2), (0, 3), (0, 4)]
    blen = _blen(edges, [1.50, 1.55, 1.60, 1.80])
    m = _make_mol(elements, edges)
    m.GetBondBetweenAtoms(0, 3).SetBondType(Chem.BondType.DOUBLE)
    fired = _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("sulfonamide"),)
    )
    assert fired == ["sulfonamide"]
    bonds = _bonds(m, 0)
    assert bonds[1] == Chem.BondType.DOUBLE
    assert bonds[2] == Chem.BondType.DOUBLE
    assert bonds[3] == Chem.BondType.SINGLE
    assert bonds[4] == Chem.BondType.SINGLE


def test_sulfonamide_needs_n():
    # an S with two terminal O's but no N (a sulfone) is not sulfonamide
    elements = ["S", "O", "O", "C"]
    edges = [(0, 1), (0, 2), (0, 3)]
    blen = _blen(edges, [1.50, 1.55, 1.80])
    m = _make_mol(elements, edges)
    assert _make_engine(elements, edges, blen).apply(
        m, set(), (_rule("sulfonamide"),)
    ) == []


# ---------------------------------------------------------------------------
# End-to-end wiring through the full geometry pipeline
# ---------------------------------------------------------------------------
def test_caffeine_carbonyls_end_to_end():
    """Caffeine (3RFM CFF) through the full geometry pipeline keeps
    both ring carbonyls double (the rescue's motivating case)."""
    with open(
        os.path.join(os.path.dirname(__file__), "fixtures", "crystal_ligands.json"),
        encoding="utf-8",
    ) as fh:
        atoms_raw = json.load(fh)["3RFM"]["CFF"][0]
    atoms = [(el, tuple(xyz)) for el, xyz in atoms_raw]
    mol, strategy = perception.perceive_bond_orders(atoms, resname=None)
    assert strategy == "geometry"
    assert mol is not None
    # [#6]=O: the ring carbons are aromatic (c=O in SMILES), and an
    # aliphatic `C` query would not match them
    assert len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]=O"))) == 2
