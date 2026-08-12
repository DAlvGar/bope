"""Bond-order perception test suite for :mod:`bope`.

* **Synthetic recoverability** - every corpus SMILES is embedded with
  ETKDG (randomSeed=42), perturbed with isotropic Gaussian coordinate
  noise, and must come back with the *exact* reference bond graph,
  the exact formula, and a successful ``AddHs`` (no over-valent atoms).
  Two noise levels are pinned: 0.0 (182 of 187 neutral molecules) and
  0.03 A bond-RMS (172 of 187) - the exclusions below are documented,
  measured properties of the current ETKDG embeds, not silent skips.
* **Crystal ground truth** - the 16 eval-dataset complexes: every
  ligand residue must reproduce the RCSB CCD formula on both the CCD
  template path (cache seeded offline) and the pure-geometry path
  (cache poisoned).
* **Charge coverage** - charged molecules degrade gracefully (never
  crash, never over-valent), and the one formally-supported charge
  case - quaternary nitrogen - recovers exactly.
* **Failure modes** - unknown HET codes, disconnected fragments,
  missing atoms, single atoms, empty input: no crash, documented
  results.
* **Tautomer sensitivity** - pyrrole-vs-pyridine N-H placement is
  checked per-atom against the reference.

The exclusions are keyed by corpus name; the corpus lives in
:mod:`bope.corpus`.  If an RDKit upgrade changes the ETKDG seed-42
geometries, re-measure the exclusion sets with the benchmark harness
before re-enabling the strict assertions.
"""

from __future__ import annotations

import json
import os
import sys

# Ensure the src directory is on the path when running without install
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pytest

import bope as perception
from bope import corpus

pytestmark = pytest.mark.skipif(
    not perception._RDKIT_AVAILABLE, reason="rdkit not installed"
)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
except ImportError:  # pragma: no cover - guarded by pytestmark above
    Chem = None  # type: ignore[assignment]

# ---------------------------------------------------------------------------
# Documented exclusions - shared with the benchmark (see corpus.py)
# ---------------------------------------------------------------------------
# The exclusion sets and the crystal ground-truth formula table live in
# :mod:`bope.corpus` so tests and the benchmark consume the same measured
# metadata.

#: Crystal-ligand atom fixtures, extracted once from the deposited PDBs:
#: {pdb_id: {het: [ [ [element, [x, y, z]], ... ], ... ]}} - one list per
#: residue, preserving the multi-residue cases (1IEP x2, 3PBL x2,
#: 5HVP ACE x1 + STA x2).
_FIXTURES = json.load(
    open(
        os.path.join(os.path.dirname(__file__), "fixtures", "crystal_ligands.json"),
        encoding="utf-8",
    )
)


def _crystal_residue_atoms(pdb_id: str, het: str) -> list[list[tuple[str, tuple]]]:
    """Per-residue (element, xyz) atom lists for one (pdb, het) pair."""
    return [
        [(el, tuple(xyz)) for el, xyz in residue]
        for residue in _FIXTURES[pdb_id][het]
    ]

#: Tautomer-sensitive molecules for the explicit-H placement test.
_TAUTOMERS: dict[str, str] = {
    "uracil": "c1c[nH]c(=O)[nH]c1=O",
    "imidazole": "c1c[nH]cn1",
    "pyridine": "c1ccncc1",
    "purine": "c1ncnc2[nH]cnc12",
    "caffeine": "Cn1cnc2c1C(=O)N(C(=O)N2C)C",
}

# ---------------------------------------------------------------------------
# Harness: SMILES -> ETKDG (seed 42) -> noise -> perceive -> exact checks
# ---------------------------------------------------------------------------

_EMBED_CACHE: dict[str, "Chem.Mol"] = {}


def _embed(smi: str) -> "Chem.Mol":
    """ETKDG embed with Hs (randomSeed=42), cached per SMILES."""
    if smi not in _EMBED_CACHE:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))  # type: ignore[attr-defined]
        AllChem.EmbedMolecule(mol, randomSeed=42)  # type: ignore[attr-defined]
        _EMBED_CACHE[smi] = mol
    return _EMBED_CACHE[smi]


def _noisy_atoms(
    molH: "Chem.Mol", bond_sigma: float, seed: int = 42
) -> list[tuple[str, tuple[float, float, float]]]:
    """Heavy-atom (element, xyz) tuples with isotropic Gaussian noise.

    The per-axis sigma is ``bond_sigma / sqrt(2)`` so that the RMS of the
    two endpoints' displacement along a bond equals ``bond_sigma`` - the
    level parameter is the bond-RMS perturbation.
    """
    sigma = bond_sigma / 2 ** 0.5
    rng = np.random.default_rng(seed)
    atoms = []
    for a in molH.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        pos = molH.GetConformer().GetAtomPosition(a.GetIdx())
        xyz = tuple(
            float(v) + float(rng.normal(0, sigma)) for v in (pos.x, pos.y, pos.z)
        )
        atoms.append((a.GetSymbol().upper(), xyz))
    return atoms


def _bond_graph(mol: "Chem.Mol") -> frozenset:
    """Exact bond graph: (low-idx, high-idx, order) triples."""
    return frozenset(
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            round(b.GetBondTypeAsDouble(), 1),
        )
        for b in mol.GetBonds()
    )


def _formula(mol: "Chem.Mol") -> str:
    return rdMolDescriptors.CalcMolFormula(mol)  # type: ignore[attr-defined]


def _addh_ok(mol: "Chem.Mol") -> bool:
    """AddHs must not raise - the over-valent-atom failure mode."""
    try:
        Chem.AddHs(mol)  # type: ignore[attr-defined]
        return True
    except Exception:  # noqa: BLE001 - RuntimeError on over-valent atoms
        return False


def _recover(
    smi: str, bond_sigma: float
) -> tuple["Chem.Mol | None", str, bool, bool, bool]:
    """Embed -> noise -> perceive, returning (mol, strategy, graph_ok,
    formula_ok, addh_ok) against the reference."""
    molH = _embed(smi)
    ref = Chem.RemoveHs(molH)  # type: ignore[attr-defined]
    mol, strategy = perception.perceive_bond_orders(
        _noisy_atoms(molH, bond_sigma), resname=None
    )
    if mol is None:
        return None, strategy, False, False, False
    return (
        mol,
        strategy,
        _bond_graph(mol) == _bond_graph(ref),
        _formula(mol) == _formula(ref),
        _addh_ok(mol),
    )


def _strict_names(bond_sigma: float) -> list[tuple[str, str]]:
    """Neutral corpus entries with exact recovery asserted at this noise
    level: 0.0 excludes the embedder-unfaithful set, 0.03 additionally the
    cutoff-boundary set.  The OpenBabel-fallback trio leaves the strict set
    when OpenBabel is unavailable (their 0-noise recovery depends on it)."""
    excluded = set(corpus.EMBED_UNFAITHFUL)
    if bond_sigma > 0.0:
        excluded |= corpus.NOISE_SENSITIVE
    if not perception._OPENBABEL_AVAILABLE:
        excluded |= corpus.OPENBABEL_FALLBACK
    return [
        (name, smi)
        for name, smi, charged in corpus.CORPUS
        if not charged and name not in excluded
    ]


# ---------------------------------------------------------------------------
# Synthetic recoverability
# ---------------------------------------------------------------------------


def test_synthetic_recovery_at_zero_noise():
    """Exact bond-graph + formula + AddHs recovery at 0 noise for the
    182 faithful neutral corpus members (187 minus the 5 documented
    embedder-unfaithful exclusions)."""
    strict = _strict_names(0.0)
    assert len(strict) >= 100  # issue #40: "~100+ molecules"
    for name, smi in strict:
        mol, strategy, gok, fok, hok = _recover(smi, 0.0)
        assert gok and fok and hok, (
            f"{name}: graph={gok} formula={fok} addh={hok} strategy={strategy}"
        )


def test_synthetic_recovery_at_03_noise():
    """Exact recovery at 0.03 A bond-RMS coordinate noise for the 172
    molecules that are robust to it (187 minus the 5 embedder-unfaithful
    and the 10 documented cutoff-boundary events; 169 when OpenBabel is
    not installed, since the sulfur-heterocycle trio recovers only via
    the OpenBabel fallback).  The perception is tuned for crystal-quality
    coordinates (correlated deviations ~0.06 A); larger noise is the
    benchmark's territory."""
    strict = _strict_names(0.03)
    assert len(strict) >= 100
    for name, smi in strict:
        mol, strategy, gok, fok, hok = _recover(smi, 0.03)
        assert gok and fok and hok, (
            f"{name}: graph={gok} formula={fok} addh={hok} strategy={strategy}"
        )


def test_recovered_molecules_never_use_distance_baseline():
    """Every 0-noise recovery must come from the geometry or OpenBabel
    tiers - never the covalent-radius distance fallback - and the in-house
    geometry path must do the bulk of the work (a regression that silently
    moved everything to OpenBabel would hide behind the exactness check)."""
    counts: dict[str, int] = {}
    for name, smi in _strict_names(0.0):
        mol, strategy, gok, fok, hok = _recover(smi, 0.0)
        if not (gok and fok and hok):
            continue  # covered by test_synthetic_recovery_at_zero_noise
        counts[strategy] = counts.get(strategy, 0) + 1
    assert counts.get("distance", 0) == 0, f"distance fallback used: {counts}"
    assert counts.get("geometry", 0) >= 170, (
        f"geometry tier degraded: {counts}"
    )


# ---------------------------------------------------------------------------
# Charge coverage
# ---------------------------------------------------------------------------


def test_charged_molecules_degrade_gracefully():
    """Charged corpus members must never crash and must never produce an
    over-valent mol (AddHs must succeed).  Exact recovery is NOT required:
    the API deliberately does not consume the ``charge`` argument, so
    formal-charge-dependent H counts differ from the reference.  The one
    exception is the quaternary nitrogen, whose +1 formal charge is forced
    by valence - asserted below as the documented charge decision."""
    for name, smi, _charged in corpus.CHARGED:
        mol, strategy = perception.perceive_bond_orders(
            _noisy_atoms(_embed(smi), 0.0), resname=None
        )
        assert mol is not None, f"{name}: perceived as None ({strategy})"
        assert _addh_ok(mol), f"{name}: AddHs failed (over-valent?)"


def test_quaternary_nitrogen_recovers_exactly():
    """The charge-decision test: tetraalkylammonium nitrogen has four
    single bonds and no double to demote - the geometry path's valence
    demotion pass assigns the +1 formal charge, which makes the perceived
    mol match the reference exactly (graph AND formula)."""
    mol, strategy, gok, fok, hok = _recover("C[N+](C)(C)C", 0.0)
    assert strategy == "geometry"
    assert gok and fok and hok


def test_charge_argument_is_ignored():
    """The API's ``charge`` argument is accepted but unused: different
    values must produce byte-identical results (documented behaviour)."""
    atoms = _noisy_atoms(_embed(_TAUTOMERS["purine"]), 0.0)
    m0, s0 = perception.perceive_bond_orders(atoms, resname=None, charge=0)
    m5, s5 = perception.perceive_bond_orders(atoms, resname=None, charge=5)
    assert m0 is not None and m5 is not None
    assert s0 == s5 == "geometry"
    assert _bond_graph(m0) == _bond_graph(m5)
    assert _formula(m0) == _formula(m5)


# ---------------------------------------------------------------------------
# Crystal ground truth on real deposited coordinates
# ---------------------------------------------------------------------------


def test_crystal_ligands_real_coordinates_ccd_path():
    """The CCD template path on real deposited coordinates: with the CCD
    cache seeded from the corpus' canonical SMILES (which ARE the CCD
    canonical SMILES), every one of the 20 ligand residues across the 16
    crystal complexes must reproduce the RCSB ground-truth formula.
    Fully offline: the cache is seeded, no network access happens."""
    seed = {name.split()[0]: smi for name, smi in corpus.CRYSTAL_LIGANDS}
    original = dict(perception._CCD_CACHE)
    perception._CCD_CACHE.update(seed)
    try:
        strategies = set()
        checked = 0
        for pdb_id, ligu in corpus.CRYSTAL_EXPECTED.items():
            for het, want in ligu.items():
                for res_i, atoms in enumerate(
                    _crystal_residue_atoms(pdb_id, het)
                ):
                    mol, strategy = perception.perceive_bond_orders(
                        atoms, resname=het
                    )
                    strategies.add(strategy)
                    checked += 1
                    assert mol is not None, (
                        f"{pdb_id}/{het} res {res_i}: perceived as None "
                        f"({strategy})"
                    )
                    got = _formula(mol)
                    assert got in want, (
                        f"{pdb_id}/{het} res {res_i}: got {got} "
                        f"want {sorted(want)} via {strategy}"
                    )
        assert checked >= 17
        # the template path must actually be exercised (the STA fragment
        # without OXT cannot match the 12-atom template and falls back)
        assert "ccd-template" in strategies, f"template never used: {strategies}"
    finally:
        perception._CCD_CACHE.clear()
        perception._CCD_CACHE.update(original)


def test_crystal_ligands_real_coordinates_geometry_path():
    """The offline geometry path on the same 20 residues, forced by
    poisoning the CCD cache: every residue must be perceived by the
    geometry strategy with the ground-truth formula - the regression that
    OpenBabel corrupts (pentavalent carbons in N-rich fused systems) is
    covered by CFF / STU / ZMA here."""
    for pdb_id, ligu in corpus.CRYSTAL_EXPECTED.items():
        for het in ligu:
            perception._CCD_CACHE[het] = None
    original = dict(perception._CCD_CACHE)
    try:
        checked = 0
        for pdb_id, ligu in corpus.CRYSTAL_EXPECTED.items():
            for het, want in ligu.items():
                for res_i, atoms in enumerate(
                    _crystal_residue_atoms(pdb_id, het)
                ):
                    mol, strategy = perception.perceive_bond_orders(
                        atoms, resname=het
                    )
                    checked += 1
                    assert strategy == "geometry", (
                        f"{pdb_id}/{het} res {res_i}: strategy {strategy}, "
                        f"expected geometry"
                    )
                    assert mol is not None
                    got = _formula(mol)
                    assert got in want, (
                        f"{pdb_id}/{het} res {res_i}: got {got} "
                        f"want {sorted(want)}"
                    )
        assert checked >= 17
    finally:
        perception._CCD_CACHE.clear()
        perception._CCD_CACHE.update(original)


# ---------------------------------------------------------------------------
# Held-out bug-family regressions (crystal100 held-out failures)
# ---------------------------------------------------------------------------
# The 27 held-out targets that drove the perception fixes (see
# benchmarks/crystal100/heldout_failures_analysis.md): every target must
# recover exactly on the geometry tier, except the five documented data
# limitations below, which must degrade gracefully instead.  The atoms are
# the deposited crystal coordinates, frozen in the fixture; the ref SMILES
# is the CCD canonical descriptor from the held-out dataset entry.

_HELDOUT_FIXTURES = json.load(
    open(
        os.path.join(
            os.path.dirname(__file__), "fixtures", "heldout_regression.json"
        ),
        encoding="utf-8",
    )
)

#: (pdb, het) targets whose deposited coordinates defeat the geometry tier.
#: 5MUY MGT: ribose refined with a C=C(O) enol and C-C at 1.296 A.  7FOZ
#: WD0: N#C at 1.172 A inside a thiazoline ring.  5ME6 M7G: charged
#: reference (the geometry tier is neutral).  6S7B KYH: ribose with a
#: ring-refined enol.  7T2X EMY: ring-chain ambiguous thioester -
#: coordinates carry both the open (CCD) and closed ring forms (C-C
#: 1.784 A closure, C-N 1.04-1.18 A, S-C 1.334 A), so the CCD descriptor
#: cannot be reproduced from them.  The CCD template path recovers all but
#: EMY, which cannot match the template graph.
_HELDOUT_LIMITATIONS = {
    "5MUY_MGT", "7FOZ_WD0", "5ME6_M7G", "6S7B_KYH", "7T2X_EMY",
}


def _heldout_metrics(mol, ref):
    """(formula, graph, exact, addh) mirroring the benchmark's _metrics:
    formula/graph on the charge-zeroed copy, exact on the raw mol,
    canonical non-isomeric SMILES throughout."""
    if mol is None:
        return False, False, False, False
    try:
        m_n = Chem.Mol(mol)  # type: ignore[attr-defined]
        for a in m_n.GetAtoms():
            a.SetFormalCharge(0)
        for a in m_n.GetAtoms():
            a.SetNoImplicit(False)
        m_n.UpdatePropertyCache(strict=False)
        r_n = Chem.Mol(ref)  # type: ignore[attr-defined]
        for a in r_n.GetAtoms():
            a.SetFormalCharge(0)
        for a in r_n.GetAtoms():
            a.SetNoImplicit(False)
        r_n.UpdatePropertyCache(strict=False)
        formula = _formula(m_n) == _formula(r_n)
        graph = Chem.MolToSmiles(  # type: ignore[attr-defined]
            m_n, isomericSmiles=False
        ) == Chem.MolToSmiles(r_n, isomericSmiles=False)  # type: ignore[attr-defined]
        exact = Chem.MolToSmiles(  # type: ignore[attr-defined]
            mol, isomericSmiles=False
        ) == Chem.MolToSmiles(ref, isomericSmiles=False)  # type: ignore[attr-defined]
        addh = _addh_ok(mol)
        return formula, graph, exact, addh
    except Exception:  # noqa: BLE001 - sanitization failure is a metric
        return False, False, False, False


def _heldout_perceive(key: str):
    """Geometry-tier perception (CCD cache poisoned) of one target."""
    d = _HELDOUT_FIXTURES[key]
    atoms = [(el, tuple(xyz)) for el, xyz in d["atoms"]]
    perception._CCD_CACHE[d["het"]] = None
    try:
        mol, strategy = perception.perceive_bond_orders(atoms, resname=None)
    finally:
        del perception._CCD_CACHE[d["het"]]
    return mol, strategy, Chem.MolFromSmiles(d["smiles"])  # type: ignore[attr-defined]


def test_heldout_regression_targets_recover_exactly():
    """Every non-limitation held-out target must recover exactly on the
    geometry tier: formula AND graph AND canonical SMILES AND AddHs.
    Covers the bug families - azole N-aryl discrimination (5KYA 6Y4),
    S-ring admission via C-S evidence (4L9Q 9TP), the spurious lactone
    C=C demotion (4L9Q 9TP), the nitro double-stacking crash (2YOH WMJ),
    and the hard-max ring rejections (7RS8 7EI, 2QG0 A94, 1WQW BT5)."""
    checked = 0
    for key, d in _HELDOUT_FIXTURES.items():
        if key in _HELDOUT_LIMITATIONS:
            continue
        mol, strategy, ref = _heldout_perceive(key)
        formula, graph, exact, addh = _heldout_metrics(mol, ref)
        assert mol is not None, f"{key}: perceived as None ({strategy})"
        assert strategy == "geometry", f"{key}: strategy {strategy}"
        assert formula and graph and exact and addh, (
            f"{key}: formula={formula} graph={graph} exact={exact} "
            f"addh={addh}"
        )
        checked += 1
    assert checked == 22


def test_heldout_documented_limitations_degrade_gracefully():
    """The five documented data limitations must never crash and must
    return a sanitizable geometry-tier mol, but are pinned NOT to recover
    exactly: if one starts recovering, move it out of _HELDOUT_LIMITATIONS
    (and re-examine its classification in heldout_failures_analysis.md)
    rather than editing the assertion."""
    for key in sorted(_HELDOUT_LIMITATIONS):
        mol, strategy, ref = _heldout_perceive(key)
        assert mol is not None, f"{key}: perceived as None ({strategy})"
        assert strategy == "geometry", f"{key}: strategy {strategy}"
        assert _addh_ok(mol), f"{key}: AddHs failed (over-valent?)"
        formula, graph, exact, addh = _heldout_metrics(mol, ref)
        assert not graph, (
            f"{key}: now recovers exactly - move out of _HELDOUT_LIMITATIONS"
        )


def test_s_heterocycle_fallback_chain():
    """The S-heterocycle fallback design (corpus.OPENBABEL_FALLBACK):
    thiazole's ETKDG embed breaks its ring (S-C 1.90 A, C-N 1.16 A) - the
    geometry tier must REFUSE the broken ring (return None) so the
    OpenBabel fallback recovers the molecule exactly; the fused
    S-heterocycles recover on the geometry tier itself.  A geometry tier
    that silently returns the length-rule ring here would regress the
    recovery and hide the fallback."""
    for name, smi, want_strategy in [
        ("thiazole", "c1cscn1", "openbabel"),
        ("benzothiazole", "c1ccc2c(c1)scn2", "geometry"),
        ("dibenzothiophene", "c1ccc2c(c1)c3ccccc3s2", "geometry"),
    ]:
        molH = _embed(smi)
        ref = Chem.RemoveHs(molH)  # type: ignore[attr-defined]
        mol, strategy = perception.perceive_bond_orders(
            _noisy_atoms(molH, 0.0), resname=None
        )
        assert strategy == want_strategy, f"{name}: {strategy}"
        assert mol is not None
        assert _bond_graph(mol) == _bond_graph(ref), f"{name}: graph mismatch"


def test_heldout_limitations_template_recovery():
    """The CCD template path recovers four of the five limitations (the
    production route for known HET codes); EMY cannot - its coordinates
    are ring-chain ambiguous and never match the template graph."""
    template_recoverable = sorted(_HELDOUT_LIMITATIONS - {"7T2X_EMY"})
    original = dict(perception._CCD_CACHE)
    try:
        for key in template_recoverable:
            d = _HELDOUT_FIXTURES[key]
            perception._CCD_CACHE[d["het"]] = d["smiles"]
            atoms = [(el, tuple(xyz)) for el, xyz in d["atoms"]]
            mol, strategy = perception.perceive_bond_orders(
                atoms, resname=d["het"]
            )
            assert mol is not None, f"{key}: template path gave None"
            assert strategy == "ccd-template", f"{key}: {strategy}"
            ref = Chem.MolFromSmiles(d["smiles"])  # type: ignore[attr-defined]
            assert _formula(mol) == _formula(ref), (
                f"{key}: template formula {_formula(mol)} != {_formula(ref)}"
            )
    finally:
        perception._CCD_CACHE.clear()
        perception._CCD_CACHE.update(original)


# ---------------------------------------------------------------------------
# Tautomer sensitivity
# ---------------------------------------------------------------------------


def test_tautomer_explicit_h_placement():
    """Pyrrole-vs-pyridine N-H placement must match the reference
    per-atom: the geometry path's Huckel judge sets the tautomer bit, and
    after AddHs every heavy atom must carry the same number of H's as the
    reference (uracil's two ring N-H's, imidazole's 1H, pyridine's 0,
    purine's 7H, caffeine's 9H-methyl-substituted ring)."""
    for name, smi in _TAUTOMERS.items():
        molH = _embed(smi)
        ref = Chem.RemoveHs(molH)  # type: ignore[attr-defined]
        mol, strategy = perception.perceive_bond_orders(
            _noisy_atoms(molH, 0.0), resname=None
        )
        assert mol is not None, f"{name}: None ({strategy})"
        assert strategy == "geometry", f"{name}: strategy {strategy}"
        refH = Chem.AddHs(ref)  # type: ignore[attr-defined]
        gotH = Chem.AddHs(mol)  # type: ignore[attr-defined]
        for i in range(refH.GetNumHeavyAtoms()):
            def h_neighbors(m, idx):
                atom = m.GetAtomWithIdx(idx)
                return sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == "H")

            assert h_neighbors(refH, i) == h_neighbors(gotH, i), (
                f"{name}: heavy atom {i} H count differs "
                f"({h_neighbors(refH, i)} vs {h_neighbors(gotH, i)})"
            )


# ---------------------------------------------------------------------------
# Failure modes
# ---------------------------------------------------------------------------


def test_unknown_het_code_falls_back_to_geometry():
    """An unknown HET code (cache poisoned, so no network) must fall
    through to the geometry tier and still recover the molecule."""
    perception._CCD_CACHE["ZZZ"] = None
    try:
        atoms = _noisy_atoms(_embed("Nc1ncnc2[nH]cnc12"), 0.0)  # adenine
        mol, strategy = perception.perceive_bond_orders(atoms, resname="ZZZ")
        assert strategy == "geometry"
        assert mol is not None
        assert _formula(mol) == "C5H5N5"
    finally:
        del perception._CCD_CACHE["ZZZ"]


def test_disconnected_fragments():
    """Benzene + ethane 20 A apart: the perception must produce a
    two-fragment sanitizable mol (bond inference is distance-based, the
    fragments simply never connect)."""
    a = _noisy_atoms(_embed("c1ccccc1"), 0.0)
    b = _noisy_atoms(_embed("CC"), 0.0)
    disconn = a + [(el, (x + 20, y + 20, z + 20)) for el, (x, y, z) in b]
    mol, strategy = perception.perceive_bond_orders(disconn, resname=None)
    assert mol is not None, f"None ({strategy})"
    assert strategy == "geometry"
    assert len(Chem.GetMolFrags(mol)) == 2  # type: ignore[attr-defined]
    assert _addh_ok(mol)


def test_missing_atom_degrades_gracefully():
    """Caffeine minus one nitrogen: the perception must not crash; it
    falls back past geometry (a broken ring) and returns a sanitizable mol
    or None - never an exception."""
    atoms = _noisy_atoms(_embed(_TAUTOMERS["caffeine"]), 0.0)
    missing = atoms[:8] + atoms[9:]
    mol, strategy = perception.perceive_bond_orders(missing, resname=None)
    assert strategy in {"geometry", "openbabel", "distance"}
    if mol is not None:
        assert _addh_ok(mol)


def test_single_atom_and_empty_input():
    """A lone atom returns a valid single-atom mol; an empty atom list
    returns (None, '') - both without raising."""
    mol, strategy = perception.perceive_bond_orders([("C", (0.0, 0.0, 0.0))])
    assert mol is not None
    mol, strategy = perception.perceive_bond_orders([])
    assert mol is None
    assert strategy == ""
