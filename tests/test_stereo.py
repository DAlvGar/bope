"""Stereochemistry perception tests for :mod:`bope.stereo`.

* **Tetrahedral centers** - (R)- and (S)-2-butanol embedded with ETKDG
  (randomSeed=42) must come back with the label the SMILES declares,
  both on a pre-built mol and through the full public pipeline
  (perceive_bond_orders -> perceive_stereochemistry).
* **E/Z double bonds** - E- and Z-2-butene must come back with STEREOE /
  STEREOZ on the double bond.
* **Purity** - the input mol is not modified; the returned mol has the
  same heavy-atom count (the H's used for assignment are stripped again,
  tags corrected).
* **Degradation** - a coordinate-less mol returns None, never raises.
"""

from __future__ import annotations

import os
import sys

# Ensure the src directory is on the path when running without install
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import pytest

import bope as perception
from bope.stereo import perceive_stereochemistry

pytestmark = pytest.mark.skipif(
    not perception._RDKIT_AVAILABLE, reason="rdkit not installed"
)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors  # noqa: F401 - availability check
except ImportError:  # pragma: no cover - guarded by pytestmark above
    Chem = None  # type: ignore[assignment]


def _embed(smi: str) -> Chem.Mol:
    """ETKDG embed with Hs (randomSeed=42) - same harness as the bond-order
    tests, so the geometries are the ones already validated there."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))  # type: ignore[attr-defined]
    AllChem.EmbedMolecule(mol, randomSeed=42)  # type: ignore[attr-defined]
    return mol


def _heavy_atoms(molH: Chem.Mol) -> list[tuple[str, tuple[float, float, float]]]:
    """(element, xyz) tuples for the heavy atoms, in mol order."""
    conf = molH.GetConformer()
    return [
        (a.GetSymbol(), tuple(conf.GetAtomPosition(a.GetIdx())))
        for a in molH.GetAtoms()
        if a.GetSymbol() != "H"
    ]


def _chiral_labels(mol: Chem.Mol) -> dict[str, str]:
    """{R/S label: count} for the declared tetrahedral centers."""
    out: dict[str, str] = {}
    for idx, label in Chem.FindMolChiralCenters(  # type: ignore[attr-defined]
        mol, useLegacyImplementation=False, includeUnassigned=True
    ):
        out[str(idx)] = label
    return out


def _stereo_bonds(mol: Chem.Mol) -> dict[tuple[int, int], str]:
    """{(begin, end): STEREOE/STEREOZ} for the labeled double bonds."""
    return {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): str(b.GetStereo())
        for b in mol.GetBonds()
        if b.GetStereo() != Chem.BondStereo.STEREONONE  # type: ignore[attr-defined]
    }


# ---------------------------------------------------------------------------
# Tetrahedral centers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("smi", "want"),
    [
        # PubChem-verified mapping (CID 444683 / ChEBI:45475):
        # (S)-2-butanol is "CC[C@H](C)O", (R)-2-butanol is
        # "CC[C@@H](C)O" - the @/@@ reading depends on the
        # substituent order as written, which is easy to get inverted.
        ("C[C@H](O)CC", "S"),
        ("C[C@@H](O)CC", "R"),
    ],
)
def test_tetrahedral_center(smi, want):
    mol = perceive_stereochemistry(_embed(smi))
    assert mol is not None
    labels = list(_chiral_labels(mol).values())
    assert labels == [want]


def test_pipeline_integration():
    """Full public path: coordinates -> bond orders -> stereo."""
    smi = "C[C@H](O)CC"
    atoms = _heavy_atoms(_embed(smi))
    mol, strategy = perception.perceive_bond_orders(atoms, resname=None)
    assert mol is not None
    assert strategy == "geometry"
    labeled = perceive_stereochemistry(mol)
    assert labeled is not None
    assert list(_chiral_labels(labeled).values()) == ["S"]


def test_isomeric_smiles_round_trip():
    """The labeled mol writes back the declared stereo."""
    smi = "C[C@H](O)CC"
    labeled = perceive_stereochemistry(_embed(smi))
    assert labeled is not None
    got = Chem.MolToSmiles(  # type: ignore[attr-defined]
        labeled, isomericSmiles=True
    )
    assert got == Chem.MolToSmiles(  # type: ignore[attr-defined]
        Chem.MolFromSmiles(smi), isomericSmiles=True  # type: ignore[attr-defined]
    )


# ---------------------------------------------------------------------------
# E/Z double bonds
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("smi", "want"),
    [
        ("C/C=C/C", "STEREOE"),
        ("C/C=C\\C", "STEREOZ"),
    ],
)
def test_double_bond_stereo(smi, want):
    labeled = perceive_stereochemistry(_embed(smi))
    assert labeled is not None
    assert list(_stereo_bonds(labeled).values()) == [want]


# ---------------------------------------------------------------------------
# Purity and degradation
# ---------------------------------------------------------------------------

def test_input_not_modified():
    mol = _embed("C[C@H](O)CC")
    before = (
        Chem.MolToSmiles(mol, isomericSmiles=True),  # type: ignore[attr-defined]
        [(a.GetIdx(), a.GetChiralTag()) for a in mol.GetAtoms()],
    )
    perceive_stereochemistry(mol)
    after = (
        Chem.MolToSmiles(mol, isomericSmiles=True),  # type: ignore[attr-defined]
        [(a.GetIdx(), a.GetChiralTag()) for a in mol.GetAtoms()],
    )
    assert after == before


def test_output_heavy_atom_count():
    mol = _embed("C[C@H](O)CC")
    labeled = perceive_stereochemistry(mol)
    assert labeled is not None
    assert labeled.GetNumAtoms() == mol.GetNumHeavyAtoms()
    assert labeled.GetNumHeavyAtoms() == mol.GetNumHeavyAtoms()


def test_no_conformer_returns_none():
    mol = Chem.MolFromSmiles("C[C@H](O)CC")  # type: ignore[attr-defined]
    assert perceive_stereochemistry(mol) is None
    assert perceive_stereochemistry(None) is None
