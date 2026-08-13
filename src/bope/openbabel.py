"""OpenBabel perception: round-trip the atom list through PDB text.

OpenBabel's PDB reader is designed for exactly this input (coordinates,
no CONECT records): it perceives connectivity, bond orders and
aromaticity, and the result is converted back to an RDKit Mol.  Kept as a
fallback for simple ligands and charged/protonated groups - on N-rich
fused heteroaromatics ``PerceiveBondOrders`` corrupts ring systems into
pentavalent carbons.
"""

from __future__ import annotations

from typing import Any

from bope._deps import _OPENBABEL_AVAILABLE, _RDKIT_AVAILABLE, Chem, _ob


def _atoms_to_pdb_text(
    elements: list[str], coords: list[tuple[float, float, float]]
) -> str:
    """Minimal PDB text (HETATM records, no CONECT) for OpenBabel input."""
    lines = []
    for i, (el, (x, y, z)) in enumerate(zip(elements, coords), start=1):
        lines.append(
            f"HETATM{i:5d} {el:>4} LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00 20.00          {el:>2}"
        )
    lines.append("END")
    return "\n".join(lines)


def perceive_bond_orders_with_openbabel(
    elements: list[str], coords: list[tuple[float, float, float]]
) -> Any | None:
    """Perceive bond orders by round-tripping the atom list through OpenBabel.

    Serialises the coordinates to PDB text, lets OpenBabel perceive
    connectivity, bond orders and aromaticity, and returns the result as an
    RDKit Mol with 3-D coordinates.  Hydrogens are dropped; RDKit callers
    add their own.

    Returns ``None`` when OpenBabel or RDKit is unavailable or perception
    fails - the caller falls back to the RDKit-only strategies.
    """
    if not (_RDKIT_AVAILABLE and _OPENBABEL_AVAILABLE):
        return None
    pdb_text = _atoms_to_pdb_text(elements, coords)
    try:
        conv = _ob.OBConversion()  # type: ignore[attr-defined]
        conv.SetInAndOutFormats("pdb", "sdf")
        obmol = _ob.OBMol()  # type: ignore[attr-defined]
        if not conv.ReadString(obmol, pdb_text):
            return None
        obmol.PerceiveBondOrders()  # type: ignore[attr-defined]
        obmol.FindRingAtomsAndBonds()  # type: ignore[attr-defined]
        obmol.DeleteHydrogens()  # type: ignore[attr-defined]
        sdf_text = conv.WriteString(obmol)
    except Exception:  # noqa: BLE001
        return None
    if not sdf_text.strip():
        return None
    try:
        mol = Chem.MolFromMolBlock(sdf_text, removeHs=True)  # type: ignore[attr-defined]
        if mol is None:
            return None
        Chem.SanitizeMol(mol)  # type: ignore[attr-defined]
        return mol
    except Exception:  # noqa: BLE001
        return None
