"""Stereochemistry perception: R/S and E/Z labels from 3-D coordinates.

Bond-order perception (:mod:`bope.geometry`) never assigns stereo, and
nothing in :mod:`bope` implicitly labels one: the perceived mol comes back
with unassigned chiral tags.  Mixing the two concerns would couple the
stereo labels to whichever bond-order strategy happened to fire, so they
are separate.  Callers that want stereo hand the perceived mol to
:func:`perceive_stereochemistry`, which adds explicit hydrogens (with
coordinates), lets RDKit's ``AssignStereochemistryFrom3D`` assign CIP R/S
labels on tetrahedral centers and E/Z labels on double bonds from the
geometry, and returns the labeled mol (hydrogens stripped again, chiral
tags preserved).

The labels are what the coordinates support: a center whose geometry is
ambiguous is left unassigned rather than guessed.  Stereochemistry cannot
rescue a wrong bond graph - a mis-perceived double bond or ring will
mislabel nearby centers, so validate the bond orders first (the strategy
string from :func:`bope.perceive_bond_orders` says how trustworthy they
are).

This mirrors what OpenBabel's ``PerceiveBondOrders`` round trip already
does as a side effect (its SDF output carries its own stereo labels); the
geometry and distance tiers get the same capability here without touching
their outputs.
"""

from __future__ import annotations

from typing import Any

from bope._deps import _RDKIT_AVAILABLE, Chem


def perceive_stereochemistry(mol: Any | None) -> Any | None:
    """Assign R/S (tetrahedral) and E/Z (double bond) labels from 3-D
    coordinates.

    Args:
        mol: an RDKit Mol with 3-D coordinates and perceived bond orders
            (e.g. the output of ``perceive_bond_orders_geometric``).  The
            input is not modified; a labeled copy is returned.

    Returns:
        A copy of *mol* with stereochemistry assigned, or ``None`` when
        RDKit is unavailable, *mol* is ``None``, or *mol* carries no
        conformer.  The returned mol has the same heavy-atom count as the
        input: explicit hydrogens are added for the assignment and removed
        afterwards, with the tetrahedral tags corrected for the removal.
    """
    if not _RDKIT_AVAILABLE or mol is None:
        return None
    try:
        if mol.GetNumConformers() == 0:
            return None
        m = Chem.Mol(mol)
        Chem.AddHs(m, addCoords=True)  # type: ignore[attr-defined]
        Chem.AssignStereochemistryFrom3D(m)  # type: ignore[attr-defined]
        return Chem.RemoveHs(m)  # type: ignore[attr-defined]
    except Exception:  # noqa: BLE001 - unkekulizable / coordinate-less mol
        return None
