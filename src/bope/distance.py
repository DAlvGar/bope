"""Distance / covalent-radius fallback perception.

Bare topology: bonds are added when the interatomic distance is shorter
than the sum of covalent radii plus 0.40 A.  All bonds start as
``SINGLE``; ``Chem.SanitizeMol`` then applies Kekulé / aromaticity
perception to upgrade bond orders where the topology allows it.  Double
bonds are also assigned directly when the distance falls below the sum of
*double-bond* covalent radii plus 0.03 A.

This is the last-resort strategy and is labelled as such: it cannot
distinguish aromatic from single bonds in fused ring systems at all.
"""

from __future__ import annotations

import math
from typing import Any

from bope._deps import Chem
from bope.helpers import _build_rwmol
from bope.tables import (
    _BOND_TOLERANCE,
    _COVALENT_RADII,
    _COVALENT_RADII_DEFAULT,
    _DOUBLE_BOND_TOLERANCE,
    _DOUBLE_COVALENT_RADII,
)


def perceive_bond_orders_distance(
    elements: list[str], coords: list[tuple[float, float, float]]
) -> Any | None:
    """Bare topology: bonds from covalent radii, orders from sanitization."""
    rwmol = _build_rwmol(elements, coords)
    n = rwmol.GetNumAtoms()
    positions = [rwmol.GetConformer().GetAtomPosition(i) for i in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            r_i = _COVALENT_RADII.get(elements[i], _COVALENT_RADII_DEFAULT)
            r_j = _COVALENT_RADII.get(elements[j], _COVALENT_RADII_DEFAULT)
            pi, pj = positions[i], positions[j]
            dist = math.sqrt(
                (pi.x - pj.x) ** 2
                + (pi.y - pj.y) ** 2
                + (pi.z - pj.z) ** 2
            )
            if dist < r_i + r_j + _BOND_TOLERANCE:
                # Use double-bond radii to distinguish C=O / C=C from C-O / C-C.
                # Atoms with no entry in _DOUBLE_COVALENT_RADII (metals, etc.)
                # always get a single bond.
                d_i = _DOUBLE_COVALENT_RADII.get(elements[i], None)
                d_j = _DOUBLE_COVALENT_RADII.get(elements[j], None)
                if (
                    d_i is not None
                    and d_j is not None
                    and dist < d_i + d_j + _DOUBLE_BOND_TOLERANCE
                ):
                    bond_type = Chem.BondType.DOUBLE  # type: ignore[attr-defined]
                else:
                    bond_type = Chem.BondType.SINGLE  # type: ignore[attr-defined]
                rwmol.AddBond(i, j, bond_type)

    try:
        mol = rwmol.GetMol()
        Chem.SanitizeMol(mol)  # type: ignore[attr-defined]
        return mol
    except Exception:  # noqa: BLE001
        return None
