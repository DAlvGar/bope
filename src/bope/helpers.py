"""Shared helpers for the perception strategies.

Distance-graph connectivity, best-fit-plane RMS and RWMol construction are
used by more than one strategy, so they live here instead of being
duplicated.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from bope._deps import Chem
from bope.tables import (
    _BOND_TOLERANCE,
    _COVALENT_RADII,
    _COVALENT_RADII_DEFAULT,
)


def _sym_pair(a: int, b: int) -> tuple[int, int]:
    return tuple(sorted((a, b)))


def _planarity_rms(coords: list[Any]) -> float:
    """RMS deviation of points from their best-fit plane (numpy SVD)."""
    pts = np.asarray(coords, dtype=float)
    if pts.shape[0] < 4:
        return 0.0
    centered = pts - pts.mean(axis=0)
    _, singular, _ = np.linalg.svd(centered)
    return float(singular[-1] / np.sqrt(pts.shape[0]))


def _distance_bond_graph(
    elem_syms: list[str], positions: list[Any]
) -> set[tuple[int, int]]:
    """Heavy-atom connectivity from interatomic distances.

    A bond exists when the distance is shorter than the sum of covalent
    radii plus :data:`_BOND_TOLERANCE`.  Returns a set of sorted ``(i, j)``
    index pairs.  On the eval-dataset structures this matches the CONECT
    records exactly.
    """
    n = len(elem_syms)
    graph: set[tuple[int, int]] = set()
    for i in range(n):
        for j in range(i + 1, n):
            r_i = _COVALENT_RADII.get(elem_syms[i], _COVALENT_RADII_DEFAULT)
            r_j = _COVALENT_RADII.get(elem_syms[j], _COVALENT_RADII_DEFAULT)
            pi, pj = positions[i], positions[j]
            dist = math.sqrt(
                (pi.x - pj.x) ** 2 + (pi.y - pj.y) ** 2 + (pi.z - pj.z) ** 2
            )
            if dist < r_i + r_j + _BOND_TOLERANCE:
                graph.add(_sym_pair(i, j))
    return graph


def _build_rwmol(
    elements: list[str], coords: list[tuple[float, float, float]]
) -> Any:
    """Build an RDKit RWMol with atoms and a 3-D conformer; no bonds added."""
    rwmol = Chem.RWMol()  # type: ignore[attr-defined]
    conf = Chem.Conformer(len(elements))  # type: ignore[attr-defined]
    for i, (el, (x, y, z)) in enumerate(zip(elements, coords)):
        rwmol.AddAtom(Chem.Atom(el.capitalize()))  # type: ignore[attr-defined]
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    rwmol.AddConformer(conf, assignId=True)
    return rwmol
