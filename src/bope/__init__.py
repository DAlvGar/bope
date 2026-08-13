"""Bond-order perception for 3-D ligand coordinates.

Takes a plain list of ``(element, xyz)`` atoms - no BioPython, no PDB
parsing - and returns an RDKit molecule with inferred connectivity and
bond orders.  Dependency-light (RDKit, numpy, stdlib only): the
strategies below are the pipeline used by the 3D-PLI-Agent project
(where a thin BioPython bridge wraps this API for protein residues).

Four strategies are attempted in order (each lives in its own module
under this package):

1. :mod:`bope.ccd` - **CCD template**: the authoritative
   bond orders, tautomer and protonation for the HET group are fetched
   from the RCSB Chemical Component Dictionary (by residue name) and
   stamped onto the crystal graph via ``AssignBondOrdersFromTemplate``.
   Exact where it applies; requires network access and unknown HET codes
   fall through.  Results are cached per HET code in
   :data:`bope.ccd._CCD_CACHE`.
2. :mod:`bope.geometry` - **geometry perception**: an
   in-house perception: planar rings whose bonds fall in aromatic
   bond-length envelopes are scored by a Hückel (4n+2) electron-count
   judge over per-atom pi assignments; the winning mask sets aromatic
   bonds and pyrrole-type N-H's, and the remaining bonds get orders from
   length thresholds with chemistry fixups (amidines, carbonyls,
   exocyclic N/S).  Reproduces the exact molecular formulas of all 16
   crystal ligands of the validation corpus (see :mod:`bope.corpus`),
   including N-rich fused
   heteroaromatics where OpenBabel's ``PerceiveBondOrders`` corrupts ring
   systems (pentavalent carbons).
3. :mod:`bope.openbabel` - **OpenBabel**: the atom list
   is serialised to PDB and read back through OpenBabel.  Kept for simple
   ligands and charged/protonated groups.
4. :mod:`bope.distance` - **distance / covalent-radius
   fallback**: bonds are added when the interatomic distance is shorter
   than the sum of covalent radii plus 0.40 A.  All bonds start as
   ``SINGLE``; ``Chem.SanitizeMol`` then applies Kekulé / aromaticity
   perception to upgrade bond orders where the topology allows it.
   Double bonds are also assigned directly when the distance falls below
   the sum of *double-bond* covalent radii plus 0.03 A.

The returned strategy string ("ccd-template", "geometry", "openbabel",
"distance" or "") tells consumers how trustworthy the perceived orders
are: the CCD template carries the authoritative tautomer/protonation,
the geometry perception is validated on the eval dataset, OpenBabel is
known to corrupt N-rich ring systems, and the distance fallback is bare
topology.

Stereochemistry is a separate concern: nothing here labels centers, and
:func:`bope.stereo.perceive_stereochemistry` assigns R/S and E/Z from
the coordinates on request.

Validation
----------
The pipeline was validated end-to-end on 16 experimental crystal
complexes (kinases, GPCRs and proteases; rigid and flexible ligands;
N-rich fused heteroaromatics such as staurosporine, ZM241385 and
caffeine):

* the exact molecular formula of every crystal ligand is reproduced,
  on both the CCD-template path and the offline geometry path
  (exercised with the CCD cache poisoned);
* no over-valent atoms: ``Chem.AddHs`` completes without
  "explicit valence ... greater than permitted" warnings on any
  ligand - the failure mode of the previous pipeline (pentavalent
  carbons) is gone;
* the distance connectivity matches the PDB CONECT records exactly;
* the strain analysis built on top of this perception reports 0-80
  kcal/mol on the crystal ligands (rigid ligands: coordinate noise;
  flexible ligands: twisted sp2 conjugation bonds in the deposited
  coordinates) instead of the 119-565 kcal/mol garbage the previous
  pipeline produced (3D-PLI-Agent issue #39).

Quantitative comparison against the existing tools (OpenBabel
``PerceiveBondOrders``, RDKit ``rdDetermineBonds``, plain distance
connectivity) lives in the benchmark at ``benchmarks/`` - see the
README's "Why not just use the existing bond-order perception?"
section.

Why a geometry fallback at all?
-------------------------------
PDB ligands arrive without bond information; the only reliable
signals are the 3-D coordinates themselves.  The in-house geometry
perception turns those coordinates into chemistry with explicit,
checkable rules - planar aromatic rings (bond-length envelopes), a
Hückel (4n+2) electron-count judge that also decides pyrrole-vs-
pyridine tautomer N-H placement, length-threshold orders with
chemistry fixups (amidines, carbonyls, exocyclic N/S), and a valence
demotion pass before sanitization.  It reproduces the exact formulas
of all 16 eval-dataset ligands offline, and the strategy string tells
the caller exactly which tier produced the result.
"""

from __future__ import annotations

from typing import Any

# Submodules and flags are re-exported intentionally: `bope.geometry`,
# `bope._OPENBABEL_AVAILABLE` and `bope._CCD_CACHE` are part of the
# documented module surface (tests and the 3D-PLI-Agent bridge use them).
from bope import (  # noqa: F401
    ccd,
    distance,
    geometry,
    helpers,
    openbabel,
    stereo,
    tables,
)
from bope._deps import (
    _OPENBABEL_AVAILABLE,  # noqa: F401
    _RDKIT_AVAILABLE,
)
from bope.ccd import _CCD_CACHE  # noqa: F401
from bope.distance import perceive_bond_orders_distance
from bope.geometry import perceive_bond_orders_geometric
from bope.helpers import _build_rwmol, _distance_bond_graph
from bope.openbabel import perceive_bond_orders_with_openbabel
from bope.stereo import perceive_stereochemistry

#: Strategy name for "no atoms" / "RDKit unavailable" / "everything failed".
NO_STRATEGY = ""

__all__ = [
    "NO_STRATEGY",
    "perceive_bond_orders",
    "perceive_stereochemistry",
]


def perceive_bond_orders(
    atoms: list[tuple[str, tuple[float, float, float]]],
    resname: str | None = None,
    charge: int = 0,
) -> tuple[Any | None, str]:
    """Perceive bond connectivity and bond orders for a 3-D atom list.

    Args:
        atoms: Sequence of ``(element, (x, y, z))`` tuples; *element* is an
            upper-case symbol (``"C"``, ``"CL"``, ``"FE"``, ...) and the
            position may be any numeric triple.
        resname: Optional 3-letter HET code (e.g. ``"CFF"``).  When given,
            the RCSB Chemical Component Dictionary template is tried first
            and its authoritative bond orders, tautomer and protonation are
            stamped onto the crystal graph.
        charge: Net formal charge of the molecule (default ``0``).
            Accepted for API compatibility with the BioPython bridge;
            no current strategy consumes it.

    Returns:
        A 2-tuple ``(mol, strategy)`` where *mol* is an RDKit Mol with 3-D
        coordinates and perceived bond orders (or ``None``), and *strategy*
        is ``"ccd-template"``, ``"geometry"``, ``"openbabel"``,
        ``"distance"``, or ``""`` when RDKit is unavailable or every
        strategy failed.
    """
    if not _RDKIT_AVAILABLE:  # pragma: no cover
        return None, NO_STRATEGY
    if not atoms:
        return None, NO_STRATEGY

    elements = [el.upper() for el, _xyz in atoms]
    coords = [tuple(float(v) for v in xyz) for _el, xyz in atoms]

    # --- Strategy 1: CCD template (authoritative orders + tautomer) ---
    mol = ccd.perceive_bond_orders_ccd(elements, coords, resname)
    if mol is not None:
        return mol, "ccd-template"

    # --- Strategy 2: geometry perception (validated on the eval dataset) ---
    positions = _build_rwmol(elements, coords).GetConformer()
    graph = _distance_bond_graph(
        elements, [positions.GetAtomPosition(i) for i in range(len(elements))]
    )
    mol, _err = perceive_bond_orders_geometric(elements, coords, graph)
    if mol is not None:
        return mol, "geometry"

    # --- Strategy 3: OpenBabel ---
    mol = perceive_bond_orders_with_openbabel(elements, coords)
    if mol is not None:
        return mol, "openbabel"

    # --- Strategy 4: distance-based connectivity + sanitization ---
    mol = perceive_bond_orders_distance(elements, coords)
    if mol is not None:
        return mol, "distance"

    return None, NO_STRATEGY
