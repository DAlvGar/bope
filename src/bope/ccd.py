"""CCD template perception: authoritative bond orders from the RCSB.

The Chemical Component Dictionary holds the authoritative bond orders,
formal charges and tautomer for every crystallographic HET group.  The
canonical SMILES is fetched from the RCSB REST API (cached per HET code for
the process lifetime) and its bond orders are stamped onto the crystal
graph with ``AssignBondOrdersFromTemplate``.
"""

from __future__ import annotations

import json
import urllib.request
from typing import Any

from bope._deps import AllChem, Chem
from bope.helpers import _build_rwmol, _distance_bond_graph

_CCD_API = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"
_CCD_TIMEOUT = 10.0
_CCD_CACHE: dict[str, str | None] = {}


def _fetch_ccd_smiles(het_code: str) -> str | None:
    """Return the canonical SMILES for a HET group from the RCSB CCD.

    The Chemical Component Dictionary holds the authoritative bond orders,
    formal charges and tautomer for every crystallographic HET group.  The
    canonical SMILES is fetched from the RCSB REST API and cached per HET
    code for the process lifetime.  Network errors and unknown HET codes
    return ``None`` (callers fall back to geometry perception), so this
    function never raises.
    """
    het = (het_code or "").strip().upper()
    if not het:
        return None
    if het in _CCD_CACHE:
        return _CCD_CACHE[het]
    try:
        with urllib.request.urlopen(  # noqa: S310 - https only
            _CCD_API.format(het), timeout=_CCD_TIMEOUT
        ) as resp:
            data = json.load(resp)
        smiles = (data.get("rcsb_chem_comp_descriptor") or {}).get("SMILES")
        _CCD_CACHE[het] = smiles
        return smiles
    except Exception:  # noqa: BLE001 - offline / unknown HET
        _CCD_CACHE[het] = None
        return None


def perceive_bond_orders_ccd(
    elements: list[str],
    coords: list[tuple[float, float, float]],
    resname: str | None,
) -> Any | None:
    """Perceive bond orders from the CCD template for the atom list.

    The CCD SMILES is fetched by the HET code and its bond orders are
    stamped onto the crystal graph (single bonds from the distance graph)
    with ``AssignBondOrdersFromTemplate``, preserving the 3-D coordinates.
    Hydrogens stay implicit; the assignment also carries the template's
    tautomer and protonation.

    Returns ``None`` when the HET code is unknown, the network is
    unavailable, the SMILES cannot be parsed, or the template graph does not
    match the crystal graph (callers fall back to geometry perception).
    """
    smiles = _fetch_ccd_smiles(resname or "")
    if not smiles:
        return None
    try:
        template = Chem.MolFromSmiles(smiles)  # type: ignore[attr-defined]
        if template is None:
            return None
        rwmol = _build_rwmol(elements, coords)
        n = rwmol.GetNumAtoms()
        positions = [rwmol.GetConformer().GetAtomPosition(i) for i in range(n)]
        # All bonds start AROMATIC: in substructure matching an aromatic
        # query atom only matches an aromatic target atom, and the template
        # is aromatic, so an all-single target never matches.
        for i, j in _distance_bond_graph(elements, positions):
            rwmol.AddBond(i, j, Chem.BondType.AROMATIC)  # type: ignore[attr-defined]
        # the function returns a NEW mol with the template's orders; the
        # input stays untouched
        ordered = AllChem.AssignBondOrdersFromTemplate(  # type: ignore[attr-defined]
            template, rwmol.GetMol()
        )
        Chem.SanitizeMol(ordered)  # type: ignore[attr-defined]
        return ordered
    except Exception:  # noqa: BLE001 - template/crystal graph mismatch
        return None
