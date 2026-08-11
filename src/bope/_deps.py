"""Optional-dependency availability flags for the perception package.

RDKit is required for every strategy; OpenBabel is required only for the
OpenBabel strategy.  Each strategy module degrades gracefully when its
dependency is missing, so the package imports without either installed.
"""

from __future__ import annotations

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdGeometry

    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    _RDKIT_AVAILABLE = False

try:
    from openbabel import openbabel as _ob

    _OPENBABEL_AVAILABLE = True
except ImportError:  # pragma: no cover
    _OPENBABEL_AVAILABLE = False
