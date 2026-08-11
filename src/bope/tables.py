"""Shared constants for bond-order perception: radii and bond-length tables.

All lengths are in Angstroms.  The tables are shared by the geometry
perception, the distance fallback and the distance bond graph used by the
CCD template path.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Covalent radii table (Alvarez 2008, doi:10.1039/b801115j)
# ---------------------------------------------------------------------------

#: Single-bond covalent radii in Angstroms, upper-case element keys.
_COVALENT_RADII: dict[str, float] = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "CL": 1.02,
    "BR": 1.20,
    "I": 1.39,
    "FE": 1.32,
    "ZN": 1.22,
    "CA": 1.76,
    "MG": 1.41,
    "NA": 1.66,
    "K": 2.03,
    "MN": 1.61,
    "CU": 1.32,
    "CO": 1.26,
    "NI": 1.24,
    "SE": 1.20,
    "SI": 1.11,
    "B": 0.84,
}

_COVALENT_RADII_DEFAULT: float = 0.77  # generic fallback (~= C)

#: Extra tolerance (A) added to the sum of covalent radii when detecting bonds.
_BOND_TOLERANCE: float = 0.40

#: Double-bond covalent radii in Angstroms.  Pairs whose summed radius
#: (plus :data:`_DOUBLE_BOND_TOLERANCE`) exceeds the observed distance are
#: assigned a ``DOUBLE`` bond in the distance-based fallback.
_DOUBLE_COVALENT_RADII: dict[str, float] = {
    "C": 0.67,
    "N": 0.60,
    "O": 0.57,
    "S": 0.94,
    "P": 1.00,
    "B": 0.78,
}

#: Tolerance (A) added to the sum of double-bond radii when deciding whether
#: to assign a ``DOUBLE`` bond in the distance-based fallback.  Kept tight
#: (0.03 A) so that long C-O single bonds (~1.34 A) are not mislabeled.
_DOUBLE_BOND_TOLERANCE: float = 0.03

# ---------------------------------------------------------------------------
# Geometry perception tables
# ---------------------------------------------------------------------------

#: (el1, el2) -> (single_max, double_max) bond lengths in Angstrom; a
#: ``None`` double_max means the pair never takes a double bond.
_BOND_ORDER_TABLE: dict[tuple[str, str], tuple[float, float | None]] = {
    ("C", "C"): (1.50, 1.38), ("C", "N"): (1.48, 1.33), ("C", "O"): (1.44, 1.30),
    ("C", "S"): (1.82, 1.72), ("C", "P"): (1.85, 1.75),
    ("C", "F"): (1.40, None), ("C", "CL"): (1.80, None), ("C", "BR"): (1.95, None),
    ("C", "I"): (2.15, None), ("N", "N"): (1.47, 1.33), ("N", "O"): (1.45, 1.28),
    ("N", "S"): (1.80, 1.70), ("N", "P"): (1.85, 1.75), ("O", "O"): (1.48, 1.28),
    ("O", "S"): (1.75, 1.55), ("O", "P"): (1.75, 1.55), ("S", "S"): (2.05, None),
    # (O,S) double_max 1.55 (not 1.50): sulfone S=O measures 1.43-1.47 A in
    # crystals but 1.50-1.53 A in ETKDG embeds (a +0.06 systematic bias),
    # while a genuine S-O single (sulfonate / sulfate ester) is 1.57+ - the
    # 0.10 A gap keeps the raised cutoff safe for both regimes.
}
#: unmistakably-short bonds: nitriles / alkynes (a crystal C=C would need to
#: be badly wrong to land here)
_TRIPLE_BOND_TABLE: dict[tuple[str, str], float] = {
    ("C", "C"): 1.26, ("C", "N"): 1.19, ("N", "N"): 1.16,
}
#: generous bounds: low-resolution crystal structures elongate aromatic
#: bonds by up to ~0.06 A, and triazole N-N's are SHORT (1.21-1.24) rather
#: than long - planarity + Huckel count still gate the chemistry
_AROMATIC_ENVELOPE: dict[tuple[str, str], tuple[float, float]] = {
    ("C", "C"): (1.28, 1.50), ("C", "N"): (1.28, 1.47), ("C", "O"): (1.28, 1.44),
    ("C", "S"): (1.66, 1.80), ("N", "N"): (1.18, 1.44), ("N", "O"): (1.28, 1.44),
    ("C", "P"): (1.75, 1.87),
}

#: a C-O bond at or below this length is a carbonyl even in crystal
#: structures, whose carbonyls often refine to 1.34-1.36 A (caffeine C2/C6
#: in 3RFM).  The pi-count C=O test uses the stricter 1.30 (the length-rule
#: double cutoff) so a phenol/enol C-O at 1.34-1.40 is not misread as a
#: carbonyl; the N-pyridone/amide discriminator and the carbonyl rescue use
#: this generous bound instead.
_CRYSTAL_CARBONYL: float = 1.40

#: a candidate aromatic ring may exceed its envelope upper bound by at most
#: ``_AROMATIC_SLACK_BOND`` per bond and ``_AROMATIC_SLACK_RING`` in total.
#: Low-resolution crystals elongate one or two aromatic bonds slightly (the
#: 07L coumarin lactone ring at 1.0 A refines its two C-O to 1.444/1.462 vs
#: the 1.44 upper, and one C-C to 1.514 vs 1.50), while a uniformly
#: elongated saturated ring (several bonds at 1.52+) must stay rejected -
#: the Huckel gate still decides the chemistry.  The per-bond cap is the
#: tighter of the two bounds: a single bond far past its envelope (the 44L
#: imidazolidinone C-C at 1.539, +0.039 over the 1.50 upper) is an sp3
#: single, not an elongated aromatic bond - such excess must be spread over
#: at least two bonds (07L's max single-bond excess is +0.022).
_AROMATIC_SLACK_BOND: float = 0.03
_AROMATIC_SLACK_RING: float = 0.06

#: when the full aromatic set fails to kekulize, candidate rings admitted
#: only through slack at or below this total excess are dropped and the
#: perception retried.  Noise-level slack (1D1's pyridinone ring at +0.003)
#: can admit a ring the chemistry rejects, and its fused neighbour makes the
#: set unkekulizable; true low-resolution aromatics sit well above the line
#: (NDP's pyridinium ring +0.012, 07L's coumarin +0.041) and survive.
_AROMATIC_SLACK_DROP: float = 0.005
_MAX_VALENCE: dict[str, int] = {"C": 4, "N": 3, "O": 2, "S": 6, "P": 5}
_HUCKEL: set[int] = {2, 6, 10, 14, 18}
