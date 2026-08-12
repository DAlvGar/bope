# bope - Bond Order Perception Engine

Perceive bond connectivity and bond orders for 3-D ligand coordinates.
Plain `(element, xyz)` atom lists in - RDKit molecule with correct bond
orders out.  No BioPython, no PDB parsing, no network required (unless
you use the CCD template strategy).

```python
from bope import perceive_bond_orders

atoms = [("C", (1.23, 4.56, 7.89)), ("N", (...)), ...]   # your ligand
mol, strategy = perceive_bond_orders(atoms, resname="CFF")

# mol      -> RDKit Mol with 3-D conformer + perceived bond orders
# strategy -> "ccd-template" | "geometry" | "openbabel" | "distance" | ""
```

This is the bond-order perception of the
[3D-PLI-Agent](https://github.com/DAlvGar/3D-PLI-Agent) project, extracted
as a standalone, dependency-light package.

## Install

```bash
pip install bope            # once published to PyPI
# or from source:
uv sync --dev --extra openbabel   # tests + OpenBabel strategy
```

Core dependencies are RDKit and numpy only.  The OpenBabel strategy is
optional (`pip install bope[openbabel]`); without it the pipeline skips
from geometry perception straight to the distance fallback.

## Usage

`perceive_bond_orders(atoms, resname=None, charge=0)` returns
`(mol, strategy)`:

- **atoms** - sequence of `(element, (x, y, z))` tuples; element is an
  upper-case symbol (`"C"`, `"CL"`, `"FE"`, ...).
- **resname** - optional 3-letter HET code.  When given, the RCSB
  Chemical Component Dictionary template is tried first.
- **charge** - accepted for API compatibility; no current strategy
  consumes it (the perception is coordinate-driven).

The strategy string tells you how trustworthy the orders are.  Four
strategies are attempted in order:

| # | strategy | what it does | when it's used |
|---|---|---|---|
| 1 | `ccd-template` | Fetch the authoritative bond orders / tautomer / protonation for the HET code from the [RCSB CCD](https://www.rcsb.org/) and stamp them onto the crystal graph | known HET code + network access (cached per code) |
| 2 | `geometry` | In-house geometric perception: planar rings in aromatic bond-length envelopes scored by a Hückel 4n+2 electron-count judge over per-atom pi assignments; length-threshold orders with chemistry fixups (amidines, carbonyls, exocyclic N/S); valence demotion before sanitization | unknown HET, offline, or template mismatch |
| 3 | `openbabel` | Atom list serialised to PDB and read back through OpenBabel `PerceiveBondOrders` | simple / charged / protonated groups |
| 4 | `distance` | Covalent-radius connectivity, then `Chem.SanitizeMol` upgrades orders where the topology allows | last resort - bare topology |

Stereochemistry is a separate, opt-in call - the bond-order API stays
pure, and stereo never alters the perceived graph:

```python
from bope import perceive_stereochemistry

labeled = perceive_stereochemistry(mol)   # R/S + E/Z from the 3-D coordinates
```

The returned molecule carries tetrahedral and double-bond labels assigned
from the coordinates (the same side-effect stereo OpenBabel emits with its
SDF output), with the same heavy-atom count as the input.

## Why not just use the existing bond-order perception?

Structures deposited in the PDB carry **coordinates, not bonds** - bond
orders and stereochemistry have to be inferred.  The obvious answer is
"call OpenBabel's `PerceiveBondOrders` or RDKit's `rdDetermineBonds` and
be done."  We measured both, plus a pure distance baseline, against bope
on identical input: 187 synthetic molecules embedded with ETKDG (noise
robustness, [`benchmarks/results.md`](benchmarks/results.md)) and 202
real crystal ligands from the RCSB across two resolution tiers, with the
RCSB Chemical Component Dictionary (CCD) canonical SMILES as ground
truth (bond orders + stereo,
[`benchmarks/crystal100/results.md`](benchmarks/crystal100/results.md),
[`results_stereo.md`](benchmarks/crystal100/results_stereo.md)).  Scope
note: every method perceives from **bare coordinates only** - the CCD is
used purely as the ground-truth reference, never as an input (the
`ccd-template` strategy is disabled in these runs; `geometry` is bope's
own perception, no templates, no network, no external models).

**Synthetic - exact recovery** (bond graph *and* formula *and*
sanitizable `AddHs`) at 0 / 0.03 Å coordinate noise:

| method | 0.00 Å | 0.03 Å |
|---|---|---|
| **bope** | **182/187 (97%)** | **173/187 (93%)** |
| OpenBabel `PerceiveBondOrders` | 136/187 (73%) | 133/187 (71%) |
| RDKit `rdDetermineBonds` | 3/187 (2%) | 3/187 (2%) |
| distance baseline | 60/187 (32%) | 51/187 (27%) |

**Real crystal ligands - bond-order recovery** (exact bond graph against
the CCD, 101 complexes per tier):

| method | 1.0-2.0 Å | 2.5-3.0 Å |
|---|---|---|
| **bope geometry** | **87/101 (86%)** | **71/101 (70%)** |
| OpenBabel `PerceiveBondOrders` | 72/101 (71%) | 63/101 (62%) |
| distance baseline | 12/101 (12%) | 11/101 (11%) |

**Real crystal ligands - stereo recovery** (of the CCD stereo-declaring
subset: 66 main-tier, 57 low-res; full-string recovery on the entries
whose perceived bond graph matches the CCD, per-center R/S precision on
every center the CCD declares):

| method | full 1.0-2.0 Å | full 2.5-3.0 Å | centers 1.0-2.0 Å | centers 2.5-3.0 Å |
|---|---|---|---|---|
| **bope geometry + RDKit** | **49/56 (88%)** | **33/41 (80%)** | **223/223 (100%)** | **155/155 (100%)** |
| OpenBabel (SDF stereo) | 43/48 (90%) | 30/38 (79%) | 189/189 (100%) | 114/115 (99%) |
| distance + RDKit | 8/8 (100%) | 7/9 (78%) | 41/41 (100%) | 28/28 (100%) |

The stereo numbers are honest: the only systematic R/S flips are at
tetrahedral phosphate P (CIP ranking flips with P-OH vs P-O-
protonation, which the crystal does not record) and the only OpenBabel
center error in ~300 comparisons is EPY (1C72).  Remaining full-string
misses are phosphate-P flips (9 entries across both tiers, every
declared center otherwise correct), extras - geometry supporting more
stereo than the deposit declares (1J5, 9XR, HUF, 8IX) - or genuine E/Z
coordinate disagreements (OLB, 6W9Z: the crystal is -173 degrees around
the C=C, genuinely E, while the CCD declares Z).  Full per-entry
classification in
[`results_stereo.md`](benchmarks/crystal100/results_stereo.md).

The mainstream options fail for structural reasons, not tuning:

- **`rdDetermineBonds` returns molecules without implicit hydrogens** -
  its formula is systematically wrong (0/187 formula matches even on
  unperturbed input), its aromatic kekulization corrupts ring systems
  (benzene becomes alternating triple bonds), and it raises
  `Final molecular charge does not match input` on many inputs.
- **OpenBabel's `PerceiveBondOrders` corrupts N-rich fused
  heteroaromatics** - staurosporine, ZM241385, caffeine-style scaffolds
  come back with pentavalent carbons and wrong formulas (24 points
  below bope on the synthetic corpus even at zero noise, and 15/8
  points behind on the two crystal tiers).  It also needs a PDB-text
  round-trip.
- The **distance baseline** never assigns aromaticity or formal charges,
  so formulas drift immediately under realistic noise.

Re-run the benchmarks yourself:

```bash
uv run python benchmarks/benchmark.py                                  # synthetic noise sweep
uv run python benchmarks/crystal100/benchmark.py                       # bond orders, main tier
uv run python benchmarks/crystal100/benchmark.py --dataset dataset_res250-300.json
uv run python benchmarks/crystal100/stereo_benchmark.py                # stereo, main tier
uv run python benchmarks/crystal100/stereo_benchmark.py --dataset dataset_res250-300.json
```

**Held out, never seen** - the tables above were measured on the set the
geometry tier was tuned against, so they are optimistic.  To get
generalization numbers the benchmark was re-run on 600 fresh ligands (300
per tier), sampled from the same RCSB universes minus every entry whose
PDB id or HET code appears in either tuning set, as 5 independent
seeded buckets of 60 per tier.  Each bucket is a simple random sample, so
the per-bucket spread is genuine sampling variation - reported as mean +/-
std (n = 5), with pooled counts alongside.  All held-out runs executed
the exact committed perception code; the protocol and per-bucket detail
are in [`results_heldout.md`](benchmarks/crystal100/results_heldout.md).

**Real crystal ligands, held out - bond-order recovery** (formula AND
graph AND AddHs, 300 ligands per tier):

| method | main tuning | main held-out | lowres tuning | lowres held-out |
|---|---|---|---|---|
| **bope geometry** | 87/101 (86%) | **223/300 (74.3% +/- 6.3)** | 71/101 (70%) | **221/300 (73.7% +/- 1.8)** |
| OpenBabel `PerceiveBondOrders` | 72/101 (71%) | 206/300 (68.7% +/- 7.9) | 63/101 (62%) | 194/300 (64.7% +/- 9.2) |
| distance baseline | 12/101 (12%) | 38/300 (12.7% +/- 2.5) | 11/101 (11%) | 32/300 (10.7% +/- 3.5) |
| RDKit `rdDetermineBonds` | 0/101 (0%) | 0/300 (0%) | 0/101 (0%) | 0/300 (0%) |

Two honest caveats from the held-out numbers.  First, the main-tier
tuning set was optimistic for geometry: 86.1% on tuning vs 74.3% held
out (OpenBabel: 71.3% vs 68.7%).  The geometry advantage over OpenBabel
shrinks from 15 points to about 6 on never-seen main-tier data, and
stays about 9 points on the low-res tier - real, but smaller than the
tuning tables suggested.  Second, the low-res tier held up (70.3% ->
73.7%): the low-res tuning set happened to be harder than its universe
average, not easier.

**Real crystal ligands, held out - stereo recovery** (of the CCD
stereo-declaring subset; full-string on entries whose perceived bond
graph matches the CCD):

| method | full main | full lowres | R/S centers main | R/S centers lowres |
|---|---|---|---|---|
| **bope geometry + RDKit** | **114/121 (94.2% +/- 2.6)** | **103/125 (82.4% +/- 9.1)** | **384/386 (99.5%)** | **395/396 (99.7%)** |
| OpenBabel (SDF stereo) | 108/121 (89.3% +/- 6.5) | 95/114 (83.3% +/- 8.1) | 403/405 (99.5%) | 392/393 (99.7%) |
| distance + RDKit | 18/22 (81.8%) | 22/22 (100%) | 96/98 (98.0%) | 108/108 (100%) |

Stereo generalizes well: full-string recovery is as good or better held
out than on the tuning set, and per-center R/S stays above 99% on both
tiers.  The only E/Z miss on the main tier (11/12) is one of the
coordinate-vs-CCD conflicts the tuning set already showed - the deposit's
geometry contradicts its declared stereo, and every method reads the
coordinates.  `rdDetermineBonds` cannot perceive stereo at all on
hydrogen-less PDB input (0 stereo-comparable entries in all 600).

## Known limitations

**Coordinate-artifact class (H4)** - entries where the deposited model
contradicts its own CCD bond orders: the measured bond length says one
order, the CCD says another.  bope reports what the coordinates say,
which is what the crystal actually contains; the CCD ground truth is
wrong at these positions (or the density is too ambiguous to refine the
bond).  Documented cases from the benchmark:

| pdb | het | bond | measured | CCD says | bope says |
|---|---|---|---|---|---|
| 9BF0 | UTP | ribose C2'-C3' | 1.338 Å | single | double |
| 1UY7 | PU4 | butyl C-C | 1.344 Å | single | double |
| 3I7E | DJR | sulfone S=O | 1.576 Å | double | single |

These are kept as documented failures, not warnings: the API returns
`(mol, strategy)` with no warning channel, and the perception is
deliberately coordinate-driven - a length-contradiction alarm would fire
precisely where the tool did its job.  The full failure tables (including
the metalloporphyrin HEM/HEC/ZNH and phosphate NDP/NAP cases, whose
metal / protonation states the coordinates do not record) are committed
with the benchmark results.

**Stereo limitations** mirror the same principle - the coordinates are
the evidence: phosphate-P R/S flips are protonation-sensitive (the
crystal does not record P-OH vs P-O-), and an E/Z "disagreement" (OLB,
6W9Z) can be the deposit's own geometry contradicting its declared
stereochemistry.  See the stereo results tables for the per-entry
classification.

## Validation

The test suite (`tests/`, offline - the CCD cache is seeded, no network)
pins, on the current RDKit / ETKDG seed-42 geometries:

- **Synthetic recoverability** - exact bond-graph + formula + `AddHs`
  recovery across the ~190-molecule corpus (aromatic and fused systems,
  N-rich heterocycles, tautomers, carbonyls, amides, nitriles, alkynes,
  halogens, strained rings, sugars, amino acids, nucleobases, drugs) at
  0 and 0.03 Å bond-RMS noise.  Documented exclusions for the five
  embedder-unfaithful and ten cutoff-boundary molecules live in
  `src/bope/corpus.py`.
- **Crystal ground truth** - all 20 ligand residues of the 16 complexes
  (extracted once into `tests/fixtures/crystal_ligands.json`) reproduce
  the RCSB CCD formula on both the CCD-template and the offline
  geometry paths, with no over-valent atoms.
- **Charged molecules** degrade gracefully; quaternary nitrogen
  recovers exactly; tautomer H-placement is checked per-atom; unknown
  HET codes, disconnected fragments, missing atoms, single atoms and
  empty input never crash.
- **Stereochemistry** - `perceive_stereochemistry` is pinned on
  PubChem-verified expectations (`C[C@H](O)CC` -> S, `C[C@@H](O)CC` ->
  R), E/Z double bonds, full-pipeline integration with the geometry
  strategy, and the no-conformer / no-RDKit degradations.

## Project structure

```
bope/
├── src/bope/          # the package
│   ├── __init__.py    # public API: perceive_bond_orders()
│   ├── ccd.py         # RCSB CCD template matching (network, cached)
│   ├── geometry.py    # in-house geometric perception (Hückel judge)
│   ├── openbabel.py   # OpenBabel PerceiveBondOrders fallback
│   ├── distance.py    # covalent-radius distance baseline
│   ├── stereo.py      # perceive_stereochemistry: R/S + E/Z from 3-D
│   ├── tables.py      # length thresholds, valences, Hückel sets
│   ├── corpus.py      # validation corpus + measured exclusions
│   └── _deps.py       # lazy-imported optional dependencies
├── tests/             # offline test suite (+ crystal fixtures)
├── benchmarks/        # synthetic noise sweep + crystal100 (bond orders
│                      # and stereo vs CCD, both resolution tiers)
└── pyproject.toml
```

## License

MIT - see [LICENSE](LICENSE).
