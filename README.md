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

## How it works

`perceive_bond_orders` tries four strategies in order - from most to
least trustworthy - and returns the first one that succeeds.  The
`strategy` string in the return value tells you which tier produced the
molecule.

| # | strategy | what it does | when it's used |
|---|---|---|---|
| 1 | `ccd-template` | Fetch the authoritative bond orders / tautomer / protonation for the HET code from the [RCSB CCD](https://www.rcsb.org/) and stamp them onto the crystal graph | known HET code + network access (cached per code) |
| 2 | `geometry` | In-house geometric perception: planar rings in aromatic bond-length envelopes scored by a Hückel 4n+2 electron-count judge over per-atom pi assignments; length-threshold orders with chemistry fixups (carbonyls, amidines, nitro, phosphate P=O, sulfonamides); valence demotion before sanitization | unknown HET, offline, or template mismatch |
| 3 | `openbabel` | Atom list serialised to PDB and read back through OpenBabel `PerceiveBondOrders` | simple / charged / protonated groups |
| 4 | `distance` | Covalent-radius connectivity, then `Chem.SanitizeMol` upgrades orders where the topology allows | last resort - bare topology |

### Extending the geometry fixup rulebook

The geometry tier's functional-group corrections (carbonyl rescue +
ester-O protection, amidine, nitro, phosphate P=O, sulfonamide) live in
`src/bope/fixups.py` as declarative `FixupRule` data evaluated by the
generic `FixupEngine` - adding a functional group is one rule entry
plus a test, no new code path.

A rule is a trigger and an action:

- **Trigger** - the center element (with optional aromaticity and
  minimum degree), named neighbor groups (`NbrGroup`: element /
  aromatic / terminal filters with `min` / `max` / `exact` count
  bounds), and gate predicates: required groups (`require`), fallback
  triggers (`require_or`), forbidden neighbor elements
  (`exclude_nbrs`), measured-length caps with gate semantics
  (`max_len` - the whole rule is skipped if any member exceeds the
  cutoff), molecule-state guards (`no_double_to`).
- **Action** - applied in order: `make_single`, then `make_double`
  (`shortest:` or `all:` of one group, with `action_len` as a filter -
  members beyond the cutoff are left untouched, unlike the gate
  semantics of `max_len`), then formal `charges`.  `only_if_single`
  keeps `make_double` from upgrading an existing double.

Worked example - the nitro rule:

```python
FixupRule(
    name="nitro",
    center="N",                      # trigger: an N with ...
    min_degree=3,                    # ... at least 3 graph neighbours
    groups=(NbrGroup("term_o", "O", terminal=True, exact=2),),
    require=("term_o",),             # ... exactly 2 terminal O's,
    max_len={"term_o": 1.45},        # both within 1.45 A (gate)
    make_double="shortest:term_o",   # the shorter N-O goes double
    make_single=("term_o_others", "non_single_non_arom"),
    charges=(("center", 1), ("term_o_others", -1)),
    note="crystal citations + threshold reasoning",  # keep it
)
```

Rules run in order, each on the molecule as the previous rule left it -
order is load-bearing.  Every rule carries a `note` with the crystal
citations behind its thresholds; keep them, they are the scientific
record.  Add a focused test in `tests/test_fixups.py` with the
`_make_mol` / `_make_engine` harness: build the smallest molecule that
exercises the trigger, assert the resulting bonds / charges and the
fired-rule names, and pin the near-miss that must NOT fire.  Then
re-run the full suite and `benchmarks/benchmark.py` - the committed
parity targets in `benchmarks/results.md` verify the rulebook against
the corpus.  Full field reference: the `FixupRule` / `NbrGroup`
docstrings in `src/bope/fixups.py`.

## Install

```bash
pip install bope
# or from source, with the optional OpenBabel extra (strategy tier 3 +
# the tests and benchmarks that compare against OpenBabel):
uv sync --dev --extra openbabel
```

Core dependencies are RDKit and numpy only.  The `openbabel` strategy is
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
python benchmarks/crystal100/run_yuelbond.py  # YuelBond head-to-head (needs the
#   yuel_bond checkout + torch env, see the script's docstring)
```

**Held out, never seen** - the tables above were measured on the set the
geometry tier was tuned against, so they are optimistic.  To get
generalization numbers the benchmark was re-run on fresh ligands (300
per tier, two generations; gen2 seed 43) sampled from the same RCSB
universes minus every entry whose PDB id or HET code appears in either
tuning set or a prior generation, as 5 independent seeded buckets of 60
per tier.  Each bucket is a simple random sample, so the per-bucket
spread is genuine sampling variation - reported as mean +/- std (n =
5), with pooled counts alongside.  All held-out runs executed the exact
committed perception code; the protocol, per-bucket detail and
verification are in
[`results_heldout_gen2.md`](benchmarks/crystal100/results_heldout_gen2.md)
(gen1 runs exist as a documented development loop, `results_heldout.md`).
The YuelBond row comes from `run_yuelbond.py`, which runs the released
model with its published weights (Zenodo record 15353365) head-to-head
on the same input - the first evaluation of any ML bond-perception
model on experimental PDB coordinates.

**Real crystal ligands, held out - bond-order recovery** (formula AND
graph AND AddHs, 300 ligands per tier):

| method | main tuning | main held-out | lowres tuning | lowres held-out |
|---|---|---|---|---|
| **bope geometry** | 89/101 (88.1%) | **231/300 (77.0% +/- 5.1)** | 74/101 (73.3%) | **227/300 (75.7% +/- 9.2)** |
| OpenBabel `PerceiveBondOrders` | 72/101 (71.3%) | 207/300 (69.0% +/- 10.1) | 63/101 (62.4%) | 194/300 (64.7% +/- 8.1) |
| distance baseline | 12/101 (11.9%) | 25/300 (8.3% +/- 3.7) | 11/101 (10.9%) | 12/300 (4.0% +/- 1.9) |
| RDKit `rdDetermineBonds` | 0/101 (0%) | 0/300 (0%) | 0/101 (0%) | 0/300 (0%) |
| YuelBond (GEOM-trained GNN, released weights) | - | 111/300 (37.0% +/- 3.6) | - | 83/300 (27.7% +/- 3.5) |

Two honest caveats from the held-out numbers.  First, the main-tier
tuning set was optimistic for geometry: 88.1% on tuning vs 77.0% held
out - an 11-point gap the paper quantifies with a 95% CI (3.2-19.0)
(OpenBabel: 71.3% vs 69.0%).  The geometry advantage over OpenBabel
shrinks from ~17 tuning points to 8.0 held-out points on the main tier
(t-CI -1.2 to 17.2: not distinguishable from zero at the bucket level)
and 11.0 points on the low-res tier (t-CI -0.4 to 22.4).  Second, the
low-res tier held up (73.3% -> 75.7%, within noise): the low-res tuning
set happened to be harder than its universe average, not easier.
Third, the ML baseline does not transfer: the same YuelBond model that
reports ~98% F1 on computed GEOM geometries recovers 37.0%/27.7% of
experimental PDB ligands, and ~40% of its outputs fail RDKit
sanitization outright.

**Real crystal ligands, held out - stereo recovery** (of the CCD
stereo-declaring subset; full-string on entries whose perceived bond
graph matches the CCD):

| method | full main | full lowres | R/S centers main | R/S centers lowres |
|---|---|---|---|---|
| **bope geometry + RDKit** | **96/99 (97.0% +/- 3.0)** | **90/107 (84.1% +/- 9.5)** | **238/239 (99.6%)** | **311/318 (97.8%)** |
| OpenBabel (SDF stereo) | 89/93 (95.7% +/- 2.8) | 91/101 (90.1% +/- 4.2) | 260/260 (100%) | 288/289 (99.7%) |
| distance + RDKit | 14/15 (93.3%) | 6/7 (85.7%) | 37/37 (100%) | 43/43 (100%) |

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
- **Smoke tests** - `tests/test_smoke.py` re-runs the geometry tier
  over both committed crystal100 fix datasets (coordinates embedded,
  fully offline) and asserts the aggregate metrics still match the
  committed sidecars (`benchmarks/crystal100/results_dataset_fix_*.json`),
  so any perception change that shifts recovery on real deposited
  coordinates fails the suite in seconds instead of requiring a full
  evaluation re-run.  The sidecar is the source of truth: a deliberate,
  validated change regenerates it with the harness and commits both
  together.  `benchmarks/crystal100/perf_smoke.py` (a script, not a
  pytest test - wall-clock asserts would be flaky) times the tier
  against the recorded 1.78 ms/perception baseline and exits nonzero
  on an order-of-magnitude regression.

## Project structure

```
bope/
├── src/bope/          # the package
│   ├── __init__.py    # public API: perceive_bond_orders()
│   ├── ccd.py         # RCSB CCD template matching (network, cached)
│   ├── geometry.py    # in-house geometric perception (Hückel judge)
│   ├── openbabel.py   # OpenBabel PerceiveBondOrders fallback
│   ├── distance.py    # covalent-radius distance baseline
│   ├── helpers.py     # shared plumbing: distance graph, plane RMS,
│   │                  # RWMol construction
│   ├── stereo.py      # perceive_stereochemistry: R/S + E/Z from 3-D
│   ├── tables.py      # length thresholds, valences, Hückel sets
│   ├── corpus.py      # validation corpus + measured exclusions
│   └── _deps.py       # lazy-imported optional dependencies
├── tests/             # offline test suite (+ crystal fixtures)
├── benchmarks/        # synthetic noise sweep + crystal100 (bond orders
│                      # and stereo vs CCD, both resolution tiers);
│                      # perf_smoke.py = wall-clock regression guard
└── pyproject.toml
```

## License

MIT - see [LICENSE](LICENSE).
