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

## Why not just use the existing bond-order perception?

Structures deposited in the PDB carry **coordinates, not bonds** - bond
orders have to be inferred.  The obvious answer is "call OpenBabel's
`PerceiveBondOrders` or RDKit's `rdDetermineBonds` and be done."  We
measured both, plus a pure distance baseline, against bope on identical
input (187 synthetic molecules embedded with ETKDG, 20 ligand residues of
16 experimental crystal complexes, and a charged corpus - full
methodology and per-metric tables in
[`benchmarks/results.md`](benchmarks/results.md)):

**Exact recovery** (bond graph *and* formula *and* sanitizable `AddHs`)
at 0 / 0.03 Å coordinate noise:

| method | 0.00 Å | 0.03 Å |
|---|---|---|
| **bope** | **182/187 (97%)** | **173/187 (93%)** |
| OpenBabel `PerceiveBondOrders` | 136/187 (73%) | 133/187 (71%) |
| RDKit `rdDetermineBonds` | 3/187 (2%) | 3/187 (2%) |
| distance baseline | 60/187 (32%) | 51/187 (27%) |

**Crystal ligands** (real deposited coordinates, formula against the RCSB
CCD ground truth): bope 20/20 (100%), OpenBabel 15/20 (75%), distance
3/20 (15%), rdDetermineBonds 0/20 (0%).

The mainstream options fail for structural reasons, not tuning:

- **`rdDetermineBonds` returns molecules without implicit hydrogens** -
  its formula is systematically wrong (0/187 formula matches even on
  unperturbed input), its aromatic kekulization corrupts ring systems
  (benzene becomes alternating triple bonds), and it raises
  `Final molecular charge does not match input` on many inputs.
- **OpenBabel's `PerceiveBondOrders` corrupts N-rich fused
  heteroaromatics** - staurosporine, ZM241385, caffeine-style scaffolds
  come back with pentavalent carbons and wrong formulas (75% on the
  crystal set; 27% of its outputs fail graph recovery even at zero
  noise).  It also needs a PDB-text round-trip.
- The **distance baseline** never assigns aromaticity or formal charges,
  so formulas drift immediately under realistic noise.

Re-run the benchmark yourself:

```bash
uv run python benchmarks/benchmark.py
```

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

## Project structure

```
bope/
├── src/bope/          # the package
│   ├── __init__.py    # public API: perceive_bond_orders()
│   ├── ccd.py         # RCSB CCD template matching (network, cached)
│   ├── geometry.py    # in-house geometric perception (Hückel judge)
│   ├── openbabel.py   # OpenBabel PerceiveBondOrders fallback
│   ├── distance.py    # covalent-radius distance baseline
│   ├── tables.py      # length thresholds, valences, Hückel sets
│   ├── corpus.py      # validation corpus + measured exclusions
│   └── _deps.py       # lazy-imported optional dependencies
├── tests/             # offline test suite (+ crystal fixtures)
├── benchmarks/        # re-runnable comparison vs OpenBabel / RDKit
└── pyproject.toml
```

## License

MIT - see [LICENSE](LICENSE).
