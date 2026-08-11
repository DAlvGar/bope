"""Benchmark: bope vs the existing bond-order perception tools.

Four methods are compared on identical input (plain ``(element, xyz)``
atom lists, no CCD templates - everything runs offline):

1. **bope** - the full pipeline (geometry perception, with OpenBabel
   and distance as internal fallbacks): ``perception.perceive_bond_orders``.
2. **OpenBabel** - ``PerceiveBondOrders`` through the package's PDB
   round-trip wrapper.
3. **RDKit rdDetermineBonds** - ``rdDetermineBonds.DetermineBonds`` on an
   RDKit mol with coordinates but no bonds.
4. **distance baseline** - covalent radii + tolerance, then sanitization
   (the package's ``perceive_bond_orders_distance``).

Three axes are measured:

* **Synthetic recoverability vs noise**: every neutral corpus molecule is
  embedded with ETKDG (randomSeed=42), perturbed with isotropic Gaussian
  noise at bond-RMS levels 0.0 / 0.03 / 0.07 / 0.14 / 0.28 / 0.5 A, and
  each method's output is checked for exact bond-graph match, formula
  correctness, sanitization success, and ``AddHs`` success (no
  over-valent atoms).  187 neutral molecules per level.
* **Charged corpus** (10 molecules, 0 noise): exact bond-graph match and
  ``AddHs`` success (formula match is not meaningful - no method consumes
  the charge).
* **Crystal ligands** (20 residues of 16 experimental crystal complexes,
  real deposited coordinates, no noise): formula reproduction against the
  RCSB CCD ground truth and ``AddHs`` success.

Results are printed and written to ``results.md`` next to this script.
Re-run with::

    uv run python benchmarks/benchmark.py
"""

from __future__ import annotations

import contextlib
import json
import os
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
)

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDetermineBonds, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

import bope as perception
from bope import corpus
from bope._deps import _ob, _RDKIT_AVAILABLE
from bope.helpers import _build_rwmol

# OpenBabel's failed-kekulization warnings (the N-rich-ring corruption this
# benchmark measures) flood stderr on every corrupted molecule; silence the
# log here - failures are recorded in the tables, not the terminal.  Level 0
# is obErrorLevel::None (some bindings expose the enum, some don't, but the
# integer is stable across OpenBabel versions).
try:
    _ob.obErrorLog.SetOutputLevel(0)  # type: ignore[attr-defined]
except Exception:  # noqa: BLE001 - API drift across OpenBabel versions
    pass

NOISE_LEVELS = [0.0, 0.03, 0.07, 0.14, 0.28, 0.5]
SEED = 42
METHODS = ["bope", "openbabel", "rdDetermineBonds", "distance"]

HERE = os.path.dirname(os.path.abspath(__file__))
# Crystal-ligand atom fixtures (same file the test suite uses):
# {pdb_id: {het: [ [ [element, [x, y, z]], ... ], ... ]}}
FIXTURES = json.load(
    open(
        os.path.join(HERE, "..", "tests", "fixtures", "crystal_ligands.json"),
        encoding="utf-8",
    )
)

# ---------------------------------------------------------------------------
# Harness (identical to tests/test_perception.py)
# ---------------------------------------------------------------------------

_EMBED_CACHE: dict[str, Chem.Mol] = {}


def embed(smi: str) -> Chem.Mol:
    if smi not in _EMBED_CACHE:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=SEED)
        _EMBED_CACHE[smi] = mol
    return _EMBED_CACHE[smi]


def noisy_atoms(
    molH: Chem.Mol, bond_sigma: float, seed: int = SEED
) -> list[tuple[str, tuple[float, float, float]]]:
    """Per-axis sigma = bond_sigma/sqrt(2) so bond_sigma is the bond-RMS."""
    sigma = bond_sigma / 2 ** 0.5
    rng = np.random.default_rng(seed)
    atoms = []
    for a in molH.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        pos = molH.GetConformer().GetAtomPosition(a.GetIdx())
        xyz = tuple(
            float(v) + float(rng.normal(0, sigma)) for v in (pos.x, pos.y, pos.z)
        )
        atoms.append((a.GetSymbol().upper(), xyz))
    return atoms


def bond_graph(mol: Chem.Mol) -> frozenset:
    return frozenset(
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            round(b.GetBondTypeAsDouble(), 1),
        )
        for b in mol.GetBonds()
    )


def formula(mol: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(mol)


def addh_ok(mol: Chem.Mol) -> bool:
    try:
        Chem.AddHs(mol)
        return True
    except Exception:  # noqa: BLE001 - RuntimeError on over-valent atoms
        return False


# ---------------------------------------------------------------------------
# The four methods
# ---------------------------------------------------------------------------


def run_bope(elements, coords):
    atoms = list(zip(elements, coords))
    mol, _strategy = perception.perceive_bond_orders(atoms, resname=None)
    return mol


@contextlib.contextmanager
def _silence_stderr():
    """OpenBabel's C++ code writes some diagnostics (e.g. SDF stereo-wedge
    notices on chiral centers without implicit H's) straight to fd 2,
    bypassing both the error log and sys.stderr.  Divert fd 2 for the
    duration - the benchmark's story is in the tables, not the terminal."""
    saved = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved, 2)
        os.close(saved)
        os.close(devnull)


def run_openbabel(elements, coords):
    return perception.perceive_bond_orders_with_openbabel(elements, coords)


def run_rd_determine_bonds(elements, coords):
    mol = _build_rwmol(elements, coords).GetMol()
    rdDetermineBonds.DetermineBonds(mol, charge=0)
    return mol


def run_distance(elements, coords):
    return perception.perceive_bond_orders_distance(elements, coords)


RUNNERS = {
    "bope": run_bope,
    "openbabel": run_openbabel,
    "rdDetermineBonds": run_rd_determine_bonds,
    "distance": run_distance,
}


def perceive(elements, coords, method):
    """One method call, never raising; None on failure.  Wrapped in the
    fd-2 diversion: bope's own OpenBabel fallback can emit the same
    C-level diagnostics the standalone OpenBabel method does."""
    with _silence_stderr():
        try:
            return RUNNERS[method](elements, coords)
        except Exception:  # noqa: BLE001 - a method may raise on bad input
            return None


def check_neutral(elements, coords, method, ref_graph, ref_formula):
    """(sanitize_ok, graph_ok, formula_ok, addh_ok) for one method call."""
    mol = perceive(elements, coords, method)
    if mol is None:
        return False, False, False, False
    try:
        Chem.SanitizeMol(mol)
        sanitize_ok = True
    except Exception:  # noqa: BLE001
        sanitize_ok = False
    return (
        sanitize_ok,
        bond_graph(mol) == ref_graph,
        formula(mol) == ref_formula,
        addh_ok(mol),
    )


def pct(ok: int, total: int) -> str:
    return f"{ok}/{total} ({100.0 * ok / total:.0f}%)"


# ---------------------------------------------------------------------------
# Axis 1: synthetic recoverability vs noise
# ---------------------------------------------------------------------------


def sweep_neutral():
    """results[method][level] = [sanitize, graph, formula, addh] counts."""
    neutral = [(n, s) for n, s, c in corpus.CORPUS if not c]
    total = len(neutral)
    results = {m: {lv: [0, 0, 0, 0] for lv in NOISE_LEVELS} for m in METHODS}
    for _name, smi in neutral:
        molH = embed(smi)
        ref = Chem.RemoveHs(molH)
        ref_graph, ref_formula = bond_graph(ref), formula(ref)
        for lv in NOISE_LEVELS:
            atoms = noisy_atoms(molH, lv)
            elements = [el for el, _xyz in atoms]
            coords = [xyz for _el, xyz in atoms]
            for m in METHODS:
                s, g, f, h = check_neutral(
                    elements, coords, m, ref_graph, ref_formula
                )
                results[m][lv][0] += s
                results[m][lv][1] += g
                results[m][lv][2] += f
                results[m][lv][3] += h
    return results, total


def sweep_charged():
    """Charged corpus at 0 noise: [graph, addh] counts (formula match is
    not meaningful - none of the methods consume the formal charge)."""
    total = len(corpus.CHARGED)
    results = {m: [0, 0] for m in METHODS}
    for _name, smi, _charged in corpus.CHARGED:
        molH = embed(smi)
        ref_graph = bond_graph(Chem.RemoveHs(molH))
        atoms = noisy_atoms(molH, 0.0)
        elements = [el for el, _xyz in atoms]
        coords = [xyz for _el, xyz in atoms]
        for m in METHODS:
            mol = perceive(elements, coords, m)
            if mol is None:
                continue
            results[m][0] += bond_graph(mol) == ref_graph
            results[m][1] += addh_ok(mol)
    return results, total


def sweep_crystal():
    """Per-method [formula_ok, addh_ok] over the crystal ligand residues."""
    total = 0
    results = {m: [0, 0] for m in METHODS}
    for pdb_id, ligu in corpus.CRYSTAL_EXPECTED.items():
        for het, want in ligu.items():
            for atoms in FIXTURES[pdb_id][het]:
                total += 1
                elements = [el for el, _xyz in atoms]
                coords = [xyz for _el, xyz in atoms]
                for m in METHODS:
                    mol = perceive(elements, coords, m)
                    if mol is None:
                        continue
                    try:
                        Chem.SanitizeMol(mol)
                    except Exception:  # noqa: BLE001 - sanitize is not the
                        pass  # metric here; formula needs it only to kekulize
                    results[m][0] += formula(mol) in want
                    results[m][1] += addh_ok(mol)
    return results, total


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def render(results, total, charged_results, charged_total, crystal_results,
           crystal_total):
    lines = []
    w = lines.append

    w("# Bond-order perception benchmark")
    w("")
    w("Auto-generated by `benchmarks/benchmark.py` - "
      "do not edit by hand.  Re-run with `uv run python "
      "benchmarks/benchmark.py`.")
    w("")
    w("Methodology: every corpus SMILES is embedded with ETKDG "
      "(randomSeed=42), perturbed with isotropic Gaussian coordinate noise "
      "(per-axis sigma = bond-RMS / sqrt(2)), and each method's output is "
      "checked for exact bond-graph match (same atom pairs with the same "
      "orders), exact molecular formula, sanitization success and `AddHs` "
      "success (no over-valent atoms).  `bope` runs its full pipeline "
      "with `resname=None` - no CCD templates, no network.")
    w("")
    w(f"Neutral corpus: {total} molecules per noise level.  The 5 "
      "documented embedder-unfaithful molecules (ZMA, thiophene, "
      "benzophenone, ketoprofen, 1,3-cyclohexadiene) are included in the "
      "denominators - their ETKDG seed-42 embeds are not faithful at zero "
      "noise, so they count against every method equally.")
    w("")

    w("## Exact recovery vs noise")
    w("")
    w("Recovery = exact bond graph AND exact formula AND `AddHs` success.")
    w("")
    w("| method | " + " | ".join(f"{lv:.2f}" for lv in NOISE_LEVELS) + " |")
    w("|" + "---|" * (len(NOISE_LEVELS) + 1))
    for m in METHODS:
        cells = []
        for lv in NOISE_LEVELS:
            s, g, f, h = results[m][lv]
            cells.append(pct(min(g, f, h), total))
        w(f"| {m} | " + " | ".join(cells) + " |")
    w("")

    w("## Metric breakdown vs noise")
    w("")
    w("| method | metric | " + " | ".join(f"{lv:.2f}" for lv in NOISE_LEVELS) + " |")
    w("|" + "---|" * (len(NOISE_LEVELS) + 2))
    metric_names = ["sanitize", "graph", "formula", "AddHs"]
    for m in METHODS:
        for mi, mname in enumerate(metric_names):
            cells = [pct(results[m][lv][mi], total) for lv in NOISE_LEVELS]
            w(f"| {m} | {mname} | " + " | ".join(cells) + " |")
    w("")

    w("## Charged corpus (0 noise)")
    w("")
    w("Exact bond-graph match and `AddHs` success on the 10 charged "
      "molecules.  Formula match is not meaningful - none of the methods "
      "consume the formal charge.")
    w("")
    w("| method | graph | AddHs |")
    w("|---" * 3 + "|")
    for m in METHODS:
        g, h = charged_results[m]
        w(f"| {m} | {pct(g, charged_total)} | {pct(h, charged_total)} |")
    w("")

    w("## Crystal ligands (real coordinates, 0 noise)")
    w("")
    w(f"Formula reproduction against the RCSB CCD ground truth and `AddHs` "
      f"success on the {crystal_total} ligand residues of the 16 eval "
      "complexes (including the N-rich fused heteroaromatics - "
      "staurosporine, ZM241385, caffeine - where OpenBabel's "
      "`PerceiveBondOrders` is known to corrupt ring systems).")
    w("")
    w("| method | formula | AddHs |")
    w("|---" * 3 + "|")
    for m in METHODS:
        f_ok, h_ok = crystal_results[m]
        w(f"| {m} | {pct(f_ok, crystal_total)} | {pct(h_ok, crystal_total)} |")
    w("")
    return "\n".join(lines)


def main() -> None:
    if not _RDKIT_AVAILABLE:
        sys.exit("RDKit not available - cannot run the benchmark")
    print("sweeping neutral corpus...")
    results, total = sweep_neutral()
    print("sweeping charged corpus...")
    charged_results, charged_total = sweep_charged()
    print("sweeping crystal ligands...")
    crystal_results, crystal_total = sweep_crystal()

    report = render(
        results, total, charged_results, charged_total,
        crystal_results, crystal_total,
    )
    out = os.path.join(HERE, "results.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(report + "\n")
    print(report)


if __name__ == "__main__":
    main()
