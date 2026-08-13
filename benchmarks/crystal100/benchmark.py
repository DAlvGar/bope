"""Crystal-set benchmark: the geometry tier on 100 unique RCSB ligands.

The question this benchmark answers: how good is the in-house GEOMETRY
perception - the main value of this package - on real deposited
coordinates, measured against the authoritative bond orders of the RCSB
Chemical Component Dictionary?

Dataset (``dataset.json`` next to this script by default, built by
``fetch_dataset.py``): 100 PDB entries, protein-only, exactly one
non-polymer entity, every ligand a distinct HET code, every ligand
residue complete (heavy-atom count matches the CCD SMILES exactly).
The resolution tier comes from the fetch parameters (``--min-res`` /
``--max-res``); each entry records its resolution.  Ground truth per
ligand: the CCD canonical SMILES (formula and bond graph).

A second tier (poorer resolution, e.g. 2.5-3.0 A) is fetched into its
own dataset file and benchmarked the same way::

    uv run python benchmarks/crystal100/fetch_dataset.py \
        --min-res 2.5 --max-res 3.0 --out dataset_res250-300.json
    uv run python benchmarks/crystal100/benchmark.py \
        --dataset dataset_res250-300.json

Keeping tiers in separate files isolates the resolution effect: the
geometry tier decides bond orders from length thresholds, so coordinate
noise at lower resolution is exactly the stress this benchmark measures.
(Note for the low-res tier: the CCD ground truth is still authoritative -
the RCSB records the chemistry - but noisy coordinates may no longer
support that chemistry faithfully; a failure there is either a
perception threshold exceeded or coordinates that underdetermine the
bond orders.)

Three methods, identical input:

1. **geometry** - the in-house perception called DIRECTLY
   (``perceive_bond_orders_geometric``): the CCD template is never
   consulted and the OpenBabel / distance fallbacks are not used.  A
   ``None`` result or a wrong formula here is a genuine geometry-tier
   failure.
2. **openbabel** - ``PerceiveBondOrders`` via the PDB round-trip.
3. **distance** - covalent-radius connectivity + sanitization.

Metrics per ligand and method:

* ``formula`` - molecular formula matches the CCD reference (formal
  charges neutralised on both sides: the perception consumes no charge
  argument, and a charged reference would otherwise fail every method
  identically)
* ``graph`` - the perceived bond graph matches the reference up to
  isomorphism and bond orders (canonical SMILES equality, charges
  neutralised, stereo-insensitive - the CCD ground truth is
  stereo-agnostic and stereo has its own benchmark): a tautomer
  difference moves the H to a different atom, so it FAILS the graph
  check while keeping the formula
* ``exact`` - canonical SMILES equality, charges and tautomer included,
  stereo-insensitive
* ``AddHs`` - no over-valent atoms
* ``recovery`` = formula AND graph AND AddHs

Tautomer detection is classified separately: ligands where the formula
matches but the graph does not are listed in the tautomer section -
the geometry tier assigns N-H placement (Hückel judge) from the
coordinates, and the coordinates can support a different tautomer than
the CCD canonical SMILES records.

Every geometry-tier failure is collected in a table (nothing is
corrected in this pass) - the hypothesis-for-fixing step consumes that
table.  Re-run with::

    uv run python benchmarks/crystal100/benchmark.py
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src")
)

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds, rdMolDescriptors

import bope as perception
from bope._deps import _ob
from bope.geometry import perceive_bond_orders_geometric
from bope.helpers import _build_rwmol, _distance_bond_graph

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))


def _env_info() -> str:
    """Software versions + commit that produced a results file - the
    reproducibility record a published benchmark must carry."""
    import subprocess

    parts = []
    try:
        import rdkit

        parts.append(f"RDKit {rdkit.__version__}")
    except Exception:  # noqa: BLE001
        parts.append("RDKit unknown")
    try:
        import openbabel.openbabel as _ob_ver

        parts.append(f"OpenBabel {_ob_ver.OBReleaseVersion()}")
    except Exception:  # noqa: BLE001
        parts.append("OpenBabel unknown")
    try:
        import importlib.metadata

        parts.append(
            f"openbabel-wheel {importlib.metadata.version('openbabel-wheel')}"
        )
    except Exception:  # noqa: BLE001
        pass
    try:
        import numpy

        parts.append(f"numpy {numpy.__version__}")
    except Exception:  # noqa: BLE001
        parts.append("numpy unknown")
    parts.append(f"Python {sys.version.split()[0]}")
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=HERE, check=True,
        ).stdout.strip()
        parts.append(f"commit {commit}")
    except Exception:  # noqa: BLE001
        parts.append("commit unknown")
    return ", ".join(parts)


def load_dataset(name: str) -> list[dict]:
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        sys.exit(f"dataset not found: {path}")
    return json.load(open(path, encoding="utf-8"))

# OpenBabel's failed-kekulization warnings on N-rich rings flood stderr;
# failures are recorded in the tables, not the terminal.
with contextlib.suppress(Exception):  # API drift across OpenBabel versions
    _ob.obErrorLog.SetOutputLevel(0)  # type: ignore[attr-defined]


@contextlib.contextmanager
def _silence_stderr():
    """OpenBabel's C++ code writes some diagnostics straight to fd 2,
    bypassing both the error log and sys.stderr."""
    saved = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved, 2)
        os.close(saved)
        os.close(devnull)


def _neutral(mol: Chem.Mol) -> Chem.Mol:
    """Copy with all formal charges zeroed (no change to the graph).

    Implicit hydrogens are recomputed for the neutral form: a carboxylate
    [O-] must count as O-H (one implicit H), not an O keeping the stale
    0-H count the charge left behind.  A quaternary N stays effectively
    [N+] (valence 4, 0 H), a pyridinium N+ becomes a neutral pyridine N.
    No re-parse: a zeroed pyridinium/quaternary state is chemically
    impossible, so the SMILES writer would emit a string the parser
    rejects.  UpdatePropertyCache alone refreshes the implicit-H counts.
    """
    m = Chem.Mol(mol)
    for a in m.GetAtoms():
        a.SetFormalCharge(0)
    for a in m.GetAtoms():
        a.SetNoImplicit(False)
    m.UpdatePropertyCache(strict=False)
    return m


def _metrics(mol: Chem.Mol | None, ref: Chem.Mol):
    """(sanitize, formula, graph, exact, addh) for one perceived mol."""
    if mol is None:
        return False, False, False, False, False
    try:
        Chem.SanitizeMol(mol)
        sanitize = True
    except Exception:  # noqa: BLE001 - sanitization failure is a metric
        sanitize = False
    try:
        formula = rdMolDescriptors.CalcMolFormula(_neutral(mol)) == (
            rdMolDescriptors.CalcMolFormula(_neutral(ref))
        )
    except Exception:  # noqa: BLE001
        formula = False
    try:
        graph = Chem.MolToSmiles(_neutral(mol), isomericSmiles=False) == (
            Chem.MolToSmiles(_neutral(ref), isomericSmiles=False)
        )
    except Exception:  # noqa: BLE001
        graph = False
    try:
        exact = Chem.MolToSmiles(mol, isomericSmiles=False) == (
            Chem.MolToSmiles(ref, isomericSmiles=False)
        )
    except Exception:  # noqa: BLE001
        exact = False
    try:
        Chem.AddHs(mol)
        addh = True
    except Exception:  # noqa: BLE001 - RuntimeError on over-valent atoms
        addh = False
    return sanitize, formula, graph, exact, addh


def run_geometry(entry):
    elements = [el for el, _xyz in entry["atoms"]]
    coords = [xyz for _el, xyz in entry["atoms"]]
    rwmol = _build_rwmol(elements, coords)
    n = rwmol.GetNumAtoms()
    positions = [rwmol.GetConformer().GetAtomPosition(i) for i in range(n)]
    graph = _distance_bond_graph(elements, positions)
    mol, err = perceive_bond_orders_geometric(elements, coords, graph)
    return mol, err


def run_openbabel(entry):
    elements = [el for el, _xyz in entry["atoms"]]
    coords = [xyz for _el, xyz in entry["atoms"]]
    with _silence_stderr():
        return perception.perceive_bond_orders_with_openbabel(elements, coords), None


def run_distance(entry):
    elements = [el for el, _xyz in entry["atoms"]]
    coords = [xyz for _el, xyz in entry["atoms"]]
    return perception.perceive_bond_orders_distance(elements, coords), None


def run_rdb(entry):
    """RDKit's rdDetermineBonds (the xyz2mol port): perceives bond
    orders from coordinates via the connect-the-dots covalent-radius
    method, then a valence/charge-constrained assignment.  Identical
    call to the synthetic benchmark (defaults, charge=0), plus a finite
    iteration bound: its default maxIterations=0 means no limit, and a
    pathological PDB ligand made the assignment loop indefinitely (the
    exception is then recorded as a failed perception)."""
    elements = [el for el, _xyz in entry["atoms"]]
    coords = [xyz for _el, xyz in entry["atoms"]]
    mol = _build_rwmol(elements, coords).GetMol()
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=0, maxIterations=5000)
        return mol, None
    except Exception as exc:  # noqa: BLE001 - a failed perception is a metric
        return None, f"{type(exc).__name__}: {exc}"


METHODS = {
    "geometry": run_geometry,
    "openbabel": run_openbabel,
    "distance": run_distance,
    "rdDetermineBonds": run_rdb,
}


def pct(ok: int, total: int) -> str:
    if total == 0:
        return "0/0 (-)"
    return f"{ok}/{total} ({100.0 * ok / total:.0f}%)"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", default="dataset.json",
                        help="dataset file next to this script "
                             "(default: dataset.json)")
    ds_name = parser.parse_args().dataset
    DATASET = load_dataset(ds_name)

    total = len(DATASET)
    res = [d.get("resolution") for d in DATASET if d.get("resolution")]
    res_band = f"{min(res):.1f}-{max(res):.1f} A" if res else "unknown"
    print(f"dataset: {total} ligands ({res_band}), unique HET codes: "
          f"{len({d['het'] for d in DATASET})}")

    counts = {m: [0, 0, 0, 0, 0] for m in METHODS}  # s, f, g, e, h
    failures = []  # geometry-only, for the fix-hypothesis pass
    tautomers = []  # formula ok, graph differs
    charges = []    # graph ok, exact differs
    rows = []

    for entry in DATASET:
        ref = Chem.MolFromSmiles(entry["smiles"])
        row = {
            "pdb": entry["pdb"], "het": entry["het"],
            "name": entry["name"], "formula": entry["formula"],
        }
        for m in METHODS:
            mol, err = METHODS[m](entry)
            s, f, g, e, h = _metrics(mol, ref)
            counts[m][0] += s
            counts[m][1] += f
            counts[m][2] += g
            counts[m][3] += e
            counts[m][4] += h
            row[f"{m}_ok"] = f and g and h
            if m == "geometry":
                row["geometry_mol"] = mol
                row["geometry_err"] = err
                if mol is None:
                    failures.append({**row, "kind": "None"})
                elif not (f and g and h):
                    failures.append({**row, "kind": "wrong"})
                    if f and not g:
                        tautomers.append(entry)
                    elif g and not e:
                        charges.append(entry)
        rows.append(row)

    print("\n== recovery (formula AND graph AND AddHs) ==")
    for m in METHODS:
        ok = sum(1 for r in rows if r[f"{m}_ok"])
        print(f"  {m:10s} {pct(ok, total)}")
    print("\n== metrics ==")
    for m in METHODS:
        s, f, g, e, h = counts[m]
        print(f"  {m:10s} sanitize={pct(s, total)} formula={pct(f, total)} "
              f"graph={pct(g, total)} exact={pct(e, total)} AddHs={pct(h, total)}")

    print(f"\ngeometry failures: {len(failures)}")
    for r in failures:
        print(f"  {r['pdb']} {r['het']} {r['name'][:30]:30s} "
              f"{r['formula']:12s} kind={r['kind']} err={r['geometry_err']}")
    print(f"tautomer-class (formula ok, graph differs): {len(tautomers)}")
    for e in tautomers:
        print(f"  {e['pdb']} {e['het']} {e['name'][:40]}")
    print(f"charge-class (graph ok, exact differs): {len(charges)}")
    for e in charges:
        print(f"  {e['pdb']} {e['het']} {e['name'][:40]}")

    # ---- render results.md ----
    stem = os.path.splitext(ds_name)[0]
    out_name = "results.md" if stem == "dataset" else f"results_{stem}.md"
    lines = []
    w = lines.append
    w(f"# Crystal-set benchmark: geometry tier on RCSB ligands ({res_band})")
    w("")
    w("Auto-generated by `benchmarks/crystal100/benchmark.py` - do not "
      "edit by hand.  Re-run with `uv run python "
      f"benchmarks/crystal100/benchmark.py --dataset {ds_name}`.  "
      f"Environment: {_env_info()}.")
    w("")
    w(f"Dataset: {total} PDB entries (protein-only, resolution band "
      f"{res_band}, one non-polymer entity each, all distinct HET "
      "codes), ligand residues extracted with complete-atom check "
      "against the CCD SMILES.  Ground truth: the RCSB Chemical "
      "Component Dictionary canonical SMILES (formula + bond graph).  "
      "`geometry` is the in-house perception called directly - the CCD "
      "template is never consulted and no fallbacks apply; a failure "
      "here is a genuine geometry-tier failure.  Metrics: `formula` "
      "(neutralised), `graph` (canonical SMILES equality, "
      "neutralised), `exact` (full canonical SMILES, charges and "
      "tautomer included), `AddHs` (no over-valent atoms), `recovery` "
      "= formula AND graph AND AddHs.")
    w("")
    w("## Summary")
    w("")
    w("| method | formula | graph | exact | AddHs | recovery |")
    w("|---|---|---|---|---|---|")
    for m in METHODS:
        s, f, g, e, h = counts[m]
        rec = sum(1 for r in rows if r[f"{m}_ok"])
        w(f"| {m} | {pct(f, total)} | {pct(g, total)} | {pct(e, total)} "
          f"| {pct(h, total)} | {pct(rec, total)} |")
    w("")
    w("## Tautomer detection (geometry)")
    w("")
    w("Ligands whose formula matches but whose bond graph does not: the "
      "geometry tier's Hückel judge placed the movable H(s) on a "
      "different atom than the CCD canonical SMILES records.  The "
      "coordinates decide - these are candidates for a legitimately "
      "perceived alternate tautomer, not graph corruption.")
    w("")
    if tautomers:
        w("| pdb | het | name | formula |")
        w("|---|---|---|---|")
        for e in tautomers:
            w(f"| {e['pdb']} | {e['het']} | {e['name']} | {e['formula']} |")
    else:
        w("(none)")
    w("")
    w("## Geometry failures (collected, not corrected)")
    w("")
    w("Every ligand the geometry tier failed to recover exactly.  This "
      "table is the input to the fix-hypothesis pass.")
    w("")
    if failures:
        w("| pdb | het | name | want formula | kind | err |")
        w("|---|---|---|---|---|---|")
        for r in failures:
            w(f"| {r['pdb']} | {r['het']} | {r['name']} | {r['formula']} "
              f"| {r['kind']} | {r['geometry_err'] or ''} |")
    else:
        w("(none)")
    w("")
    w("## Per-ligand detail")
    w("")
    w("| pdb | het | name | formula | res (A) | geom | ob | dist |")
    w("|---|---|---|---|---|---|---|---|")
    for r in rows:
        res = next((d.get("resolution") for d in DATASET
                    if d["pdb"] == r["pdb"]), None)
        res_s = f"{res:.2f}" if res else ""
        w(f"| {r['pdb']} | {r['het']} | {r['name']} | {r['formula']} "
          f"| {res_s} | {'OK' if r['geometry_ok'] else 'FAIL'} "
          f"| {'OK' if r['openbabel_ok'] else 'FAIL'} "
          f"| {'OK' if r['distance_ok'] else 'FAIL'} |")
    w("")
    # ---- machine-readable sidecar (same numbers as the markdown) ----
    # Consumed by summarize_heldout.py; keeps the aggregate step free of
    # markdown parsing.
    stats = {
        "dataset": stem,
        "dataset_file": ds_name,
        "total": total,
        "res_band": res_band,
        "methods": {
            m: {
                "sanitize": counts[m][0],
                "formula": counts[m][1],
                "graph": counts[m][2],
                "exact": counts[m][3],
                "addh": counts[m][4],
                "recovery": sum(1 for r in rows if r[f"{m}_ok"]),
            }
            for m in METHODS
        },
    }
    out_json = os.path.join(
        HERE, "results.json" if stem == "dataset" else f"results_{stem}.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=1, separators=(",", ": "))
    print(f"stats written: {out_json}")

    out = os.path.join(HERE, out_name)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"results written: {out}")


if __name__ == "__main__":
    main()
