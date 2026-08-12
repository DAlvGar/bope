"""Stereo benchmark: R/S + E/Z perception vs the CCD stereo ground truth.

The bond-order benchmark (``benchmark.py``) proves the geometry tier
recovers the CCD bond graph; this one proves it recovers the
stereochemistry the graph declares.  Ground truth per ligand: the CCD
``SMILES_stereo`` (backfilled by ``backfill_smiles_stereo.py``; only the
entries that carry one are measured - coverage is reported).  Re-fetching
PDBs is not needed: the coordinates and bond-order results are identical
to the bond-order benchmark.

Three methods, identical input (each runs its own bond-order perception
first - stereo cannot be judged on a wrong graph):

1. **geometry** - the in-house perception, then
   :func:`bope.stereo.perceive_stereochemistry` (RDKit
   ``AssignStereochemistryFrom3D`` on the perceived mol).
2. **openbabel** - the OpenBabel round trip, which assigns its own
   stereo from the same coordinates as a side effect of its SDF output.
3. **distance** - covalent-radius connectivity + sanitization, then
   the same RDKit stereo assignment.

Only the entries whose perceived bond graph matches the CCD exactly
(canonical SMILES equality, neutralised) are stereo-comparable: the
comparison maps atoms via the shared canonical skeleton.  The mapping is
verified to be a bond-for-bond isomorphism; the rare symmetric cases
where the canonical tie-break differs are reported as mapping-ambiguous
rather than silently mis-compared.

Metrics:

* coverage - entries whose CCD declares stereo (``smiles_stereo`` set)
* per-center R/S - the labels of the centers the CCD declares, compared
  atom-by-atom; unassigned-by-us centers count as errors; centers our
  geometry labels that the CCD leaves unmarked are "extra" (informational
  - the geometry supports more stereo than the deposit declares)
* E/Z - the double bonds the CCD declares, compared bond-by-bond
* full-molecule - isomeric canonical SMILES equality (secondary
  headline; fails on extras even when every declared center matches)
* **phosphate-P class** - R/S mismatches at tetrahedral P with three or
  more O neighbors.  CIP ranking at phosphate P flips with P-OH vs P-O-
  protonation, which the crystal does not record; these are reported as
  a separate class, not as stereo errors.

Re-run with::

    uv run python benchmarks/crystal100/stereo_benchmark.py
    uv run python benchmarks/crystal100/stereo_benchmark.py \
        --dataset dataset_res250-300.json
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src")
)

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

import bope as perception
from benchmark import _env_info, _neutral, load_dataset, pct, run_geometry, run_openbabel, run_distance  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

METHODS = {
    "geometry": run_geometry,
    "openbabel": run_openbabel,
    "distance": run_distance,
}
#: geometry / distance assign stereo after perception; OpenBabel's SDF
#: round trip carries its own stereo already (verified: E/Z + tetrahedral
#: tags survive into the RDKit mol).
_OPS = {
    "geometry": "rdkit-3d",
    "openbabel": "openbabel-sdf",
    "distance": "rdkit-3d",
}


def _map_atoms(
    ref: "Chem.Mol",
    per: "Chem.Mol",
    ref_centers: dict[int, str],
    per_centers: dict[int, str],
) -> dict[int, int] | None:
    """ref idx -> per idx isomorphism, maximising R/S label agreement.

    Atom maps cannot be used to align two differently-indexed molecules
    (the SMILES writer orders mapped atoms by map number), so the
    correspondence comes from graph matching.  All isomorphisms of the
    ref graph in the perceived mol are enumerated; the one with the most
    agreeing R/S labels is the stereo-consistent choice - this matters
    when two homotopic centers exist (meso compounds), where an arbitrary
    isomorphism would compare swapped centers and report phantom errors.

    Returns None when no isomorphism exists (cannot happen after the
    canonical-SMILES gate; reported as mapping-ambiguous regardless).
    """
    # The pattern must carry the ref atom indices: MolToSmarts does not
    # preserve atom indexing, so a bare pattern's atoms are ordered by the
    # Smarts string, not by the ref mol.  Atom map numbers survive the
    # Smarts round trip and recover the ref index of every pattern atom.
    m = Chem.Mol(ref)
    for a in m.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)  # type: ignore[attr-defined]
    pattern = Chem.MolFromSmarts(Chem.MolToSmarts(m))  # type: ignore[attr-defined]
    matches = per.GetSubstructMatches(pattern, uniquify=False)  # type: ignore[attr-defined]
    if not matches:
        return None

    def score(match: tuple) -> int:
        return sum(
            1
            for ri, label in ref_centers.items()
            if label != "?"
            and per_centers.get(match[pattern.GetAtomWithIdx(ri).GetAtomMapNum() - 1], "?")  # type: ignore[attr-defined]
            == label
        )

    best = max(matches, key=score)
    return {
        pattern.GetAtomWithIdx(i).GetAtomMapNum() - 1: best[i]  # type: ignore[attr-defined]
        for i in range(pattern.GetNumAtoms())
    }


def _phosphate_p(mol: "Chem.Mol", idx: int) -> bool:
    """True when atom *idx* is tetrahedral P with >= 3 O neighbors - the
    protonation-sensitive phosphate center class."""
    a = mol.GetAtomWithIdx(idx)
    if a.GetSymbol() != "P" or a.GetDegree() != 4:
        return False
    return sum(1 for n in a.GetNeighbors() if n.GetSymbol() == "O") >= 3


def _compare_stereo(
    ref: "Chem.Mol", per: "Chem.Mol"
) -> tuple[dict, dict, dict]:
    """Per-center R/S, E/Z and full-string comparison for two heavy-atom
    mols with equal plain canonical SMILES.

    Returns (centers, ez, full) where *centers* counts correct / wrong /
    unassigned / phosphate-P / extra labels, *ez* counts correct / wrong
    / extra E/Z labels, and *full* is the isomeric canonical SMILES
    equality verdict.
    """
    centers = {"correct": 0, "wrong": 0, "unassigned": 0, "phosphate_p": 0,
               "extra": 0}
    ez = {"correct": 0, "wrong": 0, "extra": 0}

    ref_centers = dict(Chem.FindMolChiralCenters(  # type: ignore[attr-defined]
        ref, useLegacyImplementation=False, includeUnassigned=True))
    per_centers = dict(Chem.FindMolChiralCenters(  # type: ignore[attr-defined]
        per, useLegacyImplementation=False, includeUnassigned=True))
    map_ref = _map_atoms(ref, per, ref_centers, per_centers)
    if map_ref is None:
        return {"mapping_ambiguous": 1, **centers}, ez, False
    for ri, label in ref_centers.items():
        if label == "?":
            continue  # CCD leaves this center unmarked - not compared
        pi = map_ref[ri]
        plabel = per_centers.get(pi, "?")
        if plabel == label:
            centers["correct"] += 1
        elif _phosphate_p(per, pi):
            centers["phosphate_p"] += 1
        elif plabel == "?":
            centers["unassigned"] += 1
        else:
            centers["wrong"] += 1
    # extra: centers we label that the CCD does not declare (including
    # potential centers the CCD leaves unmarked - the "reference declares
    # less stereo than the geometry supports" class, e.g. 1J5)
    declared = set(
        map_ref[ri] for ri, lab in ref_centers.items() if lab != "?"
    )
    for pi, plabel in per_centers.items():
        if plabel != "?" and pi not in declared:
            centers["extra"] += 1

    map_per = {pi: ri for ri, pi in map_ref.items()}  # per idx -> ref idx
    for b in ref.GetBonds():
        if b.GetStereo() == Chem.BondStereo.STEREONONE:  # type: ignore[attr-defined]
            continue
        pb = per.GetBondBetweenAtoms(map_ref[b.GetBeginAtomIdx()],
                                     map_ref[b.GetEndAtomIdx()])
        if pb is None or pb.GetStereo() != b.GetStereo():
            ez["wrong"] += 1
        else:
            ez["correct"] += 1
    for b in per.GetBonds():
        if b.GetStereo() != Chem.BondStereo.STEREONONE:  # type: ignore[attr-defined]
            pb = ref.GetBondBetweenAtoms(map_per[b.GetBeginAtomIdx()],
                                         map_per[b.GetEndAtomIdx()])
            if pb is None or pb.GetStereo() == Chem.BondStereo.STEREONONE:  # type: ignore[attr-defined]
                ez["extra"] += 1

    full = (
        Chem.MolToSmiles(per, isomericSmiles=True)  # type: ignore[attr-defined]
        == Chem.MolToSmiles(ref, isomericSmiles=True)  # type: ignore[attr-defined]
    )
    return centers, ez, full


def _method_stereo_mol(entry: dict, method: str):
    """Perceived bond-order mol with stereo labels for *method*."""
    mol, _err = METHODS[method](entry)
    if mol is None:
        return None, "no-mol"
    if _OPS[method] == "rdkit-3d":
        labeled = perception.perceive_stereochemistry(mol)
        if labeled is None:
            return None, "no-stereo"
        return labeled, None
    return mol, None  # openbabel: its own stereo rides along


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", default="dataset.json",
                        help="dataset file next to this script "
                             "(default: dataset.json)")
    ds_name = parser.parse_args().dataset
    DATASET = load_dataset(ds_name)

    # stereo-declaring entries: smiles_stereo that actually contains
    # stereo symbols (@ for tetrahedral, /\\ for E/Z) - a set field without
    # any symbol declares nothing
    stereo_entries = [
        d for d in DATASET
        if d.get("smiles_stereo") and re.search(r"[@/\\]", d["smiles_stereo"])
    ]
    total = len(DATASET)
    print(f"dataset: {total} ligands, {len(stereo_entries)} carry CCD "
          f"stereo ({100.0 * len(stereo_entries) / total:.0f}% coverage)")

    rows = []  # per-entry verdicts
    per_method = {m: {
        "no_mol": 0, "graph_mismatch": 0, "mapping_ambiguous": 0,
        "centers": {"correct": 0, "wrong": 0, "unassigned": 0,
                    "phosphate_p": 0, "extra": 0},
        "ez": {"correct": 0, "wrong": 0, "extra": 0},
        "full": 0, "comparable": 0,
    } for m in METHODS}

    for entry in stereo_entries:
        ref = Chem.MolFromSmiles(entry["smiles_stereo"])
        ref_plain = Chem.MolToSmiles(_neutral(ref), isomericSmiles=False)
        row = {"pdb": entry["pdb"], "het": entry["het"],
               "name": entry["name"], "formula": entry["formula"]}
        for m in METHODS:
            stats = per_method[m]
            mol, err = _method_stereo_mol(entry, m)
            if mol is None:
                stats["no_mol"] += 1
                row[m] = f"FAIL:{err}"
                continue
            plain = Chem.MolToSmiles(_neutral(mol), isomericSmiles=False)
            if plain != ref_plain:
                stats["graph_mismatch"] += 1
                row[m] = "FAIL:graph"
                continue
            centers, ez, full = _compare_stereo(_neutral(ref), _neutral(mol))
            if centers.get("mapping_ambiguous"):
                stats["mapping_ambiguous"] += 1
                row[m] = "FAIL:mapping"
                continue
            stats["comparable"] += 1
            for k in stats["centers"]:
                stats["centers"][k] += centers[k]
            for k in stats["ez"]:
                stats["ez"][k] += ez[k]
            stats["full"] += 1 if full else 0
            if full:
                row[m] = "OK"
            elif centers["correct"] and not (centers["wrong"] or centers["unassigned"]):
                row[m] = "centers-ok"
            else:
                bits = []
                if centers["wrong"]:
                    bits.append(f"r/s:{centers['wrong']}")
                if centers["unassigned"]:
                    bits.append(f"unassigned:{centers['unassigned']}")
                if centers["phosphate_p"]:
                    bits.append(f"po4:{centers['phosphate_p']}")
                if ez["wrong"]:
                    bits.append(f"e/z:{ez['wrong']}")
                row[m] = "FAIL:" + ",".join(bits) if bits else "FAIL"
        rows.append(row)

    # ---- console summary ----
    print("\n== stereo-comparable (bond graph matches the CCD exactly) ==")
    for m in METHODS:
        s = per_method[m]
        print(f"  {m:10s} comparable={pct(s['comparable'], len(stereo_entries))} "
              f"no-mol={s['no_mol']} graph-mismatch={s['graph_mismatch']} "
              f"mapping-ambiguous={s['mapping_ambiguous']}")
    print("\n== per-center R/S (centers the CCD declares) ==")
    for m in METHODS:
        c = per_method[m]["centers"]
        judged = c["correct"] + c["wrong"] + c["unassigned"]
        prec = f"{100.0 * c['correct'] / judged:.0f}%" if judged else "-"
        print(f"  {m:10s} correct={c['correct']} wrong={c['wrong']} "
              f"unassigned={c['unassigned']} (precision {prec})  "
              f"phosphate-P flips={c['phosphate_p']} (separate class)  "
              f"extra centers ours={c['extra']}")
    print("\n== E/Z (double bonds the CCD declares) ==")
    for m in METHODS:
        e = per_method[m]["ez"]
        judged = e["correct"] + e["wrong"]
        prec = f"{100.0 * e['correct'] / judged:.0f}%" if judged else "-"
        print(f"  {m:10s} correct={e['correct']} wrong={e['wrong']} "
              f"(precision {prec})  extra ours={e['extra']}")
    print("\n== full-molecule stereo recovery ==")
    for m in METHODS:
        s = per_method[m]
        print(f"  {m:10s} {pct(s['full'], s['comparable'])} "
              f"(of {s['comparable']} stereo-comparable)")
    # entries with a full stereo mismatch on the geometry tier
    print("\ngeometry full-mismatches:")
    for r in rows:
        if r["geometry"] != "OK":
            print(f"  {r['pdb']} {r['het']} {r['name'][:34]:34s} {r['geometry']}")

    # ---- render results_md ----
    stem = os.path.splitext(ds_name)[0]
    out_name = "results_stereo.md" if stem == "dataset" else f"results_stereo_{stem}.md"
    lines = []
    w = lines.append
    w(f"# Stereo benchmark: R/S + E/Z perception vs CCD stereo ground truth")
    w("")
    w("Auto-generated by `benchmarks/crystal100/stereo_benchmark.py` - do "
      "not edit by hand.  Re-run with `uv run python "
      f"benchmarks/crystal100/stereo_benchmark.py --dataset {ds_name}`.  "
      f"Environment: {_env_info()}.")
    w("")
    w(f"Dataset: {total} entries, {len(stereo_entries)} carry CCD "
      f"`SMILES_stereo` ({100.0 * len(stereo_entries) / total:.0f}% "
      "coverage).  Methods: `geometry` (in-house perception + RDKit "
      "stereo from 3-D), `openbabel` (its own stereo from the SDF round "
      "trip), `distance` (covalent-radius connectivity + RDKit stereo).  "
      "Only entries whose perceived bond graph matches the CCD exactly "
      "are stereo-comparable; the atom mapping uses the shared canonical "
      "skeleton and is verified bond-for-bond.  `phosphate-P flips` are "
      "R/S mismatches at tetrahedral P with 3+ O neighbors: CIP ranking "
      "flips with P-OH vs P-O- protonation, which the crystal does not "
      "record - a separate class, not a stereo error.")
    w("")
    w("## Per-center R/S (centers the CCD declares)")
    w("")
    w("| method | correct | wrong | unassigned | precision | phosphate-P flips | extra (ours) |")
    w("|---|---|---|---|---|---|---|")
    for m in METHODS:
        c = per_method[m]["centers"]
        judged = c["correct"] + c["wrong"] + c["unassigned"]
        prec = f"{100.0 * c['correct'] / judged:.0f}%" if judged else "-"
        w(f"| {m} | {c['correct']} | {c['wrong']} | {c['unassigned']} "
          f"| {prec} | {c['phosphate_p']} | {c['extra']} |")
    w("")
    w("## E/Z (double bonds the CCD declares)")
    w("")
    w("| method | correct | wrong | precision | extra (ours) |")
    w("|---|---|---|---|---|")
    for m in METHODS:
        e = per_method[m]["ez"]
        judged = e["correct"] + e["wrong"]
        prec = f"{100.0 * e['correct'] / judged:.0f}%" if judged else "-"
        w(f"| {m} | {e['correct']} | {e['wrong']} | {prec} | {e['extra']} |")
    w("")
    w("## Full-molecule stereo recovery (isomeric canonical SMILES)")
    w("")
    w("| method | recovered | comparable |")
    w("|---|---|---|")
    for m in METHODS:
        s = per_method[m]
        w(f"| {m} | {pct(s['full'], s['comparable'])} | {s['comparable']} |")
    w("")
    w("## Per-ligand detail (stereo-declaring ligands only)")
    w("")
    w("| pdb | het | name | formula | res (A) | geom | ob | dist |")
    w("|---|---|---|---|---|---|---|---|")
    for r in rows:
        res = next((d.get("resolution") for d in DATASET
                    if d["pdb"] == r["pdb"]), None)
        res_s = f"{res:.2f}" if res else ""
        w(f"| {r['pdb']} | {r['het']} | {r['name']} | {r['formula']} "
          f"| {res_s} | {r['geometry']} | {r['openbabel']} | {r['distance']} |")
    w("")
    # ---- machine-readable sidecar (same numbers as the markdown) ----
    # Consumed by summarize_heldout.py; keeps the aggregate step free of
    # markdown parsing.
    stats = {
        "dataset": stem,
        "dataset_file": ds_name,
        "total": total,
        "stereo_entries": len(stereo_entries),
        "coverage_pct": round(100.0 * len(stereo_entries) / total, 1),
        "methods": {m: per_method[m] for m in METHODS},
    }
    out_json = os.path.join(
        HERE, "results_stereo.json" if stem == "dataset"
        else f"results_stereo_{stem}.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=1, separators=(",", ": "))
    print(f"stats written: {out_json}")

    out = os.path.join(HERE, out_name)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"results written: {out}")


if __name__ == "__main__":
    main()
