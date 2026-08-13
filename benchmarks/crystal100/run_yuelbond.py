"""YuelBond head-to-head on the crystal benchmark (held-out gen2).

The published YuelBond GNN (Wang & Dokholyan, JCIM 2026; code:
bitbucket.org/dokhlab/yuel_bond) is the first ML bond-perception model to
release pretrained weights (`geom_3d.ckpt`, Zenodo record 15353365).  This
script runs it on the held-out gen2 ligands - the first evaluation of any ML
bond-perception model on experimental PDB coordinates against CCD ground
truth.

The method receives exactly the same input as every other baseline - the
heavy-atom ``(element, xyz)`` list from the dataset entry - and its output
mol is evaluated by the same ``_metrics`` as ``benchmark.py`` (formula /
graph / exact / AddHs against the CCD reference SMILES), so the numbers land
in the same tables.  Known domain mismatches are recorded, not hidden:

* the released model's atom vocabulary is C/O/N/F/S/Cl/Br/I/P; its own code
  (`atom_one_hot`) silently maps any other element to Cl.  Entries carrying
  such atoms are flagged per-entry and counted in the sidecar;
* the 3.0 A edge cutoff, the distance-only edge features and the 10-way
  bond-type argmax (SINGLE..ZERO) are the released model's, unchanged;
* the model never predicts hydrogens or charges - the neutralised
  comparison (same as all other methods) applies.

Run with the yuel_bond environment (torch + pytorch-lightning + rdkit, see
the checkpoint's own README)::

    python benchmarks/crystal100/run_yuelbond.py

writes `results_yuelbond_{stem}_k{n}.json` sidecars (same schema as
benchmark.py, method key `yuelbond`) plus per-bucket markdown.  The env and
the checkpoint sha256 are recorded in every sidecar.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from benchmark import _env_info, _metrics, load_dataset

HERE = os.path.dirname(os.path.abspath(__file__))
TIERS = ("main", "lowres")
BUCKETS = 5


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--yuelbond-dir", default="/tmp/yuel_bond",
                    help="checkout of bitbucket.org/dokhlab/yuel_bond "
                         "(default: %(default)s)")
    ap.add_argument("--checkpoint", default=None,
                    help="geom_3d.ckpt (default: <yuelbond-dir>/models/geom_3d.ckpt)")
    ap.add_argument("--prefix", default="dataset_heldout_gen2",
                    help="dataset basename (default: %(default)s)")
    ap.add_argument("--tiers", default="main,lowres")
    args = ap.parse_args()

    yb_dir = os.path.abspath(args.yuelbond_dir)
    ckpt = args.checkpoint or os.path.join(yb_dir, "models", "geom_3d.ckpt")
    if not os.path.exists(ckpt):
        sys.exit(f"checkpoint not found: {ckpt}")
    if not os.path.exists(os.path.join(yb_dir, "yuel_bond.py")):
        sys.exit(f"not a yuel_bond checkout: {yb_dir}")
    sys.path.insert(0, yb_dir)

    # Heavy imports after path setup; RDKit logs are noise at this volume.
    import torch
    from rdkit import Chem, RDLogger
    from yuel_bond import create_molecule_from_predictions

    from src import const
    from src.datasets import (
        BondDataset,
        atom_one_hot,
        collate,
        get_dataloader,
    )
    from src.lightning import YuelBond

    RDLogger.DisableLog("rdApp.*")

    model = YuelBond.load_from_checkpoint(ckpt, map_location="cpu").eval()
    print(f"model loaded: {ckpt} (sha256 {sha256(ckpt)}), "
          f"torch {torch.__version__}")

    def run_yuelbond(entry):
        """entry -> (mol, err) - identical input contract to benchmark.py."""
        elements = [el.capitalize() for el, _xyz in entry["atoms"]]
        coords = [xyz for _el, xyz in entry["atoms"]]
        one_hot = np.array([atom_one_hot(el) for el in elements])
        positions = np.array(coords, dtype=np.float64)
        raw = [{
            "name": f"{entry['pdb']}_{entry['het']}",
            "positions": positions,
            "atoms": one_hot,
            "bonds": np.array([]),
        }]
        try:
            dataset = BondDataset(raw_data=raw, device="cpu",
                                  has_bonds=False, progress_bar=False)
            dl = get_dataloader(dataset, batch_size=1, collate_fn=collate)
            data = next(iter(dl))
            with torch.no_grad():
                edge_pred = model.forward(data)
            mol = create_molecule_from_predictions(
                positions=data["positions"], one_hot=data["one_hot"],
                edge_index=data["edge_index"], edge_pred=edge_pred,
                node_mask=data["node_mask"], edge_mask=data["edge_mask"],
                name=data["name"][0])
            return mol, None
        except Exception as exc:  # noqa: BLE001 - a failed perception is a metric
            return None, f"{type(exc).__name__}: {exc}"

    env = _env_info()
    for tier in [t for t in TIERS if t in args.tiers.split(",")]:
        for k in range(1, BUCKETS + 1):
            ds_name = f"{args.prefix}_{tier}_k{k}.json"
            dataset = load_dataset(ds_name)
            total = len(dataset)
            res = [d.get("resolution") for d in dataset if d.get("resolution")]
            res_band = (f"{min(res):.1f}-{max(res):.1f} A" if res
                        else "unknown")
            counts = [0, 0, 0, 0, 0]  # s, f, g, e, h
            failures = []
            out_of_vocab = []
            rows = []
            for i, entry in enumerate(dataset):
                elements = [el.capitalize() for el, _xyz in entry["atoms"]]
                if any(el not in const.ALLOWED_ATOM_TYPES for el in elements):
                    out_of_vocab.append(entry)
                mol, err = run_yuelbond(entry)
                ref = Chem.MolFromSmiles(entry["smiles"])
                s, f, g, e, h = _metrics(mol, ref)
                counts[0] += s
                counts[1] += f
                counts[2] += g
                counts[3] += e
                counts[4] += h
                ok = f and g and h
                rows.append((entry, ok, mol, err))
                if not ok:
                    failures.append((entry, "None" if mol is None else "wrong", err))
            rec = sum(1 for _e, ok, _m, _err in rows if ok)

            # ---- render markdown ----
            stem = os.path.splitext(ds_name)[0]
            lines = []
            w = lines.append
            w(f"# YuelBond head-to-head: {tier} tier bucket k{k}")
            w("")
            w("Auto-generated by `benchmarks/crystal100/run_yuelbond.py` - "
              "do not edit by hand.  Environment: " + env + ", torch "
              f"{torch.__version__}, checkpoint sha256 {sha256(ckpt)}.")
            w("")
            w(f"Dataset: {total} PDB ligands ({stem}).  Method: YuelBond "
              "GNN (`geom_3d.ckpt`, Zenodo 15353365) on heavy-atom "
              "(element, xyz) input - identical input contract to the other "
              "baselines.  Metrics identical to benchmark.py: `recovery` = "
              "formula AND graph AND AddHs against the CCD reference.")
            w("")
            w("| method | formula | graph | exact | AddHs | recovery |")
            w("|---|---|---|---|---|---|")
            w(f"| yuelbond | {counts[1]}/{total} | {counts[2]}/{total} | "
              f"{counts[3]}/{total} | {counts[4]}/{total} | {rec}/{total} |")
            w("")
            w(f"Out-of-vocab entries (element outside the model's 9-atom "
              f"vocabulary, mapped to Cl by its own code): {len(out_of_vocab)}")
            for e in out_of_vocab:
                w(f"* {e['pdb']} {e['het']} "
                  f"{[el for el, _x in e['atoms'] if el not in const.ALLOWED_ATOM_TYPES]}")
            w("")
            w("## Failures")
            w("")
            for e, kind, err in failures:
                w(f"* {e['pdb']} {e['het']} {e['formula']:12s} kind={kind}"
                  + (f" err={err}" if err else ""))
            w("")
            out_md = os.path.join(HERE, f"results_yuelbond_{stem}.md")
            with open(out_md, "w", encoding="utf-8") as fh:
                fh.write("\n".join(lines) + "\n")

            # ---- sidecar (same schema as benchmark.py) ----
            stats = {
                "dataset": stem,
                "dataset_file": ds_name,
                "total": total,
                "res_band": res_band,
                "methods": {
                    "yuelbond": {
                        "sanitize": counts[0],
                        "formula": counts[1],
                        "graph": counts[2],
                        "exact": counts[3],
                        "addh": counts[4],
                        "recovery": rec,
                    }
                },
                "out_of_vocab": len(out_of_vocab),
                "checkpoint_sha256": sha256(ckpt),
            }
            out_json = os.path.join(HERE, f"results_yuelbond_{stem}.json")
            with open(out_json, "w", encoding="utf-8") as fh:
                json.dump(stats, fh, indent=1, separators=(",", ": "))
            print(f"{stem}: recovery {rec}/{total}, out-of-vocab "
                  f"{len(out_of_vocab)} -> {out_json}")


if __name__ == "__main__":
    main()
