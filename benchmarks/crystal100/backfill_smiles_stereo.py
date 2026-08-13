"""Backfill ``smiles_stereo`` (CCD stereo ground truth) into dataset files.

The original datasets predate the stereo benchmark and store only the
stereo-agnostic ``SMILES``.  This script adds the CCD ``SMILES_stereo``
field to every entry of every dataset file in this directory, without
re-fetching any PDB (the coordinates are unchanged).  Idempotent: entries
that already carry the field are left alone, so re-running after a partial
failure continues where it stopped.

::

    uv run python benchmarks/crystal100/backfill_smiles_stereo.py

``fetch_dataset.py`` stores the field natively for future fetches; this
script exists only for the committed datasets it predates.
"""

from __future__ import annotations

import json
import os
import time

from fetch_dataset import CCD_URL, get

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    for name in sorted(os.listdir(HERE)):
        if not (name.startswith("dataset") and name.endswith(".json")):
            continue
        path = os.path.join(HERE, name)
        with open(path, encoding="utf-8") as fh:
            ds = json.load(fh)
        missing = [d for d in ds if "smiles_stereo" not in d]
        if not missing:
            print(f"{name}: all {len(ds)} entries already carry smiles_stereo")
            continue
        for i, d in enumerate(missing, start=1):
            ccd = json.loads(get(CCD_URL.format(d["het"])))
            d["smiles_stereo"] = (
                ccd.get("rcsb_chem_comp_descriptor") or {}
            ).get("SMILES_stereo")
            if i % 20 == 0:
                print(f"{name}: {i}/{len(missing)}")
                with open(path, "w", encoding="utf-8") as fh:
                    json.dump(ds, fh, indent=1, separators=(",", ": "))
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(ds, fh, indent=1, separators=(",", ": "))
        stereo = sum(1 for d in ds if d.get("smiles_stereo"))
        print(f"{name}: {len(missing)} backfilled, {stereo}/{len(ds)} "
              "entries carry stereo")
        time.sleep(1)


if __name__ == "__main__":
    main()
