"""Sample held-out ligands for the publication benchmark.

The committed crystal100 datasets (``dataset.json``, ``dataset_res250-300.json``)
were used to tune the geometry tier: every fix in the issue pass was
informed by them.  Numbers measured on them are therefore optimistic -
the correct way to publish is a **held-out** set the tuning never saw.

Sampling design (this script):

* Universe: the same RCSB search the tuning fetch uses (protein-only,
  exactly one non-polymer entity, resolution band per tier), minus every
  entry whose **PDB id or HET code appears in either tuning dataset** -
  the exclusion bar the user chose is disjoint by entry AND chemotype.
* ``--buckets`` x ``--size`` independent samples per tier (default 5 x 60):
  the first N accepted candidates (N = buckets x size) are seeded-shuffled
  and chunked into buckets, so each bucket is a simple random sample of
  the tier universe minus exclusions.  Between-bucket spread is genuine
  sampling variation - the mean +/- std the paper reports is computed
  from these buckets by ``summarize_heldout.py``.
* Acceptance gates are identical to the tuning fetch (blacklist, CCD
  SMILES parses, 5-60 heavy atoms, complete residue, distinct HET codes
  within the pool).
* A per-tier manifest records the sampling frame (search total, excluded,
  accepted), the seeds and the bucket sizes - the numbers a paper must
  state to describe its sampling.

The held-out sets are run with the exact committed perception code - no
tuning on them (see the benchmark results for the code-freeze record).

Example::

    uv run python benchmarks/crystal100/sample_heldout.py --tier main
    uv run python benchmarks/crystal100/sample_heldout.py --tier lowres
"""

from __future__ import annotations

import argparse
import json
import os
import random
import threading
import time
from concurrent.futures import ThreadPoolExecutor

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

from fetch_dataset import (  # noqa: E402
    BLACKLIST,
    CCD_URL,
    PDB_URL,
    extract_residues,
    get,
    het_codes,
    nonpolymer_het,
    resolution_of,
    search_entries,
)

HERE = os.path.dirname(os.path.abspath(__file__))

TIERS = {
    # tier name -> (min_res, max_res) - the same bands as the tuning sets
    "main": (0.0, 2.0),
    "lowres": (2.5, 3.0),
}

TUNING = ("dataset.json", "dataset_res250-300.json")


def dataset_exclusions(files: list[str]) -> tuple[set[str], set[str]]:
    """Union of PDB ids and HET codes across the given dataset files."""
    pdbs, hets = set(), set()
    for name in files:
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            raise SystemExit(f"exclusion dataset not found: {path}")
        with open(path, encoding="utf-8") as fh:
            entries = json.load(fh)
        for entry in entries:
            pdbs.add(entry["pdb"])
            hets.add(entry["het"])
    return pdbs, hets


def tuning_exclusions() -> tuple[set[str], set[str]]:
    """Union of PDB ids and HET codes across both tuning datasets."""
    return dataset_exclusions(TUNING)


def evaluate_candidate(pdb_id: str, seen_hets: set[str]) -> dict | str:
    """One candidate -> entry dict, or a skip-reason string.

    Identical acceptance gates to ``fetch_dataset.process_candidate``,
    plus the HET exclusion: a candidate whose HET code appears in the
    tuning sets (or was already accepted into this pool) is skipped.
    """
    try:
        pdb_text = get(PDB_URL.format(pdb_id)).decode("utf-8", "replace")
        codes = het_codes(pdb_text)
        if len(codes) == 1:
            het = codes.pop()
        else:
            het = nonpolymer_het(pdb_id)
            if het is None or het not in codes:
                return "het_ambiguous"
        if het in BLACKLIST:
            return f"blacklist:{het}"
        if het in seen_hets:
            return "het_seen_or_excluded"
        ccd = json.loads(get(CCD_URL.format(het)))
        smiles = (ccd.get("rcsb_chem_comp_descriptor") or {}).get("SMILES")
        if not smiles:
            return "ccd_no_smiles"
        ref = Chem.MolFromSmiles(smiles)  # type: ignore[attr-defined]
        if ref is None:
            return "smiles_unparseable"
        n_heavy = ref.GetNumHeavyAtoms()
        if not (5 <= n_heavy <= 60):
            return "size_out_of_range"
        instances = extract_residues(pdb_text, het)
        if not instances:
            return "residue_not_found"
        atoms = next((r for r in instances if len(r) == n_heavy), None)
        if atoms is None:
            atoms = max(instances, key=len)
        if len(atoms) != n_heavy:
            return "atom_count_mismatch"
        return {
            "pdb": pdb_id,
            "het": het,
            "name": (ccd.get("name") or "").strip(),
            "smiles": Chem.MolToSmiles(ref),  # type: ignore[attr-defined]
            "smiles_stereo": (ccd.get("rcsb_chem_comp_descriptor") or {}).get(
                "SMILES_stereo"),
            "formula": rdMolDescriptors.CalcMolFormula(ref),  # type: ignore[attr-defined]
            "heavy_atoms": n_heavy,
            "resolution": resolution_of(pdb_text),
            "atoms": atoms,
        }
    except Exception as exc:  # noqa: BLE001 - network hiccups
        return f"error:{type(exc).__name__}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tier", choices=sorted(TIERS), required=True)
    parser.add_argument("--buckets", type=int, default=5)
    parser.add_argument("--size", type=int, default=60)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--label", default="",
                        help="output infix: dataset_heldout{_label}_{tier}_k*.json "
                             "(empty label keeps the default generation name)")
    parser.add_argument("--exclude-datasets", nargs="+", default=[],
                        help="extra dataset JSONs whose (pdb, het) pairs join "
                             "the exclusion - e.g. a prior held-out generation "
                             "must not re-appear in a fresh draw")
    parser.add_argument("--preflight", action="store_true",
                        help="report universe/exclusion sizes and stop before "
                             "fetching any candidate")
    args = parser.parse_args()

    min_res, max_res = TIERS[args.tier]
    need = args.buckets * args.size
    excluded_pdbs, excluded_hets = tuning_exclusions()
    if args.exclude_datasets:
        extra_pdbs, extra_hets = dataset_exclusions(args.exclude_datasets)
        excluded_pdbs |= extra_pdbs
        excluded_hets |= extra_hets
    entries = search_entries(min_res, max_res)
    candidates = [p for p in entries if p not in excluded_pdbs]
    print(f"universe: {len(entries)} candidates, {len(entries) - len(candidates)} "
          f"excluded by PDB id (tuning + {len(args.exclude_datasets)} extra "
          f"datasets), {len(candidates)} eligible")
    if len(candidates) < need:
        raise SystemExit(f"universe too small: {len(candidates)} eligible < "
                         f"{need} needed - re-decide bucket size with the user")
    if args.preflight:
        print("preflight OK - universe sufficient, no candidates fetched")
        return

    random.Random(args.seed).shuffle(candidates)
    accepted: list[dict] = []
    seen_hets = set(excluded_hets)  # HET codes already owned by the tuning set
    skipped: dict[str, int] = {}
    lock = threading.Lock()
    t0 = time.time()

    def worker(i: int) -> None:
        idx = i
        while len(accepted) < need and idx < len(candidates):
            result = evaluate_candidate(candidates[idx], seen_hets)
            idx += args.workers
            time.sleep(0.25)
            if isinstance(result, str):
                with lock:
                    skipped[result] = skipped.get(result, 0) + 1
                continue
            with lock:
                accepted.append(result)
                seen_hets.add(result["het"])
                print(f"[{len(accepted)}/{need}] {result['pdb']} "
                      f"{result['het']} {result['formula']} "
                      f"({result['resolution']} A) {result['name'][:40]}")

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(worker, range(args.workers)))

    if len(accepted) > need:
        del accepted[need:]  # worker-race overshoot, keep acceptance order
    if len(accepted) < need:
        raise SystemExit(f"only {len(accepted)} accepted of {need} needed - "
                         "re-decide bucket size with the user")

    # bucket assignment: one seeded shuffle, chunk into --buckets samples
    label = f"_{args.label}" if args.label else ""
    random.Random(args.seed + 1).shuffle(accepted)
    for k in range(args.buckets):
        bucket = accepted[k * args.size:(k + 1) * args.size]
        out = os.path.join(
            HERE, f"dataset_heldout{label}_{args.tier}_k{k + 1}.json")
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(bucket, fh, indent=1, separators=(",", ": "))
        print(f"bucket {k + 1}: {len(bucket)} entries -> {out}")

    manifest = {
        "tier": args.tier,
        "resolution_band": [min_res, max_res],
        "search_total": len(entries),
        "excluded_by_tuning_pdb": len(entries) - len(candidates),
        "excluded_datasets": args.exclude_datasets,
        "accepted": len(accepted),
        "buckets": args.buckets,
        "bucket_size": args.size,
        "seed": args.seed,
        "skipped": dict(sorted(skipped.items(), key=lambda kv: -kv[1])),
        "elapsed_s": round(time.time() - t0, 1),
    }
    out = os.path.join(HERE, f"dataset_heldout{label}_{args.tier}_manifest.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1, separators=(",", ": "))
    print(f"manifest: {out}")
    for reason, count in sorted(skipped.items(), key=lambda kv: -kv[1]):
        print(f"  skipped {count:4d}: {reason}")


if __name__ == "__main__":
    main()
