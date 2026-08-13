"""Aggregate the held-out benchmark runs into results_heldout.md.

The held-out protocol (see ``sample_heldout.py``): 5 independent buckets
of 60 ligands per resolution tier, sampled from the same RCSB universes
as the tuning sets minus every entry whose PDB id or HET code the tuning
sets own.  Each bucket is benchmarked by ``benchmark.py`` (bond orders)
and ``stereo_benchmark.py`` (stereo); each run writes a machine-readable
sidecar (``results_*.json``) next to its markdown - this script reads
the 10 bond-order sidecars + 10 stereo sidecars and aggregates:

* **pooled** counts (all 5 buckets summed),
* **mean +/- std** across the 5 buckets: each bucket is a simple random
  sample of the tier universe, so the between-bucket spread is genuine
  sampling variation - the sampling-error estimate the paper reports
  (statistics pass at write-up time can refine the CI from these).

The tuning-set numbers (sidecars of the committed ``dataset.json`` and
``dataset_res250-300.json`` runs) are read for comparison - the
tuning-vs-held-out gap is exactly the reviewer question this protocol
answers.

The YuelBond head-to-head (``run_yuelbond.py``, the released model and
weights from Zenodo record 15353365) is merged from its own sidecars
(``results_yuelbond_*.json``) into the bond-order tables when all 5
buckets of a tier are present.

Also writes the sampling frame from the manifests (search total,
exclusion counts, seeds) and verifies the held-out PDB ids and HET
codes have zero overlap with the tuning sets.

Re-run with::

    uv run python benchmarks/crystal100/summarize_heldout.py

Use ``--prefix`` and ``--out`` for labeled re-runs of the held-out sets
(e.g. after perception fixes): ``--prefix dataset_heldout_fix
--out results_heldout_fix.md``.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src")
)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from benchmark import _env_info

HERE = os.path.dirname(os.path.abspath(__file__))

TIERS = ("main", "lowres")
BUCKETS = 5
#: canonical method order - the same in benchmark.py / stereo_benchmark.py
METHODS_ORDER = ("geometry", "openbabel", "distance", "rdDetermineBonds")
TUNING_DATASETS = ("dataset.json", "dataset_res250-300.json")
TUNING_STEMS = {  # tier -> bond sidecar, stereo sidecar
    "main": ("results.json", "results_stereo.json"),
    "lowres": ("results_dataset_res250-300.json",
               "results_stereo_dataset_res250-300.json"),
}
TIER_LABEL = {"main": "1.0-2.0 A", "lowres": "2.5-3.0 A"}


def load_json(name: str):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return None
    return json.load(open(path, encoding="utf-8"))


def frac(ok: int, total: int) -> float | None:
    return 100.0 * ok / total if total else None


def mean_std(values: list) -> tuple[float | None, float | None]:
    """Mean and sample std of percentage fractions (None entries dropped)."""
    vals = [v for v in values if v is not None]
    if not vals:
        return None, None
    m = sum(vals) / len(vals)
    s = statistics.stdev(vals) if len(vals) >= 2 else None
    return m, s


def f1(v: float | None) -> str:
    return f"{v:.1f}" if v is not None else "-"


def fms(v: float | None, s: float | None) -> str:
    if v is None:
        return "-"
    return f"{v:.1f} +/- {s:.1f}" if s is not None else f"{v:.1f}"


def pct(ok: int, total: int) -> str:
    if total == 0:
        return "0/0 (-)"
    return f"{ok}/{total} ({frac(ok, total):.1f}%)"


def bond_fraction(sidecar: dict, method: str, metric: str) -> float | None:
    m = sidecar["methods"][method]
    return frac(m[metric], sidecar["total"])


def stereo_fractions(sidecar: dict, method: str) -> dict:
    """full-string, R/S precision, E/Z precision fractions for one run."""
    m = sidecar["methods"][method]
    c = m["centers"]
    judged = c["correct"] + c["wrong"] + c["unassigned"]
    ez_judged = m["ez"]["correct"] + m["ez"]["wrong"]
    return {
        "full": frac(m["full"], m["comparable"]),
        "rs": frac(c["correct"], judged),
        "ez": frac(m["ez"]["correct"], ez_judged),
        "comparable": m["comparable"],
    }


def verify_exclusion(prefix: str = "dataset_heldout") -> tuple[set, set, set, set]:
    """(held-out pdb ids, held-out het codes, overlap with tuning pdb ids,
    overlap with tuning het codes) - the overlaps must be empty for the
    protocol to hold."""
    tuning_pdbs, tuning_hets = set(), set()
    for name in TUNING_DATASETS:
        for e in load_json(name):
            tuning_pdbs.add(e["pdb"])
            tuning_hets.add(e["het"])
    ho_pdbs, ho_hets = set(), set()
    for tier in TIERS:
        for k in range(1, BUCKETS + 1):
            for e in load_json(f"{prefix}_{tier}_k{k}.json"):
                ho_pdbs.add(e["pdb"])
                ho_hets.add(e["het"])
    return ho_pdbs, ho_hets, ho_pdbs & tuning_pdbs, ho_hets & tuning_hets


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--prefix", default="dataset_heldout",
                    help="held-out dataset basename (default: %(default)s)")
    ap.add_argument("--out", default="results_heldout.md",
                    help="summary output file (default: %(default)s)")
    ap.add_argument("--prior", default=None,
                    help="earlier held-out generation prefix to verify "
                         "disjointness against (e.g. dataset_heldout when "
                         "summarizing gen2)")
    args = ap.parse_args()
    prefix = args.prefix
    env = _env_info()
    bond = {t: [load_json(f"results_{prefix}_{t}_k{k}.json")
                for k in range(1, BUCKETS + 1)] for t in TIERS}
    stereo = {t: [load_json(f"results_stereo_{prefix}_{t}_k{k}.json")
                  for k in range(1, BUCKETS + 1)] for t in TIERS}
    # YuelBond head-to-head sidecars (run_yuelbond.py) - merged into the
    # bond sidecars so the same tables render one row per method.
    yb = {t: [load_json(f"results_yuelbond_{prefix}_{t}_k{k}.json")
              for k in range(1, BUCKETS + 1)] for t in TIERS}
    have_yb = {}
    for t in TIERS:
        have_yb[t] = all(s is not None and bond[t][k] is not None
                         and s.get("total") == bond[t][k]["total"]
                         for k, s in enumerate(yb[t]))
        for k, s in enumerate(yb[t]):
            if have_yb[t] and "yuelbond" in s["methods"]:
                bond[t][k]["methods"]["yuelbond"] = s["methods"]["yuelbond"]
    manifests = {t: load_json(f"{prefix}_{t}_manifest.json")
                 for t in TIERS}
    for t in TIERS:
        if any(s is None for s in bond[t] + stereo[t]):
            raise SystemExit(f"missing held-out sidecars for tier {t} - "
                             "run benchmark.py and stereo_benchmark.py on "
                             f"all {prefix}_*_k*.json first")
    tuning_bond = {t: load_json(TUNING_STEMS[t][0]) for t in TIERS}
    tuning_stereo = {t: load_json(TUNING_STEMS[t][1]) for t in TIERS}
    if any(v is None for v in tuning_bond.values()) or \
       any(v is None for v in tuning_stereo.values()):
        print("note: tuning sidecars missing - re-run benchmark.py and "
              "stereo_benchmark.py on dataset.json / dataset_res250-300.json "
              "to include the tuning-vs-held-out comparison")

    ho_pdbs, ho_hets, overlap_pdbs, overlap_hets = verify_exclusion(prefix)
    prior_overlap_pdbs = prior_overlap_hets = None
    if args.prior:
        prior_pdbs, prior_hets = set(), set()
        for tier in TIERS:
            for k in range(1, BUCKETS + 1):
                for e in load_json(f"{args.prior}_{tier}_k{k}.json"):
                    prior_pdbs.add(e["pdb"])
                    prior_hets.add(e["het"])
        prior_overlap_pdbs = ho_pdbs & prior_pdbs
        prior_overlap_hets = ho_hets & prior_hets

    lines = []
    w = lines.append
    w("# Held-out benchmark: bond orders + stereo on 600 never-seen "
      "RCSB ligands")
    w("")
    w("Auto-generated by `benchmarks/crystal100/summarize_heldout.py` - "
      "do not edit by hand.  Re-run with `uv run python "
      "benchmarks/crystal100/summarize_heldout.py`.  "
      f"Environment: {env}.")
    w("")
    w("The tuning datasets (`dataset.json`, `dataset_res250-300.json`) "
      "drove every fix of the geometry tier, so numbers measured on them "
      "are optimistic.  The held-out sets sample the same RCSB universes "
      "(protein-only, one non-polymer entity, resolution band) minus "
      "every entry whose **PDB id or HET code** appears in either tuning "
      "dataset or any prior held-out generation (the manifest records the "
      "exact exclusion list), as 5 independent buckets of 60 per tier.  "
      "Each bucket is a simple random sample: the per-bucket spread is "
      "genuine sampling variation, reported as mean +/- std (n = 5).  "
      "All held-out runs use the exact committed perception code - no "
      "tuning on held-out results.")
    if any(have_yb.values()):
        w("YuelBond rows are the released model with its published weights "
          "(Zenodo record 15353365) run head-to-head on the same "
          "(element, xyz) input, merged from the `results_yuelbond_*.json` "
          "sidecars of `run_yuelbond.py`.")
        w("")
    w("")
    w("## Sampling frame")
    w("")
    w("| tier | resolution band | search total | excluded (PDB id) "
      "| accepted | buckets x size | seed |")
    w("|---|---|---|---|---|---|---|")
    for t in TIERS:
        mn = manifests[t]
        w(f"| {t} | {TIER_LABEL[t]} | {mn['search_total']} | "
          f"{mn['excluded_by_tuning_pdb']} | {mn['accepted']} | "
          f"{mn['buckets']} x {mn['bucket_size']} | {mn['seed']} |")
    w("")
    w("Skip reasons (top, per tier):")
    w("")
    for t in TIERS:
        skipped = sorted(manifests[t]["skipped"].items(), key=lambda kv: -kv[1])
        top = ", ".join(f"`{r}` x {n}" for r, n in skipped[:4])
        w(f"* {t}: {top}")
    w("")
    w("## Verification")
    w("")
    w(f"Held-out vs tuning overlap - PDB ids: **{len(overlap_pdbs)}**, "
      f"HET codes: **{len(overlap_hets)}** (protocol requires zero).")
    if prior_overlap_pdbs is not None:
        w(f"Held-out vs prior generation ({args.prior}) overlap - PDB ids: "
          f"**{len(prior_overlap_pdbs)}**, HET codes: "
          f"**{len(prior_overlap_hets)}** (protocol requires zero).")
    w("Per-bucket results carry the environment + commit footer of the "
      "exact code they ran.")
    w("")
    w("## Bond-order recovery (formula AND graph AND AddHs)")
    w("")
    w("| tier | method | tuning | held-out pooled | per-bucket mean +/- std (5 x 60) |")
    w("|---|---|---|---|---|")
    for t in TIERS:
        for m in METHODS_ORDER + (("yuelbond",) if have_yb[t] else ()):
            pooled_ok = sum(b["methods"][m]["recovery"] for b in bond[t])
            pooled_tot = sum(b["total"] for b in bond[t])
            per_bucket = [bond_fraction(b, m, "recovery") for b in bond[t]]
            mmean, mstd = mean_std(per_bucket)
            tb = tuning_bond[t]
            if tb and m in tb["methods"]:
                tun = (f"{tb['methods'][m]['recovery']}/{tb['total']} "
                       f"({frac(tb['methods'][m]['recovery'], tb['total']):.1f}%)")
            else:
                tun = "-"
            w(f"| {t} | {m} | {tun} | {pct(pooled_ok, pooled_tot)} "
              f"| {fms(mmean, mstd)} |")
    w("")
    w("## Bond-order metrics (held-out pooled)")
    w("")
    w("| tier | method | formula | graph | exact | AddHs |")
    w("|---|---|---|---|---|---|")
    for t in TIERS:
        for m in METHODS_ORDER + (("yuelbond",) if have_yb[t] else ()):
            tot = sum(b["total"] for b in bond[t])
            cells = [f"{sum(b['methods'][m][k] for b in bond[t])}/{tot}"
                     for k in ("formula", "graph", "exact", "addh")]
            w(f"| {t} | {m} | " + " | ".join(cells) + " |")
    w("")
    w("## Stereo: full-molecule recovery (isomeric canonical SMILES)")
    w("")
    w("| tier | method | tuning | held-out pooled | per-bucket mean +/- std |")
    w("|---|---|---|---|---|")
    for t in TIERS:
        for m in METHODS_ORDER:
            fs = [stereo_fractions(s, m) for s in stereo[t]]
            pooled_full = sum(s["methods"][m]["full"] for s in stereo[t])
            pooled_comp = sum(f["comparable"] for f in fs)
            mmean, mstd = mean_std([f["full"] for f in fs])
            ts = tuning_stereo[t]
            if ts and ts["methods"][m]["comparable"]:
                tun = pct(ts["methods"][m]["full"],
                          ts["methods"][m]["comparable"])
            else:
                tun = "-"
            w(f"| {t} | {m} | {tun} | {pct(int(pooled_full), pooled_comp)} "
              f"| {fms(mmean, mstd)} |")
    w("")
    for t in TIERS:
        entries = sum(s["stereo_entries"] for s in stereo[t])
        tot = sum(s["total"] for s in stereo[t])
        cov = [frac(s["stereo_entries"], s["total"]) for s in stereo[t]]
        cm, cs = mean_std(cov)
        w(f"Stereo coverage {t}: {entries}/{tot} ({frac(entries, tot):.1f}% "
          f"of entries, per-bucket {fms(cm, cs)})")
    w("")
    w("## Stereo: per-center R/S (centers the CCD declares)")
    w("")
    w("| tier | method | correct | wrong | unassigned | pooled precision | "
      "per-bucket precision mean +/- std |")
    w("|---|---|---|---|---|---|---|")
    for t in TIERS:
        for m in METHODS_ORDER:
            fs = [stereo_fractions(s, m) for s in stereo[t]]
            cs = [sum(s["methods"][m]["centers"]["correct"] for s in stereo[t]),
                  sum(s["methods"][m]["centers"]["wrong"] for s in stereo[t]),
                  sum(s["methods"][m]["centers"]["unassigned"] for s in stereo[t])]
            judged = sum(cs)
            mmean, mstd = mean_std([f["rs"] for f in fs])
            w(f"| {t} | {m} | {cs[0]} | {cs[1]} | {cs[2]} "
              f"| {pct(cs[0], judged)} | {fms(mmean, mstd)} |")
    w("")
    w("## Stereo: E/Z (double bonds the CCD declares)")
    w("")
    w("| tier | method | correct | wrong | pooled precision | "
      "per-bucket precision mean +/- std |")
    w("|---|---|---|---|---|---|")
    for t in TIERS:
        for m in METHODS_ORDER:
            fs = [stereo_fractions(s, m) for s in stereo[t]]
            ez = [sum(s["methods"][m]["ez"]["correct"] for s in stereo[t]),
                  sum(s["methods"][m]["ez"]["wrong"] for s in stereo[t])]
            judged = sum(ez)
            mmean, mstd = mean_std([f["ez"] for f in fs])
            w(f"| {t} | {m} | {ez[0]} | {ez[1]} | {pct(ez[0], judged)} "
              f"| {fms(mmean, mstd)} |")
    w("")
    w("## Per-bucket detail")
    w("")
    for t in TIERS:
        w(f"### {t} tier")
        w("")
        w("| bucket | bond rec (geom) | bond rec (ob) | bond rec (dist) "
          "| stereo full (geom) | R/S prec (geom) | stereo coverage |")
        w("|---|---|---|---|---|---|---|")
        for k, (b, s) in enumerate(zip(bond[t], stereo[t]), start=1):
            g = b["methods"]["geometry"]
            gf = stereo_fractions(s, "geometry")
            w(f"| k{k} | {g['recovery']}/{b['total']} "
              f"({frac(g['recovery'], b['total']):.1f}%) | "
              f"{b['methods']['openbabel']['recovery']}/{b['total']} "
              f"({frac(b['methods']['openbabel']['recovery'], b['total']):.1f}%) | "
              f"{b['methods']['distance']['recovery']}/{b['total']} "
              f"({frac(b['methods']['distance']['recovery'], b['total']):.1f}%) | "
              f"{s['methods']['geometry']['full']}/{s['methods']['geometry']['comparable']} | "
              f"{f1(gf['rs'])} | {s['stereo_entries']}/{s['total']} |")
        w("")

    out = os.path.join(HERE, args.out)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- console summary ----
    print(f"overlap with tuning: {len(overlap_pdbs)} PDB ids, "
          f"{len(overlap_hets)} HET codes (must be 0)")
    if prior_overlap_pdbs is not None:
        print(f"overlap with {args.prior}: {len(prior_overlap_pdbs)} PDB ids, "
              f"{len(prior_overlap_hets)} HET codes (must be 0)")
    print("\n== bond-order recovery ==")
    for t in TIERS:
        for m in METHODS_ORDER + (("yuelbond",) if have_yb[t] else ()):
            pooled_ok = sum(b["methods"][m]["recovery"] for b in bond[t])
            pooled_tot = sum(b["total"] for b in bond[t])
            mmean, mstd = mean_std([bond_fraction(b, m, "recovery")
                                    for b in bond[t]])
            print(f"  {t:7s} {m:10s} pooled={pct(pooled_ok, pooled_tot)} "
                  f"mean={fms(mmean, mstd)}")
    print("\n== stereo full-molecule recovery ==")
    for t in TIERS:
        for m in METHODS_ORDER:
            fs = [stereo_fractions(s, m) for s in stereo[t]]
            pooled_full = sum(s["methods"][m]["full"] for s in stereo[t])
            pooled_comp = sum(f["comparable"] for f in fs)
            mmean, mstd = mean_std([f["full"] for f in fs])
            print(f"  {t:7s} {m:10s} pooled={pct(int(pooled_full), pooled_comp)} "
                  f"mean={fms(mmean, mstd)}")
    print("\n== per-center R/S precision ==")
    for t in TIERS:
        for m in METHODS_ORDER:
            fs = [stereo_fractions(s, m) for s in stereo[t]]
            c = [sum(s["methods"][m]["centers"]["correct"] for s in stereo[t]),
                 sum(s["methods"][m]["centers"]["wrong"] for s in stereo[t]),
                 sum(s["methods"][m]["centers"]["unassigned"] for s in stereo[t])]
            mmean, mstd = mean_std([f["rs"] for f in fs])
            print(f"  {t:7s} {m:10s} pooled={pct(c[0], sum(c))} "
                  f"mean={fms(mmean, mstd)}")
    print(f"\nresults written: {out}")


if __name__ == "__main__":
    main()
