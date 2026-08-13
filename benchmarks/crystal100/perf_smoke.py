"""Perf smoke: time the geometry tier on the committed crystal100 set.

Detects order-of-magnitude performance regressions (a refactor that
reintroduces quadratic work, an accidental per-call import, ...) in
seconds, without re-running the full evaluation pipeline.

Times ``perceive_bond_orders_geometric`` over the 101-ligand
``dataset_fix_main.json`` set (bond graph pre-built, identical input
construction to ``benchmark.run_geometry``) for ``--repeats`` rounds
and compares the per-perception wall time against a recorded baseline.
A wall-clock assertion inside pytest would be a flaky test, so this is
a standalone script::

    uv run python benchmarks/crystal100/perf_smoke.py

Exit code 0: within the ceiling; 1: regression.

Baseline: 1.78 ms/perception, recorded 2026-08-13 on this machine with
20,200 perceptions (101 ligands x 200 repeats) at commit 1073e56 - the
post-refactor main, measured identically before and after the
``_GeometricPerceiver`` refactor.  The ceiling is intentionally
generous: this guards against order-of-magnitude regressions, not
measurement noise.
"""

from __future__ import annotations

import argparse
import time

from benchmark import load_dataset

from bope.geometry import perceive_bond_orders_geometric
from bope.helpers import _build_rwmol, _distance_bond_graph

BASELINE_MS = 1.78  # per perception @ 1073e56, 2026-08-13, 20,200 perceptions
DEFAULT_CEILING = 3.0  # generous: guards order-of-magnitude regressions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--repeats", type=int, default=200,
        help="times to run over the dataset (default: 200)",
    )
    parser.add_argument(
        "--ceiling", type=float, default=DEFAULT_CEILING,
        help=f"max ratio vs the {BASELINE_MS:.2f} ms baseline "
             f"(default: {DEFAULT_CEILING})",
    )
    args = parser.parse_args()

    entries = load_dataset("dataset_fix_main.json")
    jobs = []
    for entry in entries:
        elements = [el for el, _xyz in entry["atoms"]]
        coords = [xyz for _el, xyz in entry["atoms"]]
        rwmol = _build_rwmol(elements, coords)
        n = rwmol.GetNumAtoms()
        positions = [rwmol.GetConformer().GetAtomPosition(i) for i in range(n)]
        graph = _distance_bond_graph(elements, positions)
        jobs.append((elements, coords, graph))

    t0 = time.perf_counter()
    for _ in range(args.repeats):
        for elements, coords, graph in jobs:
            perceive_bond_orders_geometric(elements, coords, graph)
    dt = time.perf_counter() - t0

    total = len(jobs) * args.repeats
    ms = dt / total * 1000.0
    ratio = ms / BASELINE_MS
    print(f"molecules={len(jobs)} repeats={args.repeats} "
          f"total_perceptions={total}")
    print(f"wall={dt:.3f}s  per_perception={ms:.2f}ms  "
          f"baseline={BASELINE_MS:.2f}ms  ratio={ratio:.2f}")
    if ratio <= args.ceiling:
        print(f"OK: within the {args.ceiling:.1f}x ceiling")
        return 0
    print(f"REGRESSION: {ratio:.2f}x exceeds the {args.ceiling:.1f}x ceiling "
          f"(baseline {BASELINE_MS:.2f} ms/perception @ 1073e56)")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
