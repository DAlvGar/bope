"""Smoke tests: the committed benchmark numbers must not drift.

These tests re-run the in-house geometry tier over the two committed
crystal100 fix datasets (coordinates embedded, fully offline) and
compare the aggregate metrics against the numbers recorded in the
committed sidecars (``results_dataset_fix_*.json``).  They are the
quick divergence detector: any change to perception that shifts
recovery on real deposited coordinates fails here in seconds, without
re-running the full evaluation pipeline.

The sidecar is the source of truth.  A deliberate, validated change
regenerates it with the benchmark harness and commits both together::

    uv run python benchmarks/crystal100/benchmark.py \
        --dataset dataset_fix_main.json
    uv run python benchmarks/crystal100/benchmark.py \
        --dataset dataset_fix_res250-300.json

Only the geometry tier is checked: it is the package's own code,
deterministic given the pinned RDKit (``rdkit>=2025.9.6,<2026``).
OpenBabel / distance / rdDetermineBonds rows are dependency-version
sensitive and are not part of this smoke.
"""

from __future__ import annotations

import json
import os
import sys

# Ensure the src directory is on the path when running without install
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
# The crystal100 harness: dataset loading, run_geometry, _metrics.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "benchmarks", "crystal100"))

import benchmark as _bench
import pytest
from rdkit import Chem

_BENCH = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "crystal100")

DATASETS = ["dataset_fix_main", "dataset_fix_res250-300"]


def _sidecar_geometry_counts(stem: str) -> dict:
    path = os.path.join(_BENCH, f"results_{stem}.json")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)["methods"]["geometry"]


def _aggregate_geometry_counts(entries: list[dict]) -> dict:
    """Replicate benchmark.py's aggregate loop for the geometry tier."""
    counts = [0, 0, 0, 0, 0]  # sanitize, formula, graph, exact, addh
    recovery = 0
    for entry in entries:
        ref = Chem.MolFromSmiles(entry["smiles"])
        mol, _err = _bench.run_geometry(entry)
        s, f, g, e, h = _bench._metrics(mol, ref)
        counts[0] += s
        counts[1] += f
        counts[2] += g
        counts[3] += e
        counts[4] += h
        if f and g and h:
            recovery += 1
    return {
        "sanitize": counts[0],
        "formula": counts[1],
        "graph": counts[2],
        "exact": counts[3],
        "addh": counts[4],
        "recovery": recovery,
    }


@pytest.mark.parametrize("stem", DATASETS)
def test_geometry_tier_matches_committed_sidecar(stem: str) -> None:
    ds_file = os.path.join(_BENCH, f"{stem}.json")
    if not os.path.exists(ds_file):
        pytest.skip("benchmark datasets not present (source checkout only)")

    entries = _bench.load_dataset(f"{stem}.json")
    actual = _aggregate_geometry_counts(entries)
    expected = _sidecar_geometry_counts(stem)
    assert actual == expected, (
        f"geometry tier diverged from committed sidecar ({stem}):\n"
        f"  actual:   {actual}\n"
        f"  expected: {expected}\n"
        "If this is a deliberate perception change, regenerate the sidecars with\n"
        "  uv run python benchmarks/crystal100/benchmark.py --dataset "
        f"{stem}.json\n"
        "and commit both together."
    )
