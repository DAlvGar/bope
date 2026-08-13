"""Per-failure coordinate inspection for the held-out benchmark.

For every failure in heldout_failures.json: find the disputed bonds via
MCS (perceived vs CCD reference), measure their actual lengths from the
entry coordinates, plus full ring geometry (bond lengths + planarity RMS
via the inertia tensor) for rings holding a disputed bond.  Writes a
compact per-cluster evidence table to heldout_failures_inspect.txt.

This is the evidence base for heldout_failures_analysis.md - every
length cited there comes from this file.  Run from
benchmarks/crystal100:

    uv run python inspect_heldout_failures.py
"""

from __future__ import annotations

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src")
)

from benchmark import load_dataset, run_geometry
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

TIERS = ("main", "lowres")


def mcs_pair(per: Chem.Mol, ref: Chem.Mol):
    mcs = rdFMCS.FindMCS(
        [per, ref],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        timeout=10,
    )
    q = Chem.MolFromSmarts(mcs.smartsString)
    pm = per.GetSubstructMatch(q)
    rm = ref.GetSubstructMatch(q)
    return q, pm, rm


def fit_rms(xyzs) -> float:
    """Planarity RMS of a set of points (A), via the inertia tensor."""
    if len(xyzs) < 4:
        return 0.0
    import numpy as np

    arr = np.array(xyzs, dtype=float)
    cov = np.cov(arr, rowvar=False)
    return float(np.sqrt(max(np.min(np.linalg.eigvalsh(cov)), 0.0)))


def inspect(entry, per, ref):
    """Return multi-line report of the disputed bonds with measurements."""
    lines = []
    q, pm, rm = mcs_pair(per, ref)
    per_to_q = {p: i for i, p in enumerate(pm)}
    ref_to_q = {r: i for i, r in enumerate(rm)}
    xyz = [c for _el, c in entry["atoms"]]

    def dist(i, j):
        a, b = xyz[i], xyz[j]
        return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))

    diff_atoms = set()
    for qi in range(q.GetNumAtoms()):
        for qj in range(qi + 1, q.GetNumAtoms()):
            if not q.GetBondBetweenAtoms(qi, qj):
                continue
            pb = per.GetBondBetweenAtoms(pm[qi], pm[qj])
            rb = ref.GetBondBetweenAtoms(rm[qi], rm[qj])
            if pb is None or rb is None:
                continue
            if pb.GetBondTypeAsDouble() != rb.GetBondTypeAsDouble() or \
                    pb.GetIsAromatic() != rb.GetIsAromatic():
                i, j = pm[qi], pm[qj]
                d = dist(i, j)
                els = sorted((per.GetAtomWithIdx(i).GetSymbol(),
                              per.GetAtomWithIdx(j).GetSymbol()))
                lines.append(
                    f"    diff {''.join(els)} {rb.GetBondTypeAsDouble():.1f}"
                    f"{'a' if rb.GetIsAromatic() else ''}"
                    f"->{pb.GetBondTypeAsDouble():.1f}"
                    f"{'a' if pb.GetIsAromatic() else ''}"
                    f"  measured {d:.3f} A (atoms {i},{j})")
                diff_atoms.add(i)
                diff_atoms.add(j)
    # rings (in ref) holding a disputed atom: dump geometry
    ri = ref.GetRingInfo()
    for ring in sorted(ri.AtomRings(), key=len, reverse=True):
        ring = list(ring)
        if not any(r in ref_to_q and pm[ref_to_q[r]] in diff_atoms for r in ring):
            continue
        mapped = []
        missing = 0
        ring_lens = []
        for k in range(len(ring)):
            a = ring[k]
            b = ring[(k + 1) % len(ring)]
            if a in ref_to_q and b in ref_to_q:
                pa, pb2 = pm[ref_to_q[a]], pm[ref_to_q[b]]
                ring_lens.append(dist(pa, pb2))
            else:
                missing += 1
        for a in ring:
            if a in ref_to_q:
                mapped.append(xyz[pm[ref_to_q[a]]])
            else:
                missing += 1
        els = [ref.GetAtomWithIdx(a).GetSymbol() for a in ring]
        rms = fit_rms(mapped) if mapped else float("nan")
        lines.append(
            f"    ring {''.join(els)} n={len(ring)} rms={rms:.3f} A"
            + (f" bonds={','.join(f'{x:.3f}' for x in ring_lens)}"
               if ring_lens else "")
            + (f" [missing {missing} atoms from MCS]" if missing else ""))
    return lines, per_to_q


def main() -> None:
    with open(os.path.join(HERE, "heldout_failures.json"),
              encoding="utf-8") as fh:
        recs = json.load(fh)
    from collections import Counter
    clusters = Counter(r["sig"] for r in recs)
    out = []
    for sig, n in clusters.most_common():
        members = [r for r in recs if r["sig"] == sig]
        out.append(f"\n### {n}x [{sig}]")
        for r in members:
            entry = load_dataset(
                f"dataset_heldout_{r['tier']}_k{r['bucket']}.json",
            )
            entry = next(e for e in entry
                         if e["pdb"] == r["pdb"] and e["het"] == r["het"])
            per, _ = run_geometry(entry)
            ref = Chem.MolFromSmiles(entry["smiles"])
            if per is None or ref is None:
                out.append(f"  {r['pdb']} {r['het']} ({r['tier']} k{r['bucket']})"
                           f" ob_ok={r['ob_ok']} dist_ok={r['dist_ok']}")
                continue
            lines, _p2q = inspect(entry, per, ref)
            out.append(f"  {r['pdb']} {r['het']} ({r['tier']} k{r['bucket']})"
                       f" ob_ok={r['ob_ok']} dist_ok={r['dist_ok']}")
            out.extend(lines)
            if not lines:
                out.append("    (no bond diffs - H-only case)")
    out_path = os.path.join(HERE, "heldout_failures_inspect.txt")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"wrote {out_path} ({len(out)} lines)")


if __name__ == "__main__":
    main()
