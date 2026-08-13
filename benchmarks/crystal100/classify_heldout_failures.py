"""Mechanism classification of the held-out failures.

Tags every cluster signature with a failure mechanism using (a) the
structured diff signature, (b) whether OpenBabel / the distance baseline
recovered the entry, and (c) the coordinate evidence in
heldout_failures_inspect.txt (read by hand for the major clusters when
the report was written).  The tags are rule-based with manual
verification of representative clusters; the tally feeds
heldout_failures_analysis.md.  A few entries have multiple causes (e.g.
7FOZ WD0 is both a false-positive-aromatic ring and a documented data
limitation; 1BPE DTP is a formula-junk entry whose root cause is
triphosphate charge emission) - the rules assign the first match.

Run from benchmarks/crystal100:

    uv run python classify_heldout_failures.py
"""

from __future__ import annotations

import json
import os
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, "heldout_failures.json"),
          encoding="utf-8") as fh:
    recs = json.load(fh)

# --- mechanism rules, checked in order ------------------------------
MECHANISMS: list[tuple[str, list[str]]] = [
    ("crash", []),  # kind == "None", checked first in classify()
    ("metal/other-element", ["MoO", "OV ", "AsO", "Mo "]),
    ("formula-junk", ["+CH", "-CH", "+PH", "-PH", "+OH"]),
    ("nitro", ["NO "]),
    ("phosphate", ["OP "]),
    ("sulfone/sulfinate", ["OS ", "NS "]),
    ("false-positive-aromatic", ["1.0->1.5a", "2.0->1.5a", "3.0->1.5a"]),
    ("aromatic-rejected", ["1.5a->"]),
    ("guanidinium/amidine", ["CN 1.0->2.0", "CN 2.0->1.0"]),
    ("enol-imine/lactam", ["CO 1.0->2.0", "CO 2.0->1.0"]),
    ("C=C-length-rule", ["CC 1.0->2.0", "CC 2.0->1.0", "CC 3.0->2.0"]),
    ("H-only-tautomer", ["NH0->1", "NH1->0"]),
    ("charge/valence-edge", []),  # everything left
]


def classify(sig: str, kind: str) -> str:
    if kind == "None":
        return "crash"
    for name, needles in MECHANISMS:
        if name == "crash":
            continue
        if any(nd in sig for nd in needles):
            return name
    return "charge/valence-edge"


tally = Counter()
per_cluster: dict[str, list[tuple[str, str, str, str, bool, bool]]] = defaultdict(list)

for r in recs:
    mech = classify(r["sig"], r["kind"])
    tally[mech] += 1
    per_cluster[mech].append((r["pdb"], r["het"], r["sig"], r["kind"],
                              r["ob_ok"], r["dist_ok"]))

print(f"total failures: {len(recs)}  distinct signatures: "
      f"{len({r['sig'] for r in recs})}")
print(f"\n{'mechanism':28s} {'n':>4s} {'ob-recovers':>12s} "
      f"{'dist-recovers':>13s}")
for mech, n in tally.most_common():
    mem = per_cluster[mech]
    ob = sum(1 for m in mem if m[4])
    di = sum(1 for m in mem if m[5])
    print(f"{mech:28s} {n:4d} {ob:12d} {di:13d}")

print("\n-- examples per mechanism --")
for mech, n in tally.most_common():
    mem = per_cluster[mech]
    seen = set()
    shown = 0
    print(f"\n[{mech}] {n} failures")
    for pdb, het, sig, kind, ob, di in mem:
        if pdb in seen:
            continue
        seen.add(pdb)
        print(f"  {pdb} {het} ob={ob} dist={di} :: {sig[:110]}")
        shown += 1
        if shown >= 4:
            break
