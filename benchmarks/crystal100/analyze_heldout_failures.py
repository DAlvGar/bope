"""Per-case failure analysis of the held-out benchmark.

For every held-out entry where the geometry tier fails recovery
(formula AND graph AND AddHs), dump an atom-mapped diff of the
perceived molecule vs the CCD reference, plus whether OpenBabel and
the distance baseline recovered the same entry (a case another method
recovers proves the coordinates support the CCD chemistry, so the
failure is a perception defect, not a coordinate artifact).

Regenerates heldout_failures.json (the raw per-failure records) and
prints the cluster inventory.  Run from benchmarks/crystal100:

    uv run python analyze_heldout_failures.py
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src")
)

from benchmark import (
    _metrics,
    _neutral,
    load_dataset,
    run_distance,
    run_geometry,
    run_openbabel,
)
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

TIERS = ("main", "lowres")
BUCKETS = 5


def formula(mol: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(_neutral(mol))


def bond_diff_detail(per: Chem.Mol, ref: Chem.Mol) -> dict:
    """Atom-mapped bond-order and H differences via MCS."""
    out = {"mcs_n": 0, "bond_diffs": [], "h_diffs": [], "unmapped_per": [],
           "unmapped_ref": [], "charges": [], "note": ""}
    try:
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
        out["mcs_n"] = mcs.numAtoms
        for qi in range(q.GetNumAtoms()):
            for qj in range(qi + 1, q.GetNumAtoms()):
                if not q.GetBondBetweenAtoms(qi, qj):
                    continue
                pb = per.GetBondBetweenAtoms(pm[qi], pm[qj])
                rb = ref.GetBondBetweenAtoms(rm[qi], rm[qj])
                if pb is None or rb is None:
                    continue
                po, ro = pb.GetBondTypeAsDouble(), rb.GetBondTypeAsDouble()
                pa, ra = pb.GetIsAromatic(), rb.GetIsAromatic()
                if po != ro or pa != ra:
                    els = sorted(
                        (per.GetAtomWithIdx(pm[qi]).GetSymbol(),
                         per.GetAtomWithIdx(pm[qj]).GetSymbol())
                    )
                    out["bond_diffs"].append({
                        "els": els, "per": f"{po:.1f}{'a' if pa else ''}",
                        "ref": f"{ro:.1f}{'a' if ra else ''}",
                    })
        # per-atom H + charge diffs on the matched atoms
        for qi in range(q.GetNumAtoms()):
            pa = per.GetAtomWithIdx(pm[qi])
            ra = ref.GetAtomWithIdx(rm[qi])
            ph = pa.GetTotalNumHs()
            rh = ra.GetTotalNumHs()
            if ph != rh:
                out["h_diffs"].append({"el": pa.GetSymbol(),
                                       "per": ph, "ref": rh})
            if pa.GetFormalCharge() != ra.GetFormalCharge():
                out["charges"].append({
                    "el": pa.GetSymbol(),
                    "per": pa.GetFormalCharge(), "ref": ra.GetFormalCharge(),
                })
        per_set, ref_set = set(pm), set(rm)
        for i in range(per.GetNumAtoms()):
            if i not in per_set:
                a = per.GetAtomWithIdx(i)
                out["unmapped_per"].append(
                    f"{a.GetSymbol()}H{a.GetTotalNumHs()}")
        for i in range(ref.GetNumAtoms()):
            if i not in ref_set:
                a = ref.GetAtomWithIdx(i)
                out["unmapped_ref"].append(
                    f"{a.GetSymbol()}H{a.GetTotalNumHs()}")
    except Exception as exc:  # noqa: BLE001
        out["note"] = f"mcs failed: {exc}"
    return out


def signature(d: dict) -> str:
    """Compact cluster key for a failure."""
    parts = [d["kind"]]
    detail = d.get("detail") or {}
    if "bond_diffs" in detail:
        parts.extend(
            f"{''.join(bd['els'])} {bd['ref']}->{bd['per']}"
            for bd in sorted(detail["bond_diffs"],
                             key=lambda b: (tuple(b["els"]), b["ref"], b["per"]))
        )
    if "h_diffs" in detail:
        parts.extend(
            f"{hd['el']}H{hd['ref']}->{hd['per']}"
            for hd in sorted(detail["h_diffs"], key=lambda h: (h["el"], h["ref"]))
        )
    if "unmapped_per" in detail:
        parts.extend(f"+{u}" for u in sorted(detail["unmapped_per"]))
    if "unmapped_ref" in detail:
        parts.extend(f"-{u}" for u in sorted(detail["unmapped_ref"]))
    return " | ".join(parts)


def main() -> None:
    all_fail = []
    for tier in TIERS:
        for k in range(1, BUCKETS + 1):
            entries = load_dataset(f"dataset_heldout_{tier}_k{k}.json")
            for entry in entries:
                mol, err = run_geometry(entry)
                ref = Chem.MolFromSmiles(entry["smiles"])
                s, f, g, e, h = _metrics(mol, ref)
                ok = f and g and h
                ob_mol, _ = run_openbabel(entry)
                ob_s, ob_f, ob_g, ob_e, ob_h = _metrics(ob_mol, ref)
                ob_ok = ob_f and ob_g and ob_h
                dist_mol, _ = run_distance(entry)
                d_s, d_f, d_g, d_e, d_h = _metrics(dist_mol, ref)
                dist_ok = d_f and d_g and d_h
                if ok:
                    continue
                rec = {
                    "tier": tier, "bucket": k, "pdb": entry["pdb"],
                    "het": entry["het"], "name": entry["name"],
                    "kind": ("None" if mol is None else
                             "formula" if not f else
                             "graph" if not g else "addh"),
                    "err": err,
                    "formula_per": formula(mol) if mol is not None else None,
                    "formula_ref": formula(ref),
                    "per_smiles": (Chem.MolToSmiles(_neutral(mol))
                                   if mol is not None else None),
                    "ref_smiles": Chem.MolToSmiles(_neutral(ref)),
                    "ob_ok": ob_ok, "dist_ok": dist_ok,
                }
                if mol is not None:
                    try:
                        per_n = _neutral(mol)
                        rec["detail"] = bond_diff_detail(per_n, _neutral(ref))
                    except Exception as exc:  # noqa: BLE001
                        rec["detail"] = {"note": f"detail failed: {exc}"}
                rec["sig"] = signature(rec)
                all_fail.append(rec)
    out_path = os.path.join(HERE, "heldout_failures.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(all_fail, fh, indent=1, ensure_ascii=False)

    # ---- cluster by signature ----
    from collections import Counter
    clusters = Counter(r["sig"] for r in all_fail)
    print(f"total failures: {len(all_fail)}")
    print(f"distinct signatures: {len(clusters)}")
    print("\n== top signatures ==")
    for sig, n in clusters.most_common(40):
        ex = next(r for r in all_fail if r["sig"] == sig)
        print(f"\n{n:3d} x  [{sig}]")
        print(f"     ex: {ex['pdb']} {ex['het']} ({ex['tier']} k{ex['bucket']}) "
              f"ob_ok={ex['ob_ok']} dist_ok={ex['dist_ok']}")
        print(f"     per: {ex['per_smiles']}")
        print(f"     ref: {ex['ref_smiles']}")
    print("\n== kind counts ==")
    for kind, n in Counter(r["kind"] for r in all_fail).most_common():
        print(f"  {kind:10s} {n}")
    print("\n== also-failed-by-ob (coordinate-bound candidate) vs ob-ok ==")
    print(f"  ob also fails: {sum(1 for r in all_fail if not r['ob_ok'])}")
    print(f"  ob recovers:   {sum(1 for r in all_fail if r['ob_ok'])}")
    for tier in TIERS:
        n_ob_ok = sum(1 for r in all_fail if r["tier"] == tier and r["ob_ok"])
        n = sum(1 for r in all_fail if r["tier"] == tier)
        print(f"  {tier}: {n_ob_ok}/{n} of bope's failures recovered by OpenBabel")
    print(f"\njson: {out_path}")


if __name__ == "__main__":
    main()
