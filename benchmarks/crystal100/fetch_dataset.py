"""Fetch 100+ unique-ligand crystal structures from RCSB with CCD ground truth.

For the bope crystal-set benchmark (benchmarks/crystal100/):

* Search: protein-only PDB entries, exactly one non-polymer entity
  (single-ligand complex), resolution within [--min-res, --max-res].
  The candidate id list is fetched in full (paginated), then shuffled
  with a fixed seed so the collected set is a fair sample of the pool.
* Per entry: download the PDB, read the ligand HET code from the
  HETATM records (disambiguated via the entry's non-polymer entity when
  the file carries extra HET codes such as modified residues), fetch
  the CCD canonical SMILES for that HET code (ground-truth bond
  orders), extract the ligand residue (altloc A / highest-occupancy
  handling), require the heavy-atom count to match the SMILES exactly
  (complete, well-ordered residue).
* Uniqueness: every ligand in the dataset has a distinct HET code.
* Exclusion list: metals, ions and common crystallographic additives
  (glycerol, ethylene glycol, buffers...) - dataset curation, not
  perception filtering.
* Each entry records its resolution (parsed from the PDB header) so a
  benchmark can group tiers (e.g. <= 2.0 A vs 2.5-3.0 A).

Parallel worker pool (--workers, default 6) for the per-entry network
work; checkpointing is lock-protected and resumable - dataset.json is
rewritten after every accepted entry and .processed.txt records every
candidate seen, so a re-run continues exactly where it stopped.

Examples::

    # main tier
    uv run python benchmarks/crystal100/fetch_dataset.py
    # poorer-resolution tier, separate dataset file
    uv run python benchmarks/crystal100/fetch_dataset.py \
        --min-res 2.5 --max-res 3.0 --out dataset_res250-300.json
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

# Runs with the package installed (`uv sync`) - rdkit is a core dependency.
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
NONPOLY_URL = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{}/{}"
CCD_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"
PDB_URL = "https://files.rcsb.org/download/{}.pdb"

HERE = os.path.dirname(os.path.abspath(__file__))
PROCESSED = os.path.join(HERE, ".processed.txt")

#: Metals, ions and common crystallographic additives - never the ligand
#: a bond-order benchmark is about.
BLACKLIST = {
    # metals
    "FE", "ZN", "MN", "MG", "CA", "CO", "NI", "CU", "NA", "K", "CD", "HG",
    "PT", "AU", "AG", "AL", "LI", "RB", "CS", "SR", "BA", "CR", "MO", "V",
    "W", "SE", "PB", "TL", "GA", "GE", "IN", "SN", "YB", "LA", "CE", "EU",
    "GD", "TB", "DY", "HO", "ER", "TM", "LU", "SC", "TI", "Y", "ZR", "NB",
    "HF", "TA", "RE", "OS", "IR", "PD", "RH", "RU", "TC",
    # ions
    "CL", "BR", "I", "F", "PO4", "SO4", "NO3", "PER", "NH4", "CYN", "SCN",
    # inorganic clusters (no organic framework - 8TC3 slipped through into
    # the low-res tier as SF4 Fe4S4)
    "SF4",
    # solvents / additives / buffers
    "GOL", "EDO", "EOH", "ACT", "ACE", "PGE", "MPD", "PEG", "DMS", "DTT",
    "BME", "TRS", "HEPES", "MES", "MOPS", "BIC", "IMZ", "1PE", "2PE",
    "PG4", "PG5", "P6G", "1PG", "3PG", "FMT", "PEG400", "Urea", "UREA",
    "TAR", "CIT", "FLC", "HOH",
}

_RES_RE = re.compile(r"REMARK\s+2\s+RESOLUTION\.\s*([\d.]+)")
_ELEMENTS2 = {
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl", "Ar", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np",
    "Pu",
}


def get(url: str, retries: int = 3, timeout: float = 30.0) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            time.sleep(0.5 * (attempt + 1))
    raise last  # type: ignore[misc]


def search_entries(min_res: float, max_res: float) -> list[str]:
    """All candidate ids for the tier; paginated in full."""
    def page(start: int) -> dict:
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "range",
                        "value": {"from": min_res, "to": max_res}}},
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                        "operator": "equals", "value": 1}},
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                        "operator": "exact_match", "value": "Protein (only)"}},
                ],
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": start, "rows": 10000}},
        }
        req = urllib.request.Request(
            SEARCH_URL, data=json.dumps(query).encode(), method="POST",
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.load(resp)

    first = page(0)
    total = first["total_count"]
    ids = [h["identifier"] for h in first["result_set"]]
    for start in range(10000, total, 10000):
        data = page(start)
        ids.extend(h["identifier"] for h in data["result_set"])
    print(f"search: {total} candidates in tier, got {len(ids)}")
    return ids


def _element_from_name(name: str) -> str:
    """Element symbol when the PDB element column is blank (legacy files)."""
    n = name.lstrip("0123456789")
    if len(n) >= 2 and n[1].islower() and n[:2].capitalize() in _ELEMENTS2:
        return n[:2].capitalize()
    return n[:1]


def extract_residues(pdb_text: str, het: str) -> list[list]:
    """Heavy-atom (element, xyz) lists for every residue instance of *het*.

    A PDB file can carry several copies of the same HET group (multiple
    instances of one non-polymer entity, e.g. 4NBP has two TLA residues)
    - merging them by atom name would produce a chemically inconsistent
    atom set (4NBP: one O picked from the symmetry-related copy, at the
    other end of the cell).  Each (chain, seqnum, icode) residue is
    therefore extracted separately.

    Altloc handling within a residue: keep the highest-occupancy copy
    per atom name, with altloc rank (A > blank > B) as tiebreak - the
    same convention the parent project's loader uses.  A residue whose
    atoms duplicate after altloc resolution is dropped.
    """
    _ALT_RANK = {"A": 2, " ": 1, "B": 0}
    instances: dict[tuple[str, str, str], dict] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != het:
            continue
        key = (line[21:22], line[22:26].strip(), line[26:27])
        by_name = instances.setdefault(key, {})
        name = line[12:16].strip()
        altloc = line[16:17]
        occ = float(line[54:60]) if line[54:60].strip() else 1.0
        el = line[76:78].strip()
        if not el:
            el = _element_from_name(name)
        if el.upper() == "H":
            continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        prev = by_name.get(name)
        if prev is None:
            by_name[name] = (el, altloc, occ, *xyz)
        else:
            _p_occ, _p_rank = prev[2], _ALT_RANK.get(prev[1], 0)
            _o_rank = _ALT_RANK.get(altloc, 0)
            if occ > _p_occ or (occ == _p_occ and _o_rank > _p_rank):
                by_name[name] = (el, altloc, occ, *xyz)
    out = []
    for by_name in instances.values():
        if not by_name:
            continue
        out.append([[el, [x, y, z]] for el, _a, _o, x, y, z in by_name.values()])
    return out


def het_codes(pdb_text: str) -> set[str]:
    return {line[17:20].strip() for line in pdb_text.splitlines()
            if line.startswith("HETATM") and line[17:20].strip() != "HOH"}


def resolution_of(pdb_text: str) -> float | None:
    m = _RES_RE.search(pdb_text)
    return float(m.group(1)) if m else None


def nonpolymer_het(pdb_id: str) -> str | None:
    """HET code of the entry's non-polymer entity (entry API path)."""
    entry = json.loads(get(ENTRY_URL.format(pdb_id)))
    eids = (entry.get("rcsb_entry_container_identifiers") or {}).get(
        "non_polymer_entity_ids") or []
    if len(eids) != 1:
        return None
    ent = json.loads(get(NONPOLY_URL.format(pdb_id, eids[0])))
    het = (ent.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get(
        "chem_ref_def_id")
    return str(het).strip() if het else None


def process_candidate(pdb_id: str, args, ctx) -> None:
    """Evaluate one candidate; append to the checkpoint on accept."""
    try:
        pdb_text = get(PDB_URL.format(pdb_id)).decode("utf-8", "replace")
        codes = het_codes(pdb_text)
        if len(codes) == 1:
            het = codes.pop()
        else:
            het = nonpolymer_het(pdb_id)
            if het is None or het not in codes:
                ctx.skip("het_ambiguous")
                return
        if het in BLACKLIST:
            ctx.skip(f"blacklist:{het}")
            return
        ccd = json.loads(get(CCD_URL.format(het)))
        smiles = (ccd.get("rcsb_chem_comp_descriptor") or {}).get("SMILES")
        if not smiles:
            ctx.skip("ccd_no_smiles")
            return
        ref = Chem.MolFromSmiles(smiles)  # type: ignore[attr-defined]
        if ref is None:
            ctx.skip("smiles_unparseable")
            return
        n_heavy = ref.GetNumHeavyAtoms()
        if not (5 <= n_heavy <= 60):
            ctx.skip("size_out_of_range")
            return
        instances = extract_residues(pdb_text, het)
        if not instances:
            ctx.skip("residue_not_found")
            return
        # prefer a CCD-complete instance (one copy may be a partial
        # symmetry-related duplicate); fall back to the largest.
        atoms = next((r for r in instances if len(r) == n_heavy), None)
        if atoms is None:
            atoms = max(instances, key=len)
        if len(atoms) != n_heavy:
            ctx.skip("atom_count_mismatch")
            return
        if any(het == d["het"] for d in ctx.dataset):
            ctx.skip("duplicate_het")
            return
        formula = rdMolDescriptors.CalcMolFormula(ref)  # type: ignore[attr-defined]
        entry = {
            "pdb": pdb_id,
            "het": het,
            "name": (ccd.get("name") or "").strip(),
            "smiles": Chem.MolToSmiles(ref),  # type: ignore[attr-defined]
            # stereo ground truth for the stereo benchmark; None when the
            # CCD has no stereo (most achiral HET codes)
            "smiles_stereo": (ccd.get("rcsb_chem_comp_descriptor") or {}).get(
                "SMILES_stereo"
            ),
            "formula": formula,
            "heavy_atoms": n_heavy,
            "resolution": resolution_of(pdb_text),
            "atoms": atoms,
        }
        ctx.accept(entry)
    except Exception as exc:  # noqa: BLE001 - network hiccups
        ctx.skip(f"error:{type(exc).__name__}")


class Context:
    """Shared mutable state across workers (all access under _lock)."""

    def __init__(self, out_json: str):
        self.out_json = out_json
        self.dataset: list[dict] = []
        self.processed: set[str] = set()
        self.skipped: dict[str, int] = {}
        self.lock = threading.Lock()
        if os.path.exists(out_json):
            with open(out_json, encoding="utf-8") as fh:
                self.dataset = json.load(fh)
            self.processed = {d["pdb"] for d in self.dataset}
            print(f"resume: {len(self.dataset)} entries already collected")
        if os.path.exists(PROCESSED):
            with open(PROCESSED, encoding="utf-8") as fh:
                self.processed |= {line.strip() for line in fh if line.strip()}

    def skip(self, reason: str) -> None:
        with self.lock:
            self.skipped[reason] = self.skipped.get(reason, 0) + 1

    def accept(self, entry: dict) -> None:
        with self.lock:
            self.dataset.append(entry)
            with open(self.out_json, "w", encoding="utf-8") as fh:
                json.dump(self.dataset, fh, indent=1, separators=(",", ": "))
            with open(PROCESSED, "a", encoding="utf-8") as fh:
                fh.write(entry["pdb"] + "\n")
            print(f"[{len(self.dataset)}/{args.target}] {entry['pdb']} "
                  f"{entry['het']} {entry['formula']} "
                  f"({entry['resolution']} A) {entry['name'][:40]}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--min-res", type=float, default=0.0)
    parser.add_argument("--max-res", type=float, default=2.0)
    parser.add_argument("--target", type=int, default=100)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="dataset.json")
    global args
    args = parser.parse_args()

    out_json = os.path.join(HERE, args.out)
    ctx = Context(out_json)

    if len(ctx.dataset) >= args.target:
        print(f"already at target ({len(ctx.dataset)} entries)")
        return

    entries = search_entries(args.min_res, args.max_res)
    random.Random(args.seed).shuffle(entries)
    todo = [p for p in entries if p not in ctx.processed]
    print(f"processing {len(todo)} candidates (of {len(entries)})")

    def worker(i: int) -> None:
        idx = i
        while len(ctx.dataset) < args.target and idx < len(todo):
            process_candidate(todo[idx], args, ctx)
            idx += args.workers
            time.sleep(0.25)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(worker, range(args.workers)))

    if len(ctx.dataset) > args.target:
        # worker race: several workers can pass the target check before the
        # count updates, so the file can overshoot by a few accepts - truncate
        # to the acceptance order (the first target entries).
        del ctx.dataset[args.target:]
        with open(out_json, "w", encoding="utf-8") as fh:
            json.dump(ctx.dataset, fh, indent=1, separators=(",", ": "))
        print(f"truncated to {len(ctx.dataset)} entries (worker overshoot)")

    print(f"\ndone: {len(ctx.dataset)} entries in {out_json}")
    for reason, count in sorted(ctx.skipped.items(), key=lambda kv: -kv[1]):
        print(f"  skipped {count:4d}: {reason}")


if __name__ == "__main__":
    main()
