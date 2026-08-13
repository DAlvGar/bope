# Held-out failure analysis (131 remaining cases) - working report

Status: **analysis refreshed (2026-08-13) against the current failure
list** (`heldout_failures.json`, 131 entries = main 59 + lowres 72,
regenerated with the current pipeline on the current held-out datasets).
The 2026-08-12 analysis of the original 156 failures is superseded: the
bug families it identified are fixed and re-evaluated (see below), and
the 131 failures analyzed here are what remains after those fixes.  The
frozen benchmark numbers in `results_heldout.md` are untouched
(code-freeze discipline: fixes land in clearly-separated follow-up
commits and are never re-measured into the published numbers).  The fix
run is a NEW labeled run: `results_heldout_fix.md` (+ per-bucket
sidecars `results_dataset_heldout_fix_*_k*.json` /
`results_stereo_...`).

Question answered: of the 156 geometry-tier failures on the original
held-out run (600 entries, 5 buckets x 60 per tier), were they bugs or
per-case tuning tweaks?  Answer: a mix - **~45 genuine perception bugs
(5 families + 1 demotion-choice rule), now fixed and re-measured; the
remaining 131 split into ~12 residual defect cases in three still-open
sub-families, ~50 threshold-boundary (one shared set of global cutoffs,
each miss 0.001-0.05 A), ~55 coordinate-bound (the deposit contradicts
the CCD - the benchmark's intrinsic floor), 11 H-only tautomers, and 4
out-of-scope metals (+ 2 metal complexes counted under other
mechanisms).**

## Fix implementation and re-evaluation (2026-08-12)

All 5 bug families fixed in `src/bope` (separate commits, one per family),
plus the demotion-choice rule (8JZ7 FI7, listed under the boundary class
in the 156-case analysis) which turned out to be a genuine 6th fix: when
an atom competes between a C=C and a C=N double, demote the C=N first.
Per-family resolution:

1. **Crash family (13):** the over-valent crashes are fixed by the
   valence-demotion choices below and the nitro fixup; the kekulize crashes
   by the aromatic-set gates (5-ring exo-double block, short-bond strictness)
   plus the slack-drop retry.  Verified: 6Y4 (5KYA), 9TP (4L9Q), WMJ (3QTX)
   recover exactly.
2. **Nitro valence:** nitro fixup - an N with two terminal O's <= 1.45 A is
   +1 charged with one N=O and no H; non-O doubles are demoted first.
3. **False-positive aromatization (V0X class):** 5-ring O/S rings must carry
   at least one genuine double-length edge not incident to the O/S atom
   (gate 6); pure-carbon and pyridine-type rings stay unchecked.
4. **Phosphate P=O gap:** P with 3+ terminal O's emits one P=O on the
   shortest O plus a -1 charged O(-) - same charge-emission machinery as the
   nitro zwitterion.
5. **Azole/S-ring judge:** 5-ring N with an exo substituent is
   pyridinium-type (1 pi) when the ring holds a second heteroatom (X43);
   S/P exo bonds never block a ring (W5C); thiophene-type C-S in the ring is
   counted as conjugation evidence only when short.

Re-evaluation (new labeled run, `results_heldout_fix.md`): bond-order
recovery by the geometry tier rises from 223/300 (74.3%) to **241/300
(80.3%)** on main and 221/300 (73.7%) to **228/300 (76.0%)** on lowres
(+18 / +7 ligands).  Stereo full-string: main 114/121 (94.2%) -> **122/131
(93.1%)** with the comparable set growing 121 -> 131 (+8 exact graphs is the
fix effect), lowres 103/125 (82.4%) -> **107/129 (82.9%)**.  Per-center R/S
99.5% / 99.8%, E/Z 12/13 main.  openbabel / distance / rdDetermineBonds
columns unchanged (no code change in those tiers).  **Zero regressions** in
any bucket.  Tuning-gate numbers unchanged (89/101 main, 74/101 lowres).

**5 remaining data limitations** (kept out of the strict recovery
assertions, documented in `tests/fixtures/heldout_regression.json` +
`test_perception.py`): 5MUY MGT (ribose enol-C2'), 7FOZ WD0 (N=C 1.172 A in
thiazoline), 5ME6 M7G (charged reference, +16 H vs the CCD zwitterion),
6S7B KYH (ribose), 7T2X EMY (ring-chain ambiguous; see the correction
below).  All degrade gracefully (no crash, sane mol, wrong graph) and are
pinned as such in the regression tests.  All five are still present in the
131-failure list.

**Broken-ring refusal (thiazole fix, narrowed):** a pre-existing synthetic
test failure (thiazole at zero noise: ETKDG seed-42 refines the 5-ring to a
fragment - S-C 1.90 A, C-N 1.16 A) was fixed by refusing a molecule whose
SINGLE ring carries both an over-long edge (beyond envelope + slack) and an
under-short edge (below envelope lo - slack), so the OpenBabel fallback the
corpus reserves for S-heterocycle embeds can rebuild it.  A first version
that refused on the over-long edge alone was measured against the crystal
datasets and rejected: it fires on real low-res ligands (porphyrins HEM/HEC/
ZNH/HEB, nucleotides, 21 entries across the 12 datasets) and cost ~2
recoveries + ~14 AddHs.  Real planar rings stretch at 2.5-3.0 A resolution
but never also compress below ~1.25 A, so the both-directions condition
fits the data: zero crystal entries match it, and the final re-run
reproduces the fix-run sidecars exactly.

## What the fixes recovered (156 -> 131)

- **Crash family: 13 -> 0.**  All six kekulize + seven over-valent crashes
  are gone; the three ob-recoverable crashes (4L9Q 9TP, 5KYA 6Y4, 1GSF EAA)
  now recover exactly.
- **Nitro: 12 -> 4.**  Every textbook terminal nitro recovered (4EK8 16K,
  8BGA, 7ENE, 6SUH LVE, 3B67 B67, 2YOH were the named representatives).
  The four survivors are non-terminal sub-cases - see the residual section.
- **False-positive aromatization: 11 -> 9.**  Gate 6 recovered 5ST6 V0X and
  1WQW BT5; 8GBK HEB and 5HE2 F0S still fail but now under other mechanisms
  (formula-junk, C=C length rule).  The nine survivors are the N-only /
  pure-carbon rings the gate deliberately left unchecked - see the residual
  section.
- **Phosphate P=O: 20 -> 18.**  The short P=O cases recovered (1TPB PGH,
  2VF5 GLP, 2I22 I22, 3U4H C8R were the named representatives).  The 18
  survivors are mostly di/tri-phosphates and ester cases - see the boundary
  and coordinate sections.
- **Azole/S-ring judge:** 3QTX X43 (N-substituted thiazole), 7FOY W5C
  (triazole) and 2XDA JPS (thiophene) all recovered.
- **Demotion choice:** 8JZ7 FI7 recovered.
- **Boundary/label movement:** 8QIY VCC (pyrazole C-N floor case - the
  entry that exposed the stale failure file), 7FNL VWR (H-only), 3PD3 A3T
  (lactam-enol-imine) recovered.  The sulfone/sulfinate mechanism is now
  empty: 2O9I 444's S=O perceives correctly (sulfonamide fixup); the entry
  survives only as an aromatic-rejected ring.

## Tally (131 failures, 107 signatures)

| mechanism | n | ob recovers | dist recovers |
|---|---|---|---|
| C=C-length-rule | 27 | 12 | 0 |
| aromatic-rejected | 19 | 7 | 0 |
| guanidinium/amidine | 18 | 6 | 0 |
| phosphate | 18 | 7 | 0 |
| enol-imine/lactam | 12 | 4 | 0 |
| H-only-tautomer | 11 | 2 | 0 |
| false-positive-aromatic | 9 | 1 | 0 |
| charge/valence-edge | 5 | 1 | 0 |
| metal/other-element | 4 | 0 | 0 |
| nitro | 4 | 1 | 0 |
| formula-junk | 4 | 1 | 0 |

(Some entries have multiple causes - e.g. 7FOZ WD0 is both a
false-positive-aromatic ring and a documented data limitation; 1BPE DTP is
a formula-junk entry whose root cause is triphosphate charge emission;
2O9I 444 is an aromatic-rejected ring whose S=O the sulfonamide fixup
already corrected.  The rules assign the first match.)

## Residual defect families (still open, ~12)

### 1. N/C-only ring aromatization gate gap (the remaining false-positive-aromatic class)

Gate 6 (O/S 5-ring double-edge gate) fixed the V0X class but deliberately
left N-only and pure-carbon rings unchecked.  The nine survivors are
exactly that remainder:

- **Genuine saturated-ring aromatizations:** 1R4U OXC (hydantoin: ref
  `O=C1NC(=O)NC(C(=O)O)N1` is a saturated 5-ring, per aromatized it to a
  uracil-like ring; ring C-N 1.325-1.453, rms 0.011), 1IL3 7DG
  (pyrimidinone + saturated dihydro ring: ref `Nc1nc2c(c(=O)[nH]1)CC=N2`,
  per aromatized both), 2XUH TZ4 (ref has a tetrahydro ring
  `...c3c4c(nc5ccccc35)CCCC4`; per aromatized the four saturated C's,
  CH2->1 x4), 2BXT C2D (saturated 6-ring with C=N, CH2->1 x5), 3DCV 55E
  (dihydro-pyrimidinone; **ob recovers** - the coordinates fully support
  the saturated form, so this is a pure gate gap).
- **Deposit-favors-aromatic tension:** 7R7R AWJ (ref pyrazoline
  `N2NC=CN2`; the ring measures 1.344-1.399, planar 0.008 - the deposit
  supports per's aromatic pyrazole), 6DGL GEV (thiazolone: ref has C=O +
  C=C, per aromatized the ring; rms 0.002), 9EUT LBV (porphyrinoid double
  shuffle), 7FOZ WD0 (documented limitation - thiazoline C=N at 1.172 A).

ob recovers 1/9 (3DCV 55E).  The clean gate gap is 5-6 cases; the rest sit
where the deposit geometry is short and planar enough that the aromatic
read is defensible.

### 2. Isoxazole N-O envelope and long-nitro lowres (the 4 nitro survivors)

The terminal-nitro fixup (N with two terminal O's <= 1.45 A -> both N=O,
no H) recovered every textbook nitro.  The survivors are non-terminal:

- **Ring N-O:** 4B7P 9UN and 2AMP I12 - isoxazole rings whose N-O measures
  1.485 / 1.540 A (aromatic N-O is ~1.35; the deposit sits in the
  hydroxylamine single-bond range).  bope reads the N-O as single and the
  ring as non-aromatic, which cascades into amide->imidol and nitro
  scrambling.  4B7P 9UN is ob-recoverable - the coordinates support the
  aromatic isoxazole; a ring-N-O envelope/gate case.
- **Long nitro at lowres:** 5DV6 B4H (N-O 1.381/1.390) and 2F9B N1H
  (N-O 1.465, read as N(O)O) - the blurry nitro where the downstream
  amide/imine tautomer choices fail.

ob recovers only 9UN (1/4).

### 3. Charge emission: sulfonium and triphosphate (charge/valence-edge + formula-junk)

- **5H5F SAM:** S-adenosylmethionine sulfonium S+ (ref `C[SH](CCC...)`).
  bope reads the S+ as C=S plus a valence mess (C-S measured 1.653-1.684;
  ob recovers - the coordinates support the sulfonium).  Needs the +1
  charge-emission machinery (the same pass the P=O/nitro charges use),
  still open.
- **1BPE DTP:** ATP-like triphosphate.  bope falls into a P-O-P
  ring-closure junk (`COP2(O)(O)OP2(O)(O)OP(=O)(O)O`) instead of three
  P=O + charges; ob recovers.  The single-P P=O fix works; multi-P ester
  chains need per-center emission with ester-aware placement.
- The remaining charge/valence-edge entries are not defect cases: 8A0I CNZ
  (NH3->2), 6LW1 EX3 (CH1->0 x2) and 4B75 4VA (NH2->1) are H-only in
  nature; 8U6M VW2 is the deposit-triple coordinate case (C-C 1.257).

### 4. Hetero-double choice: thione vs enone (3NR9 NR9)

ref `O=C1N=C(NCc2cccs2)SC1=Cc1ccc2ncccc2c1` (rhodanine-type: ring C=O +
C=N + exo C=C), per demotes the ring C=N and the exo C=C and reads the
ring sulfur as a thione C=S (`O=C1NC(NCc2cccs2)=S=C1C...`).  The
demotion-choice rule (demote C=N first) does not cover the S-vs-C double
choice; ob recovers.  One case.

## Threshold-boundary class (~50): shared cutoffs, not per-case tuning

All failures here sit on one shared set of global thresholds; each miss is
0.001-0.05 A from a single cutoff.  Loosening any one fixes several cases at
once, but each threshold was set against the tuning set, so every loosen
risks new false positives elsewhere (the V0X class shows the judge is already
too permissive in places).  Treat as a threshold-sensitivity paragraph in the
paper, not a fix list.

- Envelope floors/ceilings: 5RUV W6J pyridine C=N 1.256 (0.024 below the
  1.28 C-N floor - 4 mA past the 0.02 slack; ring planar 0.009), 5U2E 837
  thiophene C-S 1.635/1.654 (0.025/0.006 below the 1.66 C-S floor - the
  1.635 edge sits 5 mA past the 0.02 slack, while the recovered 4L9Q 9TP
  thiophene rides the same slack at 1.642/1.658), 3C21 2BA (cyclic
  phosphate diester: 12-ring at rms 0.413, each P holding one terminal O,
  so the P=O fixup's 2+ terminal-O condition never fires and the 1.515
  P=O reads single; ob recovers).
- C=C length rule near the 1.38 A cutoff: demotions of elongated deposits
  (ref double, per single) - 5E1C 5K8 1.452, 3IMT IW3 1.391, 1KT7 RTL
  1.431/1.412, 3IRX UDR 1.391, 2NXX P1A 1.389, 7SRR 7LD 1.398, 3VBG 03M
  1.470, 8FEC XU0 1.442, 5HE2 F0S 1.468, 3R5M MLO 1.440/1.389, 2UXO TAC
  1.400/1.414, 3VBQ 0F5 1.410 - and the symmetric false-positive doubles
  (ref single, per double) at 1.356-1.370 (4KKO 1RE 1.356, 5P9L 7G9 1.359,
  7EJ9 DYX 1.367).  12 of the 27 are ob-recoverable - the coordinates
  support the CCD; the 1.38 A C=C cutoff sits inside the density of real
  crystal C=C lengths at lowres.
- Alkyne cutoff: 6FZM EE5 C=C 1.286 vs the 1.26 triple cutoff (ob recovers;
  1.286 is long but a crystal C=C is 1.33+, so 1.29 is a safe cutoff).
- Aromatic-rejected with ob-recovers: 3C1K T15, 5JQ8 I73, 6S9X L1W, 6OXS
  NJ4, 5SPZ QJC, 8CCC UCT, 2O9I 444 - the coordinates support aromatic;
  bope's envelope rejected the ring.  Mostly fused-ring or distorted cases
  where one bond sits just outside the envelope (2O9I 444's phenyl measures
  1.42-1.54 A).
- Amidine C=N cutoffs: 2ITW ITQ 1.452 (ob recovers), 1J0B 5PA 1.323,
  1K7W AS1 1.430, 3I4A LN5 1.420.
- Phosphate placement on single-P esters: 1CQ6 PY4 (P-O 1.490/1.430 - the
  fix's P=O choice flips; its C=O 1.304 is also at the double cutoff),
  1LZO PGA 1.505/1.493, 2CZF XMP 1.610/1.482.
- (The previous analysis' pyrazole floor case, 8QIY VCC at C-N 1.263, now
  recovers.)

## Coordinate-bound class (~55): the deposit contradicts the CCD

Verified by direct measurement - no bond of the required length exists
anywhere in the deposit.  OpenBabel fails on nearly all of these too.  This is
the benchmark's intrinsic floor: crystallographic residues where the
refined geometry disagrees with the CCD record.  Document as a limitations
section, do not tune against them.

- No C=C-length bond exists: 4Y9G MKU (prenyl C=C at 1.533), 4KY2 1W3
  1.552, 5TRG 7HJ 1.519, 7OXB 35Z 1.530, 7LID EOL 1.525, 6N8Y KFY 1.515.
- Low-res cyclohexene distortion (rings rms 0.13-0.25, no reliable short
  edge): 2QN2 0MA, 7LAI R78, 5IWD 6EV, 2UXO TAC.
- Ribose C1'-C2' shorter than a single: 6VCO APC 1.309, 7WQI 80I 1.355,
  5MUY MGT 1.296 (ribose enol; limitation).
- Tautomer coordinates (deposit refined the other tautomer): 5STO W9U
  hydroxylamine-vs-oxime (ref C-NH-OH, per C=N-O; the C-N measures 1.276 -
  genuinely oxime-length, supporting per), 8E9Y WE9 / 8E9X WEC amidine
  C=N 1.445/1.452, 7T2M EGK / 7SVN CW8 nitrile
  C=N 1.259/1.258, 3LV6 BMP uracil C-N 1.453/1.493/1.495, 7N44 06I
  pyridinone C=O 1.384-1.406, 4DPI 0N1 / 1FD7 AI1 lactam->enol-imine
  (CN 1.365 + CO 1.437 / CN 1.258 + CO 1.304), 5ZC5 09I, 1ETZ GAS,
  5GZS ARG, 7A4Q QY2, 3RIR BT5, 8QEX UFU.
- Enol-imine/lactam tautomer deposits (the 12-case mechanism): 4PWI ROA
  (1.251/1.247), 2IFC OAA (1.363/1.310), 1TYR 9CR (1.388 + distorted
  ring), 2ETK HFS 1.252, 4JIT 3ZF 1.375, 5SPM S4O 1.395, 3KOC DHI
  1.237/1.236, 2QXG K7I 1.298, 6OP0 A7A 1.395, 7OTY 1IX 1.391, 4FSH SKM
  1.242/1.233, 7VJL 7IF (nitrile 1.134 + CO 1.413).
- Deposit shows a triple the CCD calls single: 7VJL 7IF C-N 1.134,
  8U6M VW2 C-C 1.257 (bope reads these as nitrile/alkyne; the coordinates
  support it).
- Di/tri-phosphate charge placement (the deposit's P-O's are all
  ~1.48-1.52 A, so the P=O vs P-OH choice is not geometry-resolvable):
  1N1E NDE, 7Q3B 8SC, 8GJX ZNT, 9J0F NAI, 6B9U NAI, 4KSY 1SY, 9WN1 4BW,
  1MFZ GDX, 4EAQ ATM, 3LV6 BMP (uracil + phosphate), 5ME6 M7G (+16 H vs
  the CCD zwitterion; limitation).
- False-positive-aromatic entries where the deposit supports the aromatic
  read (see residual section 1): 7R7R AWJ, 6DGL GEV, 9EUT LBV, 7FOZ WD0
  (limitation).

## H-only tautomers (11)

8OV7 W3W, 2F4J VX6, 3KF4 B90, 4E20 0MY, 4ESI 0RB, 5FPO 10L, 5RKH LPZ,
5T68 77V, 6QO3 J9N (pure NH0<->NH1 swaps), plus 7T2M EGK / 7SVN CW8
(nitrile tautomers, counted here by the rules).  ob fails on 9/11.
Unresolvable without H positions - document, do not fix.  (Also H-only in
nature, counted elsewhere by the rules: 8A0I CNZ, 6LW1 EX3, 4B75 4VA.)

## Out of scope

- Metals/other elements (4): 1H9M MOO (MoO2), 6LVQ VO4, 9PTR VO4 (VO4
  vanadate), 1BEH CAC (cacodylate AsO2) - the perception targets organic
  ligands.
- Metal complexes counted under other mechanisms: 6DO0 JY1 (Rh bisphosphine
  catalyst complex - formula-junk; the per-molecule graph fragments into
  phosphine-H junk), 8GBK HEB (heme B at lowres - the Fe-porphyrin
  macrocycle with Fe-N 1.99-2.08 A; the metal coordination breaks the
  graph; formula-junk).

## Follow-up plan (open items only, ordered by leverage)

1. N/C-ring aromatic gate: exclude candidate rings containing any
   tetrahedral carbon (CH2/CH3 by the H count) from the Huckel candidate
   set, or require a genuine short edge (C-C <= 1.40 or C-N/C-O <= 1.35)
   not incident to a heteroatom.  Kills 1R4U OXC, 1IL3 7DG, 2XUH TZ4,
   2BXT C2D, 3DCV 55E (ob-provable).  Check against the tuning set - the
   pure-carbon rings were left unchecked because the judge is most
   permissive there.
2. Isoxazole N-O: treat ring N-O <= ~1.45 A as a candidate aromatic bond
   (4B7P 9UN is ob-provable; 2AMP I12 at lowres follows).
3. Charge emission: sulfonium (S with 3 C and no H -> +1, no C=S - 5H5F
   SAM, ob-provable) and multi-P phosphate chains (per-center P=O + O(-)
   with ester-aware placement - 1BPE DTP).
4. Hetero-double choice: prefer C=O over C=N over C=S when an atom competes
   (3NR9 NR9).
5. Boundary sensitivity analysis for the paper: per-threshold error curves
   (envelope floor/ceiling, C=C cutoff, slack values) with tuning-set
   counterexamples for every loosen.

Discipline: each fix in a clearly-separated commit; the held-out numbers in
`results_heldout.md` are never re-measured.  If a follow-up benchmark is ever
run after fixes, it is a NEW run with its own results file, labeled as such.

## How this was produced

- `heldout_failures.json` - per-failure records: kind (formula/graph),
  per/ref SMILES, whether OpenBabel and the distance baseline recovered the
  entry, and an MCS atom-mapped diff signature (bond-order and H diffs vs
  the CCD reference).  Generated by `analyze_heldout_failures.py`.
- `heldout_failures_inspect.txt` - the coordinate evidence: every disputed
  bond's measured length from the deposit, plus ring geometry (bond lengths +
  planarity RMS via inertia tensor).  Generated by
  `inspect_heldout_failures.py`.  Every length cited above comes from this
  file.
- `classify_heldout_failures.py` - the mechanism tally (rule-based tags,
  representative clusters hand-verified against the inspect table).
- Recovery flags are the discriminator: a case OpenBabel (or the distance
  baseline) recovers proves the coordinates support the CCD chemistry, so the
  failure is a perception defect, not a coordinate artifact.  Caveat: ob also
  recovers some genuinely ambiguous cases via zwitterionic charge guesses
  (e.g. W9U, long-nitro), so "ob recovers" alone is not proof of a bope bug.

Reproduce from `benchmarks/crystal100`:

```bash
uv run python analyze_heldout_failures.py    # regenerate heldout_failures.json
uv run python inspect_heldout_failures.py    # regenerate heldout_failures_inspect.txt
uv run python classify_heldout_failures.py   # mechanism tally
```
