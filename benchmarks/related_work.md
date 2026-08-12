# Related work: bond-order and stereo perception from 3D coordinates

Hand-written working survey (2026-08-12) for the bope paper's
related-work section.  Every claim is sourced; all citations verified
2026-08-12.

The problem: protein-bound ligands are deposited as bare coordinates -
the PDB does not require bond orders or stereo in the coordinate files,
and refinement tools routinely leave them undefined or wrong.  Anything
downstream (interaction analysis, strain calculation, docking prep,
ML featurization) must re-derive the chemistry from geometry.  The
RCSB Chemical Component Dictionary (CCD) defines the authoritative
chemistry per HET code, but coordinates are what a perception algorithm
sees.  This survey covers who else solves that problem, how, and what
they did or did not evaluate.

## 1. Classical perception: distance cutoffs + valence rules

- **Sayle 2001** [1] - "Cruft to content: Perception of bond order in
  organic molecules", Daylight MUG.  The canonical rule-based treatment:
  bond detection by element-pair distance cutoffs, then valence/charge
  constraints and aromatic ring perception.  No empirical fitting beyond
  the cutoffs.  Unrefereed (user-group proceedings) but hugely
  influential: it is the algorithm behind Open Babel's
  `PerceiveBondOrders`.
- **Baber & Hodgkin 1992** [2] - OXBRIDGE: rule-based connectivity
  assignment for CSD entries, with a dedicated test set of 91 CSD
  molecules (reported ~82% connectivity success).  The first systematic
  evaluation of connectivity perception against experimental crystal
  geometry, and the template for what a perception benchmark looks like.
- **Open Babel** [3] - `PerceiveBondOrders` (implementing Sayle's
  method); assigns bond orders on PDB round-trip.  No dedicated
  evaluation paper; accuracy inherited from Sayle's rules.  PLIP [4],
  the standard protein-ligand interaction tool, delegates bond-order
  perception to it, so PLIP's chemistry quality is Open Babel's.
- **xyz2mol / rdDetermineBonds** [5, 6] - Kim & Kim's valence
  enumeration: connectivity from covalent radii (scaled), then an
  exhaustive assignment of bond orders subject to valence/charge
  constraints, without fitted parameters.  Validated on 10,000 PubChem
  molecules (near-100% recovery of the bond-order matrix).  Practical
  caveats: needs the molecular charge, and requires hydrogens to be
  reliable - both unavailable for typical PDB ligand coordinates.
  Ported into RDKit as `rdDetermineBonds`, which this work evaluates on
  the crystal benchmark: it recovers **0/600** held-out ligands (0/202
  tuning) - without hydrogens its connect-the-dots connectivity and
  valence assignment fail before the bond-order stage, and its outputs
  cannot be sanitized (13% even on the main tier).  The finding
  mirrors the synthetic benchmark (2% at zero noise) and is structural,
  not tunable: the RDKit route assumes a hydrogen-complete input that
  PDB ligands never have.
- **Indigo** [7] - EPAM's open-source toolkit assigns ground-state bond
  orders and formal charges (energy-based Lewis-structure search; a
  temperature parameter interpolates resonance states).  No published
  accuracy benchmark against experimental coordinates.

## 2. Statistical and optimization methods from coordinates

- **Labute 2005** [8] - "On the Perception of Molecules from 3D Atomic
  Coordinates" (JCIM).  Bond orders assigned by maximum weighted
  matching on a non-bipartite graph, with weights from statistics of a
  large organic-molecule collection; assigns hybridization, bond orders
  and formal charges with or without hydrogens.  Tested on functional
  groups, heterocycles, and PDB + CSD entries.  The closest
  methodological ancestor of bope's geometry tier (statistical
  thresholds, no hydrogen requirement), and the origin of the only
  perception benchmark set the literature reuses (179 complexes, see
  Knodle below).
- **fconv** [9] - bond-order perception inside the PLANTS docking
  pipeline (Korb et al.); empirical distance rules, tuned for docked
  poses.
- **NAOMI** [10] - Urbaczek et al.: molecule perception (connectivity,
  bond orders, atom typing, stereo) for a structure-handling toolkit,
  evaluated as part of the Knodle comparison.
- **Knodle 2016** [11] - Kadukova & Grudinin: the first ML method in
  this space, an SVM trained on PDBbindCN that perceives atom types,
  hybridization and bond orders from 3D coordinates.  The benchmark
  table that defines the field's numbers: on Labute's 179-complex set,
  5-6 perception errors vs NAOMI 7, I-interpret 9, fconv 13; ~3.9%
  errors on 3,000 PDBBind complexes; ~4.5% on 332,974 Ligand Expo
  (CCD) entries.  Caveats for comparison: the Ligand Expo evaluation is
  on ideal reference geometries, not deposited coordinates; the
  PDBBind training set overlaps the evaluation domain (trained on
  PDBBind, evaluated on PDBBind); the error definition (per-structure,
  per-atom) differs from bope's per-ligand exact-recovery metric.

## 3. Machine-learning perception from coordinates (recent)

- **YuelBond 2025** [12] - GNN with edge-level prediction over pairwise
  distances, trained on the GEOM dataset (450k+ computed geometries).
  ~98% accuracy on clean coordinates, ~93% at 0.2 A noise.  Publishes
  its implementation on Bitbucket (dokhlab/yuel_bond) AND its
  pretrained weights (Zenodo record 15353365: geom_3d.ckpt etc.,
  7.7 MB each) - checked Aug 2026.  That made an independent run
  possible: this work ran geom_3d.ckpt head-to-head on the crystal
  benchmark (first evaluation of any ML perception model on
  experimental PDB coordinates against CCD ground truth).  Result:
  37.0% / 27.7% exact CCD recovery (main / low-res) vs its ~98% F1
  on clean GEOM - the GEOM-trained network does not transfer to
  deposited coordinates (its own code also silently maps any element
  outside its 9-atom vocabulary C/O/N/F/S/Cl/Br/I/P to Cl; 12 of the
  1,200 held-out ligands carry B/Fe/Se/W).  GEOM geometries are
  idealized, mostly neutral drug-like molecules - exactly the
  coordinate-quality axis the crystal benchmark splits into
  resolution tiers.
- **Uni-Bond 2026** [13] - Uni-Mol encoder + pairwise classification
  head; ~19x error reduction over baselines on GEOM scaffold split;
  zero-shot transfer to peptides/clusters.  Same evaluation domain
  caveat as YuelBond; no public code or weights found (ICML 2026
  page, checked Aug 2026).
- **CoTAR 2026** [14] - hybrid GNN/HMM reconstructing topology, formal
  charges and unpaired electrons from element types + coordinates;
  evaluated on condensed-phase MD trajectories, not crystals.
- **Knodle 2016** [11] - the only ML method designed for exactly this
  input (PDB-like coordinates, no hydrogens); its published numbers
  (5-6 errors on Labute's 179-complex set, ~3.9% on PDBBind) are the
  field's reference.  Not independently runnable here: the INRIA
  release is a research binary (no public repo, no pip package, not
  Python-3 compatible), so its numbers are cited, not re-measured.
- Position: ML perception is a genuine alternative but is trained and
  evaluated on computed, idealized geometries.  Of the 2025-26
  models, only YuelBond publishes its implementation and weights
  (Zenodo 15353365), and this work runs it head-to-head on
  experimental PDB ligand coordinates against CCD ground truth - the
  first such evaluation of any ML perception model - where it
  recovers 37.0% / 27.7% (main / low-res), 40-48 points below bope's
  geometry tier, with ~40% of its outputs unsanitizable (see
  Section 4.3 of the paper).  Uni-Bond publishes no code or weights;
  CoTAR evaluates on condensed-phase MD trajectories, not crystals;
  Knodle is a research binary (not Python-3 compatible) - their
  numbers are cited, not re-measured.  The trained models also embed
  their training-domain biases (charge assignment, metal handling,
  protonation), which is exactly what a deposit-derived ligand is
  worst at.  bope's geometry tier is parameter-light, deterministic
  and explainable - properties a publication reader can audit without
  a training corpus.

## 4. Stereo perception from 3D

- **RDKit `AssignStereochemistryFrom3D`** [15] - assigns R/S from 3D
  coordinates using the CIP rules (approximate: no explicit ranking of
  non-H substituents beyond CIP ordering); the standard programmatic
  route.  bope's stereo layer (and its stereo benchmark) use it.
- **Open Babel** [3] - assigns its own tetrahedral and E/Z stereo
  during SDF output from coordinates, without an explicit
  user-facing call.
- **Indigo** [7] - geometry-based cis/trans perception from 3D
  coordinates (`isGeomStereoBond` etc.).
- No systematic benchmark of any of these against CCD stereo ground
  truth on PDB ligands exists in the literature (see section 5: the
  nearest works validate the deposit, not the perception).

## 5. Ground truth and prior evaluations of PDB ligand chemistry

- **Chemical Component Dictionary** [16] - wwPDB's authoritative
  per-HET-code chemistry: formula, bond graph, stereo (SMILES /
  SMILES_stereo).  Successor of Ligand Depot [17].  This is the ground
  truth bope's benchmarks use.
- **ValidatorDB 2015** [18] - "validation of annotation": compares
  every deposited ligand model against its CCD reference; classifies
  chirality errors (C, metal, planar, high-order, other) and
  completeness.  At the 2014 snapshot: ~83% of validated molecules
  complete with correct chirality, <8% with chirality problems.
  Validates the deposit, not perception algorithms - but its
  error-class taxonomy is a model for reporting stereo failure modes,
  and its 83% figure bounds what deposited coordinates can support.
- **Mogul** [19] - CCDC's geometry validation: CSD-derived bond-length,
  angle and torsion distributions; flags cis-trans inconsistencies and
  strained geometry.  Validates geometry against experiment; requires
  the chemistry to be given, it does not perceive it.
- **Waibl, Liedl & Rupp 2022** [20] - documents cis-trans errors in
  deposited PDB ligands: of 24,743 HET monomers, 187 cis-only, 286
  trans-only, 46 with both isomers (case study: palmitoleic acid ~70%
  of deposits in the wrong configuration).  Direct evidence that the
  deposited ground truth is itself imperfect - the CCD stereo is the
  safer reference, and a perception benchmark against it is what bope
  provides.
- **RCSB validation reports** [21] - per-entry ligand quality scores
  (clashes, geometry, chirality) from wwPDB validation.  Aggregate
  quality statistics exist; algorithm-vs-CCD perception accuracy does
  not.

## 6. The gap bope fills

Existing perception evaluations, summarized:

| work | input | ground truth | scale | held out? | reported accuracy (their metric) |
|---|---|---|---|---|---|
| Baber & Hodgkin 1992 | CSD coords | hand-curated | 91 | no (method-developed on it) | ~82% connectivity |
| Labute 2005 | PDB + CSD coords | hand-curated | 179 complexes | no; set reused by others | no accuracy number (method paper) |
| Knodle 2016 | PDBBind coords | hand/PDB annotation | 3,000 | no - trained on PDBBind, tested on PDBBind | ~3.9% errors (PDBBind), 5-6 errors on 179-complex set |
| xyz2mol 2015 | ideal coords + H + charge | PubChem | 10,000 | no | near-100% bond-order-matrix recovery |
| YuelBond 2025 | computed coords | computed geometry | 450k | scaffold split (GEOM) | ~98% F1 (clean), 92.7% (0.2 A noise); **37.0% / 27.7% exact CCD recovery on PDB coords (run here)** |
| Uni-Bond 2026 | computed coords | computed geometry | GEOM | scaffold split (GEOM) | ~19x error reduction vs baselines; no code or weights |
| rdDetermineBonds (run here) | PDB coords, no H | CCD | 202 + 1,200 | same held-out protocol - 0/1,402 recovered | 0% exact recovery |
| **bope (this work)** | **PDB coords, no H, no charge** | **CCD (formula + graph + stereo)** | **202 tuning + 1,200 held-out (2 generations), 2 resolution tiers** | **yes - disjoint PDB ids and HET codes, 5 samples/tier** | **77.0% / 75.7% exact CCD recovery (main / low-res)** |

What no prior work does:

1. **Benchmark bond-order perception against CCD ground truth on a
   large, chemically diverse sample of experimental PDB coordinates**,
   split into resolution tiers (1.0-2.0 A vs 2.5-3.0 A) so coordinate
   noise is an explicit axis.  Knodle's PDBBind evaluation is the only
   comparably real-coordinate result, and its metric, training-set
   overlap and non-public data limit comparability.
2. **Benchmark stereo perception from coordinates at all**.  The
   closest works validate deposits (ValidatorDB, wwPDB reports) or
   document deposit errors (Waibl); nobody measures how well an
   algorithm recovers the CCD's declared stereo from coordinates,
   per-center and atom-mapped.
3. **Evaluate a template-free, parameter-light method with a
   published, reproducible protocol**: held-out sampling disjoint by
   entry and HET code, code-freeze discipline, per-bucket mean +/- std,
   environment capture in every results file, and the full pipeline
   (perception + benchmarks + datasets) in one MIT-licensed repo.

## References

1. Sayle, R. Cruft to content: perception of bond order in organic
   molecules. Daylight European MUG Meeting 2001 (unrefereed).
2. Baber, J. C.; Hodgkin, E. E. Automatic assignment of chemical
   connectivity to organic molecules in the Cambridge Structural
   Database. J. Chem. Inf. Comput. Sci. 1992, 32(5), 401-406.
   DOI 10.1021/ci00009a001.
3. O'Boyle, N. M.; Banck, M.; James, C. A.; Morley, C.; Vandermeersch,
   T.; Hutchison, G. R. Open Babel: an open chemical toolbox.
   J. Cheminform. 2011, 3, 33. DOI 10.1186/1758-2946-3-33.
4. Salentin, S.; Schreiber, S.; Haupt, V. J.; Adasme, M. F.; Schroeder,
   M. PLIP: fully automated protein-ligand interaction profiler.
   Nucleic Acids Res. 2015, 43(W1), W443-W447. DOI 10.1093/nar/gkv315.
5. Kim, Y.; Kim, W. Y. Universal structure conversion method for
   organic molecules: from atomic connectivity to three-dimensional
   geometry. Bull. Korean Chem. Soc. 2015, 36(7), 1769-1777.
   DOI 10.1002/bkcs.10334.
6. RDKit `rdDetermineBonds` (GSOC 2022 port of xyz2mol);
   https://www.rdkit.org/docs/source/rdkit.Chem.rdDetermineBonds.html
7. Indigo toolkit, EPAM; https://lifescience.opensource.epam.com/indigo/
8. Labute, P. On the perception of molecules from 3D atomic
   coordinates. J. Chem. Inf. Model. 2005, 45(2), 215-221.
   DOI 10.1021/ci049915d.
9. Korb, O.; Stutzle, T.; Exner, T. E. Empirical scoring functions for
   advanced protein-ligand docking with PLANTS. J. Chem. Inf. Model.
   2009, 49(1), 84-96. DOI 10.1021/ci800298z.
10. Urbaczek, J.; Kolodzik, A.; Fischer, J. R.; Lippert, T.; Heuser,
    S.; Groth, I.; Rarey, M. NAOMI: on the almost-trivial task of
    reading molecules from different file formats. J. Chem. Inf. Model.
    2011, 51(12), 3199-3207. DOI 10.1021/ci200324e.
11. Kadukova, M.; Grudinin, S. Knodle: a support vector machines-based
    automatic perception of organic molecules from 3D coordinates.
    J. Chem. Inf. Model. 2016, 56(8), 1410-1419.
    DOI 10.1021/acs.jcim.5b00512.
12. Wang, J.; Dokholyan, N. V. Multimodal Bond Reconstruction toward
    Generative Molecular Design. J. Chem. Inf. Model. 2026, 66(2),
    1003-1012. DOI 10.1021/acs.jcim.5c03052 (preprint: bioRxiv
    2025.05.06.652517; code: bitbucket.org/dokhlab/yuel_bond).
13. Pavlenko, S.; Maslov, P.; Alexandrovich, I. B.; Tsypin, A.;
    Telepov, A.; Ushenin, K.; Khrabrov, K.; Potapov, D.; Kadurin, A.
    Uni-Bond: Learning Chemical Bonds from Atomic Coordinates.
    In Proceedings of the 43rd International Conference on Machine
    Learning (ICML 2026), Seoul, Republic of Korea, 2026.
14. Mori, H.; Kikkawa, T.; Miyazaki, Y. CoTAR: Topology and Atomic
    State Reconstruction in Condensed Phases. arXiv:2606.27636, 2026.
15. RDKit `AssignStereochemistryFrom3D`;
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html
16. Westbrook, J. D.; Shao, C.; Feng, Z.; Zhuravleva, M.; Velankar, S.;
    Young, J. The chemical component dictionary: complete descriptions
    of constituent molecules in experimentally determined 3D
    macromolecules in the Protein Data Bank. Bioinformatics 2015,
    31(8), 1274-1278. DOI 10.1093/bioinformatics/btu789.
17. Feng, Z.; Chen, L.; Maddula, H.; Akcan, O.; Oughtred, R.; Berman,
    H. M.; Westbrook, J. Ligand Depot: a data warehouse for ligands
    bound to macromolecules. Bioinformatics 2004, 20(13), 2153-2155.
18. Sehnal, D.; Svobodova Varekova, R.; Pravda, L.; Ionescu, C.-M.;
    Geidl, S.; Horsky, V.; Jaiswal, D.; Wimmerova, M.; Koca, J.
    ValidatorDB: database of up-to-date validation results for ligands
    and non-standard residues from the Protein Data Bank.
    Nucleic Acids Res. 2015, 43(D1), D369-D375. DOI 10.1093/nar/gku1118.
19. Bruno, I. J.; Cole, J. C.; Kessler, M.; Luo, J.; Motherwell, W. D.
    S.; Purkis, L. H.; Smith, B. R.; Taylor, R.; Cooper, R. I.; Harris,
    S. E.; Orpen, A. G. Retrieval of crystallographically-derived
    molecular geometry information. J. Chem. Inf. Comput. Sci. 2004,
    44(6), 2133-2144. DOI 10.1021/ci049780b.
20. Waibl, F.; Liedl, K. R.; Rupp, B. Correcting cis-trans-
    transgressions in macromolecular structure models. FEBS J. 2022,
    289(10), 2793-2804. DOI 10.1111/febs.15884.
21. wwPDB validation reports; https://www.wwpdb.org/validation
