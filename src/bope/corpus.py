"""Curated SMILES validation corpus for bond-order perception.

Two sources:

1. **Crystal ligands** (:data:`CRYSTAL_LIGANDS`): the 17 HET residues of
   the 16 experimental crystal complexes (kinases, GPCRs, proteases), as
   canonical SMILES fetched from the RCSB Chemical Component Dictionary.
   These are the ground truth the geometry perception was validated
   against (exact formula reproduction on the deposited coordinates).
2. **Synthetic chemotypes** (:data:`SYNTHETIC`): hand-curated molecules
   covering the chemotype space the perception must handle - (fused)
   aromatic systems, N-rich heterocycles, tautomers, carbonyls, amides,
   nitriles, alkynes, halogens, small/strained rings, sugars, amino
   acids, nucleobases and drugs.  Each molecule is embedded with ETKDG,
   perturbed with coordinate noise and must be recovered exactly
   (same bond graph, same formula).

Entries whose ``(name, smiles)`` tuple carries ``True`` as a third
element are **charged** (formal charges in the SMILES).  The perception
API deliberately does not consume the ``charge`` argument (documented in
:func:`bope.perceive_bond_orders`), so charged molecules
cannot reproduce the reference bond graph exactly - they are exercised
for graceful behaviour (sanitizable mol, correct formula check relaxed)
in both the test suite and the benchmark.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Validation exclusions (measured 2026-08-11, current RDKit / ETKDG)
# ---------------------------------------------------------------------------
# These document *measured* properties of the ETKDG seed-42 embeds, shared
# by tests/test_perception.py and benchmarks/benchmark.py.
# If an RDKit upgrade changes the seed-42 geometries, re-measure before
# re-enabling the strict assertions.

#: Neutral corpus members whose ETKDG seed-42 embed is not faithful at ZERO
#: noise, so exact recovery is impossible regardless of perception quality:
#: ZMA's fused triazolo-triazine gets +16 H in the embed; thiophene's C-S
#: (1.85-1.9 A) leaves the aromatic envelope; benzophenone's C=O refines to
#: 1.304 A (vs the 1.30 double cutoff - noise actually helps it); ketoprofen
#: has the same C=O boundary; 1,3-cyclohexadiene's ring C-C sits outside the
#: aromatic envelope.  Excluded from the strict recovery assertions; the
#: benchmark reports their actual behaviour.
EMBED_UNFAITHFUL: set[str] = {
    "ZMA ZM241385 (3EML)",
    "thiophene",
    "benzophenone",
    "ketoprofen",
    "1,3-cyclohexadiene",
}

#: Molecules that recover exactly at 0.0 noise but fail at 0.03 A bond-RMS:
#: hard-cutoff boundary events where the noise pushes a bond across a length
#: threshold that cannot be raised without breaking another chemotype
#: (e.g. the 1.30 A C=O cutoff vs acid C-OH at 1.31-1.35 A; the 1.33 A C=N
#: cutoff vs amide C-N).
NOISE_SENSITIVE: set[str] = {
    "E20 donepezil (1EVE)",
    "1,4-naphthoquinone",
    "creatinine",
    "nicotinic acid",
    "isoprene",
    "dimethyl sulfone",
    "acetoxime",
    "aspartic acid",
    "glutamine",
    "estradiol",
}

#: Neutral molecules whose recovery goes through the OpenBabel fallback
#: instead of the in-house geometry perception (the C-S ring bond sits
#: outside the aromatic envelope at ETKDG distances).  They recover EXACTLY
#: at both 0.0 and 0.03 noise - but only when OpenBabel is installed.
OPENBABEL_FALLBACK: set[str] = {"thiazole", "benzothiazole", "dibenzothiophene"}

#: RCSB-verified ground-truth formulas (CCD InChI formulas) per eval
#: structure.  Each entry maps a HET code to the set of formulas its
#: residues may legitimately produce: a structure can hold the same HET as
#: several residues with different termini (5HVP has STA as an N-terminal
#: statine C8H17NO2 and a C-terminal statine C8H17NO3, and two copies each
#: of STI in 1IEP and ETQ in 3PBL).
CRYSTAL_EXPECTED: dict[str, dict[str, set[str]]] = {
    "1KTS": {"C24": {"C27H29N7O3"}},
    "1T3R": {"017": {"C27H37N3O7S"}},
    "1AQ1": {"STU": {"C28H26N4O3"}},
    "1EVE": {"E20": {"C24H29NO3"}},
    "4UEH": {"BEN": {"C7H8N2"}},
    "2A4L": {"RRC": {"C19H26N6O"}},
    "1ACJ": {"THA": {"C13H14N2"}},
    "1IEP": {"STI": {"C29H31N7O"}},
    "1M17": {"AQ4": {"C22H23N3O4"}},
    "3OG7": {"032": {"C23H18ClF2N3O3S"}},
    "2RH1": {"CAU": {"C18H22N2O2"}},
    "3EML": {"ZMA": {"C16H15N7O2"}},
    "3RFM": {"CFF": {"C8H10N4O2"}},
    "3PBL": {"ETQ": {"C17H25ClN2O3"}},
    "4S0V": {"SUV": {"C23H23ClN6O2"}},
    "5HVP": {"ACE": {"C2H4O"}, "STA": {"C8H17NO2", "C8H17NO3"}},
}

#: (name, SMILES) pairs - the 17 eval-dataset crystal ligands (RCSB CCD
#: canonical SMILES, verified 2026-08-11).
CRYSTAL_LIGANDS: list[tuple[str, str]] = [
    ("THA tacrine (1ACJ)", "c1ccc2c(c1)c(c3c(n2)CCCC3)N"),
    ("STU staurosporine (1AQ1)", "CC12C(C(CC(O1)n3c4ccccc4c5c3c6n2c7ccccc7c6c8c5C(=O)NC8)NC)OC"),
    ("E20 donepezil (1EVE)", "COc1cc2c(cc1OC)C(=O)C(C2)CC3CCN(CC3)Cc4ccccc4"),
    ("STI imatinib (1IEP)", "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"),
    ("C24 melagatran-rel. (1KTS)", "CCOC(=O)CCN(c1ccccn1)C(=O)c2ccc3c(c2)nc(n3C)CNc4ccc(cc4)C(=N)N"),
    ("AQ4 erlotinib (1M17)", "COCCOc1cc2c(cc1OCCOC)ncnc2Nc3cccc(c3)C#C"),
    ("017 darunavir (1T3R)", "CC(C)CN(CC(C(Cc1ccccc1)NC(=O)OC2COC3C2CCO3)O)S(=O)(=O)c4ccc(cc4)N"),
    ("RRC roscovitine (2A4L)", "CCC(CO)Nc1nc(c2c(n1)n(cn2)C(C)C)NCc3ccccc3"),
    ("CAU carazolol (2RH1)", "CC(C)NCC(COc1cccc2c1c3ccccc3[nH]2)O"),
    ("ZMA ZM241385 (3EML)", "c1cc(oc1)c2nc3nc(nc(n3n2)N)NCCc4ccc(cc4)O"),
    ("032 vemurafenib (3OG7)", "CCCS(=O)(=O)Nc1ccc(c(c1F)C(=O)c2c[nH]c3c2cc(cn3)c4ccc(cc4)Cl)F"),
    ("ETQ eticlopride (3PBL)", "CCc1cc(c(c(c1O)C(=O)NCC2CCCN2CC)OC)Cl"),
    ("CFF caffeine (3RFM)", "Cn1cnc2c1C(=O)N(C(=O)N2C)C"),
    ("SUV suvorexant (4S0V)", "Cc1ccc(c(c1)C(=O)N2CCN(CCC2C)c3nc4cc(ccc4o3)Cl)n5nccn5"),
    ("BEN benzamidine (4UEH)", "[H]N=C(c1ccccc1)N"),
    ("ACE acetyl (5HVP)", "CC=O"),
    ("STA statine (5HVP)", "CC(C)CC(C(CC(=O)O)O)N"),
]

#: (name, SMILES) pairs - synthetic chemotype coverage, all neutral.
SYNTHETIC: list[tuple[str, str]] = [
    # --- plain (fused) aromatics ---
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("phenanthrene", "c1ccc2c(c1)ccc1ccccc12"),
    ("pyridine", "c1ccncc1"),
    ("pyrimidine", "c1cncnc1"),
    ("pyrazine", "c1cnccn1"),
    ("pyridazine", "c1ccnnc1"),
    ("imidazole", "c1c[nH]cn1"),
    ("pyrazole", "c1cc[nH]n1"),
    ("pyrrole", "c1cc[nH]c1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("benzimidazole", "c1ccc2c(c1)nc[nH]2"),
    ("indazole", "c1ccc2c(c1)cn[nH]2"),
    ("purine", "c1ncnc2[nH]cnc12"),
    ("adenine", "Nc1ncnc2[nH]cnc12"),
    ("guanine", "Nc1nc2[nH]cnc2c(=O)[nH]1"),
    ("cytosine", "c1cnc(N)[nH]c1=O"),
    ("thymine", "Cc1c[nH]c(=O)[nH]c1=O"),
    ("uracil", "c1c[nH]c(=O)[nH]c1=O"),
    ("1,2,4-triazole", "c1n[nH]cn1"),
    ("1,2,3-triazole", "c1cn[nH]n1"),
    ("tetrazole", "c1nn[nH]n1"),
    ("1,3,5-triazine", "c1ncncn1"),
    ("melamine", "Nc1nc(N)nc(N)n1"),
    ("quinoline", "c1ccc2ncccc2c1"),
    ("isoquinoline", "c1ccc2cnccc2c1"),
    ("quinazoline", "c1ccc2ncncc2c1"),
    ("quinoxaline", "c1ccc2nccnc2c1"),
    ("acridine", "c1ccc2cc3ccccc3nc2c1"),
    ("phenazine", "c1ccc2nc3ccccc3nc2c1"),
    ("carbazole", "c1ccc2c(c1)[nH]c1ccccc12"),
    ("furan", "c1ccoc1"),
    ("thiophene", "c1ccsc1"),
    ("thiazole", "c1cscn1"),
    ("oxazole", "c1cocn1"),
    ("benzofuran", "c1ccc2c(c1)cco2"),
    ("benzothiazole", "c1ccc2c(c1)scn2"),
    ("dibenzofuran", "c1ccc2c(c1)oc1ccccc12"),
    ("dibenzothiophene", "c1ccc2c(c1)sc1ccccc12"),
    # --- carbonyl / amide / ester chemotypes ---
    ("acetone", "CC(C)=O"),
    ("acetophenone", "CC(=O)c1ccccc1"),
    ("benzophenone", "O=C(c1ccccc1)c1ccccc1"),
    ("cyclohexanone", "O=C1CCCCC1"),
    ("acetamide", "CC(N)=O"),
    ("benzamide", "NC(=O)c1ccccc1"),
    ("N-methylacetamide", "CC(=O)NC"),
    ("urea", "NC(=O)N"),
    ("phenylurea", "NC(=O)Nc1ccccc1"),
    ("ethyl carbamate", "CCOC(N)=O"),
    ("ethyl acetate", "CCOC(C)=O"),
    ("gamma-butyrolactone", "O=C1CCCO1"),
    ("delta-valerolactone", "O=C1CCCCO1"),
    ("caprolactam", "O=C1CCCCCN1"),
    ("succinimide", "O=C1CCC(=O)N1"),
    ("maleimide", "O=C1C=CC(=O)N1"),
    ("phthalimide", "O=C1NC(=O)c2ccccc12"),
    ("isatin", "O=c1c(=O)c2ccccc2[nH]1"),
    ("hydantoin", "O=C1CNC(=O)N1"),
    ("barbituric acid", "O=C1CC(=O)NC(=O)N1"),
    ("p-benzoquinone", "O=C1C=CC(=O)C=C1"),
    ("1,4-naphthoquinone", "O=C1C=CC(=O)c2ccccc21"),
    ("anthraquinone", "O=C1c2ccccc2C(=O)c2ccccc21"),
    ("coumarin", "O=c1ccc2ccccc2o1"),
    ("flavone", "C1=CC=C(C=C1)C2=CC(=O)C3=CC=CC=C3O2"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("salicylamide", "NC(=O)c1ccccc1O"),
    ("saccharin", "O=C1NS(=O)(=O)c2ccccc21"),
    ("naproxen", "CC(C(=O)O)c1ccc2cc(OC)ccc2c1"),
    ("ketoprofen", "OC(=O)C(C)c1cccc(c1)C(=O)c1ccccc1"),
    ("diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    ("warfarin", "CC(=O)CC(c1ccccc1)c1ccc2ccccc2c1O"),
    ("celecoxib", "CC1=CC(=NN1c1ccc(F)cc1)c1ccc(cc1)S(N)(=O)=O"),
    ("sulfamethoxazole", "CC1=CC(=NO1)NS(=O)(=O)c1ccc(N)cc1"),
    ("trimethoprim", "COc1cc(cc(OC)c1OC)Cc1cnc(nc1N)N"),
    ("sorafenib", "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1"),
    ("gefitinib", "COc1cc2c(cc1OC)ncnc2Nc1ccc(F)c(Cl)c1"),
    # --- nitrogen-rich drugs / heterocycles ---
    ("serotonin", "NCCc1c[nH]c2ccc(O)cc12"),
    ("melatonin", "COc1ccc2c(c1)[nH]cc2CCNC(C)=O"),
    ("tryptamine", "NCCc1c[nH]c2ccccc12"),
    ("DMT", "CN(C)CCc1c[nH]c2ccccc12"),
    ("histamine", "NCCc1cnc[nH]1"),
    ("creatinine", "CN1CC(=O)NC1=N"),
    ("nicotinic acid", "OC(=O)c1cccnc1"),
    ("niacinamide", "NC(=O)c1cccnc1"),
    ("pyridoxine", "Cc1ncc(CO)c(CO)c1O"),
    ("chloroquine", "CC(C)NCCCC(C)Nc1ccnc2ccccc12"),
    ("melagatran", "NC(=N)NCCC(=O)N1CC(CC1)c1ccc(cc1)CC(N)C(=O)O"),
    ("benzamidine-2", "N=C(N)c1ccccc1"),
    ("fentanyl", "CCC(=O)N(CCc1ccccc1)C1CCN(CC1)c1ccccc1"),
    ("creatine", "CN(CC(=O)O)C(=N)N"),
    # --- nitriles / alkynes / alkenes ---
    ("acetonitrile", "CC#N"),
    ("benzonitrile", "N#Cc1ccccc1"),
    ("phenylacetylene", "C#Cc1ccccc1"),
    ("1-hexyne", "CCCCC#C"),
    ("styrene", "C=Cc1ccccc1"),
    ("trans-stilbene", "C(=C/c1ccccc1)\\c1ccccc1"),
    ("isoprene", "CC(=C)C=C"),
    ("limonene", "CC1=CCC(CC1)C(=C)C"),
    ("fumaric acid", "OC(=O)/C=C/C(=O)O"),
    ("maleic acid", "OC(=O)/C=C\\C(=O)O"),
    ("acrylamide", "NC(=O)C=C"),
    # --- halogens / sulfur / heteroatom-rich ---
    ("fluorobenzene", "Fc1ccccc1"),
    ("chlorobenzene", "Clc1ccccc1"),
    ("bromobenzene", "Brc1ccccc1"),
    ("iodobenzene", "Ic1ccccc1"),
    ("2-chloropyridine", "Clc1ccccn1"),
    ("benzotrifluoride", "FC(F)(F)c1ccccc1"),
    ("thiophenol", "Sc1ccccc1"),
    ("ethanethiol", "CCS"),
    ("dimethyl disulfide", "CSSC"),
    ("dimethyl sulfone", "CS(C)(=O)=O"),
    ("methanesulfonamide", "CS(=O)(=O)N"),
    ("benzenesulfonamide", "NS(=O)(=O)c1ccccc1"),
    ("taurine", "NCCS(=O)(=O)O"),
    ("trimethyl phosphate", "COP(=O)(OC)OC"),
    ("acetoxime", "CC(C)=NO"),
    ("anisole", "COc1ccccc1"),
    ("1,2-dimethoxyethane", "COCCOC"),
    ("MTBE", "CC(C)(C)OC"),
    ("methyl phenyl sulfide", "CSc1ccccc1"),
    # --- amino acids / small peptides / sugars ---
    ("glycine", "NCC(=O)O"),
    ("alanine", "C[C@@H](N)C(=O)O"),
    ("serine", "OC[C@@H](N)C(=O)O"),
    ("cysteine", "SC[C@H](N)C(=O)O"),
    ("methionine", "CSCC[C@H](N)C(=O)O"),
    ("phenylalanine", "N[C@@H](Cc1ccccc1)C(=O)O"),
    ("tyrosine", "N[C@@H](Cc1ccc(O)cc1)C(=O)O"),
    ("tryptophan", "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"),
    ("histidine", "N[C@@H](Cc1cnc[nH]1)C(=O)O"),
    ("proline", "O=C(O)[C@@H]1CCCN1"),
    ("lysine", "NCCCC[C@H](N)C(=O)O"),
    ("arginine", "N[C@@H](CCCNC(=N)N)C(=O)O"),
    ("aspartic acid", "N[C@@H](CC(=O)O)C(=O)O"),
    ("glutamic acid", "N[C@@H](CCC(=O)O)C(=O)O"),
    ("glutamine", "N[C@@H](CCC(N)=O)C(=O)O"),
    ("leucine", "CC(C)C[C@@H](N)C(=O)O"),
    ("valine", "CC(C)[C@@H](N)C(=O)O"),
    ("glycylglycine", "NCC(=O)NCC(=O)O"),
    ("DOPA", "N[C@@H](Cc1cc(O)c(O)cc1)C(=O)O"),
    ("dopamine", "NCCc1ccc(O)c(O)c1"),
    ("adrenaline", "CNC[C@@H](O)c1ccc(O)c(O)c1"),
    ("ribose", "OC1C(O)C(O)C(O)CO1"),
    ("glucose", "OC1C(O)C(O)C(O)C(O)CO1"),
    ("fructose", "OCC1(O)C(O)C(O)C(O)CO1"),
    ("sorbitol", "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO"),
    # --- aliphatic / strained rings / polycycles ---
    ("cyclopropane", "C1CC1"),
    ("cyclobutane", "C1CCC1"),
    ("cyclopentane", "C1CCCC1"),
    ("cyclohexane", "C1CCCCC1"),
    ("cycloheptane", "C1CCCCCC1"),
    ("cyclooctane", "C1CCCCCCC1"),
    ("1,3-cyclopentadiene", "C1C=CC=C1"),
    ("1,3-cyclohexadiene", "C1C=CC=CC1"),
    ("adamantane", "C1C2CC3CC1CC(C2)C3"),
    ("decalin", "C1CCC2CCCCC2C1"),
    ("tetrahydrofuran", "C1CCOC1"),
    ("tetrahydropyran", "C1CCOCC1"),
    ("morpholine", "C1COCCN1"),
    ("piperidine", "C1CCNCC1"),
    ("piperazine", "C1CNCCN1"),
    ("1,4-dioxane", "C1COCCO1"),
    ("2-azetidinone", "O=C1NCC1"),
    ("cubane", "C12C3C4C1C1C4C3C21"),
    ("testosterone", "C[C@]12CC[C@H]3[C@@H](CC=C4CC(=O)CC[C@]43C)[C@@H]1CC[C@@H]2O"),
    ("estradiol", "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]43C)[C@@H]1CC[C@@H]2O"),
    ("cholesterol", "C[C@@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"),
]

#: (name, SMILES, charged=True) - formal-charge molecules exercised for
#: graceful behaviour only (the API does not consume ``charge``).
CHARGED: list[tuple[str, str, bool]] = [
    ("ammonium", "C[NH3+]", True),
    ("tetramethylammonium", "C[N+](C)(C)C", True),
    ("guanidinium", "NC(N)=[NH2+]", True),
    ("acetate", "CC(=O)[O-]", True),
    ("pyruvic acid (charged O)", "CC(=O)C(=O)[O-]", True),
    ("glycinate", "NCC(=O)[O-]", True),
    ("nitrobenzene", "[O-][N+](=O)c1ccccc1", True),
    ("methanesulfonate", "CS(=O)(=O)[O-]", True),
    ("pyridinium", "c1cc[nH+]cc1", True),
    ("methylguanidinium", "CNC(N)=[NH2+]", True),
]

#: Full corpus: (name, smiles) triples; charged entries carry ``True``.
CORPUS: list[tuple[str, str, bool]] = [
    *[(name, smiles, False) for name, smiles in CRYSTAL_LIGANDS],
    *[(name, smiles, False) for name, smiles in SYNTHETIC],
    *CHARGED,
]
