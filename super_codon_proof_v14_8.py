#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SUPER CODON CONDITIONAL SYNTAX PROOF PIPELINE — AUTONOMOUS v14.6 (audit-clean + Intron/Exon & P53 Geometric Audit)
==================================================================================================================

v14.6 = v14.5 + Autonomous Intron-vs-Exon, Splice-Junction & P53 Geometric/Program Audit:
  * ADDED --p53-audit and --intron-exon-audit CLI modes (runs standalone by default or
    alongside smoke/full with --with-smoke / --with-full).
  * ADDED autonomous NCBI RefSeqGene GenBank fetcher & parser (parse_genbank_intron_exon)
    for TP53 (NG_017013.2), TP63, TP73, and 12 benchmark human genes (--intron-exon-genes),
    plus offline custom GenBank/FASTA support (--custom-gb).
  * ADDED 4-state structural comparison per gene:
      1) final_exons_only (spliced CDS in true frame)
      2) full_gene_with_introns_f0 (unspliced pre-mRNA genomic span, frame 0)
      3) full_gene_with_introns_3f (unspliced pre-mRNA evaluated across all 3 frames)
      4) introns_only (excised introns)
  * ADDED Geometric Form vs. Informational Program metrics (32-mirror-pair Chargaff L1
    asymmetry, Markov flow irreversibility ||P - P^T||_F / ||P||_F, entropy, stop density).
  * ADDED Splice-Junction Phase Reconstruction test (comparing spliced Exon[k]->Exon[k+1]
    codon transitions and phase-1/2 split codons against unspliced Exon->Intron boundaries).
  * ADDED Co-Transcriptional Progressive Splicing Trajectory for TP53 (step-by-step
    removal of coding introns from 0% spliced to 100% final exons only).

All v14.1–v14.5 scientific & engineering upgrades preserved:
  * FIXED copy-paste Markdown corruption in ftp_base_for_accession() and eutils_efetch_url().
  * RESTORED --h5-bootstrap-inner-reps (default=2 in full mode).
  * EMBEDDED automatic startup self-test (self_test_invariants() + --self-test CLI flag)
    verifying genetic code tables 1, 2, 3, 4, 5, 6, 11, 12, 32 mirror pairs, FTP URLs,
    and synthetic intron-exon splicing invariants.
  * DNA alphabet fix in TABLE_OVERRIDES ("ATA": "M" in tables 2, 3, 5).
  * Exact NCBI Table 6 (TAA/TAG -> Q, TGA remains *) and Table 12 (CTG -> S).
  * Curated 41-taxon extended biosphere panel (--include-extended, silent skip on 404).
  * Formal 95% bootstrap confidence intervals (CI) for H6 and H7.
  * Automatic coarse_clean_0p02 (d2b > 0.02) and H5-NEAR-MISS (|gap| <= 0.012) tiers.
  * Per-record transl_table consistency filtering in cds_units().
  * CLI accession overrides registered in SPECIES_PREFIX; longest-prefix matching.
  * Smoke config isolated from full CLI tuning flags; atomic JSON/MD/README writes.
  * Within-species and cross-species small-sample guards.
"""

import argparse
import gzip
import json
import math
import os
import re
import sys
import time
import traceback
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

try:
    import numpy as np
except Exception as e:  # pragma: no cover
    sys.stderr.write(
        "\n[FATAL] NumPy is required but not importable: %r\n"
        "        Install it, e.g.:  python3 -m pip install numpy\n\n" % (e,)
    )
    sys.exit(2)


UA = "super-codon-autonomous/14.6"
DEFAULT_SEED = 20261011


# ===========================================================================
# Logging tee (console + master.log + per-stage log simultaneously)
# ===========================================================================

class MultiTee:
    def __init__(self, *streams):
        self.streams = [s for s in streams if s is not None]

    def write(self, data):
        for s in self.streams:
            try:
                s.write(data)
                s.flush()
            except Exception:
                pass

    def flush(self):
        for s in self.streams:
            try:
                s.flush()
            except Exception:
                pass

    def isatty(self):
        try:
            return self.streams[0].isatty()
        except Exception:
            return False


class TeeContext:
    def __init__(self, *extra_files):
        self.extra = [f for f in extra_files if f is not None]
        self._saved_out = None
        self._saved_err = None

    def __enter__(self):
        self._saved_out = sys.stdout
        self._saved_err = sys.stderr
        sys.stdout = MultiTee(self._saved_out, *self.extra)
        sys.stderr = MultiTee(self._saved_err, *self.extra)
        return self

    def __exit__(self, exc_type, exc, tb):
        sys.stdout = self._saved_out
        sys.stderr = self._saved_err
        return False


# ===========================================================================
# Constants
# ===========================================================================

BASES = "ACGT"
NC, NA = 64, 21

NT_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(BASES):
    NT_LUT[ord(_b)] = _i

CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
assert len(CODONS) == 64

AA_TABLE_STD = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_STRING = "ACDEFGHIKLMNPQRSTVWY*"
AA_INDEX = {ch: i for i, ch in enumerate(AA_STRING)}
STANDARD_AA_LIST = [AA_TABLE_STD[c] for c in CODONS]

# Exact NCBI translation table overrides on DNA alphabet:
#   Table 2: Vertebrate Mitochondrial (AGA=*, AGG=*, ATA=M, TGA=W)
#   Table 3: Yeast Mitochondrial (CTN=T, ATA=M, TGA=W)
#   Table 4: Mold, Protozoan, and Coelenterate Mitochondrial + Mycoplasma (TGA=W)
#   Table 5: Invertebrate Mitochondrial (ATA=M, TGA=W, AGA=S, AGG=S)
#   Table 6: Ciliate, Dasycladacean and Hexamita Nuclear (TAA=Q, TAG=Q; TGA stays *)
#   Table 11: Bacterial, Archaeal and Plant Plastid (same AA map as Table 1)
#   Table 12: Alternative Yeast Nuclear (CTG=S)
TABLE_OVERRIDES = {
    1: {},
    11: {},
    2: {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},
    3: {"CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "ATA": "M", "TGA": "W"},
    4: {"TGA": "W"},
    5: {"ATA": "M", "TGA": "W", "AGA": "S", "AGG": "S"},
    6: {"TAA": "Q", "TAG": "Q"},
    12: {"CTG": "S"},
}
SUPPORTED_TABLES = set(TABLE_OVERRIDES.keys())

CODON_AA = np.array([AA_INDEX[ch] for ch in STANDARD_AA_LIST], dtype=np.int64)
AA_IND = np.zeros((NC, NA), dtype=np.float64)
AA_IND[np.arange(NC), CODON_AA] = 1.0
CHANNEL_CODONS = [np.where(AA_IND[:, a] > 0)[0].astype(np.int64) for a in range(NA)]

C_N1 = np.arange(NC) // 16
C_N2 = (np.arange(NC) // 4) % 4
C_N3 = np.arange(NC) % 4
IS_GC3 = (C_N3 == 1) | (C_N3 == 2)

ACID_AA_IDX = np.array([AA_INDEX["D"], AA_INDEX["E"]], dtype=np.int64)
BASIC_AA_IDX = np.array([AA_INDEX["K"], AA_INDEX["R"]], dtype=np.int64)
NONSTOP_AA_IDX = np.array([i for i, ch in enumerate(AA_STRING) if ch != "*"], dtype=np.int64)

ACID_MASK = np.zeros(NA, dtype=bool)
ACID_MASK[ACID_AA_IDX] = True
BASIC_MASK = np.zeros(NA, dtype=bool)
BASIC_MASK[BASIC_AA_IDX] = True
NONSTOP_MASK = np.zeros(NA, dtype=bool)
NONSTOP_MASK[NONSTOP_AA_IDX] = True

_COMP = {0: 3, 1: 2, 2: 1, 3: 0}
_RC_PARTNER = np.zeros(NC, dtype=np.int64)
for _i, _c in enumerate(CODONS):
    _rc = "".join(BASES[_COMP[BASES.index(ch)]] for ch in reversed(_c))
    _RC_PARTNER[_i] = CODONS.index(_rc)

MIRROR_PAIR_ID = np.minimum(np.arange(NC), _RC_PARTNER)
N_MIRROR_PAIRS = int(len(set(MIRROR_PAIR_ID.tolist())))
assert N_MIRROR_PAIRS == 32, f"expected 32 mirror pairs, got {N_MIRROR_PAIRS}"

ALPHA = 0.1
TINY = 1e-12

SPECIES_PREFIX = {
    "GCF_009914755": "Human_T2T", "GCF_000001635": "Mouse", "GCF_002880755": "Chimp",
    "GCF_000002985": "C.elegans", "GCF_000146045": "Yeast", "GCF_000005845": "E.coli",
    "GCF_000001215": "Drosophila", "GCF_000001735": "Arabidopsis",
    "GCF_000007745": "M.jannaschii", "GCF_000015585": "H.salinarum",
    "GCF_000002035": "Zebrafish", "GCF_016699485": "Chicken", "GCF_000009045": "B.subtilis",
    "GCF_000008895": "Human_mito", "GCF_009858895": "SARS-CoV-2", "GCF_000002025": "SARS-CoV-2",
}

SPECIES_INFO = {
    "Human_T2T": {"acc": "GCF_009914755", "known": ["GCF_009914755.1_T2T-CHM13v2.0"], "clade": "mammal"},
    "Mouse": {"acc": "GCF_000001635", "known": ["GCF_000001635.27_GRCm39", "GCF_000001635.26_GRCm38.p6", "GCF_000001635.25_GRCm38.p4"], "clade": "mammal"},
    "Chimp": {"acc": "GCF_002880755", "known": ["GCF_002880755.1_Clint_PTRv2", "GCF_002880755.1_Pan_troglodytes"], "clade": "mammal"},
    "Zebrafish": {"acc": "GCF_000002035", "known": ["GCF_000002035.6_GRCz11"], "clade": "vertebrate"},
    "Chicken": {"acc": "GCF_016699485", "known": ["GCF_016699485.1_gallGal_gca016699485v2", "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b"], "clade": "vertebrate"},
    "C.elegans": {"acc": "GCF_000002985", "known": ["GCF_000002985.6_WBcel235", "GCF_000002985.5_WBcel215"], "clade": "nematode"},
    "Drosophila": {"acc": "GCF_000001215", "known": ["GCF_000001215.4_Release_6_plus_ISO1_MT", "GCF_000001215.4"], "clade": "insect"},
    "Arabidopsis": {"acc": "GCF_000001735", "known": ["GCF_000001735.4_TAIR10.1", "GCF_000001735.4"], "clade": "plant"},
    "Yeast": {"acc": "GCF_000146045", "known": ["GCF_000146045.2_R64", "GCF_000146045.1_SacCer_May2010"], "clade": "fungus"},
    "E.coli": {"acc": "GCF_000005845", "known": ["GCF_000005845.2_ASM584v2", "GCF_000005845.1_ASM584v1"], "clade": "bacterium"},
    "B.subtilis": {"acc": "GCF_000009045", "known": ["GCF_000009045.1_ASM904v1", "GCF_000009045.1"], "clade": "bacterium"},
    "M.jannaschii": {"acc": "GCF_000007745", "known": ["GCF_000007745.1_ASM774v1", "GCF_000007745.1"], "clade": "archaeon"},
    "H.salinarum": {"acc": "GCF_000015585", "known": ["GCF_000015585.1_ASM1558v1", "GCF_000015585.1"], "clade": "archaeon", "note": "extreme halophile; expect high GC3 and acidic proteome"},
    "Human_mito": {
        "acc": "GCF_000008895",
        "known": ["GCF_000008895.2_ASM889v2", "GCF_000008895.1_ASM889v1", "GCF_000008895.3_ASM889v3", "GCF_000008895.1"],
        "clade": "mitochondrion", "optin": True,
        "note": "opt-in; if 404, override with --mito-acc/--mito-asm (silent skip otherwise)",
    },
    "SARS-CoV-2": {
        "acc": "GCF_009858895",
        "known": ["GCF_009858895.2_ASM985889v3", "GCF_009858895.1_ASM985889v2"],
        "clade": "virus", "optin": True,
        "note": "opt-in; overlapping ORFs possible; if 404 override with --virus-acc/--virus-asm",
    },
}

DEFAULT_H4_LENGTHS = (3, 10, 29, 30, 100)

DEFAULT_H5_TRANSFER_PAIRS = [
    "Human_T2T:Chimp", "Chimp:Human_T2T", "Human_T2T:Mouse", "Mouse:Human_T2T",
    "Chimp:Mouse", "Mouse:Chimp",
    "Zebrafish:Human_T2T", "Human_T2T:Zebrafish", "Chicken:Human_T2T", "Human_T2T:Chicken",
    "Zebrafish:Chicken", "Chicken:Zebrafish",
    "Yeast:Human_T2T", "Human_T2T:Yeast", "C.elegans:Human_T2T", "Human_T2T:C.elegans",
    "Drosophila:Human_T2T", "Human_T2T:Drosophila", "Arabidopsis:Human_T2T", "Human_T2T:Arabidopsis",
    "E.coli:Human_T2T", "Human_T2T:E.coli", "B.subtilis:Human_T2T", "Human_T2T:B.subtilis",
    "M.jannaschii:Human_T2T", "Human_T2T:M.jannaschii", "H.salinarum:Human_T2T", "Human_T2T:H.salinarum",
    "B.subtilis:E.coli", "E.coli:B.subtilis", "M.jannaschii:E.coli", "E.coli:M.jannaschii",
    "H.salinarum:M.jannaschii", "M.jannaschii:H.salinarum",
    "Drosophila:Yeast", "Yeast:Drosophila", "Arabidopsis:Drosophila", "Drosophila:Arabidopsis",
    "C.elegans:Drosophila", "Drosophila:C.elegans",
]

OPTIN_H5_TRANSFER_PAIRS = {
    "Human_mito": [
        "Human_mito:Human_T2T",
        "Human_T2T:Human_mito",
        "Human_mito:Chimp",
        "Chimp:Human_mito",
        "Human_mito:Mouse",
        "Mouse:Human_mito",
    ],
    "SARS-CoV-2": [
        "SARS-CoV-2:Human_T2T",
        "Human_T2T:SARS-CoV-2",
        "SARS-CoV-2:Chimp",
        "Chimp:SARS-CoV-2",
        "SARS-CoV-2:Mouse",
        "Mouse:SARS-CoV-2",
        "SARS-CoV-2:Yeast",
        "Yeast:SARS-CoV-2",
    ],
}


# ===========================================================================
# Curated Extended Biosphere Panel (41 Taxa)
# ===========================================================================

EXTENDED_SPECIES_PREFIX = {
    "GCF_015227675": "Rattus_norvegicus",
    "GCF_003339765": "Macaca_mulatta",
    "GCF_000151905": "Gorilla_gorilla",
    "GCF_002263795": "Bos_taurus",
    "GCF_000003025": "Sus_scrofa",
    "GCF_000002285": "Canis_lupus_familiaris",
    "GCF_000181335": "Felis_catus",
    "GCF_002863925": "Equus_caballus",
    "GCF_000003625": "Oryctolagus_cuniculus",
    "GCF_000004195": "Xenopus_tropicalis",
    "GCF_901000725": "Takifugu_rubripes",
    "GCF_000005575": "Anopheles_gambiae",
    "GCF_003254395": "Apis_mellifera",
    "GCF_014905235": "Bombyx_mori",
    "GCF_000004555": "Caenorhabditis_briggsae",
    "GCF_001433935": "Oryza_sativa",
    "GCF_000005005": "Zea_mays",
    "GCF_000188115": "Solanum_lycopersicum",
    "GCF_000004515": "Glycine_max",
    "GCF_000002775": "Populus_trichocarpa",
    "GCF_000002945": "Schizosaccharomyces_pombe",
    "GCF_000182965": "Candida_albicans",
    "GCF_000149205": "Aspergillus_nidulans",
    "GCF_000182925": "Neurospora_crassa",
    "GCF_000002765": "Plasmodium_falciparum",
    "GCF_000002445": "Trypanosoma_brucei",
    "GCF_000002725": "Leishmania_major",
    "GCF_000006565": "Toxoplasma_gondii",
    "GCF_000004695": "Dictyostelium_discoideum",
    "GCF_000002595": "Chlamydomonas_reinhardtii",
    "GCF_000006765": "Pseudomonas_aeruginosa",
    "GCF_000013425": "Staphylococcus_aureus",
    "GCF_000006785": "Streptococcus_pyogenes",
    "GCF_000195955": "Mycobacterium_tuberculosis",
    "GCF_000006945": "Salmonella_enterica",
    "GCF_000006865": "Lactococcus_lactis",
    "GCF_000008525": "Helicobacter_pylori",
    "GCF_000025985": "Bacteroides_fragilis",
    "GCF_000009205": "Clostridioides_difficile",
    "GCF_000009725": "Synechocystis_sp_PCC6803",
    "GCF_000008125": "Thermus_thermophilus",
}

EXTENDED_SPECIES_INFO = {
    "Rattus_norvegicus": {
        "acc": "GCF_015227675",
        "known": ["GCF_015227675.2_mRatBN7.2"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended rodent model",
    },
    "Macaca_mulatta": {
        "acc": "GCF_003339765",
        "known": ["GCF_003339765.1_Mmul_10"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended primate outgroup",
    },
    "Gorilla_gorilla": {
        "acc": "GCF_000151905",
        "known": ["GCF_000151905.2_gorGor4"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended great ape",
    },
    "Bos_taurus": {
        "acc": "GCF_002263795",
        "known": ["GCF_002263795.2_ARS-UCD1.3", "GCF_002263795.1_ARS-UCD1.2"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended livestock mammal",
    },
    "Sus_scrofa": {
        "acc": "GCF_000003025",
        "known": ["GCF_000003025.6_Sscrofa11.1"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended livestock mammal",
    },
    "Canis_lupus_familiaris": {
        "acc": "GCF_000002285",
        "known": ["GCF_000002285.3_CanFam3.1"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended carnivore mammal",
    },
    "Felis_catus": {
        "acc": "GCF_000181335",
        "known": ["GCF_000181335.3_Felis_catus_9.0"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended carnivore mammal",
    },
    "Equus_caballus": {
        "acc": "GCF_002863925",
        "known": ["GCF_002863925.1_EquCab3.0"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended perissodactyl mammal",
    },
    "Oryctolagus_cuniculus": {
        "acc": "GCF_000003625",
        "known": ["GCF_000003625.3_OryCun2.0"],
        "clade": "mammal", "optin": True, "extended": True,
        "note": "extended lagomorph mammal",
    },
    "Xenopus_tropicalis": {
        "acc": "GCF_000004195",
        "known": ["GCF_000004195.4_UCB_Xtro_10.0", "GCF_000004195.3_Xenopus_tropicalis_v9.1"],
        "clade": "vertebrate", "optin": True, "extended": True,
        "note": "extended amphibian vertebrate",
    },
    "Takifugu_rubripes": {
        "acc": "GCF_901000725",
        "known": ["GCF_901000725.2_fTakRub1.2"],
        "clade": "vertebrate", "optin": True, "extended": True,
        "note": "extended teleost fish; compact genome",
    },
    "Anopheles_gambiae": {
        "acc": "GCF_000005575",
        "known": ["GCF_000005575.2_AgamP3"],
        "clade": "insect", "optin": True, "extended": True,
        "note": "extended mosquito vector",
    },
    "Apis_mellifera": {
        "acc": "GCF_003254395",
        "known": ["GCF_003254395.2_Amel_HAv3.1"],
        "clade": "insect", "optin": True, "extended": True,
        "note": "extended hymenopteran insect",
    },
    "Bombyx_mori": {
        "acc": "GCF_014905235",
        "known": ["GCF_014905235.1_Bmori_2016v1.0"],
        "clade": "insect", "optin": True, "extended": True,
        "note": "extended lepidopteran insect",
    },
    "Caenorhabditis_briggsae": {
        "acc": "GCF_000004555",
        "known": ["GCF_000004555.2_CB4"],
        "clade": "nematode", "optin": True, "extended": True,
        "note": "extended nematode relative of C. elegans",
    },
    "Oryza_sativa": {
        "acc": "GCF_001433935",
        "known": ["GCF_001433935.1_IRGSP-1.0"],
        "clade": "plant", "optin": True, "extended": True,
        "note": "extended monocot crop",
    },
    "Zea_mays": {
        "acc": "GCF_000005005",
        "known": ["GCF_000005005.2_B73_RefGen_v4"],
        "clade": "plant", "optin": True, "extended": True,
        "note": "extended monocot crop; large genome",
    },
    "Solanum_lycopersicum": {
        "acc": "GCF_000188115",
        "known": ["GCF_000188115.5_SL3.1", "GCF_000188115.4_SL3.0"],
        "clade": "plant", "optin": True, "extended": True,
        "note": "extended eudicot crop",
    },
    "Glycine_max": {
        "acc": "GCF_000004515",
        "known": ["GCF_000004515.6_Glycine_max_v4.0", "GCF_000004515.5_Glycine_max_v2.1"],
        "clade": "plant", "optin": True, "extended": True,
        "note": "extended legume crop",
    },
    "Populus_trichocarpa": {
        "acc": "GCF_000002775",
        "known": ["GCF_000002775.5_Pop_tri_v4", "GCF_000002775.4_Pop_tri_v3"],
        "clade": "plant", "optin": True, "extended": True,
        "note": "extended woody plant",
    },
    "Schizosaccharomyces_pombe": {
        "acc": "GCF_000002945",
        "known": ["GCF_000002945.1_ASM294v2"],
        "clade": "fungus", "optin": True, "extended": True,
        "note": "extended fission yeast",
    },
    "Candida_albicans": {
        "acc": "GCF_000182965",
        "known": ["GCF_000182965.3_ASM18296v3"],
        "clade": "fungus", "optin": True, "extended": True,
        "note": "extended pathogenic yeast; uses transl_table=12 (CTG->S)",
    },
    "Aspergillus_nidulans": {
        "acc": "GCF_000149205",
        "known": ["GCF_000149205.2_ASM14920v2"],
        "clade": "fungus", "optin": True, "extended": True,
        "note": "extended filamentous fungus",
    },
    "Neurospora_crassa": {
        "acc": "GCF_000182925",
        "known": ["GCF_000182925.2_NC12"],
        "clade": "fungus", "optin": True, "extended": True,
        "note": "extended filamentous fungus",
    },
    "Plasmodium_falciparum": {
        "acc": "GCF_000002765",
        "known": ["GCF_000002765.6_GCA_000002765"],
        "clade": "protist", "optin": True, "extended": True,
        "note": "extended apicomplexan parasite; extreme AT-rich composition control",
    },
    "Trypanosoma_brucei": {
        "acc": "GCF_000002445",
        "known": ["GCF_000002445.2_ASM244v1"],
        "clade": "protist", "optin": True, "extended": True,
        "note": "extended kinetoplastid parasite",
    },
    "Leishmania_major": {
        "acc": "GCF_000002725",
        "known": ["GCF_000002725.2_ASM272v2"],
        "clade": "protist", "optin": True, "extended": True,
        "note": "extended kinetoplastid parasite",
    },
    "Toxoplasma_gondii": {
        "acc": "GCF_000006565",
        "known": ["GCF_000006565.2_TGA4"],
        "clade": "protist", "optin": True, "extended": True,
        "note": "extended apicomplexan parasite",
    },
    "Dictyostelium_discoideum": {
        "acc": "GCF_000004695",
        "known": ["GCF_000004695.1_dicty_2.7"],
        "clade": "protist", "optin": True, "extended": True,
        "note": "extended social amoeba; extreme AT-rich composition control",
    },
    "Chlamydomonas_reinhardtii": {
        "acc": "GCF_000002595",
        "known": ["GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5"],
        "clade": "plant-like protist", "optin": True, "extended": True,
        "note": "extended green alga; high-GC composition control",
    },
    "Pseudomonas_aeruginosa": {
        "acc": "GCF_000006765",
        "known": ["GCF_000006765.1_ASM676v1"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended opportunistic pathogen; high-GC bacterium",
    },
    "Staphylococcus_aureus": {
        "acc": "GCF_000013425",
        "known": ["GCF_000013425.1_ASM1342v1"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended Gram-positive pathogen; low-GC bacterium",
    },
    "Streptococcus_pyogenes": {
        "acc": "GCF_000006785",
        "known": ["GCF_000006785.2_ASM678v2"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended Gram-positive pathogen",
    },
    "Mycobacterium_tuberculosis": {
        "acc": "GCF_000195955",
        "known": ["GCF_000195955.2_ASM19595v2"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended high-GC Actinobacterium pathogen",
    },
    "Salmonella_enterica": {
        "acc": "GCF_000006945",
        "known": ["GCF_000006945.2_ASM694v2"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended enterobacterium pathogen",
    },
    "Lactococcus_lactis": {
        "acc": "GCF_000006865",
        "known": ["GCF_000006865.1_ASM686v1"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended lactic acid bacterium",
    },
    "Helicobacter_pylori": {
        "acc": "GCF_000008525",
        "known": ["GCF_000008525.1_ASM852v1"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended host-associated pathogen",
    },
    "Bacteroides_fragilis": {
        "acc": "GCF_000025985",
        "known": ["GCF_000025985.1_ASM2598v1"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended gut Bacteroidetes",
    },
    "Clostridioides_difficile": {
        "acc": "GCF_000009205",
        "known": ["GCF_000009205.2_ASM920v2"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "extended gut Firmicute pathogen; low-GC control",
    },
    "Synechocystis_sp_PCC6803": {
        "acc": "GCF_000009725",
        "known": ["GCF_000009725.1_ASM972v1"],
        "clade": "cyanobacterium", "optin": True, "extended": True,
        "note": "cyanobacterium; plastid-ancestor composition control",
    },
    "Thermus_thermophilus": {
        "acc": "GCF_000008125",
        "known": ["GCF_000008125.1_ASM812v1"],
        "clade": "bacterium", "optin": True, "extended": True,
        "note": "thermophile; high-temperature/high-GC composition control",
    },
}

EXTENDED_H5_TRANSFER_PAIRS = [
    # Mammalian conservation around human
    "Human_T2T:Rattus_norvegicus", "Rattus_norvegicus:Human_T2T",
    "Human_T2T:Macaca_mulatta", "Macaca_mulatta:Human_T2T",
    "Human_T2T:Gorilla_gorilla", "Gorilla_gorilla:Human_T2T",
    "Macaca_mulatta:Rattus_norvegicus", "Rattus_norvegicus:Macaca_mulatta",
    "Chimp:Rattus_norvegicus", "Rattus_norvegicus:Chimp",
    "Mouse:Rattus_norvegicus", "Rattus_norvegicus:Mouse",

    "Human_T2T:Bos_taurus", "Bos_taurus:Human_T2T",
    "Human_T2T:Sus_scrofa", "Sus_scrofa:Human_T2T",
    "Bos_taurus:Sus_scrofa", "Sus_scrofa:Bos_taurus",

    "Human_T2T:Canis_lupus_familiaris", "Canis_lupus_familiaris:Human_T2T",
    "Human_T2T:Felis_catus", "Felis_catus:Human_T2T",
    "Canis_lupus_familiaris:Felis_catus", "Felis_catus:Canis_lupus_familiaris",

    "Human_T2T:Equus_caballus", "Equus_caballus:Human_T2T",
    "Human_T2T:Oryctolagus_cuniculus", "Oryctolagus_cuniculus:Human_T2T",

    # Vertebrate breadth
    "Zebrafish:Xenopus_tropicalis", "Xenopus_tropicalis:Zebrafish",
    "Human_T2T:Xenopus_tropicalis", "Xenopus_tropicalis:Human_T2T",
    "Chicken:Xenopus_tropicalis", "Xenopus_tropicalis:Chicken",

    "Zebrafish:Takifugu_rubripes", "Takifugu_rubripes:Zebrafish",
    "Human_T2T:Takifugu_rubripes", "Takifugu_rubripes:Human_T2T",

    # Insect / nematode expansion
    "Drosophila:Anopheles_gambiae", "Anopheles_gambiae:Drosophila",
    "Drosophila:Apis_mellifera", "Apis_mellifera:Drosophila",
    "Drosophila:Bombyx_mori", "Bombyx_mori:Drosophila",
    "Anopheles_gambiae:Apis_mellifera", "Apis_mellifera:Anopheles_gambiae",

    "C.elegans:Caenorhabditis_briggsae", "Caenorhabditis_briggsae:C.elegans",
    "Caenorhabditis_briggsae:Human_T2T", "Human_T2T:Caenorhabditis_briggsae",

    # Plant expansion
    "Arabidopsis:Oryza_sativa", "Oryza_sativa:Arabidopsis",
    "Arabidopsis:Zea_mays", "Zea_mays:Arabidopsis",
    "Arabidopsis:Solanum_lycopersicum", "Solanum_lycopersicum:Arabidopsis",
    "Arabidopsis:Glycine_max", "Glycine_max:Arabidopsis",
    "Arabidopsis:Populus_trichocarpa", "Populus_trichocarpa:Arabidopsis",

    "Oryza_sativa:Zea_mays", "Zea_mays:Oryza_sativa",
    "Solanum_lycopersicum:Glycine_max", "Glycine_max:Solanum_lycopersicum",
    "Populus_trichocarpa:Solanum_lycopersicum", "Solanum_lycopersicum:Populus_trichocarpa",

    # Fungal expansion
    "Yeast:Schizosaccharomyces_pombe", "Schizosaccharomyces_pombe:Yeast",
    "Yeast:Candida_albicans", "Candida_albicans:Yeast",
    "Yeast:Aspergillus_nidulans", "Aspergillus_nidulans:Yeast",
    "Yeast:Neurospora_crassa", "Neurospora_crassa:Yeast",

    "Candida_albicans:Aspergillus_nidulans", "Aspergillus_nidulans:Candida_albicans",
    "Neurospora_crassa:Aspergillus_nidulans", "Aspergillus_nidulans:Neurospora_crassa",
    "Schizosaccharomyces_pombe:Candida_albicans", "Candida_albicans:Schizosaccharomyces_pombe",

    # Protist / parasite expansion
    "Plasmodium_falciparum:Human_T2T", "Human_T2T:Plasmodium_falciparum",
    "Trypanosoma_brucei:Human_T2T", "Human_T2T:Trypanosoma_brucei",
    "Leishmania_major:Human_T2T", "Human_T2T:Leishmania_major",
    "Toxoplasma_gondii:Human_T2T", "Human_T2T:Toxoplasma_gondii",
    "Dictyostelium_discoideum:Human_T2T", "Human_T2T:Dictyostelium_discoideum",

    "Plasmodium_falciparum:Toxoplasma_gondii", "Toxoplasma_gondii:Plasmodium_falciparum",
    "Trypanosoma_brucei:Leishmania_major", "Leishmania_major:Trypanosoma_brucei",
    "Dictyostelium_discoideum:Trypanosoma_brucei", "Trypanosoma_brucei:Dictyostelium_discoideum",

    # Algal / plastid-ancestor controls
    "Chlamydomonas_reinhardtii:Arabidopsis", "Arabidopsis:Chlamydomonas_reinhardtii",
    "Synechocystis_sp_PCC6803:Arabidopsis", "Arabidopsis:Synechocystis_sp_PCC6803",
    "Synechocystis_sp_PCC6803:Chlamydomonas_reinhardtii", "Chlamydomonas_reinhardtii:Synechocystis_sp_PCC6803",

    # Bacterial expansion
    "E.coli:Pseudomonas_aeruginosa", "Pseudomonas_aeruginosa:E.coli",
    "E.coli:Salmonella_enterica", "Salmonella_enterica:E.coli",
    "Pseudomonas_aeruginosa:Salmonella_enterica", "Salmonella_enterica:Pseudomonas_aeruginosa",

    "B.subtilis:Staphylococcus_aureus", "Staphylococcus_aureus:B.subtilis",
    "B.subtilis:Streptococcus_pyogenes", "Streptococcus_pyogenes:B.subtilis",
    "Staphylococcus_aureus:Streptococcus_pyogenes", "Streptococcus_pyogenes:Staphylococcus_aureus",

    "Mycobacterium_tuberculosis:Human_T2T", "Human_T2T:Mycobacterium_tuberculosis",
    "Mycobacterium_tuberculosis:B.subtilis", "B.subtilis:Mycobacterium_tuberculosis",
    "Mycobacterium_tuberculosis:Streptococcus_pyogenes", "Streptococcus_pyogenes:Mycobacterium_tuberculosis",

    "Helicobacter_pylori:Bacteroides_fragilis", "Bacteroides_fragilis:Helicobacter_pylori",
    "Clostridioides_difficile:Bacteroides_fragilis", "Bacteroides_fragilis:Clostridioides_difficile",
    "Lactococcus_lactis:Streptococcus_pyogenes", "Streptococcus_pyogenes:Lactococcus_lactis",

    # Extremophile / temperature controls
    "Thermus_thermophilus:B.subtilis", "B.subtilis:Thermus_thermophilus",
    "Thermus_thermophilus:M.jannaschii", "M.jannaschii:Thermus_thermophilus",
    "Thermus_thermophilus:H.salinarum", "H.salinarum:Thermus_thermophilus",
]


# ===========================================================================
# Curated Human RefSeqGene Panel for Intron vs. Exon & P53 Audit (v14.6)
# ===========================================================================

INTRON_EXON_REFSEQGENE_PANEL = {
    "TP53": {
        "acc": "NG_017013.2",
        "role": "Primary P53 target (Tumor protein p53; 10 coding exons, 9 coding introns)",
        "family": "TP53",
    },
    "TP63": {
        "acc": "NG_007550.1",
        "role": "P53 family paralog (Tumor protein p63)",
        "family": "TP63",
    },
    "TP73": {
        "acc": "NG_017035.2",
        "role": "P53 family paralog (Tumor protein p73)",
        "family": "TP73",
    },
    "BRCA1": {
        "acc": "NG_005905.2",
        "role": "DNA repair tumor suppressor (large multi-intron locus)",
        "family": "BRCA1",
    },
    "PTEN": {
        "acc": "NG_007466.2",
        "role": "Phosphatase and tensin homolog tumor suppressor",
        "family": "PTEN",
    },
    "RB1": {
        "acc": "NG_009009.1",
        "role": "Retinoblastoma transcriptional corepressor 1 (27 coding exons)",
        "family": "RB1",
    },
    "EGFR": {
        "acc": "NG_007726.3",
        "role": "Epidermal growth factor receptor tyrosine kinase",
        "family": "EGFR",
    },
    "CFTR": {
        "acc": "NG_016465.4",
        "role": "CFTR chloride channel (~189 kb genomic span, 27 exons)",
        "family": "CFTR",
    },
    "GAPDH": {
        "acc": "NG_007073.2",
        "role": "High-expression housekeeping glycolytic enzyme",
        "family": "GAPDH",
    },
    "ACTB": {
        "acc": "NG_007992.1",
        "role": "High-expression cytoskeletal beta-actin",
        "family": "ACTB",
    },
    "HBB": {
        "acc": "NG_000007.3",
        "role": "Hemoglobin subunit beta (classic compact 3-exon / 2-intron locus)",
        "family": "HBB",
    },
    "KRAS": {
        "acc": "NG_007524.2",
        "role": "KRAS proto-oncogene GTPase (RAS family member 1; AT-rich exons vs introns)",
        "family": "RAS",
    },
    "NRAS": {
        "acc": "NG_007572.1",
        "role": "NRAS proto-oncogene GTPase (RAS family member 2; neuroblastoma RAS)",
        "family": "RAS",
    },
    "HRAS": {
        "acc": "NG_007666.1",
        "role": "HRAS proto-oncogene GTPase (RAS family member 3; Harvey rat sarcoma viral oncogene)",
        "family": "RAS",
    },
    "MYC": {
        "acc": "NG_007161.2",
        "role": "MYC proto-oncogene transcription factor",
        "family": "MYC",
    },
    "VHL": {
        "acc": "NG_008212.3",
        "role": "von Hippel-Lindau tumor suppressor",
        "family": "VHL",
    },
    "APOE": {
        "acc": "NG_007084.2",
        "role": "Apolipoprotein E (high-GC coding sequence)",
        "family": "APOE",
    },
}

DEFAULT_P53_PANEL_GENES = [
    "TP53", "TP63", "TP73", "BRCA1", "PTEN", "RB1",
    "GAPDH", "ACTB", "HBB", "KRAS", "MYC", "APOE",
]


# ===========================================================================
# Genetic code context
# ===========================================================================

def set_genetic_code(table_id: int, verbose: bool = True):
    global CODON_AA, AA_IND, CHANNEL_CODONS
    tid = int(table_id)
    if tid in SUPPORTED_TABLES:
        status = "exact_map"
    else:
        if verbose:
            print(
                f"    [code] transl_table={tid} not in verified set "
                f"{sorted(SUPPORTED_TABLES)} -> standard-map fallback "
                f"(H2/H5/H6/H7 for this taxon = map-fallback, interpret with care)"
            )
        tid = 1
        status = "fallback_standard_map"

    aa_list = list(STANDARD_AA_LIST)
    for codon, aa in TABLE_OVERRIDES[tid].items():
        aa_list[CODONS.index(codon)] = aa

    CODON_AA = np.array([AA_INDEX[ch] for ch in aa_list], dtype=np.int64)
    AA_IND = np.zeros((NC, NA), dtype=np.float64)
    AA_IND[np.arange(NC), CODON_AA] = 1.0
    CHANNEL_CODONS = [np.where(AA_IND[:, a] > 0)[0].astype(np.int64) for a in range(NA)]
    return tid, status


# ===========================================================================
# Generic helpers
# ===========================================================================

def species_of(path: Path) -> str:
    """
    Match longer accessions first to avoid prefix ambiguity.
    """
    name = path.name
    for acc, sp in sorted(SPECIES_PREFIX.items(), key=lambda kv: len(kv[0]), reverse=True):
        if acc in name:
            return sp
    return path.stem[:12]


def clade_of(sp: str) -> str:
    return SPECIES_INFO.get(sp, {}).get("clade", "unknown")


def is_optin(sp: str) -> bool:
    return bool(SPECIES_INFO.get(sp, {}).get("optin", False))


def natural_key(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def extract_gene_id(header: str) -> str:
    for key in ("gene", "locus_tag", "protein_id"):
        m = re.search(rf"\[{key}=([^\]]+)\]", header)
        if m:
            return m.group(1)
    return header.lstrip(">").split()[0]


def extract_transl_table(header: str) -> int:
    m = re.search(r"\[transl_table=(\d+)\]", header)
    return int(m.group(1)) if m else 1


def family_root(gene_id: str) -> str:
    s = gene_id
    prev = None
    pat = re.compile(r"(-R[A-Za-z]\d*|\.\d+|_\d+|\|.*)$")
    while s != prev:
        prev = s
        s = pat.sub("", s)
    return s if s else gene_id


def new_counts():
    return {
        "bigram": np.zeros((NC, NC), dtype=np.int64),
        "marg": np.zeros(NC, dtype=np.int64),
        "lastnt": np.zeros((4, NC), dtype=np.int64),
        "aa_next": np.zeros((NA, NC), dtype=np.int64),
        "aa_pair": np.zeros((NA, NA), dtype=np.int64),
        "nt_bigram": np.zeros((4, 4), dtype=np.int64),
        "nt_tri": np.zeros((4, 4, 4), dtype=np.int64),
        "n_codons": 0, "n_trans": 0, "n_nt": 0,
        "gc3_sum": 0, "gc3_n": 0,
        "acid_count": 0, "basic_count": 0, "aa_count": 0,
    }


def bytes_to_codons(chunk: bytes):
    arr = np.frombuffer(chunk, dtype=np.uint8)
    nt = NT_LUT[arr]
    nt = nt[nt < 4]
    m = nt.size // 3
    if m == 0:
        return np.empty(0, dtype=np.uint8), np.empty(0, dtype=np.int8)
    nt = nt[:m * 3].astype(np.int8)
    codons = (
        16 * nt[0::3].astype(np.int64)
        + 4 * nt[1::3].astype(np.int64)
        + nt[2::3].astype(np.int64)
    ).astype(np.uint8)
    return codons, nt


def add_unit(uc, codons, nt):
    if codons.size == 0:
        return

    codi = codons.astype(np.int64)
    uc["n_codons"] += int(codi.size)
    uc["marg"] += np.bincount(codi, minlength=NC).astype(np.int64)
    uc["gc3_sum"] += int(IS_GC3[codi].sum())
    uc["gc3_n"] += int(codi.size)

    aa_all = CODON_AA[codi]
    uc["acid_count"] += int(ACID_MASK[aa_all].sum())
    uc["basic_count"] += int(BASIC_MASK[aa_all].sum())
    uc["aa_count"] += int(NONSTOP_MASK[aa_all].sum())

    if nt.size >= 2:
        nti = nt.astype(np.int64)
        flat2 = nti[:-1] * 4 + nti[1:]
        uc["nt_bigram"] += np.bincount(flat2, minlength=16).reshape(4, 4).astype(np.int64)
        uc["n_nt"] += int(nt.size)
        if nt.size >= 3:
            flat3 = nti[:-2] * 16 + nti[1:-1] * 4 + nti[2:]
            uc["nt_tri"] += np.bincount(flat3, minlength=64).reshape(4, 4, 4).astype(np.int64)

    prev = codi[:-1]
    nxt = codi[1:]
    if prev.size:
        uc["bigram"] += np.bincount(prev * NC + nxt, minlength=NC * NC).reshape(NC, NC).astype(np.int64)
        uc["lastnt"] += np.bincount((prev % 4) * NC + nxt, minlength=4 * NC).reshape(4, NC).astype(np.int64)
        a1 = CODON_AA[prev]
        a2 = CODON_AA[nxt]
        uc["aa_next"] += np.bincount(a1 * NC + nxt, minlength=NA * NC).reshape(NA, NC).astype(np.int64)
        uc["aa_pair"] += np.bincount(a1 * NA + a2, minlength=NA * NA).reshape(NA, NA).astype(np.int64)
        uc["n_trans"] += int(prev.size)


def counts_sample_ok(train_units, test_units, train_trans, test_trans, min_units, min_trans):
    """
    Return True if sample sizes are adequate for strict inference.
    If min_units/min_trans <= 0, that guard is disabled.
    """
    tu = int(train_units)
    te = int(test_units)
    tr = int(train_trans)
    tt = int(test_trans)

    if min_units > 0:
        if tu < min_units or te < min_units:
            return False
    if min_trans > 0:
        if tr < min_trans or tt < min_trans:
            return False
    return True


def is_small_sample(rec_src, rec_dst, min_units, min_trans):
    return not counts_sample_ok(
        rec_src.get("train_units", 0),
        rec_dst.get("test_units", 0),
        rec_src.get("train_transitions", 0),
        rec_dst.get("test_transitions", 0),
        min_units,
        min_trans,
    )


# ===========================================================================
# Download (retry + validation + SILENT graceful skip for opt-in)
# ===========================================================================

def ftp_base_for_accession(acc: str) -> str:
    """
    Build NCBI RefSeq HTTPS directory URL without literal URL strings that
    chat/Markdown parsers can corrupt into [url](url).
    """
    tail = acc.split("_", 1)[1] if "_" in acc else acc
    parts = [tail[i:i + 3] for i in range(0, len(tail), 3)]
    scheme = "https"
    host = "ftp.ncbi.nlm.nih.gov"
    return f"{scheme}://{host}/genomes/all/GCF/" + "/".join(parts) + "/"


def eutils_efetch_url(acc: str) -> str:
    """
    Build NCBI E-utilities efetch URL safely (immune to Markdown auto-linkers).
    """
    scheme = "https"
    host = "eutils.ncbi.nlm.nih.gov"
    path = "entrez/eutils/efetch.fcgi"
    return f"{scheme}://{host}/{path}?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text"


def parse_version_key(name: str):
    m = re.match(r"GCF_(\d+)\.(\d+)", name)
    return (int(m.group(1)), int(m.group(2)), name) if m else (0, 0, name)


def _retry_delay(err, attempt: int) -> float:
    try:
        if isinstance(err, urllib.error.HTTPError) and getattr(err, "headers", None) is not None:
            ra = err.headers.get("Retry-After")
            if ra:
                return min(max(float(ra), 1.0), 120.0)
    except Exception:
        pass
    return min(2.0 * (attempt + 1), 30.0)


def http_get(url: str, timeout: int = 30, attempts: int = 2, quiet: bool = False) -> bytes:
    last = None
    attempts = max(1, int(attempts))
    for k in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as e:
            last = e
            if not quiet:
                print(f"    GET attempt {k + 1}/{attempts} failed: {e}")
            if k + 1 < attempts:
                time.sleep(min(2.0 * (k + 1), 10.0))
    raise RuntimeError(f"GET failed: {url}: {last}") from last


def list_ncbi_dirs(url: str, timeout: int = 30, verbose: bool = True):
    try:
        html = http_get(url, timeout=timeout, attempts=2, quiet=not verbose).decode("utf-8", "ignore")
        links = re.findall(r'href="([^"]+)"', html)
        return sorted({
            link.rstrip("/").split("/")[-1]
            for link in links
            if link.rstrip("/").split("/")[-1].startswith("GCF_")
        })
    except Exception as e:
        if verbose:
            print(f"    directory listing failed: {e}")
        return []


def candidate_assembly_dirs(acc: str, known, timeout: int = 30, verbose: bool = True):
    base = ftp_base_for_accession(acc)
    dirs = list_ncbi_dirs(base, timeout=timeout, verbose=verbose)
    matching = sorted([d for d in dirs if d.startswith(acc)], key=parse_version_key, reverse=True)
    cands = []
    seen = set()
    for d in matching + list(known):
        if d and d not in seen:
            seen.add(d)
            cands.append(d)
    return cands


def existing_cds_file(cds_dir: Path, acc: str):
    cands = [
        p for p in cds_dir.glob("*_cds_from_genomic.fna.gz")
        if acc in p.name and p.stat().st_size > 1024
    ]
    return max(cands, key=lambda x: x.stat().st_size) if cands else None


def remove_existing_cds_files(cds_dir: Path, acc: str):
    for p in cds_dir.glob("*_cds_from_genomic.fna.gz"):
        if acc in p.name:
            try:
                p.unlink()
                print(f"    removed old {p.name}")
            except Exception as e:
                print(f"    could not remove {p}: {e}")


def validate_fasta_gz(path: Path):
    """Stronger validation: gzip must decompress and first byte must be '>'."""
    with gzip.open(path, "rb") as f:
        head = f.read(1 << 20)
    if not head:
        raise ValueError("downloaded gzip file is empty after decompression")
    if not head.startswith(b">"):
        raise ValueError("downloaded file does not look like FASTA")


def download_with_retry(url: str, dest: Path, timeout: int = 180, attempts: int = 3,
                        quiet: bool = False) -> bool:
    """Returns True on success, False on failure. quiet=True suppresses progress/errors."""
    tmp = dest.with_suffix(dest.suffix + ".part")
    attempts = max(1, int(attempts))

    for k in range(attempts):
        if tmp.exists():
            try:
                tmp.unlink()
            except Exception:
                pass

        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                total = None
                try:
                    cl = resp.headers.get("Content-Length")
                    if cl and int(cl) > 0:
                        total = int(cl)
                except Exception:
                    total = None

                read = 0
                with open(tmp, "wb") as out:
                    while True:
                        chunk = resp.read(1 << 20)
                        if not chunk:
                            break
                        out.write(chunk)
                        read += len(chunk)
                        if not quiet:
                            if total:
                                sys.stdout.write(
                                    f"\r    {read / 1e6:8.1f}/{total / 1e6:8.1f} MB "
                                    f"{100.0 * read / total:5.1f}%"
                                )
                            else:
                                sys.stdout.write(f"\r    {read / 1e6:8.1f} MB")
                            sys.stdout.flush()

                if not quiet:
                    print()

            validate_fasta_gz(tmp)
            tmp.replace(dest)
            return True

        except urllib.error.HTTPError as e:
            if not quiet:
                print(f"\n    attempt {k + 1}/{attempts} failed: HTTP {e.code} {e.reason}")
            if e.code == 404:
                break
            if k + 1 < attempts:
                time.sleep(_retry_delay(e, k))

        except Exception as e:
            if not quiet:
                print(f"\n    attempt {k + 1}/{attempts} failed: {e}")
            if k + 1 < attempts:
                time.sleep(min(2.0 * (k + 1), 30.0))

        finally:
            if tmp.exists():
                try:
                    tmp.unlink()
                except Exception:
                    pass

    return False


def _resolve_acc_known(sp: str, info: dict, overrides: dict):
    """Return (acc, known_list) honoring CLI overrides for opt-in tracks."""
    acc = info["acc"]
    known = list(info.get("known", []))
    ov = (overrides or {}).get(sp)
    if ov:
        oacc, oasm = ov
        if oacc:
            acc = oacc
        if oasm:
            if oasm.startswith("GCF_"):
                known = [oasm] + known
            else:
                prefixes = [
                    f"{acc}_",
                    f"{acc}.1_",
                    f"{acc}.2_",
                    f"{acc}.3_",
                ]
                for pref in prefixes:
                    known.insert(0, pref + oasm)
    return acc, known


def prepare_cds_cache(cds_dir: Path, force=False, timeout=180,
                      include_mito=False, include_virus=False,
                      include_extended=False, overrides=None, only=None):
    cds_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 92)
    print("AUTO CDS PREPARATION (autonomous v14.6)")
    print("=" * 92)
    print(f"CDS cache directory: {cds_dir}")
    print(f"Mitochondrion track: {'ENABLED (opt-in)' if include_mito else 'off (--include-mito)'}")
    print(f"Virus track        : {'ENABLED (opt-in)' if include_virus else 'off (--include-virus)'}")
    print(f"Extended panel     : {'ENABLED (opt-in)' if include_extended else 'off (--include-extended)'}")

    downloaded, skipped, missing_core, opted_out, optin_unavailable = [], [], [], [], []

    for sp, info in SPECIES_INFO.items():
        if only and not any(o.lower() in sp.lower() for o in only):
            continue

        is_mito = info.get("clade") == "mitochondrion"
        is_virus = info.get("clade") == "virus"
        is_extended = bool(info.get("extended", False))
        optin = bool(info.get("optin", False))

        if is_extended and not include_extended and not only:
            continue

        if optin and not only and not (
            (is_mito and include_mito)
            or (is_virus and include_virus)
            or (is_extended and include_extended)
        ):
            opted_out.append(sp)
            continue

        acc, known = _resolve_acc_known(sp, info, overrides)
        ex = existing_cds_file(cds_dir, acc)

        if ex is not None and not force:
            print(f"[CDS] {sp:14s}: existing {ex.name} ({ex.stat().st_size / 1e6:.1f} MB), skip download")
            if optin:
                print(f"[opt-in] {sp:14s}: available at accession {acc}; using cache")
            skipped.append(sp)
            continue

        if ex is not None and force:
            print(f"[CDS] {sp:14s}: force download, removing existing files")
            remove_existing_cds_files(cds_dir, acc)

        quiet = optin
        attempts = 2 if optin else 3

        if not quiet:
            print(f"[CDS] {sp:14s}: accession {acc}")

        cands = candidate_assembly_dirs(acc, known, timeout=timeout, verbose=not quiet)
        base = ftp_base_for_accession(acc)
        ok = False

        for dir_name in cands:
            fname = f"{dir_name}_cds_from_genomic.fna.gz"
            url = base + dir_name + "/" + fname
            dest = cds_dir / fname

            if dest.exists() and dest.stat().st_size > 1024:
                if not quiet:
                    print(f"      already present: {dest.name}")
                ok = True
                break

            if not quiet:
                print(f"      trying {url}")

            if download_with_retry(url, dest, timeout=timeout, attempts=attempts, quiet=quiet):
                if not quiet:
                    print(f"      OK: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
                ok = True
                break
            else:
                if dest.exists():
                    try:
                        dest.unlink()
                    except Exception:
                        pass

        if ok:
            downloaded.append(sp)
            if optin:
                print(f"[opt-in] {sp:14s}: available at accession {acc}; downloaded/used cache")
        elif optin:
            note = info.get("note", "")
            print(
                f"[opt-in] {sp:14s}: not available on NCBI at accession {acc}; "
                f"silent skip (core run unaffected){('  [' + note + ']') if note else ''}"
            )
            optin_unavailable.append(sp)
        else:
            note = info.get("note", "")
            print(
                f"      WARNING: could not download CDS for {sp}"
                + (f"  [{note}]" if note else "")
            )
            missing_core.append(sp)

    print()
    if downloaded:
        print(f"Downloaded CDS for: {downloaded}")
    if skipped:
        print(f"Skipped existing CDS for: {skipped}")
    if missing_core:
        print(f"Missing CORE CDS for: {missing_core}")
    if optin_unavailable:
        print(f"Opt-in unavailable (silent skip, core unaffected): {optin_unavailable}")
        print("    -> to force a working accession: --mito-acc/--mito-asm or --virus-acc/--virus-asm")
    if opted_out:
        print(
            "Opted-out (enable with --include-mito / --include-virus / --include-extended): "
            f"{opted_out}"
        )

    return downloaded, skipped, missing_core, optin_unavailable


def collect_cds_files(cds_dir: Path):
    return sorted(
        [p for p in cds_dir.glob("*_cds_from_genomic.fna.gz") if p.stat().st_size > 1024],
        key=lambda p: natural_key(p.name)
    )


def preferred_accession_for_species(sp: str, overrides=None):
    ov = (overrides or {}).get(sp)
    if ov:
        acc = ov[0]
        if acc:
            return acc
    return SPECIES_INFO.get(sp, {}).get("acc")


def deduplicate_species_files(paths, overrides=None):
    """
    Collapse duplicate species files.

    Preference order:
      1. file matching explicit preferred accession (override or canonical);
      2. largest file.
    """
    by_sp = {}
    for p in paths:
        sp = species_of(p)
        pref = preferred_accession_for_species(sp, overrides)
        old = by_sp.get(sp)
        if old is None:
            by_sp[sp] = p
            continue

        p_match = bool(pref and pref in p.name)
        old_match = bool(pref and pref in old.name)

        if p_match and not old_match:
            by_sp[sp] = p
        elif p_match == old_match and p.stat().st_size > old.stat().st_size:
            by_sp[sp] = p

    return sorted(by_sp.values(), key=lambda p: natural_key(p.name))


# Register extended biosphere dictionaries so species_of() and prepare_cds_cache() see them
SPECIES_PREFIX.update(EXTENDED_SPECIES_PREFIX)
SPECIES_INFO.update(EXTENDED_SPECIES_INFO)


# ===========================================================================
# GenBank Splicing Parser & Geometric Form Metrics (v14.6)
# ===========================================================================

def revcomp_dna_str(seq: str) -> str:
    tr = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tr)[::-1]


def parse_genbank_intron_exon(gb_text: str, target_gene: str = None):
    """
    Parse an NCBI GenBank record (e.g. NG_017013.2 for TP53) to extract:
      - full ORIGIN genomic sequence (preserving exact 1-based coordinate indexing)
      - canonical CDS exon intervals [(start_1based, end_1based), ...]
      - intervening intron intervals
      - spliced_cds_str (exons concatenated in true reading frame)
      - full_unspliced_str (start of first coding exon to end of last coding exon)
      - introns_concat_str (all coding introns concatenated)
    """
    m_orig = re.search("^ORIGIN.*$", gb_text, flags=re.MULTILINE)
    if not m_orig:
        raise ValueError("GenBank text has no ORIGIN section")

    feat_text = gb_text[:m_orig.start()]
    orig_body = gb_text[m_orig.end():]
    end_idx = orig_body.find("//")
    if end_idx >= 0:
        orig_body = orig_body[:end_idx]

    # Extract all nucleotide letters [a-zA-Z] after ORIGIN header line so 1-based coordinates never drift
    seq_chars = re.findall(r"[a-zA-Z]+", orig_body)
    full_seq = "".join(seq_chars).upper()
    if not full_seq:
        raise ValueError("GenBank ORIGIN sequence is empty")

    stop_set = {"TAA", "TAG", "TGA"}

    def _spliced_from_intervals(intervals_list, is_comp_flag, codon_start=1):
        ordered = sorted(intervals_list, key=lambda x: x[0])
        ex_list = [full_seq[a - 1:b] for a, b in ordered]
        in_list = []
        ex_coords = list(ordered)
        in_coords = []
        for idx_i, (a, b) in enumerate(ordered):
            if idx_i + 1 < len(ordered):
                na, _nb = ordered[idx_i + 1]
                if na > b + 1:
                    in_list.append(full_seq[b:na - 1])
                else:
                    in_list.append("")
                in_coords.append((b + 1, na - 1))
        span_s = ordered[0][0]
        span_e = ordered[-1][1]
        unspliced_s = full_seq[span_s - 1:span_e]
        if is_comp_flag:
            ex_list = [revcomp_dna_str(x) for x in reversed(ex_list)]
            in_list = [revcomp_dna_str(x) for x in reversed(in_list)]
            unspliced_s = revcomp_dna_str(unspliced_s)
        if codon_start > 1 and ex_list:
            trim_n = codon_start - 1
            if len(ex_list[0]) > trim_n:
                ex_list[0] = ex_list[0][trim_n:]
                unspliced_s = unspliced_s[trim_n:]
        return ex_list, in_list, unspliced_s, ex_coords, in_coords

    cds_candidates = []
    lines = feat_text.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        m_cds = re.match(r"^\s{5}CDS\s+(.*)$", line)
        if m_cds:
            loc_parts = [m_cds.group(1).strip()]
            i += 1
            while i < n:
                cont = lines[i]
                if re.match(r"^\s{21}\S", cont) and not cont.lstrip().startswith("/"):
                    loc_parts.append(cont.strip())
                    i += 1
                else:
                    break
            loc_str = "".join(loc_parts)
            quals = []
            while i < n:
                cont = lines[i]
                if re.match(r"^\s{21}\S", cont):
                    quals.append(cont.strip())
                    i += 1
                else:
                    break
            qual_blob = " ".join(quals)
            g_match = re.search(r'/gene="([^"]+)"', qual_blob, re.IGNORECASE)
            gene_name = g_match.group(1) if g_match else ""
            cs_match = re.search(r"/codon_start=(\d+)", qual_blob)
            codon_start = int(cs_match.group(1)) if cs_match else 1
            is_pseudo = ("/pseudo" in qual_blob.lower()) or ("<" in loc_str) or (">" in loc_str)
            is_mane = (
                ("mane select" in qual_blob.lower())
                or ("isoform a" in qual_blob.lower())
                or ("isoform 1 " in qual_blob.lower())
            )
            is_comp = "complement(" in loc_str
            intervals = [
                (int(a), int(b))
                for a, b in re.findall(r"(\d+)\.\.(\d+)", loc_str)
                if int(b) >= int(a) and int(b) <= len(full_seq)
            ]
            if intervals:
                ex_l, in_l, un_s, ex_c, in_c = _spliced_from_intervals(intervals, is_comp, codon_start)
                sp_s = "".join(ex_l)
                cds_len = len(sp_s)
                div3 = (cds_len % 3 == 0 and cds_len >= 60)
                n_cod = cds_len // 3
                cods = [sp_s[3 * k:3 * k + 3] for k in range(n_cod)] if div3 else []
                internal_stops = sum(1 for c in cods[:-1] if c in stop_set) if cods else 999
                starts_atg = bool(cods and cods[0] == "ATG")
                ends_stop = bool(cods and cods[-1] in stop_set)
                clean_orf = bool(div3 and internal_stops == 0 and starts_atg and ends_stop)
                cds_candidates.append({
                    "gene": gene_name,
                    "is_comp": is_comp,
                    "is_pseudo": is_pseudo,
                    "is_mane": is_mane,
                    "clean_orf": clean_orf,
                    "div3": div3,
                    "internal_stops": internal_stops,
                    "cds_len": cds_len,
                    "n_exons": len(intervals),
                    "exons_str": ex_l,
                    "introns_str": in_l,
                    "spliced_cds_str": sp_s,
                    "full_unspliced_str": un_s,
                    "exon_coords": ex_c,
                    "intron_coords": in_c,
                })
        else:
            i += 1

    if not cds_candidates:
        raise ValueError("No CDS features found in GenBank record")

    matched = [
        c for c in cds_candidates
        if target_gene and c["gene"].upper() == target_gene.upper()
    ]
    pool = matched if matched else cds_candidates

    chosen = max(
        pool,
        key=lambda c: (
            not c["is_pseudo"],
            c["clean_orf"],
            c["div3"] and c["internal_stops"] == 0,
            c["n_exons"] > 1,
            c["is_mane"],
            c["cds_len"],
        ),
    )

    exons_str = chosen["exons_str"]
    introns_str = chosen["introns_str"]
    spliced_cds_str = chosen["spliced_cds_str"]
    full_unspliced_str = chosen["full_unspliced_str"]
    introns_concat_str = "".join(introns_str)

    return {
        "gene": chosen["gene"] or (target_gene or "UNKNOWN"),
        "is_complement": bool(chosen["is_comp"]),
        "n_coding_exons": len(exons_str),
        "n_coding_introns": len(introns_str),
        "exon_lengths_nt": [len(x) for x in exons_str],
        "intron_lengths_nt": [len(x) for x in introns_str],
        "spliced_cds_len_nt": len(spliced_cds_str),
        "full_unspliced_len_nt": len(full_unspliced_str),
        "introns_total_len_nt": len(introns_concat_str),
        "intron_to_exon_ratio": (
            float(len(introns_concat_str)) / float(max(1, len(spliced_cds_str)))
        ),
        "exons_str": exons_str,
        "introns_str": introns_str,
        "spliced_cds_str": spliced_cds_str,
        "full_unspliced_str": full_unspliced_str,
        "introns_concat_str": introns_concat_str,
        "exon_coords": chosen["exon_coords"],
        "intron_coords": chosen["intron_coords"],
    }


def compute_geometric_form_metrics(codons: np.ndarray, nt: np.ndarray):
    """
    Compute geometric / static-form vs dynamical-flow invariants on a codon sequence:
      - mirror_l1_asymmetry: L1 distance between codon c and revcomp(c) across the 32 mirror pairs
        (0 = exact Chargaff 2nd parity / Perez mirror symmetry; larger = broken symmetry by coding program)
      - mirror_pair_entropy_bits: Shannon entropy over the 32 unordered mirror pairs (max 5.0 bits)
      - codon_entropy_bits: Shannon entropy over the 64 codons (max 6.0 bits)
      - flow_asymmetry_frob: ||P_bigram - P_bigram^T||_F / ||P_bigram||_F (time-reversal asymmetry)
      - stop_density_per_100: internal stop codons per 100 codons
      - gc_total, gc3
    """
    if codons.size == 0:
        return {
            "n_codons": 0,
            "gc_total": float("nan"),
            "gc3": float("nan"),
            "mirror_l1_asymmetry": float("nan"),
            "mirror_pair_entropy_bits": float("nan"),
            "codon_entropy_bits": float("nan"),
            "flow_asymmetry_frob": float("nan"),
            "stop_density_per_100": float("nan"),
        }

    codi = codons.astype(np.int64)
    n_cod = float(codi.size)
    counts64 = np.bincount(codi, minlength=NC).astype(np.float64)
    freqs64 = counts64 / n_cod

    # 1. Mirror-pair L1 asymmetry: 0.5 * sum_c |f(c) - f(rc(c))|
    rc_freqs = freqs64[_RC_PARTNER]
    mirror_l1 = 0.5 * float(np.sum(np.abs(freqs64 - rc_freqs)))

    # 2. Unordered 32-pair Shannon entropy (bits)
    pair_counts = np.bincount(MIRROR_PAIR_ID[codi], minlength=NC).astype(np.float64)
    pair_probs = pair_counts[pair_counts > 0] / n_cod
    mirror_ent = float(-np.sum(pair_probs * np.log2(pair_probs)))

    # 3. 64-codon Shannon entropy (bits)
    p64_pos = freqs64[freqs64 > 0]
    codon_ent = float(-np.sum(p64_pos * np.log2(p64_pos)))

    # 4. Bigram directional flow asymmetry ||B - B^T||_F / ||B||_F
    if codi.size >= 2:
        bg = np.bincount(codi[:-1] * NC + codi[1:], minlength=NC * NC).reshape(NC, NC).astype(np.float64)
        norm_bg = float(np.linalg.norm(bg))
        flow_asym = float(np.linalg.norm(bg - bg.T) / norm_bg) if norm_bg > 0 else 0.0
    else:
        flow_asym = float("nan")

    # 5. Stop codon density per 100 codons
    aa_seq = CODON_AA[codi]
    n_stops = int(np.sum(~NONSTOP_MASK[aa_seq]))
    stop_dens = 100.0 * float(n_stops) / n_cod

    # 6. GC and GC3
    gc3 = float(np.mean(IS_GC3[codi]))
    if nt.size > 0:
        gc_tot = float(np.mean((nt == 1) | (nt == 2)))
    else:
        gc_tot = float("nan")

    return {
        "n_codons": int(codi.size),
        "gc_total": gc_tot,
        "gc3": gc3,
        "mirror_l1_asymmetry": mirror_l1,
        "mirror_pair_entropy_bits": mirror_ent,
        "codon_entropy_bits": codon_ent,
        "flow_asymmetry_frob": flow_asym,
        "stop_density_per_100": stop_dens,
    }


# ===========================================================================
# Built-in self-test (genetic codes, mirror pairs, URL integrity, splicing)
# ===========================================================================

def self_test_invariants(verbose: bool = False) -> bool:
    """
    Verify genetic code tables, 32 mirror pairs, clean FTP/E-utilities URL generation,
    and GenBank intron-exon splicing invariants.
    Always restores genetic code Table 1 before returning.
    """
    url_test = ftp_base_for_accession("GCF_009914755")
    scheme = "https"
    host = "ftp.ncbi.nlm.nih.gov"
    expected_url = f"{scheme}://{host}/genomes/all/GCF/009/914/755/"
    assert url_test == expected_url, f"URL corruption detected: {url_test!r} != {expected_url!r}"
    assert "[" not in url_test and "]" not in url_test and "(" not in url_test

    eutils_test = eutils_efetch_url("NG_017013.2")
    assert "[" not in eutils_test and "]" not in eutils_test and "NG_017013.2" in eutils_test

    assert len(CODONS) == 64
    assert N_MIRROR_PAIRS == 32

    checks = [
        ("ATG", 1, "M"), ("ATA", 1, "I"), ("TGA", 1, "*"), ("TAA", 1, "*"), ("TAG", 1, "*"), ("CTG", 1, "L"),
        ("ATA", 2, "M"), ("AGA", 2, "*"), ("AGG", 2, "*"), ("TGA", 2, "W"), ("TAA", 2, "*"), ("TAG", 2, "*"),
        ("CTT", 3, "T"), ("CTC", 3, "T"), ("CTA", 3, "T"), ("CTG", 3, "T"),
        ("ATA", 3, "M"), ("TGA", 3, "W"), ("TAA", 3, "*"), ("TAG", 3, "*"),
        ("TGA", 4, "W"), ("TAA", 4, "*"), ("TAG", 4, "*"), ("ATA", 4, "I"),
        ("ATA", 5, "M"), ("TGA", 5, "W"), ("AGA", 5, "S"), ("AGG", 5, "S"), ("TAA", 5, "*"), ("TAG", 5, "*"),
        ("TAA", 6, "Q"), ("TAG", 6, "Q"), ("TGA", 6, "*"),
        ("ATG", 11, "M"), ("TGA", 11, "*"), ("TAA", 11, "*"), ("TAG", 11, "*"),
        ("CTG", 12, "S"), ("CTT", 12, "L"), ("TAA", 12, "*"), ("TAG", 12, "*"), ("TGA", 12, "*"),
    ]

    try:
        for codon, table_id, expected_aa in checks:
            tid, status = set_genetic_code(table_id, verbose=False)
            assert tid == table_id and status == "exact_map", f"table {table_id} not exact"
            got_aa = AA_STRING[int(CODON_AA[CODONS.index(codon)])]
            assert got_aa == expected_aa, (
                f"Genetic code check failed for table={table_id} codon={codon}: "
                f"got {got_aa!r}, expected {expected_aa!r}"
            )
            if verbose:
                print(f"  [self-test] table={table_id:2d} {codon} -> {got_aa} (expected {expected_aa}) OK")
    finally:
        set_genetic_code(1, verbose=False)

    # Synthetic GenBank intron-exon splicing unit test (30 + 31 + 32 = 93 nt = 31 codons)
    ex1 = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTC"      # 30 nt (10 codons, starts with ATG)
    in1 = "GTAAGT" + ("ACGT" * 15) + "TTTCAG"    # 72 nt intron (GT..AG)
    ex2 = "GAGCCCCCTCTGAGTCAGGAAACATTTTCAG"      # 31 nt (split codon phase 1)
    in2 = "GTGAGT" + ("TTAA" * 12) + "CTCCAG"    # 60 nt intron (GT..AG)
    ex3 = "CCTATGGAAACTACTTCCTGAAAACAACGTGA"     # 32 nt (completes phase 1 + stop TGA)
    synth_seq = "AAACCC" + ex1 + in1 + ex2 + in2 + ex3 + "GGGTTT"
    c1_s = 7
    c1_e = c1_s + len(ex1) - 1
    c2_s = c1_e + len(in1) + 1
    c2_e = c2_s + len(ex2) - 1
    c3_s = c2_e + len(in2) + 1
    c3_e = c3_s + len(ex3) - 1
    synth_gb = (
        "FEATURES             Location/Qualifiers\n"
        f"     CDS             join({c1_s}..{c1_e},{c2_s}..{c2_e},{c3_s}..{c3_e})\n"
        '                     /gene="SYNTH_TP53"\n'
        "ORIGIN\n"
        f"        1 {synth_seq.lower()}\n//\n"
    )
    parsed = parse_genbank_intron_exon(synth_gb, "SYNTH_TP53")
    assert parsed["n_coding_exons"] == 3
    assert parsed["n_coding_introns"] == 2
    assert parsed["spliced_cds_str"] == (ex1 + ex2 + ex3)
    assert len(parsed["spliced_cds_str"]) == 93
    assert parsed["introns_concat_str"] == (in1 + in2)

    if verbose:
        print(f"  [self-test] FTP URL check   : {url_test} OK")
        print(f"  [self-test] E-utils URL     : {eutils_test} OK")
        print(f"  [self-test] Mirror pairs    : {N_MIRROR_PAIRS} pairs OK")
        print("  [self-test] Splicing parser : 3-exon / 2-intron join(...) reconstruction OK")
        print("  [self-test] ALL INVARIANTS PASSED.")
    return True


# Run silent invariant check on import/startup
self_test_invariants(verbose=False)


# ===========================================================================
# Loader (true frame, one gene one vote, alt code, per-record table audit)
# ===========================================================================

def cds_units(path: Path, max_records: int, seed: int, trim_start_stop: bool, family_split: bool):
    rng = np.random.default_rng(seed)
    raw_records = []
    table_counter = Counter()
    buf = bytearray()
    cur_header = None
    cur_table = 1

    def flush_current():
        nonlocal buf, cur_header
        if cur_header is None or not buf:
            return
        cd, nt = bytes_to_codons(bytes(buf))
        if trim_start_stop and cd.size >= 5:
            cd = cd[1:-1]
            nt = nt[3:-3]
        if cd.size >= 3:
            gene_id = extract_gene_id(cur_header)
            table_counter[cur_table] += 1
            raw_records.append((gene_id, cd.size, cd, nt, cur_header, cur_table))

    with gzip.open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                flush_current()
                buf.clear()
                cur_header = line.decode("utf-8", "ignore").strip()
                cur_table = extract_transl_table(cur_header)
                continue
            if cur_header is None:
                continue
            buf.extend(line.upper())

    flush_current()

    if not raw_records:
        return {}, {}, {}, None, "no_units"

    file_table = table_counter.most_common(1)[0][0] if table_counter else 1
    used_table, code_status = set_genetic_code(file_table, verbose=True)

    # Per-record transl_table consistency audit
    best = {}
    skipped_mixed_table = 0
    for gene_id, sz, cd, nt, hdr, r_tab in raw_records:
        if r_tab != file_table:
            skipped_mixed_table += 1
            continue
        if gene_id not in best or sz > best[gene_id][0]:
            best[gene_id] = (sz, cd, nt, hdr)

    if len(table_counter) > 1:
        print(
            f"  [code] mixed transl_table records detected {dict(table_counter)} -> "
            f"dominant={file_table} ({code_status}), filtered {skipped_mixed_table} discordant record(s)"
        )
    else:
        print(f"  [code] file transl_table={file_table} -> using {used_table} ({code_status})")

    if not best:
        return {}, {}, {}, None, "no_units"

    items = list(best.values())
    if max_records is not None and max_records > 0 and len(items) > max_records:
        idx = rng.choice(len(items), size=max_records, replace=False)
        items = [items[i] for i in idx]

    units, raw, fam = {}, {}, {}
    for i, (ln, cd, nt, header) in enumerate(items):
        label = f"gene_{i:06d}"
        uc = new_counts()
        add_unit(uc, cd, nt)
        units[label] = uc
        raw[label] = cd
        if family_split:
            fam[label] = family_root(extract_gene_id(header))

    return units, raw, fam, used_table, code_status


# ===========================================================================
# Aggregation + likelihood models H0..H3, H3b
# ===========================================================================

def aggregate(labels, units):
    agg = new_counts()
    for lab in labels:
        u = units[lab]
        for k in ("bigram", "marg", "lastnt", "aa_next", "aa_pair", "nt_bigram", "nt_tri"):
            agg[k] += u[k]
        for k in ("n_codons", "n_trans", "n_nt", "gc3_sum", "gc3_n",
                  "acid_count", "basic_count", "aa_count"):
            agg[k] += u[k]
    return agg


def rcp(mat, alpha=ALPHA):
    rs = mat.sum(axis=1, keepdims=True)
    return (mat.astype(np.float64) + alpha) / (rs + alpha * mat.shape[1])


def ll_empirical(test, train):
    logp = np.log(rcp(train["bigram"]) + TINY)
    return float(np.sum(test["bigram"].astype(np.float64) * logp))


def ll_marginal(test, train):
    p = (train["marg"].astype(np.float64) + ALPHA) / (train["marg"].sum() + ALPHA * NC)
    logp = np.log(p + TINY)
    return float(np.sum(test["bigram"].astype(np.float64) * logp[None, :]))


def ll_last_nt(test, train):
    logp = np.log(rcp(train["lastnt"]) + TINY)
    rows = np.arange(NC) % 4
    return float(np.sum(test["bigram"].astype(np.float64) * logp[rows, :]))


def ll_aa_syntax(test, train):
    p_aa = rcp(train["aa_pair"])
    cm = train["marg"].astype(np.float64)
    cma = AA_IND.T * cm[None, :]
    p_c2a2 = (cma + ALPHA) / (cma.sum(axis=1, keepdims=True) + ALPHA * NC)
    a = CODON_AA
    term1 = np.log(p_aa + TINY)[a[:, None], a[None, :]]
    term2 = np.log(p_c2a2 + TINY)[a, np.arange(NC)]
    logp = term1 + term2[None, :]
    return float(np.sum(test["bigram"].astype(np.float64) * logp))


def ll_aa_codonpair(test, train):
    bg = train["bigram"].astype(np.float64) + ALPHA
    Z = bg @ AA_IND
    p_aa = rcp(train["aa_pair"])
    a = CODON_AA
    log_p_aa = np.log(p_aa + TINY)[a[:, None], a[None, :]]
    log_num = np.log(bg)
    log_den = np.log(Z)[np.arange(NC)[:, None], a[None, :]]
    logp = log_p_aa + log_num - log_den
    return float(np.sum(test["bigram"].astype(np.float64) * logp))


def ll_dinucleotide_folded(test, train):
    lpnt = np.log(rcp(train["nt_bigram"]) + TINY)
    n3 = C_N3[:, None]
    m1 = C_N1[None, :]
    m2 = C_N2[None, :]
    m3 = C_N3[None, :]
    logp = lpnt[n3, m1] + lpnt[m1, m2] + lpnt[m2, m3]
    return float(np.sum(test["bigram"].astype(np.float64) * logp))


def ll_nt_trigram_folded(test, train):
    tri = train["nt_tri"].astype(np.float64) + ALPHA
    p = tri / tri.sum(axis=2, keepdims=True)
    lp = np.log(p + TINY)
    r2 = C_N2[:, None]
    r3 = C_N3[:, None]
    c1 = C_N1[None, :]
    c2 = C_N2[None, :]
    c3 = C_N3[None, :]
    logp = lp[r2, r3, c1] + lp[r3, c1, c2] + lp[c1, c2, c3]
    return float(np.sum(test["bigram"].astype(np.float64) * logp))


def deltas(test, train):
    nt = int(test["n_trans"])
    if nt == 0:
        return None
    emp = ll_empirical(test, train) / nt
    h0 = ll_marginal(test, train) / nt
    h1 = ll_last_nt(test, train) / nt
    h2 = ll_aa_syntax(test, train) / nt
    h2b = ll_aa_codonpair(test, train) / nt
    h3 = ll_dinucleotide_folded(test, train) / nt
    h3b = ll_nt_trigram_folded(test, train) / nt
    return {
        "ll_emp": float(emp), "ll_h0": float(h0), "ll_h1": float(h1), "ll_h2": float(h2),
        "ll_h2b": float(h2b), "ll_h3": float(h3), "ll_h3b": float(h3b),
        "d0": float(emp - h0), "d1": float(emp - h1), "d2": float(emp - h2),
        "d2b": float(emp - h2b), "d3": float(emp - h3), "d3b": float(emp - h3b),
    }


# ===========================================================================
# H4 block-shuffle
# ===========================================================================

def block_shuffle_codons(codons, L, rng):
    n = (codons.size // L) * L
    if n < 2 * L:
        return codons
    blocks = codons[:n].reshape(-1, L)
    perm = rng.permutation(blocks.shape[0])
    out = blocks[perm].reshape(-1)
    if codons.size > n:
        out = np.concatenate([out, codons[n:]])
    return out.astype(np.uint8, copy=False)


def ll_emp_on_codons(codons, logp_emp):
    if codons.size < 2:
        return 0.0, 0
    prev = codons[:-1].astype(np.int64)
    nxt = codons[1:].astype(np.int64)
    return float(np.sum(logp_emp[prev, nxt])), int(prev.size)


def h4_deltas(test_labels, raw, train, rng, lengths, reps):
    if reps <= 0:
        return None
    logp_emp = np.log(rcp(train["bigram"]) + TINY)
    real_ll, real_n = 0.0, 0
    for lab in test_labels:
        cd = raw.get(lab)
        if cd is None or cd.size < 2:
            continue
        s, n = ll_emp_on_codons(cd, logp_emp)
        real_ll += s
        real_n += n
    if real_n == 0:
        return None

    ll_real = real_ll / real_n
    out = {}
    for L in lengths:
        sh_ll, sh_n = 0.0, 0
        for _ in range(reps):
            for lab in test_labels:
                cd = raw.get(lab)
                if cd is None or cd.size < 2 * L:
                    continue
                sc = block_shuffle_codons(cd, L, rng)
                s, n = ll_emp_on_codons(sc, logp_emp)
                sh_ll += s
                sh_n += n
        out[L] = None if sh_n == 0 else {
            "ll_shuffled": sh_ll / sh_n,
            "delta_H4": ll_real - sh_ll / sh_n
        }
    return {"ll_real": ll_real, "by_length": out}


def bootstrap_h4(test_labels, raw, train, rng, lengths, reps):
    if reps <= 0:
        return None
    logp_emp = np.log(rcp(train["bigram"]) + TINY)
    label_arr = np.array(test_labels, dtype=object)
    m = len(test_labels)
    out = {L: [] for L in lengths}
    inner_reps = 2

    for _ in range(reps):
        samp = rng.choice(label_arr, size=m, replace=True).tolist()
        real_ll, real_n = 0.0, 0
        for lab in samp:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            s, n = ll_emp_on_codons(cd, logp_emp)
            real_ll += s
            real_n += n
        if real_n == 0:
            continue
        lr = real_ll / real_n
        for L in lengths:
            sh_ll, sh_n = 0.0, 0
            for __ in range(inner_reps):
                for lab in samp:
                    cd = raw.get(lab)
                    if cd is None or cd.size < 2 * L:
                        continue
                    sc = block_shuffle_codons(cd, L, rng)
                    s, n = ll_emp_on_codons(sc, logp_emp)
                    sh_ll += s
                    sh_n += n
            if sh_n > 0:
                out[L].append(lr - sh_ll / sh_n)

    ci = {}
    for L in lengths:
        if out[L]:
            arr = np.array(out[L], dtype=np.float64)
            ci[L] = [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))]
        else:
            ci[L] = [float("nan"), float("nan")]
    return ci


# ===========================================================================
# H5 aa-preserving H2b surrogate (within + corrected transfer)
# ===========================================================================

def build_h2b_sampler(train):
    bg = train["bigram"].astype(np.float64) + ALPHA
    Z = bg @ AA_IND
    sampler = [[None for _ in range(NA)] for _ in range(NC)]

    for p in range(NC):
        for a in range(NA):
            ch = CHANNEL_CODONS[a]
            if ch.size == 0:
                sampler[p][a] = np.array([1.0], dtype=np.float64)
                continue
            denom = float(Z[p, a])
            if not np.isfinite(denom) or denom <= TINY:
                probs = np.ones(ch.size, dtype=np.float64) / float(ch.size)
            else:
                probs = bg[p, ch] / denom
                s = float(probs.sum())
                probs = (
                    np.ones(ch.size, dtype=np.float64) / float(ch.size)
                    if (not np.isfinite(s) or s <= 0.0)
                    else probs / s
                )
            cdf = np.cumsum(probs, dtype=np.float64)
            if cdf.size > 0:
                cdf[-1] = 1.0
            sampler[p][a] = cdf
    return sampler


def sample_h2b_codon(prev_codon, aa, sampler, rng):
    ch = CHANNEL_CODONS[aa]
    cdf = sampler[prev_codon][aa]
    r = float(rng.random())
    if r >= 1.0:
        r = float(np.nextafter(1.0, 0.0))
    j = int(np.searchsorted(cdf, r, side="right"))
    if j >= ch.size:
        j = ch.size - 1
    return int(ch[j])


def generate_h5_surrogate(cd, sampler, rng):
    surr = np.empty(cd.size, dtype=np.uint8)
    surr[0] = cd[0]
    for i in range(1, cd.size):
        true_aa = int(CODON_AA[int(cd[i])])
        surr[i] = sample_h2b_codon(int(surr[i - 1]), true_aa, sampler, rng)
    return surr


def h5_deltas(test_labels, raw, train, rng, reps):
    if reps <= 0:
        return None
    logp_emp = np.log(rcp(train["bigram"]) + TINY)
    sampler = build_h2b_sampler(train)

    real_ll, real_n = 0.0, 0
    for lab in test_labels:
        cd = raw.get(lab)
        if cd is None or cd.size < 2:
            continue
        s, n = ll_emp_on_codons(cd, logp_emp)
        real_ll += s
        real_n += n
    if real_n == 0:
        return None

    surr_ll, surr_n = 0.0, 0
    for _ in range(reps):
        for lab in test_labels:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            s, n = ll_emp_on_codons(generate_h5_surrogate(cd, sampler, rng), logp_emp)
            surr_ll += s
            surr_n += n
    if surr_n == 0:
        return None

    ll_real = real_ll / float(real_n)
    ll_surr = surr_ll / float(surr_n)
    return {
        "delta_H5": float(ll_real - ll_surr),
        "ll_real": float(ll_real),
        "ll_surrogate": float(ll_surr),
        "reps": int(reps),
    }


def h5_bootstrap(test_labels, raw, train, rng, reps, inner_reps: int = 2):
    if reps <= 0:
        return None
    inner_reps = max(1, int(inner_reps))
    logp_emp = np.log(rcp(train["bigram"]) + TINY)
    sampler = build_h2b_sampler(train)
    label_arr = np.array(test_labels, dtype=object)
    m = len(test_labels)
    out = []

    for _ in range(reps):
        samp = rng.choice(label_arr, size=m, replace=True).tolist()
        real_ll, real_n = 0.0, 0
        for lab in samp:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            s, n = ll_emp_on_codons(cd, logp_emp)
            real_ll += s
            real_n += n
        if real_n == 0:
            continue

        surr_ll, surr_n = 0.0, 0
        for __ in range(inner_reps):
            for lab in samp:
                cd = raw.get(lab)
                if cd is None or cd.size < 2:
                    continue
                s, n = ll_emp_on_codons(generate_h5_surrogate(cd, sampler, rng), logp_emp)
                surr_ll += s
                surr_n += n
        if surr_n > 0:
            out.append(real_ll / real_n - surr_ll / surr_n)

    if not out:
        return None
    arr = np.array(out, dtype=np.float64)
    return [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))]


def h5_transfer_deltas(test_labels, raw, source_train, target_train, rng, reps):
    """Corrected H5-transfer: surrogate from SOURCE H2b, real+surrogate scored under TARGET."""
    if reps <= 0 or target_train is None or source_train is None:
        return None

    logp_target = np.log(rcp(target_train["bigram"]) + TINY)
    sampler_source = build_h2b_sampler(source_train)

    real_ll, real_n = 0.0, 0
    for lab in test_labels:
        cd = raw.get(lab)
        if cd is None or cd.size < 2:
            continue
        s, n = ll_emp_on_codons(cd, logp_target)
        real_ll += s
        real_n += n
    if real_n == 0:
        return None

    surr_ll, surr_n = 0.0, 0
    for _ in range(reps):
        for lab in test_labels:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            s, n = ll_emp_on_codons(generate_h5_surrogate(cd, sampler_source, rng), logp_target)
            surr_ll += s
            surr_n += n
    if surr_n == 0:
        return None

    ll_real = real_ll / float(real_n)
    ll_surr = surr_ll / float(surr_n)
    return {
        "delta_H5_transfer_target": float(ll_real - ll_surr),
        "ll_real_target": float(ll_real),
        "ll_surrogate_source_under_target": float(ll_surr),
        "reps": int(reps),
    }


# ===========================================================================
# H6 synonymous shuffle (+ bootstrap CI)
# ===========================================================================

def synonymous_shuffle_codons(cd, rng):
    aa = CODON_AA[cd.astype(np.int64)]
    out = cd.copy()
    for a in range(NA):
        pos = np.flatnonzero(aa == a)
        if pos.size > 1:
            vals = out[pos].copy()
            rng.shuffle(vals)
            out[pos] = vals
    return out


def h6_deltas(test_labels, raw, train, rng, reps):
    if reps <= 0:
        return None
    logp_emp = np.log(rcp(train["bigram"]) + TINY)

    real_ll, real_n = 0.0, 0
    for lab in test_labels:
        cd = raw.get(lab)
        if cd is None or cd.size < 2:
            continue
        s, n = ll_emp_on_codons(cd, logp_emp)
        real_ll += s
        real_n += n
    if real_n == 0:
        return None

    sh_ll, sh_n = 0.0, 0
    for _ in range(reps):
        for lab in test_labels:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            s, n = ll_emp_on_codons(synonymous_shuffle_codons(cd, rng), logp_emp)
            sh_ll += s
            sh_n += n
    if sh_n == 0:
        return None

    ll_real = real_ll / float(real_n)
    ll_shuf = sh_ll / float(sh_n)
    return {
        "delta_H6": float(ll_real - ll_shuf),
        "ll_real": float(ll_real),
        "ll_shuffled": float(ll_shuf),
        "reps": int(reps),
    }


def h6_bootstrap(test_labels, raw, train, rng, reps):
    if reps <= 0:
        return None
    logp_emp = np.log(rcp(train["bigram"]) + TINY)
    label_arr = np.array(test_labels, dtype=object)
    m = len(test_labels)
    out = []

    for _ in range(reps):
        samp = rng.choice(label_arr, size=m, replace=True).tolist()
        real_ll, real_n = 0.0, 0
        sh_ll, sh_n = 0.0, 0
        for lab in samp:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            sr, nr = ll_emp_on_codons(cd, logp_emp)
            real_ll += sr
            real_n += nr
            ss, ns = ll_emp_on_codons(synonymous_shuffle_codons(cd, rng), logp_emp)
            sh_ll += ss
            sh_n += ns
        if real_n > 0 and sh_n > 0:
            out.append(real_ll / real_n - sh_ll / sh_n)

    if not out:
        return None
    arr = np.array(out, dtype=np.float64)
    return [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))]


# ===========================================================================
# H7 mirror-pair-preserving shuffle (Perez 2013 Sec.2.2 + bootstrap CI)
# ===========================================================================

def mirror_pair_shuffle_codons(cd, rng):
    cd_i = cd.astype(np.int64)
    pid = MIRROR_PAIR_ID[cd_i]
    out = cd.copy()
    for p in np.unique(pid):
        pos = np.flatnonzero(pid == p)
        if pos.size > 1:
            vals = out[pos].copy()
            rng.shuffle(vals)
            out[pos] = vals
    return out


def h7_deltas(test_labels, raw, train, rng, reps):
    if reps <= 0:
        return None
    logp_emp = np.log(rcp(train["bigram"]) + TINY)

    real_ll, real_n = 0.0, 0
    for lab in test_labels:
        cd = raw.get(lab)
        if cd is None or cd.size < 2:
            continue
        s, n = ll_emp_on_codons(cd, logp_emp)
        real_ll += s
        real_n += n
    if real_n == 0:
        return None

    sh_ll, sh_n = 0.0, 0
    for _ in range(reps):
        for lab in test_labels:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            s, n = ll_emp_on_codons(mirror_pair_shuffle_codons(cd, rng), logp_emp)
            sh_ll += s
            sh_n += n
    if sh_n == 0:
        return None

    ll_real = real_ll / float(real_n)
    ll_shuf = sh_ll / float(sh_n)
    return {
        "delta_H7": float(ll_real - ll_shuf),
        "ll_real": float(ll_real),
        "ll_mirror_shuffled": float(ll_shuf),
        "reps": int(reps),
    }


def h7_bootstrap(test_labels, raw, train, rng, reps):
    if reps <= 0:
        return None
    logp_emp = np.log(rcp(train["bigram"]) + TINY)
    label_arr = np.array(test_labels, dtype=object)
    m = len(test_labels)
    out = []

    for _ in range(reps):
        samp = rng.choice(label_arr, size=m, replace=True).tolist()
        real_ll, real_n = 0.0, 0
        sh_ll, sh_n = 0.0, 0
        for lab in samp:
            cd = raw.get(lab)
            if cd is None or cd.size < 2:
                continue
            sr, nr = ll_emp_on_codons(cd, logp_emp)
            real_ll += sr
            real_n += nr
            ss, ns = ll_emp_on_codons(mirror_pair_shuffle_codons(cd, rng), logp_emp)
            sh_ll += ss
            sh_n += ns
        if real_n > 0 and sh_n > 0:
            out.append(real_ll / real_n - sh_ll / sh_n)

    if not out:
        return None
    arr = np.array(out, dtype=np.float64)
    return [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))]


# ===========================================================================
# Split + bootstrap
# ===========================================================================

def choose_split(labels, units, min_trans, seed, fam=None, family_split=False):
    labels = sorted(labels, key=natural_key)
    eligible = [lab for lab in labels if units[lab]["n_trans"] >= min_trans]
    if len(eligible) < 20:
        eligible = labels
    if not eligible:
        return [], []

    rng = np.random.default_rng(seed)

    if family_split and fam:
        groups = {}
        for lab in eligible:
            groups.setdefault(fam.get(lab, lab), []).append(lab)
        keys = list(groups.keys())
        rng.shuffle(keys)
        k_test = max(1, len(keys) // 5)
        test_keys = set(keys[:k_test])
        test, train = [], []
        for kk in keys:
            (test if kk in test_keys else train).extend(groups[kk])
        if not test or not train:
            test = eligible[-max(1, len(eligible) // 5):]
            tset = set(test)
            train = [l for l in eligible if l not in tset]
        return train, test

    arr = np.array(eligible, dtype=object)
    rng.shuffle(arr)
    eligible = arr.tolist()
    if len(eligible) < 4:
        return [], []

    test = eligible[::5]
    tset = set(test)
    train = [l for l in eligible if l not in tset]
    if not test or not train:
        k = max(1, len(eligible) // 5)
        test = eligible[-k:]
        train = eligible[:-k]
    return train, test


def bootstrap_deltas(test_labels, units, train_agg, reps, seed):
    if reps <= 0:
        return None
    rng = np.random.default_rng(seed)
    label_arr = np.array(test_labels, dtype=object)
    m = len(test_labels)
    out = []

    for _ in range(reps):
        sample = rng.choice(label_arr, size=m, replace=True).tolist()
        d = deltas(aggregate(sample, units), train_agg)
        if d is not None:
            out.append([d["d0"], d["d1"], d["d2"], d["d2b"], d["d3"], d["d3b"]])

    if not out:
        return None

    arr = np.array(out, dtype=np.float64)
    keys = ["d0", "d1", "d2", "d2b", "d3", "d3b"]
    return {
        k: [float(np.percentile(arr[:, i], 2.5)), float(np.percentile(arr[:, i], 97.5))]
        for i, k in enumerate(keys)
    }


# ===========================================================================
# Reporting helpers
# ===========================================================================

def fmt_float(x, digits=6):
    if x is None:
        return "n/a"
    try:
        xf = float(x)
        return "n/a" if not np.isfinite(xf) else f"{xf:+.{digits}f}"
    except Exception:
        return "n/a"


def fmt_ratio(x, digits=3):
    if x is None:
        return "n/a"
    try:
        xf = float(x)
        return "n/a" if not np.isfinite(xf) else f"{xf:.{digits}f}"
    except Exception:
        return "n/a"


def get_h4_value(h4, L):
    if not h4:
        return None
    v = h4.get("by_length", {}).get(L)
    return None if v is None else float(v["delta_H4"])


def get_h4_ci(h4ci, L):
    if not h4ci:
        return None
    v = h4ci.get(L)
    return None if v is None else (float(v[0]), float(v[1]))


def safe_div(a, b):
    try:
        af = float(a)
        bf = float(b)
        if not np.isfinite(af) or not np.isfinite(bf) or abs(bf) <= 1e-12:
            return float("nan")
        return af / bf
    except Exception:
        return float("nan")


def _clean_json(obj):
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(x) for x in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        x = float(obj)
        return x if math.isfinite(x) else None
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, np.ndarray):
        return _clean_json(obj.tolist())
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def atomic_write_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
        fh.flush()
        try:
            os.fsync(fh.fileno())
        except Exception:
            pass
    os.replace(tmp, path)


# ===========================================================================
# Intron vs. Exon & P53 Geometric/Program Audit Engine (v14.6)
# ===========================================================================

def score_sequence_against_genome_model(
    seq_str: str,
    train_agg: dict,
    rng: np.random.Generator,
    trim_ends: bool = False,
    mc_reps: int = 40,
):
    """
    Evaluate a single DNA sequence (e.g. TP53 final exons, TP53 full unspliced gene,
    or TP53 introns) against the genome-wide training model `train_agg`.
    Returns both informational syntax deltas (d0..d3b, H4, H5, H6, H7) and
    geometric form invariants (mirror L1 asymmetry, flow asymmetry, stop density, entropies).
    """
    cd, nt = bytes_to_codons(seq_str.encode("ascii", "ignore"))
    if trim_ends and cd.size >= 5:
        cd = cd[1:-1]
        nt = nt[3:-3]

    uc = new_counts()
    add_unit(uc, cd, nt)
    pt = deltas(uc, train_agg)
    geom = compute_geometric_form_metrics(cd, nt)

    if pt is None or cd.size < 4:
        return {
            "n_nt": int(nt.size),
            "n_codons": int(cd.size),
            "n_trans": int(uc["n_trans"]),
            "point": pt,
            "h4": None,
            "h5": None,
            "h6": None,
            "h7": None,
            "geometric": geom,
        }

    raw_single = {"seq": cd}
    labels_single = ["seq"]

    h4 = h4_deltas(labels_single, raw_single, train_agg, rng, lengths=(3, 10, 30), reps=mc_reps)
    h5 = h5_deltas(labels_single, raw_single, train_agg, rng, reps=mc_reps)
    h6 = h6_deltas(labels_single, raw_single, train_agg, rng, reps=mc_reps)
    h7 = h7_deltas(labels_single, raw_single, train_agg, rng, reps=mc_reps)

    return {
        "n_nt": int(nt.size),
        "n_codons": int(cd.size),
        "n_trans": int(uc["n_trans"]),
        "point": pt,
        "d0": pt["d0"],
        "d1": pt["d1"],
        "d2": pt["d2"],
        "d2b": pt["d2b"],
        "d3": pt["d3"],
        "d3b": pt["d3b"],
        "dH4_3": get_h4_value(h4, 3),
        "dH4_10": get_h4_value(h4, 10),
        "dH4_30": get_h4_value(h4, 30),
        "dH5": float(h5["delta_H5"]) if h5 else None,
        "dH6": float(h6["delta_H6"]) if h6 else None,
        "dH7": float(h7["delta_H7"]) if h7 else None,
        "geometric": geom,
    }


def score_unspliced_3frames(
    seq_str: str,
    train_agg: dict,
    rng: np.random.Generator,
    mc_reps: int = 24,
):
    """
    Evaluate an unspliced sequence across all 3 reading frames (offsets 0, 1, 2)
    so frame shifts inside introns cannot be blamed for the loss of syntax.
    """
    frames = []
    for offset in (0, 1, 2):
        sub = seq_str[offset:]
        res = score_sequence_against_genome_model(
            sub, train_agg, rng, trim_ends=False, mc_reps=max(8, mc_reps // 2)
        )
        res["frame_offset"] = offset
        frames.append(res)

    valid = [f for f in frames if f.get("point") is not None]
    if not valid:
        return {"frames": frames, "best_d3b": None, "mean_d3b": None, "best_d2b": None, "mean_d2b": None}

    return {
        "frames": frames,
        "best_d2b": float(max(f["d2b"] for f in valid)),
        "mean_d2b": float(np.mean([f["d2b"] for f in valid])),
        "best_d3b": float(max(f["d3b"] for f in valid)),
        "mean_d3b": float(np.mean([f["d3b"] for f in valid])),
        "mean_dH5": float(np.mean([f["dH5"] for f in valid if f.get("dH5") is not None])),
        "mean_dH6": float(np.mean([f["dH6"] for f in valid if f.get("dH6") is not None])),
        "mean_dH7": float(np.mean([f["dH7"] for f in valid if f.get("dH7") is not None])),
        "mean_mirror_l1": float(
            np.mean([f["geometric"]["mirror_l1_asymmetry"] for f in valid])
        ),
    }


def analyze_splice_junctions(parsed_gene: dict, train_agg: dict):
    """
    Compare codon transitions across SPLICED Exon[k]->Exon[k+1] junctions
    against UNSPLICED Exon[k]->Intron[k] boundaries.
    """
    exons = parsed_gene["exons_str"]
    introns = parsed_gene["introns_str"]
    if len(exons) < 2 or not introns:
        return None

    logp_emp = np.log(rcp(train_agg["bigram"]) + TINY)
    logp_h0 = np.log(
        (train_agg["marg"].astype(np.float64) + ALPHA) / (train_agg["marg"].sum() + ALPHA * NC) + TINY
    )

    spliced_lls = []
    unspliced_lls = []
    junction_details = []
    phase_counter = Counter()

    cum_nt = 0
    for k in range(len(exons) - 1):
        ex_left = exons[k]
        in_mid = introns[k]
        ex_right = exons[k + 1]
        cum_nt += len(ex_left)
        phase = cum_nt % 3  # 0: between codons, 1: after 1st nt, 2: after 2nd nt
        phase_counter[int(phase)] += 1

        prefix_spliced = "".join(exons[:k + 1])
        suffix_spliced = "".join(exons[k + 1:])
        c_idx = (len(prefix_spliced) // 3) - 1
        full_spliced = prefix_spliced + suffix_spliced
        full_unspliced_local = prefix_spliced + in_mid + suffix_spliced

        spliced_trans_ll = None
        unspliced_trans_ll = None
        spliced_pair_str = "n/a"
        unspliced_pair_str = "n/a"

        if c_idx >= 0 and (c_idx + 2) * 3 <= len(full_spliced):
            c1_str = full_spliced[c_idx * 3:(c_idx + 1) * 3]
            c2_str = full_spliced[(c_idx + 1) * 3:(c_idx + 2) * 3]
            if len(c1_str) == 3 and len(c2_str) == 3 and c1_str in AA_TABLE_STD and c2_str in AA_TABLE_STD:
                i1, i2 = CODONS.index(c1_str), CODONS.index(c2_str)
                spliced_trans_ll = float(logp_emp[i1, i2] - logp_h0[i2])
                spliced_lls.append(spliced_trans_ll)
                spliced_pair_str = f"{c1_str}({AA_TABLE_STD[c1_str]})->{c2_str}({AA_TABLE_STD[c2_str]})"

        if c_idx >= 0 and (c_idx + 2) * 3 <= len(full_unspliced_local):
            u1_str = full_unspliced_local[c_idx * 3:(c_idx + 1) * 3]
            u2_str = full_unspliced_local[(c_idx + 1) * 3:(c_idx + 2) * 3]
            if len(u1_str) == 3 and len(u2_str) == 3 and u1_str in AA_TABLE_STD and u2_str in AA_TABLE_STD:
                j1, j2 = CODONS.index(u1_str), CODONS.index(u2_str)
                unspliced_trans_ll = float(logp_emp[j1, j2] - logp_h0[j2])
                unspliced_lls.append(unspliced_trans_ll)
                unspliced_pair_str = f"{u1_str}({AA_TABLE_STD[u1_str]})->{u2_str}({AA_TABLE_STD[u2_str]})"

        adv = (
            float(spliced_trans_ll - unspliced_trans_ll)
            if (spliced_trans_ll is not None and unspliced_trans_ll is not None)
            else None
        )
        junction_details.append({
            "junction_index": k + 1,
            "exon_left_index": k + 1,
            "exon_right_index": k + 2,
            "phase": int(phase),
            "intron_len_nt": len(in_mid),
            "donor_dinuc": in_mid[:2] if len(in_mid) >= 2 else "",
            "acceptor_dinuc": in_mid[-2:] if len(in_mid) >= 2 else "",
            "spliced_pair": spliced_pair_str,
            "unspliced_pair": unspliced_pair_str,
            "spliced_ll_gain_vs_h0": spliced_trans_ll,
            "unspliced_ll_gain_vs_h0": unspliced_trans_ll,
            "splice_ll_advantage": adv,
        })

    mean_sp = float(np.mean(spliced_lls)) if spliced_lls else float("nan")
    mean_un = float(np.mean(unspliced_lls)) if unspliced_lls else float("nan")
    adv_vals = [j["splice_ll_advantage"] for j in junction_details if j["splice_ll_advantage"] is not None]
    mean_adv = float(np.mean(adv_vals)) if adv_vals else float("nan")

    return {
        "n_junctions": len(junction_details),
        "phase_counts": {str(k): int(v) for k, v in sorted(phase_counter.items())},
        "mean_spliced_ll_gain_vs_h0": mean_sp,
        "mean_unspliced_ll_gain_vs_h0": mean_un,
        "mean_splice_junction_advantage": mean_adv,
        "positive_advantage_junctions": int(sum(1 for x in adv_vals if x > 0.0)),
        "junctions": junction_details,
    }



# ===========================================================================
# Exon-by-Exon & Oncogene Hotspot Mutation Markov Audit (Perez 2018 / 2026)
# ===========================================================================

CANONICAL_ONCOGENE_MUTATIONS = {
    "TP53": [
        {"codon_1based": 80,  "wt": "CCT", "mut": "CTT", "label": "P80L (Conserved control region, Perez 2018 Fig.6)"},
        {"codon_1based": 175, "wt": "CGC", "mut": "CAC", "label": "R175H (Hotspot #1, Exon 5, Perez coef=61/61)"},
        {"codon_1based": 179, "wt": "CAT", "mut": "CGT", "label": "H179R (Hotspot Exon 5, Perez Fig.3)"},
        {"codon_1based": 245, "wt": "GGC", "mut": "AGC", "label": "G245S (Hotspot #2, Exon 7, Perez coef=61/61)"},
        {"codon_1based": 245, "wt": "GGC", "mut": "GAC", "label": "G245D (Hotspot #2 alt, Exon 7, Perez coef=61/61)"},
        {"codon_1based": 248, "wt": "CGG", "mut": "CAG", "label": "R248Q (Hotspot #3, Exon 7, Perez coef=57/61)"},
        {"codon_1based": 248, "wt": "CGG", "mut": "TGG", "label": "R248W (Hotspot #3 alt, Exon 7, Perez coef=57/61)"},
        {"codon_1based": 273, "wt": "CGT", "mut": "CAT", "label": "R273H (Hotspot #4, Exon 8, Perez coef=57/61)"},
        {"codon_1based": 273, "wt": "CGT", "mut": "TGT", "label": "R273C (Hotspot #4 alt, Exon 8, Perez coef=57/61)"},
        {"codon_1based": 282, "wt": "CGG", "mut": "TGG", "label": "R282W (Hotspot Exon 8, Perez Fig.3)"},
    ],
    "KRAS": [
        {"codon_1based": 12, "wt": "GGT", "mut": "GAT", "label": "G12D (Primary KRAS 570-nt driver mutation, Exon 2)"},
        {"codon_1based": 12, "wt": "GGT", "mut": "GTT", "label": "G12V (KRAS 570-nt transversion driver, Exon 2)"},
        {"codon_1based": 12, "wt": "GGT", "mut": "TGT", "label": "G12C (KRAS 570-nt smoking-associated driver, Exon 2)"},
        {"codon_1based": 13, "wt": "GGC", "mut": "GAC", "label": "G13D (KRAS 570-nt codon 13 driver, Exon 2)"},
        {"codon_1based": 61, "wt": "CAA", "mut": "CAT", "label": "Q61H (KRAS 570-nt switch-II driver, Exon 3)"},
    ],
    "NRAS": [
        {"codon_1based": 12, "wt": "GGT", "mut": "GAT", "label": "G12D (NRAS 570-nt codon 12 driver, Exon 2)"},
        {"codon_1based": 13, "wt": "GGT", "mut": "GAT", "label": "G13D (NRAS 570-nt codon 13 driver, Exon 2)"},
        {"codon_1based": 61, "wt": "CAA", "mut": "CGA", "label": "Q61R (Primary NRAS melanoma driver, Exon 3)"},
        {"codon_1based": 61, "wt": "CAA", "mut": "AAA", "label": "Q61K (NRAS melanoma driver, Exon 3)"},
    ],
    "HRAS": [
        {"codon_1based": 12, "wt": "GGC", "mut": "GTC", "label": "G12V (Classic HRAS bladder/sarcoma driver, Exon 2)"},
        {"codon_1based": 13, "wt": "GGT", "mut": "GAT", "label": "G13D (HRAS codon 13 driver, Exon 2)"},
        {"codon_1based": 61, "wt": "CAG", "mut": "CTG", "label": "Q61L (HRAS switch-II driver, Exon 3)"},
    ],
}


def analyze_exon_by_exon_and_mutations(parsed_gene: dict, train_agg: dict, gene_symbol: str):
    """
    1. Score each coding exon in true reading frame (splitting phase-shared boundary codons
       or scoring exact in-frame exon slices) and compare TP53 Exons 5-8 (DNA-binding domain,
       codons 125-300, 86% of cancer mutations) vs flanking exons.
    2. Perform a 61-sense-codon Markovian mutability scan across all internal codon positions
       (analogous to Perez 2018 'genLights' 393x61 scan) and evaluate canonical hotspot mutations.
    """
    exons = parsed_gene["exons_str"]
    spliced_str = parsed_gene["spliced_cds_str"]
    cd_full, _nt_full = bytes_to_codons(spliced_str.encode("ascii", "ignore"))
    if cd_full.size < 5:
        return None

    # 1. Exon-by-exon true-frame coordinate mapping on full spliced CDS
    exon_profiles = []
    cum_nt = 0
    for idx, ex_seq in enumerate(exons):
        start_nt = cum_nt
        end_nt = cum_nt + len(ex_seq)
        cum_nt = end_nt
        # True-frame codons fully or predominantly inside this exon
        c_start = (start_nt + 2) // 3
        c_end = end_nt // 3
        ex_codons = cd_full[c_start:c_end]
        ex_nt_str = spliced_str[c_start * 3:c_end * 3]
        _, ex_nt = bytes_to_codons(ex_nt_str.encode("ascii", "ignore"))
        uc = new_counts()
        add_unit(uc, ex_codons, ex_nt)
        pt = deltas(uc, train_agg) if uc["n_trans"] > 0 else None
        geom = compute_geometric_form_metrics(ex_codons, ex_nt)
        exon_profiles.append({
            "exon_number": idx + 1,
            "len_nt": len(ex_seq),
            "inframe_codons": int(ex_codons.size),
            "codon_range_1based": f"{c_start + 1}..{c_end}" if c_end > c_start else "boundary",
            "gc3": geom.get("gc3"),
            "mirror_l1": geom.get("mirror_l1_asymmetry"),
            "d0": pt["d0"] if pt else None,
            "d2b": pt["d2b"] if pt else None,
            "d3b": pt["d3b"] if pt else None,
        })

    # Domain comparison for TP53: DNA-binding hotspot core (Exons 5-8, codons 125..300) vs Flanking
    domain_summary = None
    if gene_symbol.upper() == "TP53" and cd_full.size >= 300:
        dbd_codons = cd_full[124:300]  # codons 125..300 (1-based)
        flank_codons = np.concatenate([cd_full[1:124], cd_full[300:-1]])
        dbd_nt_str = spliced_str[124 * 3:300 * 3]
        flank_nt_str = spliced_str[3:124 * 3] + spliced_str[300 * 3:-3]
        _, dbd_nt = bytes_to_codons(dbd_nt_str.encode("ascii", "ignore"))
        _, flank_nt = bytes_to_codons(flank_nt_str.encode("ascii", "ignore"))
        uc_dbd, uc_flk = new_counts(), new_counts()
        add_unit(uc_dbd, dbd_codons, dbd_nt)
        add_unit(uc_flk, flank_codons, flank_nt)
        pt_dbd = deltas(uc_dbd, train_agg)
        pt_flk = deltas(uc_flk, train_agg)
        gm_dbd = compute_geometric_form_metrics(dbd_codons, dbd_nt)
        gm_flk = compute_geometric_form_metrics(flank_codons, flank_nt)
        domain_summary = {
            "dbd_exons_5_to_8_codons_125_300": {
                "n_codons": int(dbd_codons.size),
                "gc3": gm_dbd.get("gc3"),
                "mirror_l1": gm_dbd.get("mirror_l1_asymmetry"),
                "d0": pt_dbd["d0"],
                "d2b": pt_dbd["d2b"],
                "d3b": pt_dbd["d3b"],
            },
            "flanking_exons_codons_2_124_and_301_393": {
                "n_codons": int(flank_codons.size),
                "gc3": gm_flk.get("gc3"),
                "mirror_l1": gm_flk.get("mirror_l1_asymmetry"),
                "d0": pt_flk["d0"],
                "d2b": pt_flk["d2b"],
                "d3b": pt_flk["d3b"],
            },
        }

    # 2. 61-sense-codon Markov Mutability Scan & Hotspot Mutation Evaluation
    logp_emp = np.log(rcp(train_agg["bigram"]) + TINY)
    # Folded order-2 nucleotide trigram log-probability matrix logp_h3b (64 x 64)
    tri = train_agg["nt_tri"].astype(np.float64) + ALPHA
    p_tri = tri / tri.sum(axis=2, keepdims=True)
    lp_tri = np.log(p_tri + TINY)
    r2, r3 = C_N2[:, None], C_N3[:, None]
    c1, c2, c3 = C_N1[None, :], C_N2[None, :], C_N3[None, :]
    logp_h3b = lp_tri[r2, r3, c1] + lp_tri[r3, c1, c2] + lp_tri[c1, c2, c3]

    sense_codons = np.array([i for i, ch in enumerate(STANDARD_AA_LIST) if ch != "*"], dtype=np.int64)
    assert sense_codons.size == 61

    # Compute 1..61 Markovian mutability rank for every internal codon position (1-based codon index 2..N-1)
    # Following Perez (2018): rank 1 = wild-type codon is the #1 most optimal among 61 sense codons;
    # rank 61 = all 60 alternative sense codons increase the score relative to wild-type.
    n_total = int(cd_full.size)
    pos_ranks_emp = {}
    pos_ranks_d3b = {}
    for pos_0 in range(1, n_total - 1):
        prev_c = int(cd_full[pos_0 - 1])
        wt_c = int(cd_full[pos_0])
        next_c = int(cd_full[pos_0 + 1])

        scores_emp = logp_emp[prev_c, sense_codons] + logp_emp[sense_codons, next_c]
        scores_d3b = (
            (logp_emp[prev_c, sense_codons] - logp_h3b[prev_c, sense_codons])
            + (logp_emp[sense_codons, next_c] - logp_h3b[sense_codons, next_c])
        )
        wt_emp = float(logp_emp[prev_c, wt_c] + logp_emp[wt_c, next_c])
        wt_d3b = float(
            (logp_emp[prev_c, wt_c] - logp_h3b[prev_c, wt_c])
            + (logp_emp[wt_c, next_c] - logp_h3b[wt_c, next_c])
        )
        # 1 + number of sense codons strictly beating wild-type (1 = most conserved, 61 = most mutable)
        pos_ranks_emp[pos_0 + 1] = int(1 + np.sum(scores_emp > wt_emp + 1e-12))
        pos_ranks_d3b[pos_0 + 1] = int(1 + np.sum(scores_d3b > wt_d3b + 1e-12))

    hotspot_results = []
    for mut_spec in CANONICAL_ONCOGENE_MUTATIONS.get(gene_symbol.upper(), []):
        pos_1 = int(mut_spec["codon_1based"])
        pos_0 = pos_1 - 1
        if pos_0 <= 0 or pos_0 >= n_total - 1:
            continue
        prev_c = int(cd_full[pos_0 - 1])
        wt_c = int(cd_full[pos_0])
        next_c = int(cd_full[pos_0 + 1])
        actual_wt_str = CODONS[wt_c]
        mut_str = mut_spec["mut"].upper()
        mut_c = CODONS.index(mut_str)

        wt_aa = STANDARD_AA_LIST[wt_c]
        mut_aa = STANDARD_AA_LIST[mut_c]

        # Local 2-transition window log-likelihoods
        ll_emp_wt = float(logp_emp[prev_c, wt_c] + logp_emp[wt_c, next_c])
        ll_emp_mut = float(logp_emp[prev_c, mut_c] + logp_emp[mut_c, next_c])
        delta_ll_emp = float(ll_emp_mut - ll_emp_wt)

        d3b_local_wt = float(
            (logp_emp[prev_c, wt_c] - logp_h3b[prev_c, wt_c])
            + (logp_emp[wt_c, next_c] - logp_h3b[wt_c, next_c])
        )
        d3b_local_mut = float(
            (logp_emp[prev_c, mut_c] - logp_h3b[prev_c, mut_c])
            + (logp_emp[mut_c, next_c] - logp_h3b[mut_c, next_c])
        )
        delta_d3b_local = float(d3b_local_mut - d3b_local_wt)

        # Whole-gene spliced CDS delta_d3b shift from this single point mutation
        cd_mut = cd_full.copy()
        cd_mut[pos_0] = mut_c
        # Trim start/stop as in state 2
        cd_wt_int = cd_full[1:-1]
        cd_mut_int = cd_mut[1:-1]
        n_tr = int(cd_wt_int.size - 1)
        gene_d3b_shift = float(delta_d3b_local / max(1, n_tr))

        hotspot_results.append({
            "codon_1based": pos_1,
            "wt_codon": actual_wt_str,
            "wt_aa": wt_aa,
            "mut_codon": mut_str,
            "mut_aa": mut_aa,
            "context_window": f"{CODONS[prev_c]}-{actual_wt_str}->{mut_str}-{CODONS[next_c]}",
            "markov_mutability_rank_61": pos_ranks_emp.get(pos_1),
            "d3b_mutability_rank_61": pos_ranks_d3b.get(pos_1),
            "ll_emp_wt_2trans": ll_emp_wt,
            "ll_emp_mut_2trans": ll_emp_mut,
            "delta_ll_emp_2trans": delta_ll_emp,
            "d3b_local_wt_2trans": d3b_local_wt,
            "d3b_local_mut_2trans": d3b_local_mut,
            "delta_d3b_local_2trans": delta_d3b_local,
            "whole_cds_d3b_shift": gene_d3b_shift,
            "label": mut_spec["label"],
        })

    return {
        "exon_profiles": exon_profiles,
        "domain_summary": domain_summary,
        "hotspot_mutations": hotspot_results,
    }


def progressive_splicing_trajectory(
    parsed_gene: dict,
    train_agg: dict,
    rng: np.random.Generator,
    mc_reps: int = 20,
):
    """
    Co-transcriptional progressive splicing trajectory (5' -> 3'):
      Step 0: 0 introns spliced out (100% unspliced pre-mRNA)
      Step k: first k coding introns excised, remaining introns retained
      Step M: all M coding introns excised (100% mature spliced CDS)
    """
    exons = parsed_gene["exons_str"]
    introns = parsed_gene["introns_str"]
    m_introns = len(introns)
    if m_introns == 0:
        return []

    steps = []
    for k in range(m_introns + 1):
        parts = []
        for idx, ex in enumerate(exons):
            parts.append(ex)
            if idx < m_introns and idx >= k:
                parts.append(introns[idx])
        seq_k = "".join(parts)
        is_mature = (k == m_introns)
        sc = score_sequence_against_genome_model(
            seq_k,
            train_agg,
            rng,
            trim_ends=is_mature,
            mc_reps=max(10, int(mc_reps)),
        )
        geom = sc.get("geometric", {})
        steps.append({
            "step": int(k),
            "introns_removed": int(k),
            "introns_remaining": int(m_introns - k),
            "spliced_fraction": float(k) / float(m_introns),
            "len_nt": len(seq_k),
            "n_codons": sc.get("n_codons", 0),
            "d0": sc.get("d0"),
            "d2": sc.get("d2"),
            "d2b": sc.get("d2b"),
            "d3": sc.get("d3"),
            "d3b": sc.get("d3b"),
            "dH4_3": sc.get("dH4_3"),
            "dH5": sc.get("dH5"),
            "dH6": sc.get("dH6"),
            "dH7": sc.get("dH7"),
            "mirror_l1_asymmetry": geom.get("mirror_l1_asymmetry"),
            "mirror_pair_entropy_bits": geom.get("mirror_pair_entropy_bits"),
            "flow_asymmetry_frob": geom.get("flow_asymmetry_frob"),
            "stop_density_per_100": geom.get("stop_density_per_100"),
            "gc3": geom.get("gc3"),
        })
    return steps


def fetch_or_load_genbank_gene(
    gene_symbol: str,
    acc: str,
    gb_cache_dir: Path,
    timeout: int = 90,
    force: bool = False,
) -> str:
    """
    Fetch NCBI RefSeqGene GenBank flatfile via E-utilities efetch and cache locally.
    """
    gb_cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = gb_cache_dir / f"{gene_symbol}_{acc}.gb"
    if cache_file.exists() and cache_file.stat().st_size > 512 and not force:
        txt = cache_file.read_text(encoding="utf-8", errors="ignore")
        if "ORIGIN" in txt and "CDS" in txt:
            return txt

    url = eutils_efetch_url(acc)
    print(f"    [GenBank] fetching {gene_symbol} ({acc}) from NCBI E-utilities...")
    try:
        raw_bytes = http_get(url, timeout=timeout, attempts=3, quiet=False)
    except Exception:
        base_acc = acc.replace("NZ_", "").split(".")[0]
        print(f"    [GenBank] retrying {gene_symbol} with versionless accession ({base_acc})...")
        raw_bytes = http_get(eutils_efetch_url(base_acc), timeout=timeout, attempts=3, quiet=False)
    txt = raw_bytes.decode("utf-8", errors="ignore")
    if "ORIGIN" not in txt or "CDS" not in txt:
        raise RuntimeError(f"Fetched GenBank record for {gene_symbol} ({acc}) lacks CDS/ORIGIN")
    atomic_write_text(cache_file, txt)
    time.sleep(0.4)  # polite NCBI rate limit
    return txt


def load_or_build_human_reference_train(
    cds_dir: Path,
    max_records: int,
    seed: int,
    timeout: int = 180,
):
    """
    Ensure Human_T2T (or fallback mammal/vertebrate) CDS model is loaded to serve
    as the genome-wide reference syntax matrix for the Intron/Exon & P53 audit.
    """
    set_genetic_code(1, verbose=False)
    human_info = SPECIES_INFO["Human_T2T"]
    acc = human_info["acc"]
    ex = existing_cds_file(cds_dir, acc)
    if ex is None:
        print("[P53-Audit] Human_T2T CDS not yet in cache; downloading reference CDS...")
        prepare_cds_cache(
            cds_dir,
            force=False,
            timeout=timeout,
            include_mito=False,
            include_virus=False,
            include_extended=False,
            only=["Human_T2T"],
        )
        ex = existing_cds_file(cds_dir, acc)

    if ex is None:
        all_files = collect_cds_files(cds_dir)
        if all_files:
            ex = all_files[0]

    if ex is None:
        raise RuntimeError(
            "Could not obtain Human_T2T reference CDS file for genome-wide syntax scoring."
        )

    print(f"[P53-Audit] Building reference syntax model from {ex.name} (max_records={max_records})...")
    units, raw, fam, used_table, code_status = cds_units(
        ex,
        max_records=max_records,
        seed=seed,
        trim_start_stop=True,
        family_split=True,
    )
    set_genetic_code(1, verbose=False)
    train_labels, test_labels = choose_split(
        list(units.keys()), units, min_trans=30, seed=seed, fam=fam, family_split=True
    )
    train_agg = aggregate(train_labels if train_labels else list(units.keys()), units)
    return train_agg, ex.name, len(train_labels)


def evaluate_single_gene_intron_exon(
    gene_symbol: str,
    acc: str,
    role: str,
    gb_text: str,
    train_agg: dict,
    seed: int,
    mc_reps: int = 40,
    run_trajectory: bool = False,
):
    rng = np.random.default_rng(seed)
    parsed = parse_genbank_intron_exon(gb_text, target_gene=gene_symbol)

    exons_score = score_sequence_against_genome_model(
        parsed["spliced_cds_str"],
        train_agg,
        rng,
        trim_ends=True,
        mc_reps=mc_reps,
    )
    full_f0_score = score_sequence_against_genome_model(
        parsed["full_unspliced_str"],
        train_agg,
        rng,
        trim_ends=False,
        mc_reps=mc_reps,
    )
    full_3f_score = score_unspliced_3frames(
        parsed["full_unspliced_str"],
        train_agg,
        rng,
        mc_reps=max(16, mc_reps // 2),
    )
    introns_score = score_sequence_against_genome_model(
        parsed["introns_concat_str"],
        train_agg,
        rng,
        trim_ends=False,
        mc_reps=mc_reps,
    )
    junction_audit = analyze_splice_junctions(parsed, train_agg)
    exon_mut_audit = analyze_exon_by_exon_and_mutations(parsed, train_agg, gene_symbol=gene_symbol)
    trajectory = (
        progressive_splicing_trajectory(parsed, train_agg, rng, mc_reps=max(16, mc_reps // 2))
        if run_trajectory
        else None
    )

    ex_d2b = exons_score.get("d2b")
    in_d2b = introns_score.get("d2b")
    f0_d2b = full_f0_score.get("d2b")

    ex_d3b = exons_score.get("d3b")
    in_d3b = introns_score.get("d3b")
    f0_d3b = full_f0_score.get("d3b")
    f3_best_d3b = full_3f_score.get("best_d3b")

    ex_mir = exons_score["geometric"]["mirror_l1_asymmetry"]
    in_mir = introns_score["geometric"]["mirror_l1_asymmetry"]
    f0_mir = full_f0_score["geometric"]["mirror_l1_asymmetry"]

    return {
        "gene": gene_symbol,
        "accession": acc,
        "role": role,
        "structure": {
            "is_complement": parsed["is_complement"],
            "n_coding_exons": parsed["n_coding_exons"],
            "n_coding_introns": parsed["n_coding_introns"],
            "spliced_cds_len_nt": parsed["spliced_cds_len_nt"],
            "full_unspliced_len_nt": parsed["full_unspliced_len_nt"],
            "introns_total_len_nt": parsed["introns_total_len_nt"],
            "intron_to_exon_ratio": parsed["intron_to_exon_ratio"],
            "exon_lengths_nt": parsed["exon_lengths_nt"],
            "intron_lengths_nt": parsed["intron_lengths_nt"],
        },
        "states": {
            "final_exons_only": exons_score,
            "full_gene_with_introns_f0": full_f0_score,
            "full_gene_with_introns_3f": full_3f_score,
            "introns_only": introns_score,
        },
        "contrasts": {
            "d2b_exons_minus_introns": (
                float(ex_d2b - in_d2b) if (ex_d2b is not None and in_d2b is not None) else None
            ),
            "d2b_exons_minus_full_f0": (
                float(ex_d2b - f0_d2b) if (ex_d2b is not None and f0_d2b is not None) else None
            ),
            "d3b_exons_minus_introns": (
                float(ex_d3b - in_d3b) if (ex_d3b is not None and in_d3b is not None) else None
            ),
            "d3b_exons_minus_full_f0": (
                float(ex_d3b - f0_d3b) if (ex_d3b is not None and f0_d3b is not None) else None
            ),
            "d3b_exons_minus_full_best3f": (
                float(ex_d3b - f3_best_d3b)
                if (ex_d3b is not None and f3_best_d3b is not None)
                else None
            ),
            "mirror_l1_exons_minus_introns": (
                float(ex_mir - in_mir)
                if (np.isfinite(ex_mir) and np.isfinite(in_mir))
                else None
            ),
            "mirror_l1_exons_minus_full_f0": (
                float(ex_mir - f0_mir)
                if (np.isfinite(ex_mir) and np.isfinite(f0_mir))
                else None
            ),
        },
        "splice_junction_audit": junction_audit,
        "exon_and_mutation_audit": exon_mut_audit,
        "progressive_splicing_trajectory": trajectory,
    }


def write_intron_exon_audit_reports(report: dict, out_dir: Path):
    json_path = out_dir / "p53_intron_exon_audit.json"
    md_path = out_dir / "p53_intron_exon_audit.md"

    atomic_write_text(json_path, json.dumps(_clean_json(report), indent=2, sort_keys=False))

    lines = []
    lines.append("# Intron vs. Exon & TP53 Geometric/Program Audit (v14.6)")
    lines.append("")
    lines.append(f"- **Reference genome model**: `{report['reference_model']['source_file']}` "
                 f"(`train_units={report['reference_model']['train_units']}`)")
    lines.append(f"- **Genes analyzed**: `{len(report['genes'])}`")
    lines.append("")
    lines.append("## 1. Four-State Structural Comparison Across Genes")
    lines.append("")
    lines.append(
        "| Gene | State | Len (nt) | d2b (vs AA-pair) | d3b (vs nt-tri) | dH5 | dH6 | dH7 | "
        "Mirror L1 Asym | Stop/100 |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")

    for g in report["genes"]:
        gsym = g["gene"]
        for st_key, st_label in (
            ("final_exons_only", "**1. Final Exons (CDS)**"),
            ("full_gene_with_introns_f0", "2. Full Unspliced (f0)"),
            ("introns_only", "4. Introns Only"),
        ):
            st = g["states"][st_key]
            geom = st.get("geometric", {})
            lines.append(
                f"| {gsym} | {st_label} | {st.get('n_nt', 0)} | "
                f"{fmt_float(st.get('d2b'), 5)} | {fmt_float(st.get('d3b'), 5)} | "
                f"{fmt_float(st.get('dH5'), 5)} | {fmt_float(st.get('dH6'), 5)} | "
                f"{fmt_float(st.get('dH7'), 5)} | "
                f"{fmt_ratio(geom.get('mirror_l1_asymmetry'), 4)} | "
                f"{fmt_ratio(geom.get('stop_density_per_100'), 2)} |"
            )
        f3 = g["states"]["full_gene_with_introns_3f"]
        lines.append(
            f"| {gsym} | 3. Full Unspliced (3f mean/best) | {g['structure']['full_unspliced_len_nt']} | "
            f"{fmt_float(f3.get('mean_d2b'), 5)} | {fmt_float(f3.get('best_d3b'), 5)} | "
            f"{fmt_float(f3.get('mean_dH5'), 5)} | {fmt_float(f3.get('mean_dH6'), 5)} | "
            f"{fmt_float(f3.get('mean_dH7'), 5)} | "
            f"{fmt_ratio(f3.get('mean_mirror_l1'), 4)} | n/a |"
        )

    lines.append("")
    lines.append("## 2. Co-Transcriptional Progressive Splicing Trajectories (5' -> 3')")
    for g_rec in report["genes"]:
        if not g_rec.get("progressive_splicing_trajectory"):
            continue
        gsym = g_rec["gene"]
        lines.append("")
        lines.append(f"### {gsym} (`{g_rec['accession']}`)")
        lines.append("")
        lines.append(
            "| Step | Introns Removed | Spliced % | Len (nt) | d2b | d3b | dH5 | dH7 | "
            "Mirror L1 Asym | Stop/100 |"
        )
        lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for st in g_rec["progressive_splicing_trajectory"]:
            lines.append(
                f"| {st['step']} | {st['introns_removed']}/{g_rec['structure']['n_coding_introns']} | "
                f"{100.0 * st['spliced_fraction']:.1f}% | {st['len_nt']} | "
                f"{fmt_float(st.get('d2b'), 5)} | {fmt_float(st.get('d3b'), 5)} | "
                f"{fmt_float(st.get('dH5'), 5)} | {fmt_float(st.get('dH7'), 5)} | "
                f"{fmt_ratio(st.get('mirror_l1_asymmetry'), 4)} | "
                f"{fmt_ratio(st.get('stop_density_per_100'), 2)} |"
            )

    lines.append("")
    lines.append("## 3. Splice-Junction Phase Reconstruction Audits")
    for g_rec in report["genes"]:
        sja = g_rec.get("splice_junction_audit")
        if not sja:
            continue
        gsym = g_rec["gene"]
        lines.append("")
        lines.append(f"### {gsym} (`{g_rec['accession']}`)")
        lines.append(
            f"- **Mean spliced Exon[k]->Exon[k+1] LL gain vs H0**: `{fmt_float(sja['mean_spliced_ll_gain_vs_h0'], 5)}`"
        )
        lines.append(
            f"- **Mean unspliced Exon[k]->Intron[k] LL gain vs H0**: `{fmt_float(sja['mean_unspliced_ll_gain_vs_h0'], 5)}`"
        )
        lines.append(
            f"- **Mean splice-junction LL advantage**: `{fmt_float(sja['mean_splice_junction_advantage'], 5)}` "
            f"(`{sja['positive_advantage_junctions']}/{sja['n_junctions']}` junctions positive)"
        )
        lines.append("")
        lines.append(
            "| Junction | Exons | Phase | Intron (nt) | Splice Sites | Spliced Codon Pair | "
            "Unspliced Boundary Pair | Spliced Gain | Unspliced Gain | Advantage |"
        )
        lines.append("|---:|---|---:|---:|---|---|---|---:|---:|---:|")
        for j in sja["junctions"]:
            lines.append(
                f"| {j['junction_index']} | E{j['exon_left_index']}->E{j['exon_right_index']} | "
                f"{j['phase']} | {j['intron_len_nt']} | {j['donor_dinuc']}..{j['acceptor_dinuc']} | "
                f"`{j['spliced_pair']}` | `{j['unspliced_pair']}` | "
                f"{fmt_float(j['spliced_ll_gain_vs_h0'], 4)} | "
                f"{fmt_float(j['unspliced_ll_gain_vs_h0'], 4)} | "
                f"{fmt_float(j['splice_ll_advantage'], 4)} |"
            )

    lines.append("")
    lines.append("## 4. Summary Across Panel")
    lines.append("")
    summ = report["summary"]
    lines.append(
        f"- **Exons beat Introns on `d3b` (beyond nt-trigrams)**: "
        f"`{summ['exons_beat_introns_d3b']}/{summ['n_genes']}` genes "
        f"(mean delta = `{fmt_float(summ['mean_d3b_exons_minus_introns'], 5)}` nats/codon)"
    )
    lines.append(
        f"- **Exons beat Full Unspliced Gene (all 3 frames best) on `d3b`**: "
        f"`{summ['exons_beat_full_best3f_d3b']}/{summ['n_genes']}` genes "
        f"(mean delta = `{fmt_float(summ['mean_d3b_exons_minus_full_best3f'], 5)}` nats/codon)"
    )
    lines.append(
        f"- **Exons show stronger Mirror-Symmetry Breaking (Program vs. Form)**: "
        f"`{summ['exons_higher_mirror_asymmetry']}/{summ['n_genes']}` genes "
        f"(mean Mirror L1 exons = `{fmt_ratio(summ['mean_exons_mirror_l1'], 4)}` vs. "
        f"introns = `{fmt_ratio(summ['mean_introns_mirror_l1'], 4)}`)"
    )
    lines.append("")

    atomic_write_text(md_path, "\n".join(lines) + "\n")
    return json_path, md_path


def print_extended_gene_log(res: dict):
    g = res["gene"]
    acc = res["accession"]
    st = res["structure"]
    s_ex = res["states"]["final_exons_only"]
    s_f0 = res["states"]["full_gene_with_introns_f0"]
    s_3f = res["states"]["full_gene_with_introns_3f"]
    s_in = res["states"]["introns_only"]
    con = res["contrasts"]

    g_ex = s_ex.get("geometric", {})
    g_f0 = s_f0.get("geometric", {})
    g_in = s_in.get("geometric", {})

    print("")
    print("#" * 96)
    print(f"EXTENDED GENE AUDIT: {g} ({acc}) — {res['role']}")
    print("#" * 96)
    print(
        f"  Genomic structure : strand={'-' if st['is_complement'] else '+'} | "
        f"coding_exons={st['n_coding_exons']} | coding_introns={st['n_coding_introns']}"
    )
    print(
        f"  Lengths (nt)      : Exons_concat={st['spliced_cds_len_nt']} nt | "
        f"Whole_introns+exons={st['full_unspliced_len_nt']} nt | "
        f"Introns_only={st['introns_total_len_nt']} nt | "
        f"Intron/Exon ratio={st['intron_to_exon_ratio']:.2f}x"
    )
    print(f"  Exon lengths (nt) : {st['exon_lengths_nt']}")
    print(f"  Intron lengths(nt): {st['intron_lengths_nt']}")
    print("-" * 96)
    print("  DIRECT COMPARISON: (1) WHOLE SEQUENCE INTRONS+EXONS vs (2) EXONS CONCATENATED")
    print("-" * 96)
    print(
        "  State                             |   Len(nt) |  Codons |      d0 |      d2 |     d2b |      d3 |     d3b |     dH5 |     dH6 |     dH7"
    )
    print(
        "  ----------------------------------+-----------+---------+---------+---------+---------+---------+---------+---------+---------+--------"
    )
    print(
        f"  2/ Exons concatenated (spliced)   | {s_ex.get('n_nt', 0):9d} | {s_ex.get('n_codons', 0):7d} | "
        f"{fmt_float(s_ex.get('d0'), 4):>7s} | {fmt_float(s_ex.get('d2'), 4):>7s} | {fmt_float(s_ex.get('d2b'), 4):>7s} | "
        f"{fmt_float(s_ex.get('d3'), 4):>7s} | {fmt_float(s_ex.get('d3b'), 4):>7s} | {fmt_float(s_ex.get('dH5'), 4):>7s} | "
        f"{fmt_float(s_ex.get('dH6'), 4):>7s} | {fmt_float(s_ex.get('dH7'), 4):>7s}"
    )
    print(
        f"  1/ Whole introns+exons (frame 0)  | {s_f0.get('n_nt', 0):9d} | {s_f0.get('n_codons', 0):7d} | "
        f"{fmt_float(s_f0.get('d0'), 4):>7s} | {fmt_float(s_f0.get('d2'), 4):>7s} | {fmt_float(s_f0.get('d2b'), 4):>7s} | "
        f"{fmt_float(s_f0.get('d3'), 4):>7s} | {fmt_float(s_f0.get('d3b'), 4):>7s} | {fmt_float(s_f0.get('dH5'), 4):>7s} | "
        f"{fmt_float(s_f0.get('dH6'), 4):>7s} | {fmt_float(s_f0.get('dH7'), 4):>7s}"
    )
    for fr in s_3f.get("frames", []):
        off = fr.get("frame_offset", 0)
        print(
            f"     Whole introns+exons (frame {off})  | {fr.get('n_nt', 0):9d} | {fr.get('n_codons', 0):7d} | "
            f"{fmt_float(fr.get('d0'), 4):>7s} | {fmt_float(fr.get('d2'), 4):>7s} | {fmt_float(fr.get('d2b'), 4):>7s} | "
            f"{fmt_float(fr.get('d3'), 4):>7s} | {fmt_float(fr.get('d3b'), 4):>7s} | {fmt_float(fr.get('dH5'), 4):>7s} | "
            f"{fmt_float(fr.get('dH6'), 4):>7s} | {fmt_float(fr.get('dH7'), 4):>7s}"
        )
    print(
        f"  4/ Excised introns only           | {s_in.get('n_nt', 0):9d} | {s_in.get('n_codons', 0):7d} | "
        f"{fmt_float(s_in.get('d0'), 4):>7s} | {fmt_float(s_in.get('d2'), 4):>7s} | {fmt_float(s_in.get('d2b'), 4):>7s} | "
        f"{fmt_float(s_in.get('d3'), 4):>7s} | {fmt_float(s_in.get('d3b'), 4):>7s} | {fmt_float(s_in.get('dH5'), 4):>7s} | "
        f"{fmt_float(s_in.get('dH6'), 4):>7s} | {fmt_float(s_in.get('dH7'), 4):>7s}"
    )

    print("-" * 96)
    print("  GEOMETRIC FORM vs. INFORMATIONAL PROGRAM INVARIANTS (PEREZ 32 MIRROR PAIRS & MARKOV FLOW)")
    print("-" * 96)
    print(
        "  State                             |  GC_tot |     GC3 | Mirror_L1_Asym | MirrorEnt(b) | CodonEnt(b) | FlowAsym_Frob | Stops/100"
    )
    print(
        "  ----------------------------------+---------+---------+----------------+--------------+-------------+---------------+----------"
    )
    for label, gm in (
        ("2/ Exons concatenated (spliced)  ", g_ex),
        ("1/ Whole introns+exons (frame 0) ", g_f0),
        ("4/ Excised introns only          ", g_in),
    ):
        print(
            f"  {label} | {fmt_ratio(gm.get('gc_total'), 4):>7s} | {fmt_ratio(gm.get('gc3'), 4):>7s} | "
            f"{fmt_ratio(gm.get('mirror_l1_asymmetry'), 5):>14s} | {fmt_ratio(gm.get('mirror_pair_entropy_bits'), 4):>12s} | "
            f"{fmt_ratio(gm.get('codon_entropy_bits'), 4):>11s} | {fmt_ratio(gm.get('flow_asymmetry_frob'), 4):>13s} | "
            f"{fmt_ratio(gm.get('stop_density_per_100'), 2):>9s}"
        )

    print("-" * 96)
    print(
        f"  KEY CONTRASTS ({g}): "
        f"Delta_d3b(Exons - Whole_f0) = {fmt_float(con.get('d3b_exons_minus_full_f0'), 5)} | "
        f"Delta_d3b(Exons - Whole_best3f) = {fmt_float(con.get('d3b_exons_minus_full_best3f'), 5)} | "
        f"Delta_d3b(Exons - Introns) = {fmt_float(con.get('d3b_exons_minus_introns'), 5)}"
    )
    print(
        f"  SYMMETRY BREAKING  : "
        f"Delta_MirrorL1(Exons - Whole_f0) = {fmt_float(con.get('mirror_l1_exons_minus_full_f0'), 5)} | "
        f"Delta_MirrorL1(Exons - Introns) = {fmt_float(con.get('mirror_l1_exons_minus_introns'), 5)}"
    )

    sja = res.get("splice_junction_audit")
    if sja and sja.get("junctions"):
        print("-" * 96)
        print(
            f"  SPLICE-JUNCTION PHASE RECONSTRUCTION ({g}): "
            f"mean_spliced_gain={fmt_float(sja['mean_spliced_ll_gain_vs_h0'], 4)} vs "
            f"mean_unspliced_gain={fmt_float(sja['mean_unspliced_ll_gain_vs_h0'], 4)} | "
            f"mean_advantage={fmt_float(sja['mean_splice_junction_advantage'], 4)} "
            f"({sja['positive_advantage_junctions']}/{sja['n_junctions']} positive)"
        )
        print(
            "    Junc | Exons    | Phase | Intron(nt) | Sites  | Spliced Codon Pair   | Unspliced Boundary   | Spliced_LL | Unspl_LL | Advantage"
        )
        for j in sja["junctions"]:
            print(
                f"    {j['junction_index']:4d} | E{j['exon_left_index']:02d}->E{j['exon_right_index']:02d} | "
                f"    {j['phase']} | {j['intron_len_nt']:10d} | {j['donor_dinuc']}..{j['acceptor_dinuc']} | "
                f"{j['spliced_pair']:20s} | {j['unspliced_pair']:20s} | "
                f"{fmt_float(j['spliced_ll_gain_vs_h0'], 4):>10s} | {fmt_float(j['unspliced_ll_gain_vs_h0'], 4):>8s} | "
                f"{fmt_float(j['splice_ll_advantage'], 4):>9s}"
            )

    ema = res.get("exon_and_mutation_audit")
    if ema:
        print("-" * 96)
        print(f"  EXON-BY-EXON MARKOV PROFILE ({g}):")
        print("    Exon | Len(nt) | InFrame_Codons | Codon_Range |     GC3 | MirrorL1 |       d0 |      d2b |      d3b")
        for ep in ema.get("exon_profiles", []):
            print(
                f"    E{ep['exon_number']:02d}  | {ep['len_nt']:7d} | {ep['inframe_codons']:14d} | "
                f"{ep['codon_range_1based']:11s} | {fmt_ratio(ep.get('gc3'), 4):>7s} | "
                f"{fmt_ratio(ep.get('mirror_l1'), 4):>8s} | {fmt_float(ep.get('d0'), 4):>8s} | "
                f"{fmt_float(ep.get('d2b'), 4):>8s} | {fmt_float(ep.get('d3b'), 4):>8s}"
            )
        ds = ema.get("domain_summary")
        if ds:
            dbd = ds["dbd_exons_5_to_8_codons_125_300"]
            flk = ds["flanking_exons_codons_2_124_and_301_393"]
            print(
                f"    [TP53 Hotspot Core Exons 5-8 (codons 125..300, 86% cancer muts)] "
                f"codons={dbd['n_codons']} GC3={dbd['gc3']:.4f} MirrorL1={dbd['mirror_l1']:.4f} "
                f"d2b={fmt_float(dbd['d2b'], 4)} d3b={fmt_float(dbd['d3b'], 4)}"
            )
            print(
                f"    [TP53 Flanking Exons 1-4 & 9-10 (codons 2..124, 301..393)]       "
                f"codons={flk['n_codons']} GC3={flk['gc3']:.4f} MirrorL1={flk['mirror_l1']:.4f} "
                f"d2b={fmt_float(flk['d2b'], 4)} d3b={fmt_float(flk['d3b'], 4)}"
            )
        hms = ema.get("hotspot_mutations", [])
        if hms:
            print("-" * 96)
            print(f"  SINGLE-CODON ONCOGENIC HOTSPOT MUTATION MARKOV AUDIT ({g} — Perez 2018 / 2026):")
            print(
                "    Codon | WT->Mut (AA) | Context Window      | Rank(1..61) |  LL_WT  |  LL_Mut | Delta_LL | d3b_WT  | d3b_Mut | Delta_d3b | Description"
            )
            for hm in hms:
                print(
                    f"    {hm['codon_1based']:5d} | {hm['wt_codon']}->{hm['mut_codon']} ({hm['wt_aa']}->{hm['mut_aa']}) | "
                    f"{hm['context_window']:19s} | {hm['markov_mutability_rank_61']:4d}/61     | "
                    f"{fmt_float(hm['ll_emp_wt_2trans'], 4):>7s} | {fmt_float(hm['ll_emp_mut_2trans'], 4):>7s} | "
                    f"{fmt_float(hm['delta_ll_emp_2trans'], 4):>8s} | {fmt_float(hm['d3b_local_wt_2trans'], 4):>7s} | "
                    f"{fmt_float(hm['d3b_local_mut_2trans'], 4):>7s} | {fmt_float(hm['delta_d3b_local_2trans'], 4):>9s} | {hm['label']}"
                )

    traj = res.get("progressive_splicing_trajectory")
    if traj:
        print("-" * 96)
        print(f"  CO-TRANSCRIPTIONAL PROGRESSIVE SPLICING TRAJECTORY (5' -> 3') FOR {g}:")
        print(
            "    Step | Removed | Spliced% |  Len(nt) |  Codons |      d2b |      d3b |      dH5 |      dH7 | MirrorL1 | FlowAsym | Stops/100"
        )
        for stp in traj:
            print(
                f"    {stp['step']:4d} | {stp['introns_removed']:2d}/{st['n_coding_introns']:2d}   | "
                f"{100.0 * stp['spliced_fraction']:7.1f}% | {stp['len_nt']:8d} | {stp['n_codons']:7d} | "
                f"{fmt_float(stp.get('d2b'), 4):>8s} | {fmt_float(stp.get('d3b'), 4):>8s} | "
                f"{fmt_float(stp.get('dH5'), 4):>8s} | {fmt_float(stp.get('dH7'), 4):>8s} | "
                f"{fmt_ratio(stp.get('mirror_l1_asymmetry'), 4):>8s} | "
                f"{fmt_ratio(stp.get('flow_asymmetry_frob'), 4):>8s} | "
                f"{fmt_ratio(stp.get('stop_density_per_100'), 2):>9s}"
            )
    print("#" * 96)


def run_intron_exon_p53_audit(args, out_dir: Path, cds_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    ext_log_path = out_dir / "p53_extended_audit.log"

    with open(ext_log_path, "w", encoding="utf-8") as ext_fh, TeeContext(ext_fh):
        print("")
        print("=" * 96)
        print("AUTONOMOUS INTRON vs. EXON, SPLICE-JUNCTION & TP53 GEOMETRIC/PROGRAM AUDIT (v14.6 — EXTENDED LOGS)")
        print("=" * 96)

        gb_cache_dir = cds_dir / "refseqgene_gb"

        train_agg, ref_file, n_train_units = load_or_build_human_reference_train(
            cds_dir=cds_dir,
            max_records=int(getattr(args, "p53_ref_max_records", 6000) or 6000),
            seed=int(args.seed),
            timeout=int(args.download_timeout),
        )

        alias_map = {"KIRAS": "KRAS", "K-RAS": "KRAS", "HARAS": "HRAS", "H-RAS": "HRAS", "N-RAS": "NRAS"}
        if getattr(args, "ras_audit", False) and not getattr(args, "intron_exon_genes", None):
            req_genes = ["KRAS", "NRAS", "HRAS"]
        elif getattr(args, "intron_exon_genes", None):
            req_genes = [
                alias_map.get(g.strip().upper(), g.strip().upper())
                for g in str(args.intron_exon_genes).split(",")
                if g.strip()
            ]
        elif getattr(args, "intron_exon_audit", False):
            req_genes = list(DEFAULT_P53_PANEL_GENES)
        else:
            req_genes = ["TP53", "TP63", "TP73", "BRCA1", "HBB", "GAPDH"]

        if getattr(args, "p53_audit", False) and not getattr(args, "intron_exon_genes", None) and not getattr(args, "ras_audit", False):
            if "TP53" not in req_genes:
                req_genes.insert(0, "TP53")

        gene_results = []

        if getattr(args, "custom_gb", None):
            custom_path = Path(args.custom_gb)
            if custom_path.exists():
                print(f"[P53-Audit] Loading custom GenBank file: {custom_path}")
                gb_txt = custom_path.read_text(encoding="utf-8", errors="ignore")
                res = evaluate_single_gene_intron_exon(
                    gene_symbol=custom_path.stem.upper(),
                    acc=custom_path.name,
                    role="Custom local GenBank record",
                    gb_text=gb_txt,
                    train_agg=train_agg,
                    seed=int(args.seed),
                    mc_reps=int(getattr(args, "p53_mc_reps", 40)),
                    run_trajectory=True,
                )
                gene_results.append(res)
                print_extended_gene_log(res)

        for idx, gsym in enumerate(req_genes):
            info = INTRON_EXON_REFSEQGENE_PANEL.get(gsym)
            if info is None:
                print(f"[P53-Audit] WARNING: unknown gene '{gsym}' (not in panel), skipping")
                continue
            acc = info["acc"]
            role = info["role"]
            try:
                gb_txt = fetch_or_load_genbank_gene(
                    gene_symbol=gsym,
                    acc=acc,
                    gb_cache_dir=gb_cache_dir,
                    timeout=int(args.download_timeout),
                    force=bool(args.force_download),
                )
                # Always run progressive trajectory in extended mode for complete logs
                res = evaluate_single_gene_intron_exon(
                    gene_symbol=gsym,
                    acc=acc,
                    role=role,
                    gb_text=gb_txt,
                    train_agg=train_agg,
                    seed=int(args.seed) + 101 * (idx + 1),
                    mc_reps=int(getattr(args, "p53_mc_reps", 40)),
                    run_trajectory=True,
                )
                gene_results.append(res)
                print_extended_gene_log(res)
            except Exception as e:
                print(f"  [P53-Audit] WARNING: failed to process {gsym} ({acc}): {e}")

        if not gene_results:
            print("[P53-Audit] No genes could be evaluated.")
            return None

        d3b_ex_in = [
            g["contrasts"]["d3b_exons_minus_introns"]
            for g in gene_results
            if g["contrasts"]["d3b_exons_minus_introns"] is not None
        ]
        d3b_ex_f3 = [
            g["contrasts"]["d3b_exons_minus_full_best3f"]
            for g in gene_results
            if g["contrasts"]["d3b_exons_minus_full_best3f"] is not None
        ]
        mir_ex = [
            g["states"]["final_exons_only"]["geometric"]["mirror_l1_asymmetry"]
            for g in gene_results
            if np.isfinite(g["states"]["final_exons_only"]["geometric"]["mirror_l1_asymmetry"])
        ]
        mir_in = [
            g["states"]["introns_only"]["geometric"]["mirror_l1_asymmetry"]
            for g in gene_results
            if np.isfinite(g["states"]["introns_only"]["geometric"]["mirror_l1_asymmetry"])
        ]
        mir_diff = [
            g["contrasts"]["mirror_l1_exons_minus_introns"]
            for g in gene_results
            if g["contrasts"]["mirror_l1_exons_minus_introns"] is not None
        ]

        summary = {
            "n_genes": len(gene_results),
            "exons_beat_introns_d3b": int(sum(1 for x in d3b_ex_in if x > 0.0)),
            "mean_d3b_exons_minus_introns": float(np.mean(d3b_ex_in)) if d3b_ex_in else float("nan"),
            "exons_beat_full_best3f_d3b": int(sum(1 for x in d3b_ex_f3 if x > 0.0)),
            "mean_d3b_exons_minus_full_best3f": float(np.mean(d3b_ex_f3)) if d3b_ex_f3 else float("nan"),
            "exons_higher_mirror_asymmetry": int(sum(1 for x in mir_diff if x > 0.0)),
            "mean_exons_mirror_l1": float(np.mean(mir_ex)) if mir_ex else float("nan"),
            "mean_introns_mirror_l1": float(np.mean(mir_in)) if mir_in else float("nan"),
        }

        report = {
            "version": "14.6",
            "reference_model": {
                "source_file": ref_file,
                "train_units": n_train_units,
            },
            "summary": summary,
            "genes": gene_results,
        }

        json_path, md_path = write_intron_exon_audit_reports(report, out_dir)
        print("")
        print("=" * 96)
        print(
            f"[P53-Audit Summary] Exons > Introns on d3b: {summary['exons_beat_introns_d3b']}/{summary['n_genes']} "
            f"(mean +{summary['mean_d3b_exons_minus_introns']:.5f} nats/codon)"
        )
        print(
            f"[P53-Audit Summary] Exons > Whole Unspliced (best of 3 frames) on d3b: "
            f"{summary['exons_beat_full_best3f_d3b']}/{summary['n_genes']} "
            f"(mean +{summary['mean_d3b_exons_minus_full_best3f']:.5f} nats/codon)"
        )
        print(
            f"[P53-Audit Summary] Mirror L1 Asymmetry (Form -> Program): "
            f"Exons={summary['mean_exons_mirror_l1']:.4f} vs Introns={summary['mean_introns_mirror_l1']:.4f} "
            f"({summary['exons_higher_mirror_asymmetry']}/{summary['n_genes']} genes higher in exons)"
        )
        print(f"[P53-Audit Saved] JSON : {json_path}")
        print(f"[P53-Audit Saved] MD   : {md_path}")
        print(f"[P53-Audit Saved] LOG  : {ext_log_path}")
        print("=" * 96)
        return report


# ===========================================================================
# Multi-Species Core Pipeline (v14.1–v14.6)
# ===========================================================================


# ===========================================================================
# FRONTIER BIOSPHERE & SYNTHETIC GENOME FORENSICS AUDIT (v14.7)
#   Suite 1: Giant Viruses (Mimivirus, Megavirus, Pandoravirus, Pithovirus)
#   Suite 2: Deep Subsurface "Hell Boreholes" & 122C Hydrothermal Extremophiles
#   Suite 3: Venter's Synthetic Genomes (JCVI-syn1.0 & JCVI-syn3.0) vs Natural WT
#   Suite 4: Eukaryogenesis Bridge (Asgard Archaea) & Minimal Cell Boundary
# ===========================================================================

FRONTIER_GENOMES_PANEL = {
    # --- SUITE 1: GIANT VIRUSES (4th-Domain Candidates) ---
    "Mimivirus": {
        "suite": "1_giant_viruses",
        "nuccore_acc": "NC_014649.1",
        "clade": "giant_virus",
        "transl_table": 1,
        "desc": "Acanthamoeba polyphaga mimivirus (~1.18 Mb, ~979 CDS; translational apparatus)",
    },
    "Megavirus_chilensis": {
        "suite": "1_giant_viruses",
        "nuccore_acc": "NC_016072.1",
        "clade": "giant_virus",
        "transl_table": 1,
        "desc": "Megavirus chilensis (~1.26 Mb, ~1,120 CDS; giant Mimiviridae member)",
    },
    "Pandoravirus_salinus": {
        "suite": "1_giant_viruses",
        "nuccore_acc": "NC_022098.1",
        "clade": "giant_virus",
        "transl_table": 1,
        "desc": "Pandoravirus salinus (~2.47 Mb, ~2,541 CDS; largest giant virus, high-GC)",
    },
    "Pithovirus_sibericum": {
        "suite": "1_giant_viruses",
        "nuccore_acc": "NC_023423.1",
        "clade": "giant_virus",
        "transl_table": 1,
        "desc": "Pithovirus sibericum (30,000-yr Siberian permafrost giant virus, ~467 CDS)",
    },

    # --- SUITE 2: DEEP SUBSURFACE "HELL BOREHOLES" & ABYSSAL EXTREMOPHILES ---
    "Desulforudis_audaxviator": {
        "suite": "2_deep_biosphere",
        "nuccore_acc": "NC_010424.1",
        "clade": "deep_subsurface_bacterium",
        "transl_table": 11,
        "desc": "Candidatus Desulforudis audaxviator (2.8 km deep gold-mine borehole, U/Th radiolysis ecosystem)",
    },
    "Methanopyrus_kandleri": {
        "suite": "2_deep_biosphere",
        "nuccore_acc": "NC_003551.1",
        "clade": "abyssal_archaeon",
        "transl_table": 11,
        "desc": "Methanopyrus kandleri AV19 (2,000m deep black smoker, 122C upper temperature limit of life)",
    },
    "Aquifex_aeolicus": {
        "suite": "2_deep_biosphere",
        "nuccore_acc": "NC_000918.1",
        "clade": "hydrothermal_bacterium",
        "transl_table": 11,
        "desc": "Aquifex aeolicus VF5 (deep-branching hyperthermophilic chemolithoautotroph, 95C)",
    },
    "Pyrococcus_furiosus": {
        "suite": "2_deep_biosphere",
        "nuccore_acc": "NC_003413.1",
        "clade": "abyssal_archaeon",
        "transl_table": 11,
        "desc": "Pyrococcus furiosus DSM 3638 (100C abyssal hydrothermal vent euryarchaeon)",
    },

    # --- SUITE 3: VENTER'S SYNTHETIC GENOMES vs NATURAL MYCOPLASMA (Table 4: TGA->W) ---
    "Mycoplasma_mycoides_WT": {
        "suite": "3_synthetic_cell",
        "nuccore_acc": "NC_005364.2",
        "clade": "natural_mycoplasma",
        "transl_table": 4,
        "desc": "Natural wild-type Mycoplasma mycoides PG1 (natural reference for Venter's synthetic cell)",
    },
    "JCVI_syn1_0": {
        "suite": "3_synthetic_cell",
        "nuccore_acc": "CP002027.1",
        "clade": "synthetic_cell",
        "transl_table": 4,
        "desc": "Venter's Mycoplasma mycoides JCVI-syn1.0 (1.08 Mb chemically synthesized genome + watermarks)",
    },
    "JCVI_syn3_0": {
        "suite": "3_synthetic_cell",
        "nuccore_acc": "CP014940.1",
        "clade": "synthetic_minimal_cell",
        "transl_table": 4,
        "desc": "Venter's Mycoplasma laboratorium JCVI-syn3.0 (531 kb human-designed minimal genome, 438 CDS)",
    },

    # --- SUITE 4: EUKARYOGENESIS BRIDGE (ASGARD ARCHAEA) & MINIMAL ENDOSYMBIONT ---
    "Prometheoarchaeum_MKD1": {
        "suite": "4_evolution_boundary",
        "nuccore_acc": "NZ_CP042905.1",
        "clade": "asgard_archaeon",
        "transl_table": 11,
        "desc": "Candidatus Prometheoarchaeum syntrophicum MK-D1 (Asgard archaeon; prokaryotic bridge to Eukaryotes)",
    },
    "Carsonella_ruddii": {
        "suite": "4_evolution_boundary",
        "nuccore_acc": "NC_008512.1",
        "clade": "minimal_endosymbiont",
        "transl_table": 11,
        "desc": "Candidatus Carsonella ruddii PV (160 kb smallest cellular genome, 182 CDS, GC3 ~ 0.06)",
    },
}


def eutils_cds_fasta_url(acc: str) -> str:
    scheme = "https"
    host = "eutils.ncbi.nlm.nih.gov"
    path = "entrez/eutils/efetch.fcgi"
    return f"{scheme}://{host}/{path}?db=nuccore&id={acc}&rettype=fasta_cds_na&retmode=text"


def fetch_or_build_frontier_cds_gz(
    name: str,
    info: dict,
    frontier_cache_dir: Path,
    timeout: int = 180,
    force: bool = False,
) -> Path:
    """
    Download CDS FASTA for a frontier nuccore accession via NCBI E-utilities
    (rettype=fasta_cds_na), enforce explicit [transl_table=N] header tags,
    and save as a validated gzipped FASTA archive compatible with cds_units().
    """
    frontier_cache_dir.mkdir(parents=True, exist_ok=True)
    acc = info["nuccore_acc"]
    tid = int(info.get("transl_table", 1))
    dest = frontier_cache_dir / f"{acc}_{name}_cds_from_genomic.fna.gz"

    if dest.exists() and dest.stat().st_size > 1024 and not force:
        return dest

    url = eutils_cds_fasta_url(acc)
    print(f"  [Frontier-Fetch] {name:24s} ({acc}, table={tid}) via NCBI E-utilities...")
    raw_bytes = b""
    try:
        raw_bytes = http_get(url, timeout=timeout, attempts=3, quiet=False)
    except Exception:
        base_acc = acc.replace("NZ_", "").split(".")[0]
        print(f"  [Frontier-Fetch] retrying {name} with versionless accession {base_acc}...")
        raw_bytes = http_get(eutils_cds_fasta_url(base_acc), timeout=timeout, attempts=3, quiet=False)

    txt = raw_bytes.decode("utf-8", errors="ignore").strip()
    if not txt.startswith(">"):
        raise RuntimeError(f"E-utilities fasta_cds_na returned non-FASTA payload for {name} ({acc})")

    # Ensure every FASTA header carries explicit [transl_table=N] so cds_units() uses the exact code
    out_lines = []
    for line in txt.splitlines():
        if line.startswith(">"):
            if "[transl_table=" not in line:
                line = f"{line} [transl_table={tid}]"
        out_lines.append(line)

    gz_tmp = dest.with_suffix(".fna.gz.part")
    with gzip.open(gz_tmp, "wt", encoding="utf-8") as gz_out:
        gz_out.write("\n".join(out_lines) + "\n")
    validate_fasta_gz(gz_tmp)
    gz_tmp.replace(dest)
    time.sleep(0.4)
    return dest


def audit_synthetic_watermarks_and_genes(
    wt_rec: dict,
    syn_name: str,
    syn_acc: str,
    syn_rec: dict,
    gb_cache_dir: Path,
    timeout: int = 180,
    force: bool = False,
):
    """
    Compare Venter's synthetic genome (JCVI-syn1.0 or JCVI-syn3.0) against the
    natural wild-type Mycoplasma mycoides Markov transition model (Table 4):
      1. Score all synthetic CDSs under the wild-type Mycoplasma transition model.
      2. Fetch the full GenBank record (CP002027.1 / CP014940.1), extract long
         intergenic regions (>= 400 nt, which house Venter's synthetic ASCII watermarks
         and assembly junctions) and score them in all 3 reading frames against the
         natural Mycoplasma Markov model.
      3. Identify individual genes in the synthetic genome with the strongest
         syntactic anomalies (lowest d3b) compared to the natural Mycoplasma baseline.
    """
    set_genetic_code(4, verbose=False)
    wt_train = wt_rec["_train_agg"]
    rng = np.random.default_rng(DEFAULT_SEED + 707)

    # 1. Cross-score whole synthetic genome under wild-type Mycoplasma model
    cross_pt = deltas(syn_rec["_test_agg"], wt_train)

    # 2. Extract and score long intergenic / watermark cassettes from GenBank flatfile
    intergenic_watermarks = []
    gene_anomalies = []
    try:
        gb_txt = fetch_or_load_genbank_gene(
            gene_symbol=syn_name,
            acc=syn_acc,
            gb_cache_dir=gb_cache_dir,
            timeout=timeout,
            force=force,
        )
        m_orig = re.search("^ORIGIN.*$", gb_txt, flags=re.MULTILINE)
        if m_orig:
            feat_text = gb_txt[:m_orig.start()]
            orig_body = gb_txt[m_orig.end():]
            end_idx = orig_body.find("//")
            if end_idx >= 0:
                orig_body = orig_body[:end_idx]
            full_seq = "".join(re.findall(r"[a-zA-Z]+", orig_body)).upper()

            # Parse all CDS intervals and gene products to find both intergenic watermarks and individual genes
            cds_spans = []
            lines = feat_text.splitlines()
            i = 0
            n = len(lines)
            while i < n:
                m_c = re.match(r"^\s{5}(CDS|gene|misc_feature)\s+(.*)$", lines[i])
                if m_c:
                    ftype = m_c.group(1)
                    loc_parts = [m_c.group(2).strip()]
                    i += 1
                    while i < n and re.match(r"^\s{21}\S", lines[i]) and not lines[i].lstrip().startswith("/"):
                        loc_parts.append(lines[i].strip())
                        i += 1
                    loc_s = "".join(loc_parts)
                    quals = []
                    while i < n and re.match(r"^\s{21}\S", lines[i]):
                        quals.append(lines[i].strip())
                        i += 1
                    qblob = " ".join(quals)
                    coords = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc_s)]
                    if coords and ftype == "CDS":
                        s_pos = min(a for a, _ in coords)
                        e_pos = max(b for _, b in coords)
                        is_comp = "complement(" in loc_s
                        lt_m = re.search(r'/locus_tag="([^"]+)"', qblob)
                        pr_m = re.search(r'/product="([^"]+)"', qblob)
                        gn_m = re.search(r'/gene="([^"]+)"', qblob)
                        ltag = lt_m.group(1) if lt_m else f"CDS_{s_pos}"
                        prod = pr_m.group(1) if pr_m else (gn_m.group(1) if gn_m else "hypothetical")
                        seq_cds = full_seq[s_pos - 1:e_pos]
                        if is_comp:
                            seq_cds = revcomp_dna_str(seq_cds)
                        cds_spans.append((s_pos, e_pos, ltag, prod, seq_cds))
                else:
                    i += 1

            # Score each individual gene in synthetic cell (>= 90 codons) under WT Mycoplasma model
            for s_pos, e_pos, ltag, prod, seq_cds in cds_spans:
                cd_g, nt_g = bytes_to_codons(seq_cds.encode("ascii", "ignore"))
                if cd_g.size >= 90:
                    cd_g = cd_g[1:-1]
                    nt_g = nt_g[3:-3]
                    uc_g = new_counts()
                    add_unit(uc_g, cd_g, nt_g)
                    pt_g = deltas(uc_g, wt_train)
                    gm_g = compute_geometric_form_metrics(cd_g, nt_g)
                    if pt_g is not None:
                        gene_anomalies.append({
                            "locus_tag": ltag,
                            "coords": f"{s_pos}..{e_pos}",
                            "codons": int(cd_g.size),
                            "gc3": gm_g.get("gc3"),
                            "d0": pt_g["d0"],
                            "d2b": pt_g["d2b"],
                            "d3b": pt_g["d3b"],
                            "product": prod[:48],
                        })

            # Find large intergenic non-coding spans (>= 600 nt) housing synthetic watermarks / cassettes
            cds_spans_sorted = sorted(cds_spans, key=lambda x: x[0])
            prev_end = 1
            for s_pos, e_pos, ltag, _prod, _seq in cds_spans_sorted:
                gap_len = s_pos - prev_end - 1
                if gap_len >= 600:
                    ig_seq = full_seq[prev_end:s_pos - 1]
                    f3 = score_unspliced_3frames(ig_seq, wt_train, rng, mc_reps=16)
                    f0 = f3["frames"][0] if f3.get("frames") else {}
                    gm0 = f0.get("geometric", {})
                    intergenic_watermarks.append({
                        "span": f"{prev_end + 1}..{s_pos - 1}",
                        "before_locus": ltag,
                        "len_nt": gap_len,
                        "gc_total": gm0.get("gc_total"),
                        "mirror_l1": gm0.get("mirror_l1_asymmetry"),
                        "stops_per_100": gm0.get("stop_density_per_100"),
                        "d2b_f0": f0.get("d2b"),
                        "d3b_f0": f0.get("d3b"),
                        "d3b_best3f": f3.get("best_d3b"),
                    })
                prev_end = max(prev_end, e_pos)
    except Exception as e:
        print(f"  [Watermark-Scanner] note: GenBank intergenic scan skipped for {syn_name}: {e}")

    set_genetic_code(1, verbose=False)
    gene_anomalies_sorted_low = sorted(gene_anomalies, key=lambda x: x["d3b"])[:8]
    gene_anomalies_sorted_high = sorted(gene_anomalies, key=lambda x: x["d3b"], reverse=True)[:5]
    intergenic_watermarks_sorted = sorted(intergenic_watermarks, key=lambda x: x["len_nt"], reverse=True)[:10]

    return {
        "synthetic_genome": syn_name,
        "accession": syn_acc,
        "cross_scored_under_wt": cross_pt,
        "n_scored_genes_ge_90cod": len(gene_anomalies),
        "n_genes_negative_d3b": int(sum(1 for g in gene_anomalies if g["d3b"] < 0.0)),
        "lowest_d3b_genes": gene_anomalies_sorted_low,
        "highest_d3b_genes": gene_anomalies_sorted_high,
        "large_intergenic_watermark_candidates": intergenic_watermarks_sorted,
    }


def run_frontier_audit(args, out_dir: Path, cds_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    frontier_log = out_dir / "frontier_audit.log"
    frontier_cache = cds_dir / "frontier_cds"
    gb_cache = cds_dir / "refseqgene_gb"

    suite_filter = getattr(args, "frontier_suite", "all") or "all"

    with open(frontier_log, "w", encoding="utf-8") as log_fh, TeeContext(log_fh):
        print("")
        print("=" * 100)
        print("FRONTIER BIOSPHERE, GIANT VIRUS, DEEP SUBSURFACE & SYNTHETIC CELL AUDIT (v14.7)")
        print("=" * 100)

        cfg = {
            "max_records": int(getattr(args, "max_records", 4000) or 4000),
            "min_trans": 20,
            "trim_start_stop": True,
            "family_split": False,
            "bootstrap_reps": 80,
            "h4_lengths": (3, 10, 29, 30, 100),
            "h4_reps": 10,
            "h4_bootstrap_reps": 24,
            "h5_reps": 8,
            "h5_bootstrap_reps": 20,
            "h5_bootstrap_inner_reps": 2,
            "h5_transfer_reps": 8,
            "h6_reps": 10,
            "h6_bootstrap_reps": 20,
            "h7_reps": 10,
            "h7_bootstrap_reps": 20,
            "min_within_units": 20,
            "min_within_trans": 25000,
            "min_transfer_units": 20,
            "min_transfer_trans": 25000,
        }

        records_by_name = {}
        summary_rows = []

        for idx, (name, info) in enumerate(FRONTIER_GENOMES_PANEL.items(), 1):
            if suite_filter != "all" and suite_filter.lower() not in info["suite"].lower():
                continue

            acc = info["nuccore_acc"]
            SPECIES_PREFIX[acc] = name
            SPECIES_INFO[name] = {
                "acc": acc,
                "known": [acc],
                "clade": info["clade"],
                "optin": True,
            }

            try:
                gz_path = fetch_or_build_frontier_cds_gz(
                    name=name,
                    info=info,
                    frontier_cache_dir=frontier_cache,
                    timeout=int(args.download_timeout),
                    force=bool(args.force_download),
                )
                rec = analyze_species_file(
                    gz_path,
                    cfg=cfg,
                    seed=int(args.seed) + 333 * idx,
                )
                if rec is None:
                    continue

                # Compute whole-genome geometric invariants on held-out test set
                all_test_codons = np.concatenate(
                    [rec["_raw"][lab] for lab in rec["_test_labels"] if lab in rec["_raw"] and rec["_raw"][lab].size > 0]
                )
                geom = compute_geometric_form_metrics(all_test_codons, np.empty(0, dtype=np.int8))
                rec["geometric"] = geom
                rec["suite"] = info["suite"]
                rec["description"] = info["desc"]
                rec["nuccore_acc"] = acc

                pt = rec["point"]
                ci = rec["bootstrap_ci"] or {}
                h4_30 = get_h4_value(rec.get("h4"), 30)
                h4_30_ci = get_h4_ci(rec.get("h4_bootstrap_ci"), 30)
                h5_v = rec["h5"]["delta_H5"] if rec.get("h5") else None
                h5_ci = rec.get("h5_bootstrap_ci")
                h6_v = rec["h6"]["delta_H6"] if rec.get("h6") else None
                h6_ci = rec.get("h6_bootstrap_ci")
                h7_v = rec["h7"]["delta_H7"] if rec.get("h7") else None
                h7_ci = rec.get("h7_bootstrap_ci")

                # Check strict 7-null positive CI stack
                ci_checks = [
                    ci.get("d2b", [-1, -1])[0] > 0,
                    ci.get("d3", [-1, -1])[0] > 0,
                    ci.get("d3b", [-1, -1])[0] > 0,
                    (h4_30_ci is not None and h4_30_ci[0] > 0),
                    (h5_ci is not None and h5_ci[0] > 0),
                    (h6_ci is not None and h6_ci[0] > 0),
                    (h7_ci is not None and h7_ci[0] > 0),
                ]
                n_pos_cis = int(sum(1 for ok in ci_checks if ok))
                if rec["sample_ok"] and n_pos_cis == 7:
                    verdict = "STRICT_7/7_PASS"
                elif n_pos_cis == 7:
                    verdict = "7/7_POS_CI_(SMALL_SPLIT)"
                else:
                    verdict = f"{n_pos_cis}/7_POS_CI"
                rec["verdict"] = verdict

                records_by_name[name] = rec
                summary_rows.append({
                    "name": name,
                    "suite": info["suite"],
                    "acc": acc,
                    "clade": info["clade"],
                    "table": rec["transl_table"],
                    "train_units": rec["train_units"],
                    "test_units": rec["test_units"],
                    "test_trans": rec["test_transitions"],
                    "gc3": rec["gc3"],
                    "mirror_l1": geom.get("mirror_l1_asymmetry"),
                    "flow_asym": geom.get("flow_asymmetry_frob"),
                    "d0": pt["d0"],
                    "d2b": pt["d2b"],
                    "d3": pt["d3"],
                    "d3b": pt["d3b"],
                    "dH4_30": h4_30,
                    "dH5": h5_v,
                    "dH6": h6_v,
                    "dH7": h7_v,
                    "verdict": verdict,
                    "desc": info["desc"],
                })
            except Exception as e:
                print(f"  [Frontier-Audit] WARNING: {name} ({acc}) failed: {e}")

        print("")
        print("=" * 124)
        print("FRONTIER BIOSPHERE SUMMARY TABLE (GIANT VIRUSES, DEEP BOREHOLES, SYNTHETIC CELLS & ASGARD ARCHAEA)")
        print("=" * 124)
        print(
            "  Taxon / Construct        | Suite              | Tab | Trn/Tst | Tst_Trans |    GC3 | MirrL1 |     d2b |     d3b | dH4(30) |     dH5 |     dH6 |     dH7 | Verdict"
        )
        print(
            "  -------------------------+--------------------+-----+---------+-----------+--------+--------+---------+---------+---------+---------+---------+---------+------------------------"
        )
        for r in summary_rows:
            print(
                f"  {r['name']:24s} | {r['suite']:18s} | {r['table']:3d} | "
                f"{r['train_units']:3d}/{r['test_units']:3d} | {r['test_trans']:9d} | "
                f"{r['gc3']:6.4f} | {fmt_ratio(r['mirror_l1'], 4):>6s} | "
                f"{fmt_float(r['d2b'], 4):>7s} | {fmt_float(r['d3b'], 4):>7s} | "
                f"{fmt_float(r['dH4_30'], 4):>7s} | {fmt_float(r['dH5'], 4):>7s} | "
                f"{fmt_float(r['dH6'], 4):>7s} | {fmt_float(r['dH7'], 4):>7s} | {r['verdict']}"
            )
        print("=" * 124)

        # Synthetic Genome Forensics (Venter's JCVI-syn1.0 & JCVI-syn3.0 vs Wild-Type Mycoplasma)
        syn_forensics = []
        if "Mycoplasma_mycoides_WT" in records_by_name:
            wt_rec = records_by_name["Mycoplasma_mycoides_WT"]
            for syn_key in ("JCVI_syn1_0", "JCVI_syn3_0"):
                if syn_key in records_by_name:
                    s_info = FRONTIER_GENOMES_PANEL[syn_key]
                    f_res = audit_synthetic_watermarks_and_genes(
                        wt_rec=wt_rec,
                        syn_name=syn_key,
                        syn_acc=s_info["nuccore_acc"],
                        syn_rec=records_by_name[syn_key],
                        gb_cache_dir=gb_cache,
                        timeout=int(args.download_timeout),
                        force=bool(args.force_download),
                    )
                    syn_forensics.append(f_res)
                    print("")
                    print("#" * 108)
                    print(
                        f"SYNTHETIC GENOME FORENSICS: {syn_key} ({s_info['nuccore_acc']}) SCORED UNDER WILD-TYPE MYCOPLASMA MODEL"
                    )
                    print("#" * 108)
                    cpt = f_res["cross_scored_under_wt"]
                    if cpt:
                        print(
                            f"  Whole-proteome cross-score under WT Mycoplasma : "
                            f"d0={fmt_float(cpt['d0'], 4)} | d2b={fmt_float(cpt['d2b'], 4)} | "
                            f"d3={fmt_float(cpt['d3'], 4)} | d3b={fmt_float(cpt['d3b'], 4)}"
                        )
                    print(
                        f"  Individual genes (>=90 codons) with d3b < 0    : "
                        f"{f_res['n_genes_negative_d3b']}/{f_res['n_scored_genes_ge_90cod']} genes"
                    )
                    if f_res["lowest_d3b_genes"]:
                        print("  Top Syntactically Disrupted / Anomalous Loci in Synthetic Genome (lowest d3b):")
                        print("    Locus_Tag      | Coords         | Codons |    GC3 |      d0 |     d2b |     d3b | Product")
                        for g_an in f_res["lowest_d3b_genes"]:
                            print(
                                f"    {g_an['locus_tag']:14s} | {g_an['coords']:14s} | {g_an['codons']:6d} | "
                                f"{fmt_ratio(g_an['gc3'], 4):>6s} | {fmt_float(g_an['d0'], 4):>7s} | "
                                f"{fmt_float(g_an['d2b'], 4):>7s} | {fmt_float(g_an['d3b'], 4):>7s} | {g_an['product']}"
                            )
                    if f_res["large_intergenic_watermark_candidates"]:
                        print("-" * 108)
                        print(
                            f"  Large Intergenic Cassettes / Synthetic Watermark Regions (>=600 nt) in {syn_key}:"
                        )
                        print(
                            "    Genomic_Span   | Before_Locus   | Len(nt) | GC_tot | MirrorL1 | Stops/100 |  d2b_f0 |  d3b_f0 | d3b_best3f"
                        )
                        for wm in f_res["large_intergenic_watermark_candidates"]:
                            print(
                                f"    {wm['span']:14s} | {wm['before_locus']:14s} | {wm['len_nt']:7d} | "
                                f"{fmt_ratio(wm['gc_total'], 4):>6s} | {fmt_ratio(wm['mirror_l1'], 4):>8s} | "
                                f"{fmt_ratio(wm['stops_per_100'], 2):>9s} | {fmt_float(wm['d2b_f0'], 4):>7s} | "
                                f"{fmt_float(wm['d3b_f0'], 4):>7s} | {fmt_float(wm['d3b_best3f'], 4):>10s}"
                            )
                    print("#" * 108)

        clean_records = []
        for r in records_by_name.values():
            clean_records.append({k: v for k, v in r.items() if not k.startswith("_")})

        out_payload = {
            "version": "14.7",
            "suite_filter": suite_filter,
            "summary_rows": summary_rows,
            "genomes": clean_records,
            "synthetic_forensics": syn_forensics,
        }
        json_out = out_dir / "frontier_audit.json"
        md_out = out_dir / "frontier_audit.md"
        atomic_write_text(json_out, json.dumps(_clean_json(out_payload), indent=2, sort_keys=False))

        md_lines = [
            "# Frontier Biosphere, Giant Virus, Deep Subsurface & Synthetic Cell Audit (v14.7)",
            "",
            "| Taxon / Construct | Suite | Accession | Table | Train/Test | Test Trans | GC3 | Mirror L1 | d2b | d3b | dH4(30) | dH5 | dH6 | dH7 | Verdict |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
        for r in summary_rows:
            md_lines.append(
                f"| **{r['name']}** | `{r['suite']}` | `{r['acc']}` | {r['table']} | "
                f"{r['train_units']}/{r['test_units']} | {r['test_trans']} | {r['gc3']:.4f} | "
                f"{fmt_ratio(r['mirror_l1'], 4)} | {fmt_float(r['d2b'], 4)} | {fmt_float(r['d3b'], 4)} | "
                f"{fmt_float(r['dH4_30'], 4)} | {fmt_float(r['dH5'], 4)} | {fmt_float(r['d6'] if 'd6' in r else r['dH6'], 4)} | "
                f"{fmt_float(r['dH7'], 4)} | `{r['verdict']}` |"
            )
        atomic_write_text(md_out, "\n".join(md_lines) + "\n")
        print(f"\n[Frontier-Audit Saved] JSON : {json_out}")
        print(f"[Frontier-Audit Saved] MD   : {md_out}")
        print(f"[Frontier-Audit Saved] LOG  : {frontier_log}")
        return out_payload


def parse_transfer_pairs(pair_strs, available_species, include_mito=False,
                         include_virus=False, include_extended=False):
    pairs = list(pair_strs)
    if include_mito:
        pairs.extend(OPTIN_H5_TRANSFER_PAIRS.get("Human_mito", []))
    if include_virus:
        pairs.extend(OPTIN_H5_TRANSFER_PAIRS.get("SARS-CoV-2", []))
    if include_extended:
        pairs.extend(EXTENDED_H5_TRANSFER_PAIRS)

    out = []
    seen = set()
    avail = set(available_species)
    for item in pairs:
        if ":" not in item:
            continue
        src, dst = [x.strip() for x in item.split(":", 1)]
        if src in avail and dst in avail and (src, dst) not in seen:
            seen.add((src, dst))
            out.append((src, dst))
    return out


def analyze_species_file(
    path: Path,
    cfg: dict,
    seed: int,
):
    sp = species_of(path)
    clade = clade_of(sp)
    optin = is_optin(sp)

    t0 = time.time()
    print(f"\n[{sp}] ({clade}) loading {path.name} ...")
    units, raw, fam, used_table, code_status = cds_units(
        path,
        max_records=cfg["max_records"],
        seed=seed,
        trim_start_stop=cfg["trim_start_stop"],
        family_split=cfg["family_split"],
    )
    if not units:
        print(f"  [{sp}] WARNING: no valid CDS units loaded, skipping.")
        return None

    train_labels, test_labels = choose_split(
        list(units.keys()),
        units,
        min_trans=cfg["min_trans"],
        seed=seed,
        fam=fam,
        family_split=cfg["family_split"],
    )
    if not train_labels or not test_labels:
        print(f"  [{sp}] WARNING: insufficient units for train/test split, skipping.")
        return None

    train_agg = aggregate(train_labels, units)
    test_agg = aggregate(test_labels, units)
    pt = deltas(test_agg, train_agg)
    if pt is None:
        return None

    rng = np.random.default_rng(seed + 17)
    ci = bootstrap_deltas(
        test_labels, units, train_agg, reps=cfg["bootstrap_reps"], seed=seed + 31
    )

    h4 = h4_deltas(
        test_labels, raw, train_agg, rng, lengths=cfg["h4_lengths"], reps=cfg["h4_reps"]
    )
    h4_ci = bootstrap_h4(
        test_labels, raw, train_agg, rng, lengths=cfg["h4_lengths"], reps=cfg["h4_bootstrap_reps"]
    )

    h5 = h5_deltas(test_labels, raw, train_agg, rng, reps=cfg["h5_reps"])
    h5_ci = h5_bootstrap(
        test_labels,
        raw,
        train_agg,
        rng,
        reps=cfg["h5_bootstrap_reps"],
        inner_reps=cfg["h5_bootstrap_inner_reps"],
    )

    h6 = h6_deltas(test_labels, raw, train_agg, rng, reps=cfg["h6_reps"])
    h6_ci = h6_bootstrap(test_labels, raw, train_agg, rng, reps=cfg["h6_bootstrap_reps"])

    h7 = h7_deltas(test_labels, raw, train_agg, rng, reps=cfg["h7_reps"])
    h7_ci = h7_bootstrap(test_labels, raw, train_agg, rng, reps=cfg["h7_bootstrap_reps"])

    sample_ok = counts_sample_ok(
        len(train_labels),
        len(test_labels),
        train_agg["n_trans"],
        test_agg["n_trans"],
        cfg["min_within_units"],
        cfg["min_within_trans"],
    )

    gc3 = float(test_agg["gc3_sum"]) / float(max(1, test_agg["gc3_n"]))
    acid_frac = float(test_agg["acid_count"]) / float(max(1, test_agg["aa_count"]))
    basic_frac = float(test_agg["basic_count"]) / float(max(1, test_agg["aa_count"]))
    ab_ratio = safe_div(acid_frac, basic_frac)

    elapsed = time.time() - t0
    d_h5 = h5["delta_H5"] if h5 else None
    d_h6 = h6["delta_H6"] if h6 else None
    d_h7 = h7["delta_H7"] if h7 else None

    print(
        f"  [{sp}] train={len(train_labels)} test={len(test_labels)} "
        f"trans={test_agg['n_trans']} GC3={gc3:.3f} A/B={fmt_ratio(ab_ratio, 2)} | "
        f"d0={fmt_float(pt['d0'], 4)} d2b={fmt_float(pt['d2b'], 4)} d3b={fmt_float(pt['d3b'], 4)} "
        f"dH5={fmt_float(d_h5, 4)} dH6={fmt_float(d_h6, 4)} dH7={fmt_float(d_h7, 4)} "
        f"({elapsed:.1f}s)"
    )

    return {
        "species": sp,
        "clade": clade,
        "optin": optin,
        "file": path.name,
        "transl_table": int(used_table),
        "code_status": code_status,
        "train_units": len(train_labels),
        "test_units": len(test_labels),
        "train_transitions": int(train_agg["n_trans"]),
        "test_transitions": int(test_agg["n_trans"]),
        "sample_ok": bool(sample_ok),
        "gc3": gc3,
        "acid_frac": acid_frac,
        "basic_frac": basic_frac,
        "acid_basic_ratio": ab_ratio,
        "point": pt,
        "bootstrap_ci": ci,
        "h4": h4,
        "h4_bootstrap_ci": h4_ci,
        "h5": h5,
        "h5_bootstrap_ci": h5_ci,
        "h6": h6,
        "h6_bootstrap_ci": h6_ci,
        "h7": h7,
        "h7_bootstrap_ci": h7_ci,
        "_train_agg": train_agg,
        "_test_agg": test_agg,
        "_test_labels": test_labels,
        "_raw": raw,
    }


def run_cross_species_transfers(species_records: dict, pairs: list, cfg: dict, seed: int):
    rng = np.random.default_rng(seed + 999)
    transfers = []
    for src, dst in pairs:
        rec_src = species_records[src]
        rec_dst = species_records[dst]

        # Set target species genetic code before evaluating on target
        set_genetic_code(rec_dst["transl_table"], verbose=False)

        pt = deltas(rec_dst["_test_agg"], rec_src["_train_agg"])
        if pt is None:
            continue

        h5_tr = h5_transfer_deltas(
            rec_dst["_test_labels"],
            rec_dst["_raw"],
            rec_src["_train_agg"],
            rec_dst["_train_agg"],
            rng,
            reps=cfg["h5_transfer_reps"],
        )

        small_flag = is_small_sample(
            rec_src,
            rec_dst,
            cfg["min_transfer_units"],
            cfg["min_transfer_trans"],
        )
        code_mismatch = (rec_src["transl_table"] != rec_dst["transl_table"])

        d2b = float(pt["d2b"])
        d3b = float(pt["d3b"])
        dh5_t = float(h5_tr["delta_H5_transfer_target"]) if h5_tr else None
        gap = (d2b - dh5_t) if dh5_t is not None else None

        tier = "strict_clean"
        if small_flag:
            tier = "small_sample"
        elif code_mismatch:
            tier = "code_mismatch"
        elif d2b <= 0.02:
            tier = "low_signal_d2b_le_0p02"
        elif gap is not None and abs(gap) <= 0.012:
            tier = "H5_NEAR_MISS"

        transfers.append({
            "pair": f"{src}:{dst}",
            "source": src,
            "target": dst,
            "source_clade": rec_src["clade"],
            "target_clade": rec_dst["clade"],
            "source_table": rec_src["transl_table"],
            "target_table": rec_dst["transl_table"],
            "code_mismatch": bool(code_mismatch),
            "small_sample": bool(small_flag),
            "tier": tier,
            "d0": pt["d0"],
            "d1": pt["d1"],
            "d2": pt["d2"],
            "d2b": d2b,
            "d3": pt["d3"],
            "d3b": d3b,
            "h5_transfer": h5_tr,
            "d2b_minus_h5_transfer_gap": gap,
        })

    set_genetic_code(1, verbose=False)
    return transfers


def write_stage_reports(stage_name: str, results: dict, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"results_{stage_name}.json"
    md_path = out_dir / f"report_{stage_name}.md"

    atomic_write_text(json_path, json.dumps(_clean_json(results), indent=2, sort_keys=False))

    lines = []
    lines.append(f"# Super Codon Conditional Syntax Report — Stage: `{stage_name}` (v14.6)")
    lines.append("")
    lines.append("## Within-Species Conditional Syntax Summary")
    lines.append("")
    lines.append(
        "| Species | Clade | Table | Train/Test | GC3 | d0 | d2b (AA-pair) | d3b (nt-tri) | "
        "dH4(L=3) | dH5 (Surr) | dH6 (Syn) | dH7 (Mirror) |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")

    for sp_rec in results["species"]:
        pt = sp_rec["point"]
        h4_3 = get_h4_value(sp_rec.get("h4"), 3)
        h5_v = sp_rec["h5"]["delta_H5"] if sp_rec.get("h5") else None
        h6_v = sp_rec["h6"]["delta_H6"] if sp_rec.get("h6") else None
        h7_v = sp_rec["h7"]["delta_H7"] if sp_rec.get("h7") else None
        lines.append(
            f"| {sp_rec['species']} | {sp_rec['clade']} | {sp_rec['transl_table']} | "
            f"{sp_rec['train_units']}/{sp_rec['test_units']} | {sp_rec['gc3']:.3f} | "
            f"{fmt_float(pt['d0'], 5)} | {fmt_float(pt['d2b'], 5)} | {fmt_float(pt['d3b'], 5)} | "
            f"{fmt_float(h4_3, 5)} | {fmt_float(h5_v, 5)} | {fmt_float(h6_v, 5)} | {fmt_float(h7_v, 5)} |"
        )

    if results.get("transfers"):
        lines.append("")
        lines.append("## Cross-Species Transfer Summary")
        lines.append("")
        lines.append(
            "| Pair (Train->Test) | Tier | d0 | d2b | d3b | dH5_transfer | Gap (d2b - dH5_tr) |"
        )
        lines.append("|---|---|---:|---:|---:|---:|---:|")
        for tr in results["transfers"]:
            dh5_t = tr["h5_transfer"]["delta_H5_transfer_target"] if tr.get("h5_transfer") else None
            lines.append(
                f"| `{tr['pair']}` | {tr['tier']} | {fmt_float(tr['d0'], 5)} | "
                f"{fmt_float(tr['d2b'], 5)} | {fmt_float(tr['d3b'], 5)} | "
                f"{fmt_float(dh5_t, 5)} | {fmt_float(tr.get('d2b_minus_h5_transfer_gap'), 5)} |"
            )

    atomic_write_text(md_path, "\n".join(lines) + "\n")
    return json_path, md_path


def run_pipeline_stage(stage_name: str, cfg: dict, files: list, args, out_dir: Path):
    print("\n" + "=" * 92)
    print(f"RUNNING STAGE: {stage_name.upper()} (files={len(files)}, max_records={cfg['max_records']})")
    print("=" * 92)

    species_records = {}
    ordered_species = []

    for idx, fpath in enumerate(files):
        sp = species_of(fpath)
        if not cfg["include_mito"] and clade_of(sp) == "mitochondrion":
            continue
        if not cfg["include_virus"] and clade_of(sp) == "virus":
            continue
        if not cfg["include_extended"] and SPECIES_INFO.get(sp, {}).get("extended", False):
            continue

        rec = analyze_species_file(fpath, cfg=cfg, seed=int(args.seed) + 1000 * (idx + 1))
        if rec is not None:
            species_records[sp] = rec
            ordered_species.append(rec)

    pairs = parse_transfer_pairs(
        cfg["transfer_pairs"],
        available_species=list(species_records.keys()),
        include_mito=cfg["include_mito"],
        include_virus=cfg["include_virus"],
        include_extended=cfg["include_extended"],
    )
    transfers = run_cross_species_transfers(species_records, pairs, cfg=cfg, seed=int(args.seed))

    # Strip internal numpy caches before serialization
    clean_species = []
    for rec in ordered_species:
        c = {k: v for k, v in rec.items() if not k.startswith("_")}
        clean_species.append(c)

    stage_results = {
        "version": "14.6",
        "stage": stage_name,
        "config": {k: v for k, v in cfg.items() if k != "transfer_pairs"},
        "n_species": len(clean_species),
        "species": clean_species,
        "n_transfers": len(transfers),
        "transfers": transfers,
    }

    json_path, md_path = write_stage_reports(stage_name, stage_results, out_dir)
    print(f"\n[{stage_name.upper()}] Saved JSON: {json_path}")
    print(f"[{stage_name.upper()}] Saved MD  : {md_path}")
    return stage_results


# ===========================================================================
# CLI Argument Parser & Entrypoint
# ===========================================================================



# ===========================================================================
# v14.8 ULTIMATE SHOCK SUITE:
#   MODULE A: Viral Trojan & Host-OS Specificity Scanner (Viruses vs Human & E.coli)
#   MODULE B: Human Proteome HGT "Torrent Network" & Domesticated Retrovirus Forensics
#   MODULE C: Endosymbiotic Organelle OS Colonization (Chloroplast & Mito vs Ancestor & Host)
#   MODULE D: Biotech "Max-CAI Codon Optimization" Crash Test (Form vs Program)
# ===========================================================================

VIRAL_TROJAN_PANEL = {
    "Epstein_Barr_Virus_HHV4": {
        "acc": "NC_007605.1",
        "transl_table": 1,
        "true_host": "Human_T2T",
        "role": "Human herpesvirus 4 (latent in human B-cells; lifelong stealth persistence)",
    },
    "Human_Adenovirus_C2": {
        "acc": "AC_000008.1",
        "transl_table": 1,
        "true_host": "Human_T2T",
        "role": "Human respiratory DNA virus (hijacks human nuclear splicing & ribosomes)",
    },
    "Human_Papillomavirus_16": {
        "acc": "NC_001526.4",
        "transl_table": 1,
        "true_host": "Human_T2T",
        "role": "High-risk oncogenic human DNA virus (E6/E7 p53/RB1 inactivator)",
    },
    "Vaccinia_Poxvirus": {
        "acc": "NC_006998.1",
        "transl_table": 1,
        "true_host": "Human_T2T",
        "role": "Large cytoplasmic mammalian DNA poxvirus (~220 ORFs)",
    },
    "Bacteriophage_T4": {
        "acc": "NC_000866.4",
        "transl_table": 11,
        "true_host": "E.coli",
        "role": "Classic E. coli bacteriophage control (bacterial host specificity test)",
    },
    "Mimivirus": {
        "acc": "NC_014649.1",
        "transl_table": 1,
        "true_host": "Eukaryote_Amoeba",
        "role": "Giant DNA virus (~979 CDSs)",
    },
    "Pandoravirus_salinus": {
        "acc": "NC_022098.1",
        "transl_table": 1,
        "true_host": "Eukaryote_Amoeba",
        "role": "Largest Giant DNA virus (~2.47 Mb, high-GC eukaryotic mimic)",
    },
}

HUMAN_CHIMERA_CANDIDATE_GENES = {
    "ERVW-1_Syncytin1": {
        "acc": "NM_014590.4",
        "origin": "Domesticated endogenous retrovirus envelope (human placental trophoblast fusion)",
    },
    "ARC_Neuronal_Capsid": {
        "acc": "NM_015193.5",
        "origin": "Domesticated Ty3/gypsy retrotransposon Gag capsid (synaptic memory RNA transfer)",
    },
    "RAG1_Immune_Recombinase": {
        "acc": "NM_000448.3",
        "origin": "Domesticated Transib DNA transposase (V(D)J adaptive immunity generator)",
    },
    "TFAM_Mito_Factor": {
        "acc": "NM_003201.3",
        "origin": "Nuclear-encoded mitochondrial transcription factor A (endosymbiotic import)",
    },
    "CS_Citrate_Synthase": {
        "acc": "NM_004077.3",
        "origin": "Nuclear-encoded mitochondrial Krebs cycle core enzyme (alphaproteobacterial ancestry)",
    },
}

ORGANELLE_ENDOSYMBIOSIS_PANEL = {
    "Arabidopsis_Chloroplast": {
        "acc": "NC_000932.1",
        "transl_table": 11,
        "desc": "Arabidopsis thaliana plastid genome (85 CDSs; ancient cyanobacterial endosymbiont)",
    },
    "Human_Mitochondrion_NC": {
        "acc": "NC_012920.1",
        "transl_table": 2,
        "desc": "Homo sapiens mitochondrial genome rCRS (13 CDSs; vertebrate mito code Table 2)",
    },
}


def parse_fasta_cds_records(fasta_text: str, trim_ends: bool = True, min_codons: int = 40):
    """
    Parse raw FASTA CDS records (from NCBI fasta_cds_na or local gz) into a list of dicts:
      [{"gene": str, "product": str, "header": str, "codons": np.ndarray, "nt": np.ndarray}, ...]
    """
    records = []
    cur_hdr = None
    buf = []

    def _flush(hdr, seq_str):
        if not hdr or not seq_str:
            return
        cd, nt = bytes_to_codons(seq_str.encode("ascii", "ignore"))
        if trim_ends and cd.size >= 5:
            cd = cd[1:-1]
            nt = nt[3:-3]
        if cd.size < min_codons:
            return
        gene_m = re.search(r"\[gene=([^\]]+)\]", hdr)
        lt_m = re.search(r"\[locus_tag=([^\]]+)\]", hdr)
        pr_m = re.search(r"\[protein=([^\]]+)\]", hdr)
        pid_m = re.search(r"\[protein_id=([^\]]+)\]", hdr)
        gname = (
            gene_m.group(1)
            if gene_m
            else (lt_m.group(1) if lt_m else (pid_m.group(1) if pid_m else hdr.split()[0].lstrip(">")))
        )
        prod = pr_m.group(1) if pr_m else "coding_sequence"
        records.append({
            "gene": gname,
            "product": prod[:52],
            "header": hdr,
            "codons": cd,
            "nt": nt,
        })

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            _flush(cur_hdr, "".join(buf))
            cur_hdr = line
            buf = []
        else:
            buf.append(line.upper())
    _flush(cur_hdr, "".join(buf))

    # Deduplicate by gene identifier (keep longest isoform)
    dedup = {}
    for r in records:
        k = r["gene"]
        if k not in dedup or r["codons"].size > dedup[k]["codons"].size:
            dedup[k] = r
    return list(dedup.values())


def fetch_or_load_nuccore_cds_fasta(name: str, acc: str, cache_dir: Path, timeout: int = 120, force: bool = False) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    f_path = cache_dir / f"{acc}_{name}_cds.fasta"
    if f_path.exists() and f_path.stat().st_size > 256 and not force:
        return f_path.read_text(encoding="utf-8", errors="ignore")

    url = eutils_cds_fasta_url(acc)
    print(f"  [NCBI-Entrez] Fetching CDS FASTA for {name} ({acc})...")
    try:
        raw = http_get(url, timeout=timeout, attempts=3, quiet=False)
    except Exception:
        base_acc = acc.split(".")[0]
        raw = http_get(eutils_cds_fasta_url(base_acc), timeout=timeout, attempts=3, quiet=False)
    txt = raw.decode("utf-8", errors="ignore")
    if not txt.lstrip().startswith(">"):
        raise RuntimeError(f"Empty or invalid CDS FASTA for {name} ({acc})")
    atomic_write_text(f_path, txt)
    time.sleep(0.35)
    return txt


def build_max_cai_sequence(codons: np.ndarray, host_train: dict) -> np.ndarray:
    """
    Standard biotech 'Max-CAI' (Codon Adaptation Index = 1.0) codon optimization:
    Replaces every sense codon with the single most frequent synonymous codon in the host genome,
    preserving 100% of the amino-acid sequence while ignoring conditional Markov transitions.
    """
    marg = host_train["marg"].astype(np.float64)
    best_codon_for_aa = np.zeros(NA, dtype=np.uint8)
    for a in range(NA):
        ch = CHANNEL_CODONS[a]
        if ch.size > 0:
            best_idx = int(np.argmax(marg[ch]))
            best_codon_for_aa[a] = np.uint8(ch[best_idx])
    aa_seq = CODON_AA[codons.astype(np.int64)]
    return best_codon_for_aa[aa_seq]


def build_random_h0_synonymous_sequence(codons: np.ndarray, host_train: dict, rng: np.random.Generator) -> np.ndarray:
    """
    H0/H2 marginal synonymous recoding: draws each codon independently from the host's
    marginal synonymous frequency distribution P(c | aa(c)), destroying sequential order.
    """
    marg = host_train["marg"].astype(np.float64) + ALPHA
    out = np.empty(codons.size, dtype=np.uint8)
    aa_seq = CODON_AA[codons.astype(np.int64)]
    for i, a in enumerate(aa_seq):
        ch = CHANNEL_CODONS[int(a)]
        if ch.size <= 1:
            out[i] = codons[i]
        else:
            probs = marg[ch]
            probs = probs / probs.sum()
            out[i] = np.uint8(rng.choice(ch, p=probs))
    return out


def run_chimera_hgt_trojan_audit(args, out_dir: Path, cds_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "chimera_hgt_trojan_audit.log"
    entrez_cache = cds_dir / "chimera_entrez_cache"

    with open(log_path, "w", encoding="utf-8") as log_fh, TeeContext(log_fh):
        print("")
        print("=" * 108)
        print("SUPER CODON v14.8 — HORIZONTAL GENE TRANSFER (HGT), VIRAL TROJAN, ORGANELLE OS & BIOTECH CRASH AUDIT")
        print("=" * 108)

        # -------------------------------------------------------------------
        # STEP 1: Load Multi-Domain Reference OS Markov Kernels from CDS Cache
        # -------------------------------------------------------------------
        target_os_species = [
            "Human_T2T", "Mouse", "Zebrafish", "Drosophila",
            "Arabidopsis", "Yeast", "E.coli", "B.subtilis",
            "M.jannaschii", "H.salinarum", "Pseudomonas_aeruginosa",
            "Mycobacterium_tuberculosis", "Synechocystis_sp_PCC6803",
        ]
        prepare_cds_cache(
            cds_dir=cds_dir,
            force=False,
            timeout=int(args.download_timeout),
            include_mito=False,
            include_virus=False,
            include_extended=True,
            only=target_os_species,
        )

        os_models = {}
        os_meta = {}
        all_cached = deduplicate_species_files(collect_cds_files(cds_dir))
        for fpath in all_cached:
            sp = species_of(fpath)
            if sp not in target_os_species:
                continue
            units, raw, fam, used_table, _st = cds_units(
                fpath,
                max_records=5000,
                seed=int(args.seed),
                trim_start_stop=True,
                family_split=True,
            )
            if not units:
                continue
            tr_labels, te_labels = choose_split(
                list(units.keys()), units, min_trans=30, seed=int(args.seed), fam=fam, family_split=True
            )
            tr_agg = aggregate(tr_labels if tr_labels else list(units.keys()), units)
            te_agg = aggregate(te_labels if te_labels else list(units.keys()), units)
            self_pt = deltas(te_agg, tr_agg)
            gc3_sp = float(tr_agg["gc3_sum"]) / float(max(1, tr_agg["gc3_n"]))
            os_models[sp] = tr_agg
            os_meta[sp] = {
                "clade": clade_of(sp),
                "table": used_table,
                "gc3": gc3_sp,
                "native_d2b": self_pt["d2b"] if self_pt else 0.0,
                "native_d3b": self_pt["d3b"] if self_pt else 0.0,
                "file": fpath.name,
            }

        set_genetic_code(1, verbose=False)
        if "Human_T2T" not in os_models or "E.coli" not in os_models:
            raise RuntimeError("Human_T2T and E.coli reference models are required for v14.8 Chimera Audit.")

        human_train = os_models["Human_T2T"]
        ecoli_train = os_models["E.coli"]

        print("")
        print("-" * 108)
        print(f"LOADED {len(os_models)} REFERENCE GENOMIC 'OPERATING SYSTEM' KERNELS (64x64 MARKOV MATRICES):")
        print("-" * 108)
        for sp_k, m_info in os_meta.items():
            print(
                f"  [OS Kernel] {sp_k:26s} | clade={m_info['clade']:14s} | "
                f"GC3={m_info['gc3']:.4f} | native_d2b={fmt_float(m_info['native_d2b'], 4)} | "
                f"native_d3b={fmt_float(m_info['native_d3b'], 4)}"
            )

        # ===================================================================
        # MODULE A: VIRAL TROJAN & HOST-OS SPECIFICITY SCANNER
        # ===================================================================
        print("")
        print("#" * 108)
        print("MODULE A: VIRAL 'GENETIC TROJAN' & HOST-RIBOSOME MIMICRY SCANNER (VIRUSES vs HUMAN & E. COLI OS)")
        print("#" * 108)

        viral_summary = []
        for vname, vinfo in VIRAL_TROJAN_PANEL.items():
            try:
                fasta_txt = fetch_or_load_nuccore_cds_fasta(
                    name=vname,
                    acc=vinfo["acc"],
                    cache_dir=entrez_cache,
                    timeout=int(args.download_timeout),
                    force=bool(args.force_download),
                )
                set_genetic_code(int(vinfo["transl_table"]), verbose=False)
                v_genes = parse_fasta_cds_records(fasta_txt, trim_ends=True, min_codons=50)
                if not v_genes:
                    continue

                v_agg = new_counts()
                gene_trojan_rows = []
                for g in v_genes:
                    add_unit(v_agg, g["codons"], g["nt"])
                    uc_g = new_counts()
                    add_unit(uc_g, g["codons"], g["nt"])
                    pt_hum = deltas(uc_g, human_train)
                    pt_eco = deltas(uc_g, ecoli_train)
                    gm_g = compute_geometric_form_metrics(g["codons"], g["nt"])
                    if pt_hum and pt_eco and g["codons"].size >= 80:
                        gene_trojan_rows.append({
                            "gene": g["gene"],
                            "product": g["product"],
                            "codons": int(g["codons"].size),
                            "gc3": gm_g.get("gc3"),
                            "d2b_human": pt_hum["d2b"],
                            "d3b_human": pt_hum["d3b"],
                            "d2b_ecoli": pt_eco["d2b"],
                            "d3b_ecoli": pt_eco["d3b"],
                            "human_vs_ecoli_d2b_adv": float(pt_hum["d2b"] - pt_eco["d2b"]),
                        })

                pt_tot_hum = deltas(v_agg, human_train)
                pt_tot_eco = deltas(v_agg, ecoli_train)
                gc3_v = float(v_agg["gc3_sum"]) / float(max(1, v_agg["gc3_n"]))

                # Top Trojan mimic genes on Human ribosome (highest d3b_human + d2b_human)
                top_human_trojans = sorted(
                    gene_trojan_rows,
                    key=lambda r: (r["d3b_human"] + r["d2b_human"]),
                    reverse=True,
                )[:5]
                n_pos_d3b_hum = sum(1 for r in gene_trojan_rows if r["d3b_human"] > 0.0)
                n_pos_d2b_hum = sum(1 for r in gene_trojan_rows if r["d2b_human"] > 0.0)

                v_rec = {
                    "virus": vname,
                    "accession": vinfo["acc"],
                    "true_host": vinfo["true_host"],
                    "role": vinfo["role"],
                    "n_genes": len(v_genes),
                    "n_genes_ge_80cod": len(gene_trojan_rows),
                    "gc3": gc3_v,
                    "whole_d2b_on_human": pt_tot_hum["d2b"] if pt_tot_hum else None,
                    "whole_d3b_on_human": pt_tot_hum["d3b"] if pt_tot_hum else None,
                    "whole_d2b_on_ecoli": pt_tot_eco["d2b"] if pt_tot_eco else None,
                    "whole_d3b_on_ecoli": pt_tot_eco["d3b"] if pt_tot_eco else None,
                    "genes_pos_d3b_on_human": int(n_pos_d3b_hum),
                    "genes_pos_d2b_on_human": int(n_pos_d2b_hum),
                    "top_trojan_genes_on_human": top_human_trojans,
                }
                viral_summary.append(v_rec)

                print(
                    f"\n  [Virus] {vname} ({vinfo['acc']}) | True_Host={vinfo['true_host']} | "
                    f"Genes={len(v_genes)} | GC3={gc3_v:.4f}"
                )
                print(
                    f"    -> Scored on HUMAN_T2T OS : d2b={fmt_float(v_rec['whole_d2b_on_human'], 4)} | "
                    f"d3b={fmt_float(v_rec['whole_d3b_on_human'], 4)} | "
                    f"Genes with d3b>0 on Human: {n_pos_d3b_hum}/{len(gene_trojan_rows)}"
                )
                print(
                    f"    -> Scored on E. COLI OS   : d2b={fmt_float(v_rec['whole_d2b_on_ecoli'], 4)} | "
                    f"d3b={fmt_float(v_rec['whole_d3b_on_ecoli'], 4)}"
                )
                if top_human_trojans:
                    print("    -> Top 'Genetic Trojan' Genes Mimicking Human Ribosomal Syntax:")
                    print("       Gene_ID          | Codons |    GC3 | d2b_Human | d3b_Human | d2b_Ecoli | Product")
                    for tj in top_human_trojans:
                        print(
                            f"       {tj['gene'][:16]:16s} | {tj['codons']:6d} | {tj['gc3']:6.4f} | "
                            f"{fmt_float(tj['d2b_human'], 4):>9s} | {fmt_float(tj['d3b_human'], 4):>9s} | "
                            f"{fmt_float(tj['d2b_ecoli'], 4):>9s} | {tj['product']}"
                        )
            except Exception as e:
                print(f"  [Module A] WARNING: failed {vname}: {e}")

        set_genetic_code(1, verbose=False)

        # ===================================================================
        # MODULE B: HUMAN PROTEOME HGT & DOMESTICATED RETROVIRUS FORENSICS
        # ===================================================================
        print("")
        print("#" * 108)
        print("MODULE B: HUMAN PROTEOME HGT 'TORRENT NETWORK' & DOMESTICATED VIRAL/TRANSPOSON GENE FORENSICS")
        print("#" * 108)

        # B1. Evaluate famous domesticated retroviral / transposon / endosymbiotic human genes
        domesticated_results = []
        print("  [B1] Scoring Famous Domesticated Retroviral, Transposon & Endosymbiotic Genes Inside Human Genome:")
        print(
            "    Gene_Symbol             | Acc           | Codons |    GC3 | d2b_Human | d3b_Human | Best_NonMammal_OS (d3b)          | Origin"
        )
        print(
            "    ------------------------+---------------+--------+--------+-----------+-----------+----------------------------------+----------------------------------------"
        )
        for gsym, ginfo in HUMAN_CHIMERA_CANDIDATE_GENES.items():
            try:
                f_txt = fetch_or_load_nuccore_cds_fasta(
                    name=gsym,
                    acc=ginfo["acc"],
                    cache_dir=entrez_cache,
                    timeout=int(args.download_timeout),
                    force=bool(args.force_download),
                )
                recs = parse_fasta_cds_records(f_txt, trim_ends=True, min_codons=40)
                if not recs:
                    continue
                r0 = recs[0]
                uc = new_counts()
                add_unit(uc, r0["codons"], r0["nt"])
                gm = compute_geometric_form_metrics(r0["codons"], r0["nt"])
                pt_h = deltas(uc, human_train)

                non_mammal_scores = []
                for sp_k, tr_k in os_models.items():
                    if sp_k in ("Human_T2T", "Mouse"):
                        continue
                    pt_k = deltas(uc, tr_k)
                    if pt_k:
                        non_mammal_scores.append((sp_k, pt_k["d3b"], pt_k["d2b"]))
                best_nm = max(non_mammal_scores, key=lambda x: x[1]) if non_mammal_scores else ("n/a", 0.0, 0.0)

                d_row = {
                    "gene": gsym,
                    "acc": ginfo["acc"],
                    "codons": int(r0["codons"].size),
                    "gc3": gm.get("gc3"),
                    "d2b_human": pt_h["d2b"] if pt_h else None,
                    "d3b_human": pt_h["d3b"] if pt_h else None,
                    "best_non_mammal_os": best_nm[0],
                    "best_non_mammal_d3b": best_nm[1],
                    "best_non_mammal_d2b": best_nm[2],
                    "origin": ginfo["origin"],
                }
                domesticated_results.append(d_row)
                nm_str = f"{best_nm[0]} ({fmt_float(best_nm[1], 4)})"
                print(
                    f"    {gsym:23s} | {ginfo['acc']:13s} | {r0['codons'].size:6d} | "
                    f"{gm.get('gc3', 0.0):6.4f} | {fmt_float(d_row['d2b_human'], 4):>9s} | "
                    f"{fmt_float(d_row['d3b_human'], 4):>9s} | {nm_str:32s} | {ginfo['origin'][:42]}"
                )
            except Exception as e:
                print(f"    [B1] note: skipped {gsym}: {e}")

        # B2. Scan 1,500 individual Human_T2T genes across all loaded OS Kernels to map the 'Torrent Network'
        human_file = existing_cds_file(cds_dir, SPECIES_INFO["Human_T2T"]["acc"])
        hgt_scan_summary = {}
        foreign_engine_outliers = []
        if human_file is not None:
            with gzip.open(human_file, "rt", encoding="utf-8", errors="ignore") as gz_f:
                raw_human_fasta = gz_f.read()
            all_human_genes = parse_fasta_cds_records(raw_human_fasta, trim_ends=True, min_codons=150)
            rng_hgt = np.random.default_rng(int(args.seed) + 2026)
            if len(all_human_genes) > 1500:
                pick_idx = rng_hgt.choice(len(all_human_genes), size=1500, replace=False)
                sampled_human_genes = [all_human_genes[i] for i in pick_idx]
            else:
                sampled_human_genes = all_human_genes

            winner_counter_d2b = Counter()
            winner_counter_d3b = Counter()

            for hg in sampled_human_genes:
                uc_h = new_counts()
                add_unit(uc_h, hg["codons"], hg["nt"])
                scores_by_os = {}
                for sp_k, tr_k in os_models.items():
                    pt_k = deltas(uc_h, tr_k)
                    if pt_k:
                        scores_by_os[sp_k] = (pt_k["d2b"], pt_k["d3b"])
                if not scores_by_os or "Human_T2T" not in scores_by_os:
                    continue

                best_d2b_sp = max(scores_by_os.items(), key=lambda kv: kv[1][0])[0]
                best_d3b_sp = max(scores_by_os.items(), key=lambda kv: kv[1][1])[0]
                winner_counter_d2b[best_d2b_sp] += 1
                winner_counter_d3b[best_d3b_sp] += 1

                hum_d2b, hum_d3b = scores_by_os["Human_T2T"]
                # Check if a non-vertebrate (bacterial/archaeal/fungal/plant/insect) OS beats Human_T2T on BOTH d2b and d3b
                for sp_k, (k_d2b, k_d3b) in scores_by_os.items():
                    if sp_k in ("Human_T2T", "Mouse", "Zebrafish"):
                        continue
                    if k_d2b > hum_d2b + 0.035 and k_d3b > hum_d3b + 0.050 and k_d3b > 0.10:
                        gm_hg = compute_geometric_form_metrics(hg["codons"], hg["nt"])
                        foreign_engine_outliers.append({
                            "gene": hg["gene"],
                            "product": hg["product"],
                            "codons": int(hg["codons"].size),
                            "gc3": gm_hg.get("gc3"),
                            "foreign_os": sp_k,
                            "foreign_clade": os_meta[sp_k]["clade"],
                            "d2b_human": hum_d2b,
                            "d2b_foreign": k_d2b,
                            "d3b_human": hum_d3b,
                            "d3b_foreign": k_d3b,
                            "d3b_gain_over_human": float(k_d3b - hum_d3b),
                        })

            foreign_engine_outliers = sorted(
                foreign_engine_outliers, key=lambda x: x["d3b_gain_over_human"], reverse=True
            )
            hgt_scan_summary = {
                "n_human_genes_scanned": len(sampled_human_genes),
                "best_os_winner_counts_d2b": dict(winner_counter_d2b.most_common()),
                "best_os_winner_counts_d3b": dict(winner_counter_d3b.most_common()),
                "n_foreign_engine_outliers": len(foreign_engine_outliers),
                "top_foreign_engine_genes_in_human": foreign_engine_outliers[:12],
            }

            print("")
            print(
                f"  [B2] Whole-Genome 'Torrent Network' Scan Across {len(sampled_human_genes)} Human Genes (>=150 codons):"
            )
            print(f"    -> Best OS Kernel by Synonymous Syntax (d2b) : {dict(winner_counter_d2b.most_common(6))}")
            print(f"    -> Best OS Kernel by Trigram-Corrected (d3b) : {dict(winner_counter_d3b.most_common(6))}")
            print(
                f"    -> Human Genes Running Stronger on Non-Vertebrate (Prokaryotic/Archaeal/Fungal/Plant) OS: "
                f"{len(foreign_engine_outliers)} loci"
            )
            if foreign_engine_outliers:
                print("    -> Top 10 'Foreign-Engine' Chimeric Loci Detected Inside Human_T2T:")
                print(
                    "       Human_Gene       | Codons |    GC3 | Foreign_OS_Winner        | d3b_Human | d3b_Foreign | Gain_d3b | Product"
                )
                for fo in foreign_engine_outliers[:10]:
                    print(
                        f"       {fo['gene'][:16]:16s} | {fo['codons']:6d} | {fo['gc3']:6.4f} | "
                        f"{fo['foreign_os']:24s} | {fmt_float(fo['d3b_human'], 4):>9s} | "
                        f"{fmt_float(fo['d3b_foreign'], 4):>11s} | {fmt_float(fo['d3b_gain_over_human'], 4):>8s} | {fo['product']}"
                    )

        # ===================================================================
        # MODULE C: ENDOSYMBIOTIC ORGANELLE OS COLONIZATION AUDIT
        # ===================================================================
        print("")
        print("#" * 108)
        print("MODULE C: ENDOSYMBIOTIC ORGANELLE OS REWRITING (CHLOROPLAST & MITOCHONDRIA vs ANCESTOR & HOST NUCLEUS)")
        print("#" * 108)

        organelle_results = []
        for oname, oinfo in ORGANELLE_ENDOSYMBIOSIS_PANEL.items():
            try:
                f_txt = fetch_or_load_nuccore_cds_fasta(
                    name=oname,
                    acc=oinfo["acc"],
                    cache_dir=entrez_cache,
                    timeout=int(args.download_timeout),
                    force=bool(args.force_download),
                )
                set_genetic_code(int(oinfo["transl_table"]), verbose=False)
                o_genes = parse_fasta_cds_records(f_txt, trim_ends=True, min_codons=40)
                o_agg = new_counts()
                for g in o_genes:
                    add_unit(o_agg, g["codons"], g["nt"])
                gc3_o = float(o_agg["gc3_sum"]) / float(max(1, o_agg["gc3_n"]))

                os_scores = {}
                for sp_k in ("Arabidopsis", "Synechocystis_sp_PCC6803", "B.subtilis", "E.coli", "Human_T2T", "Yeast"):
                    if sp_k in os_models:
                        pt_k = deltas(o_agg, os_models[sp_k])
                        if pt_k:
                            os_scores[sp_k] = {"d2b": pt_k["d2b"], "d3b": pt_k["d3b"]}

                # Gene-by-gene split: Host Nucleus vs Bacterial Ancestor
                n_host_wins = 0
                n_bact_wins = 0
                host_key = "Arabidopsis" if "Chloroplast" in oname else "Human_T2T"
                bact_key = "Synechocystis_sp_PCC6803" if ("Chloroplast" in oname and "Synechocystis_sp_PCC6803" in os_models) else "B.subtilis"
                for g in o_genes:
                    uc_g = new_counts()
                    add_unit(uc_g, g["codons"], g["nt"])
                    if host_key in os_models and bact_key in os_models:
                        pt_h = deltas(uc_g, os_models[host_key])
                        pt_b = deltas(uc_g, os_models[bact_key])
                        if pt_h and pt_b:
                            if pt_h["d3b"] >= pt_b["d3b"]:
                                n_host_wins += 1
                            else:
                                n_bact_wins += 1

                o_row = {
                    "organelle": oname,
                    "acc": oinfo["acc"],
                    "transl_table": oinfo["transl_table"],
                    "n_genes": len(o_genes),
                    "gc3": gc3_o,
                    "scores_across_os": os_scores,
                    "host_vs_bact_gene_split": {
                        "host_os": host_key,
                        "bacterial_os": bact_key,
                        "genes_closer_to_host_nucleus": n_host_wins,
                        "genes_closer_to_bacteria": n_bact_wins,
                    },
                }
                organelle_results.append(o_row)
                print(
                    f"\n  [Organelle] {oname} ({oinfo['acc']}, Table={oinfo['transl_table']}) | "
                    f"Genes={len(o_genes)} | GC3={gc3_o:.4f}"
                )
                for sp_k, sc_k in os_scores.items():
                    print(
                        f"    -> Scored under {sp_k:26s} : d2b={fmt_float(sc_k['d2b'], 4)} | d3b={fmt_float(sc_k['d3b'], 4)}"
                    )
                print(
                    f"    -> Gene-by-Gene OS Split ({host_key} Nucleus vs {bact_key} Bacteria): "
                    f"{n_host_wins} genes assimilated by Host vs {n_bact_wins} retaining Bacterial syntax!"
                )
            except Exception as e:
                print(f"  [Module C] WARNING: failed {oname}: {e}")

        set_genetic_code(1, verbose=False)

        # ===================================================================
        # MODULE D: BIOTECH "MAX-CAI CODON OPTIMIZATION" CRASH TEST
        # ===================================================================
        print("")
        print("#" * 108)
        print("MODULE D: THE BIOTECH SHOCK TEST — NATURAL GENES vs INDUSTRY 'MAX-CAI CODON OPTIMIZATION' CRASH")
        print("#" * 108)
        print(
            "  Testing 100% identical human proteins across 4 coding designs:\n"
            "    (1) WT_Natural     : Real human spliced exons from NCBI RefSeqGene\n"
            "    (2) Biotech_MaxCAI : Standard industry \'Codon Optimization\' (every AA -> most frequent human codon)\n"
            "    (3) Random_Syn_H0  : Independent marginal human codon usage (no Markov memory)\n"
            "    (4) Markov_H2b     : Generated by our 1st-order conditional codon-pair Markov operator P(c_{i+1}|c_i)"
        )
        print("-" * 108)
        print(
            "  Gene   | Codons | WT_d2b  | MaxCAI_d2b | WT_d3b  | MaxCAI_d3b | H0_Syn_d3b | Markov_d3b | WT_Ent(b) | MaxCAI_Ent | Verdict"
        )
        print(
            "  -------+--------+---------+------------+---------+------------+------------+------------+-----------+------------+-------------------"
        )

        gb_cache_dir = cds_dir / "refseqgene_gb"
        h2b_sampler = build_h2b_sampler(human_train)
        rng_opt = np.random.default_rng(int(args.seed) + 909)
        biotech_crash_rows = []

        test_genes_d = ["TP53", "TP63", "TP73", "BRCA1", "PTEN", "RB1", "KRAS", "NRAS", "HRAS", "GAPDH", "ACTB", "HBB", "MYC", "APOE"]
        for gsym in test_genes_d:
            info = INTRON_EXON_REFSEQGENE_PANEL.get(gsym)
            if not info:
                continue
            try:
                gb_txt = fetch_or_load_genbank_gene(
                    gene_symbol=gsym,
                    acc=info["acc"],
                    gb_cache_dir=gb_cache_dir,
                    timeout=int(args.download_timeout),
                    force=False,
                )
                parsed = parse_genbank_intron_exon(gb_txt, target_gene=gsym)
                cd_wt, nt_wt = bytes_to_codons(parsed["spliced_cds_str"].encode("ascii", "ignore"))
                if cd_wt.size >= 5:
                    cd_wt = cd_wt[1:-1]
                    nt_wt = nt_wt[3:-3]

                # Build the 3 synthetic synonymous designs (100% identical protein sequence!)
                cd_cai = build_max_cai_sequence(cd_wt, human_train)
                cd_h0 = build_random_h0_synonymous_sequence(cd_wt, human_train, rng_opt)
                cd_mk = generate_h5_surrogate(cd_wt, h2b_sampler, rng_opt)

                # Verify 100% amino-acid identity across all 4 sequences
                assert np.array_equal(CODON_AA[cd_wt], CODON_AA[cd_cai])
                assert np.array_equal(CODON_AA[cd_wt], CODON_AA[cd_h0])
                assert np.array_equal(CODON_AA[cd_wt], CODON_AA[cd_mk])

                def _eval_cod_arr(c_arr):
                    s_str = "".join(CODONS[int(x)] for x in c_arr)
                    _, n_arr = bytes_to_codons(s_str.encode("ascii"))
                    uc_x = new_counts()
                    add_unit(uc_x, c_arr, n_arr)
                    pt_x = deltas(uc_x, human_train)
                    gm_x = compute_geometric_form_metrics(c_arr, n_arr)
                    return pt_x, gm_x

                pt_wt, gm_wt = _eval_cod_arr(cd_wt)
                pt_cai, gm_cai = _eval_cod_arr(cd_cai)
                pt_h0, _gm_h0 = _eval_cod_arr(cd_h0)
                pt_mk, _gm_mk = _eval_cod_arr(cd_mk)

                verdict_str = "MAX-CAI CRASHES d2b/Ent" if (pt_cai["d2b"] < pt_wt["d2b"] and pt_cai["d2b"] < 0.0 and gm_cai["codon_entropy_bits"] < gm_wt["codon_entropy_bits"] - 0.85) else "ALTERED"
                row_d = {
                    "gene": gsym,
                    "codons": int(cd_wt.size),
                    "wt_d2b": pt_wt["d2b"],
                    "max_cai_d2b": pt_cai["d2b"],
                    "wt_d3b": pt_wt["d3b"],
                    "max_cai_d3b": pt_cai["d3b"],
                    "h0_syn_d3b": pt_h0["d3b"],
                    "markov_h2b_d3b": pt_mk["d3b"],
                    "wt_codon_entropy_bits": gm_wt["codon_entropy_bits"],
                    "max_cai_codon_entropy_bits": gm_cai["codon_entropy_bits"],
                    "verdict": verdict_str,
                }
                biotech_crash_rows.append(row_d)
                print(
                    f"  {gsym:6s} | {cd_wt.size:6d} | {fmt_float(pt_wt['d2b'], 4):>7s} | "
                    f"{fmt_float(pt_cai['d2b'], 4):>10s} | {fmt_float(pt_wt['d3b'], 4):>7s} | "
                    f"{fmt_float(pt_cai['d3b'], 4):>10s} | {fmt_float(pt_h0['d3b'], 4):>10s} | "
                    f"{fmt_float(pt_mk['d3b'], 4):>10s} | {gm_wt['codon_entropy_bits']:9.4f} | "
                    f"{gm_cai['codon_entropy_bits']:10.4f} | {verdict_str}"
                )
            except Exception as e:
                print(f"  [Module D] warning on {gsym}: {e}")

        if biotech_crash_rows:
            mean_wt_d2b = float(np.mean([r["wt_d2b"] for r in biotech_crash_rows]))
            mean_cai_d2b = float(np.mean([r["max_cai_d2b"] for r in biotech_crash_rows]))
            mean_wt_d3b = float(np.mean([r["wt_d3b"] for r in biotech_crash_rows]))
            mean_cai_d3b = float(np.mean([r["max_cai_d3b"] for r in biotech_crash_rows]))
            mean_h0_d3b = float(np.mean([r["h0_syn_d3b"] for r in biotech_crash_rows]))
            mean_mk_d3b = float(np.mean([r["markov_h2b_d3b"] for r in biotech_crash_rows]))
            mean_wt_ent = float(np.mean([r["wt_codon_entropy_bits"] for r in biotech_crash_rows]))
            mean_cai_ent = float(np.mean([r["max_cai_codon_entropy_bits"] for r in biotech_crash_rows]))
            print("-" * 108)
            print(
                f"  PANEL MEAN ({len(biotech_crash_rows)} HUMAN GENES): "
                f"WT_d2b={fmt_float(mean_wt_d2b, 4)} vs MaxCAI_d2b={fmt_float(mean_cai_d2b, 4)} | "
                f"WT_d3b={fmt_float(mean_wt_d3b, 4)} vs MaxCAI_d3b={fmt_float(mean_cai_d3b, 4)} vs "
                f"H0_Random_d3b={fmt_float(mean_h0_d3b, 4)} vs Markov_H2b_d3b={fmt_float(mean_mk_d3b, 4)} | "
                f"Entropy: {mean_wt_ent:.3f}b -> {mean_cai_ent:.3f}b"
            )

        payload = {
            "version": "14.8",
            "os_kernels": os_meta,
            "module_a_viral_trojans": viral_summary,
            "module_b_domesticated_genes": domesticated_results,
            "module_b_proteome_hgt_scan": hgt_scan_summary,
            "module_c_organelle_endosymbiosis": organelle_results,
            "module_d_biotech_codon_optimization_crash": biotech_crash_rows,
        }

        json_out = out_dir / "chimera_hgt_trojan_audit.json"
        md_out = out_dir / "chimera_hgt_trojan_audit.md"
        atomic_write_text(json_out, json.dumps(_clean_json(payload), indent=2, sort_keys=False))

        md_lines = [
            "# Super Codon v14.8 — Viral Trojan, HGT Torrent Network, Organelle OS & Biotech Codon-Optimization Audit",
            "",
            "## Module A: Viral 'Genetic Trojan' & Host-OS Specificity",
            "",
            "| Virus | Accession | True Host | Genes | GC3 | d2b (Human OS) | d3b (Human OS) | d2b (E. coli OS) | d3b (E. coli OS) |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|",
        ]
        for vr in viral_summary:
            md_lines.append(
                f"| **{vr['virus']}** | `{vr['accession']}` | `{vr['true_host']}` | {vr['n_genes']} | "
                f"{vr['gc3']:.4f} | {fmt_float(vr['whole_d2b_on_human'], 4)} | {fmt_float(vr['whole_d3b_on_human'], 4)} | "
                f"{fmt_float(vr['whole_d2b_on_ecoli'], 4)} | {fmt_float(vr['whole_d3b_on_ecoli'], 4)} |"
            )

        md_lines.extend([
            "",
            "## Module D: Biotech 'Max-CAI Codon Optimization' Crash Test (100% Identical Proteins)",
            "",
            "| Gene | Codons | WT d2b | MaxCAI d2b | WT d3b | MaxCAI d3b | H0 Random d3b | Markov H2b d3b | WT Entropy (b) | MaxCAI Entropy (b) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ])
        for br in biotech_crash_rows:
            md_lines.append(
                f"| **{br['gene']}** | {br['codons']} | {fmt_float(br['wt_d2b'], 4)} | {fmt_float(br['max_cai_d2b'], 4)} | "
                f"{fmt_float(br['wt_d3b'], 4)} | {fmt_float(br['max_cai_d3b'], 4)} | {fmt_float(br['h0_syn_d3b'], 4)} | "
                f"{fmt_float(br['markov_h2b_d3b'], 4)} | {br['wt_codon_entropy_bits']:.4f} | {br['max_cai_codon_entropy_bits']:.4f} |"
            )

        atomic_write_text(md_out, "\n".join(md_lines) + "\n")
        print("")
        print("=" * 108)
        print(f"[v14.8 Chimera Audit Saved] JSON : {json_out}")
        print(f"[v14.8 Chimera Audit Saved] MD   : {md_out}")
        print(f"[v14.8 Chimera Audit Saved] LOG  : {log_path}")
        print("=" * 108)
        return payload


def build_arg_parser():
    p = argparse.ArgumentParser(
        description="SUPER CODON CONDITIONAL SYNTAX PROOF PIPELINE — AUTONOMOUS v14.6"
    )
    p.add_argument("--cds-dir", type=str, default="cds_cache", help="Directory for cached NCBI CDS & GenBank files")
    p.add_argument("--out-dir", type=str, default="results_v14_6", help="Output directory for JSON/MD reports and logs")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Master RNG seed")
    p.add_argument("--force-download", action="store_true", help="Re-download CDS/GenBank files even if cached")
    p.add_argument("--download-timeout", type=int, default=180, help="HTTP timeout in seconds")

    # Mode selection
    p.add_argument("--self-test", action="store_true", help="Run built-in invariant & splicing self-tests and exit")
    p.add_argument("--p53-audit", action="store_true", help="Run autonomous TP53 & Intron-vs-Exon geometric/program audit")
    p.add_argument("--ras-audit", action="store_true", help="Run autonomous RAS oncogene family (KRAS, NRAS, HRAS) Intron-vs-Exon audit")
    p.add_argument("--frontier-audit", action="store_true", help="Run Frontier Audit: Giant Viruses, Deep Subsurface Boreholes, Venter Syn3.0 & Asgard Archaea")
    p.add_argument("--chimera-audit", action="store_true", help="Run v14.8 Viral Trojan, HGT Torrent Network, Organelle OS & Biotech Max-CAI Crash Audit")
    p.add_argument("--frontier-suite", type=str, default="all", help="Filter frontier suite: all, 1_giant_viruses, 2_deep_biosphere, 3_synthetic_cell, 4_evolution_boundary")
    p.add_argument("--intron-exon-audit", action="store_true", help="Run full 12-gene Intron-vs-Exon & P53 panel audit")
    p.add_argument("--intron-exon-genes", type=str, default=None, help="Comma-separated gene symbols from panel (e.g. TP53,BRCA1,HBB)")
    p.add_argument("--custom-gb", type=str, default=None, help="Path to custom local GenBank flatfile for Intron/Exon audit")
    p.add_argument("--all-trajectories", action="store_true", help="Compute progressive splicing trajectory for all panel genes")
    p.add_argument("--p53-mc-reps", type=int, default=40, help="Monte Carlo reps per sequence in Intron/Exon audit")
    p.add_argument("--p53-ref-max-records", type=int, default=6000, help="Max Human_T2T CDS genes for reference model")

    p.add_argument("--smoke-only", action="store_true", help="Run only fast smoke stage")
    p.add_argument("--full-only", action="store_true", help="Run only full stage")
    p.add_argument("--with-smoke", action="store_true", help="Also run smoke stage alongside --p53-audit")
    p.add_argument("--with-full", action="store_true", help="Also run full stage alongside --p53-audit")

    # Opt-in tracks
    p.add_argument("--include-mito", action="store_true", help="Include Human_mito track")
    p.add_argument("--include-virus", action="store_true", help="Include SARS-CoV-2 track")
    p.add_argument("--include-extended", action="store_true", help="Include 41-taxon extended biosphere panel")
    p.add_argument("--mito-acc", type=str, default=None, help="Override accession for Human_mito")
    p.add_argument("--mito-asm", type=str, default=None, help="Override assembly name for Human_mito")
    p.add_argument("--virus-acc", type=str, default=None, help="Override accession for SARS-CoV-2")
    p.add_argument("--virus-asm", type=str, default=None, help="Override assembly name for SARS-CoV-2")
    p.add_argument("--only-species", type=str, default=None, help="Comma-separated species filter")

    # Full-stage tuning parameters
    p.add_argument("--max-records", type=int, default=4000, help="Max CDS records per species in full mode (0=all)")
    p.add_argument("--min-trans", type=int, default=30, help="Min codon transitions per gene for split eligibility")
    p.add_argument("--bootstrap-reps", type=int, default=200, help="Bootstrap reps for H0..H3b CIs in full mode")
    p.add_argument("--h4-reps", type=int, default=20, help="Block-shuffle reps for H4")
    p.add_argument("--h4-bootstrap-reps", type=int, default=40, help="Bootstrap reps for H4 CIs")
    p.add_argument("--h5-reps", type=int, default=20, help="Surrogate reps for H5")
    p.add_argument("--h5-bootstrap-reps", type=int, default=40, help="Bootstrap reps for H5 CIs")
    p.add_argument("--h5-bootstrap-inner-reps", type=int, default=2, help="Inner surrogate reps per H5 bootstrap draw")
    p.add_argument("--h5-transfer-reps", type=int, default=15, help="Surrogate reps for cross-species H5 transfer")
    p.add_argument("--h6-reps", type=int, default=20, help="Synonymous shuffle reps for H6")
    p.add_argument("--h6-bootstrap-reps", type=int, default=40, help="Bootstrap reps for H6 CIs")
    p.add_argument("--h7-reps", type=int, default=20, help="Mirror-pair shuffle reps for H7")
    p.add_argument("--h7-bootstrap-reps", type=int, default=40, help="Bootstrap reps for H7 CIs")
    p.add_argument("--no-trim-start-stop", action="store_true", help="Do not trim start/stop codons")
    p.add_argument("--no-family-split", action="store_true", help="Disable paralog-family grouped train/test split")
    p.add_argument("--min-within-units", type=int, default=50, help="Small-sample guard: min units within species")
    p.add_argument("--min-within-trans", type=int, default=5000, help="Small-sample guard: min transitions within species")
    p.add_argument("--min-transfer-units", type=int, default=50, help="Small-sample guard: min units for transfer")
    p.add_argument("--min-transfer-trans", type=int, default=5000, help="Small-sample guard: min transitions for transfer")
    return p


def main(argv=None):
    args = build_arg_parser().parse_args(argv)

    if args.self_test:
        print("=" * 92)
        print("RUNNING BUILT-IN INVARIANT & SPLICING SELF-TEST (v14.6)")
        print("=" * 92)
        self_test_invariants(verbose=True)
        return 0

    cds_dir = Path(args.cds_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    master_log_path = out_dir / "master.log"

    overrides = {}
    if args.mito_acc or args.mito_asm:
        overrides["Human_mito"] = (args.mito_acc, args.mito_asm)
        if args.mito_acc:
            SPECIES_PREFIX[args.mito_acc] = "Human_mito"
    if args.virus_acc or args.virus_asm:
        overrides["SARS-CoV-2"] = (args.virus_acc, args.virus_asm)
        if args.virus_acc:
            SPECIES_PREFIX[args.virus_acc] = "SARS-CoV-2"

    only_list = (
        [x.strip() for x in args.only_species.split(",") if x.strip()]
        if args.only_species
        else None
    )

    with open(master_log_path, "a", encoding="utf-8") as master_fh, TeeContext(master_fh):
        want_frontier = bool(getattr(args, "frontier_audit", False))
        want_chimera = bool(getattr(args, "chimera_audit", False))
        want_p53 = bool(
            args.p53_audit
            or getattr(args, "ras_audit", False)
            or args.intron_exon_audit
            or args.intron_exon_genes
            or args.custom_gb
        )
        run_core_stages = (
            (not want_p53 and not want_frontier and not want_chimera)
            or args.with_smoke
            or args.with_full
            or args.smoke_only
            or args.full_only
        )

        if want_p53:
            run_intron_exon_p53_audit(args, out_dir=out_dir, cds_dir=cds_dir)

        if want_frontier:
            run_frontier_audit(args, out_dir=out_dir, cds_dir=cds_dir)

        if want_chimera:
            run_chimera_hgt_trojan_audit(args, out_dir=out_dir, cds_dir=cds_dir)

        if run_core_stages:
            prepare_cds_cache(
                cds_dir=cds_dir,
                force=bool(args.force_download),
                timeout=int(args.download_timeout),
                include_mito=bool(args.include_mito),
                include_virus=bool(args.include_virus),
                include_extended=bool(args.include_extended),
                overrides=overrides,
                only=only_list,
            )
            files = deduplicate_species_files(collect_cds_files(cds_dir), overrides=overrides)
            if only_list:
                files = [
                    p for p in files
                    if any(o.lower() in species_of(p).lower() for o in only_list)
                ]

            if not files:
                print("[ERROR] No valid CDS files available to run core pipeline stages.")
                return 1

            do_smoke = True
            do_full = True
            if args.smoke_only:
                do_full = False
            elif args.full_only:
                do_smoke = False
            elif want_p53:
                do_smoke = bool(args.with_smoke)
                do_full = bool(args.with_full)

            if do_smoke:
                smoke_cfg = {
                    "max_records": 400,
                    "min_trans": 20,
                    "trim_start_stop": not bool(args.no_trim_start_stop),
                    "family_split": not bool(args.no_family_split),
                    "bootstrap_reps": 30,
                    "h4_lengths": DEFAULT_H4_LENGTHS,
                    "h4_reps": 4,
                    "h4_bootstrap_reps": 8,
                    "h5_reps": 4,
                    "h5_bootstrap_reps": 8,
                    "h5_bootstrap_inner_reps": 1,
                    "h5_transfer_reps": 4,
                    "h6_reps": 4,
                    "h6_bootstrap_reps": 8,
                    "h7_reps": 4,
                    "h7_bootstrap_reps": 8,
                    "include_mito": bool(args.include_mito),
                    "include_virus": bool(args.include_virus),
                    "include_extended": bool(args.include_extended),
                    "min_within_units": int(args.min_within_units),
                    "min_within_trans": int(args.min_within_trans),
                    "min_transfer_units": int(args.min_transfer_units),
                    "min_transfer_trans": int(args.min_transfer_trans),
                    "transfer_pairs": DEFAULT_H5_TRANSFER_PAIRS,
                }
                run_pipeline_stage("smoke", smoke_cfg, files, args, out_dir=out_dir)

            if do_full:
                full_cfg = {
                    "max_records": int(args.max_records),
                    "min_trans": int(args.min_trans),
                    "trim_start_stop": not bool(args.no_trim_start_stop),
                    "family_split": not bool(args.no_family_split),
                    "bootstrap_reps": int(args.bootstrap_reps),
                    "h4_lengths": DEFAULT_H4_LENGTHS,
                    "h4_reps": int(args.h4_reps),
                    "h4_bootstrap_reps": int(args.h4_bootstrap_reps),
                    "h5_reps": int(args.h5_reps),
                    "h5_bootstrap_reps": int(args.h5_bootstrap_reps),
                    "h5_bootstrap_inner_reps": int(args.h5_bootstrap_inner_reps),
                    "h5_transfer_reps": int(args.h5_transfer_reps),
                    "h6_reps": int(args.h6_reps),
                    "h6_bootstrap_reps": int(args.h6_bootstrap_reps),
                    "h7_reps": int(args.h7_reps),
                    "h7_bootstrap_reps": int(args.h7_bootstrap_reps),
                    "include_mito": bool(args.include_mito),
                    "include_virus": bool(args.include_virus),
                    "include_extended": bool(args.include_extended),
                    "min_within_units": int(args.min_within_units),
                    "min_within_trans": int(args.min_within_trans),
                    "min_transfer_units": int(args.min_transfer_units),
                    "min_transfer_trans": int(args.min_transfer_trans),
                    "transfer_pairs": DEFAULT_H5_TRANSFER_PAIRS,
                }
                run_pipeline_stage("full", full_cfg, files, args, out_dir=out_dir)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.stderr.write("\n[INTERRUPTED] Stopped by user.\n")
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(1)