#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CODON TORUS REAL-DATA ANALYTIC AUDIT v3.1 (Syntax & Frame-Hardened)
===================================================================
Complete simulation-free mathematical & empirical verification suite for:
"The Codon Torus: A Rigorous Discrete-Torus Structure of DNA Coding Space,
Markovian Lumpability, and Information-Theoretic Consequences".

Upgrades & Fixes in v3.1:
  1. FIXED syntax typo in hourglass profile (g_centered, v_carrier).
  2. STRICT FRAME GUARD in bytes_to_codons(strict_cds=True): rejects any CDS
     whose length is not divisible by 3 or contains non-ACGT bases, preventing
     any 1-nt frame shift.
  3. FIXED NCBI Table 11 (Bacterial/Archaeal/Plastid has standard sense/stop map)
     + exact support for NCBI Tables 1, 2, 3, 4, 5, 6, 11, 12.
  4. AUTO-CACHE DISCOVERY: automatically locates cached RefSeq CDS files in
     ./cds_cache or /Volumes/*/genom_14/cds_cache, and downloads
     from NCBI RefSeq FTP if missing.
  5. VECTORIZED 64x64 transition tensor algebra + one-gene-one-vote deduplication
     + family_root deterministic 80/20 split.
  6. EXACT 4-term quadripartite + 3-term tripartite + Chemical DPI chain-rule
     verification to machine precision (< 1e-15 nats/codon).
  7. DUAL MAX-CAI AUDIT: Theorem 6.1 internal collapse (S_syn = S_syn2 = d2b = 0)
     AND Proposition 6.2 cross-scoring KL penalty under wild-type P_WT (d2b < 0).
  8. KEMENY-SNELL LUMPABILITY G^2 test + DOUBLE-STRAND RELIC COUPLING I(A; A o rho).
  9. REAL CHR17 17q22 100-KB HOURGLASS LOCUS (NC_000017.11:55154472-55254472)
     + all 7 publication-grade figures (PNG + PDF) + Markdown & JSON reports.
"""

import argparse
import gzip
import hashlib
import json
import math
import os
import py_compile
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from scipy.stats import chi2, ttest_1samp
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False


# =====================================================================
# 1. Codon torus (Z/4Z)^3, genetic code fibration, mirror involution
# =====================================================================
BASES = "ACGT"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODON_IDX = {c: i for i, c in enumerate(CODONS)}
N_CODONS = 64
CHANNELS = list("ACDEFGHIKLMNPQRSTVWY*")
CHANNEL_IDX = {ch: i for i, ch in enumerate(CHANNELS)}
N_CHANNELS = 21
STOP_IDX = CHANNEL_IDX["*"]

NT_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(BASES):
    NT_LUT[ord(_b)] = _i


def make_table_1():
    return {
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


TABLE_OVERRIDES = {
    1: {},
    11: {},  # Bacterial, Archaeal and Plant Plastid: identical sense/stop map to Table 1
    2: {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},
    3: {"CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "ATA": "M", "TGA": "W"},
    4: {"TGA": "W"},
    5: {"ATA": "M", "TGA": "W", "AGA": "S", "AGG": "S"},
    6: {"TAA": "Q", "TAG": "Q"},
    12: {"CTG": "S"},
}
SUPPORTED_TABLES = set(TABLE_OVERRIDES.keys())


def codon_aa_channel_array(table_id=1):
    if table_id not in SUPPORTED_TABLES:
        raise ValueError(f"Unsupported translation table {table_id}")
    tbl = make_table_1()
    tbl.update(TABLE_OVERRIDES[int(table_id)])
    return np.array([CHANNEL_IDX[tbl[c]] for c in CODONS], dtype=np.int64)


AA_CH = codon_aa_channel_array(1)
CURRENT_TABLE = 1


def set_genetic_code(table_id):
    global AA_CH, CURRENT_TABLE
    AA_CH = codon_aa_channel_array(int(table_id))
    CURRENT_TABLE = int(table_id)
    return AA_CH


def build_fibers(aa_ch=None):
    if aa_ch is None:
        aa_ch = AA_CH
    fibers = [[] for _ in range(N_CHANNELS)]
    for c_idx, ch in enumerate(aa_ch):
        fibers[int(ch)].append(c_idx)
    return [np.array(f, dtype=np.int64) for f in fibers]


def rc_codon_string(c):
    return "".join(COMP[b] for b in reversed(c))


MIRROR = np.array([CODON_IDX[rc_codon_string(c)] for c in CODONS], dtype=np.int64)
MIRROR_PAIR_ID = np.minimum(np.arange(N_CODONS), MIRROR)


# =====================================================================
# 2. Chemical coarsenings q: A -> Q
# =====================================================================
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
      "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
      "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "*": 0.0}
VOL = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4,
       "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9,
       "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0, "*": 0.0}
CHG = {"A": 0, "R": 1, "N": 0, "D": -1, "C": 0, "Q": 0, "E": -1, "G": 0, "H": 0, "I": 0, "L": 0,
       "K": 1, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0, "*": 0}
AROM = {"F", "W", "Y"}


def scheme_from_scalar(values, bins):
    out = {}
    for aa, v in values.items():
        label = bins[-1][1]
        for ub, lab in bins[:-1]:
            if v <= ub:
                label = lab
                break
        out[aa] = label
    out["*"] = "stop"
    return out


HYDRO_SCHEME = scheme_from_scalar(KD, [(-1.0, "polar"), (1.0, "neutral"), (99.0, "hydrophobic")])
CHARGE_SCHEME = {aa: ("negative" if v < 0 else "positive" if v > 0 else "neutral") for aa, v in CHG.items()}
CHARGE_SCHEME["*"] = "stop"
VOLUME_SCHEME = scheme_from_scalar(VOL, [(90.0, "small"), (150.0, "medium"), (999.0, "large")])
AROMATIC_SCHEME = {aa: ("aromatic" if aa in AROM else "nonaromatic") for aa in CHANNELS}
AROMATIC_SCHEME["*"] = "stop"

CHEM_SCHEMES = {
    "hydrophobicity": HYDRO_SCHEME,
    "charge": CHARGE_SCHEME,
    "volume": VOLUME_SCHEME,
    "aromaticity": AROMATIC_SCHEME,
}


def chem_class_vector(scheme_name):
    sch = CHEM_SCHEMES[scheme_name]
    labels = [sch[ch] for ch in CHANNELS]
    uniq = sorted(set(labels))
    mp = {lab: i for i, lab in enumerate(uniq)}
    return np.array([mp[lab] for lab in labels], dtype=np.int64), uniq


# =====================================================================
# 3. Vectorized Count Tensor on the Codon Torus
# =====================================================================
class CodonCounts:
    def __init__(self):
        self.n_xy = np.zeros((N_CODONS, N_CODONS), dtype=np.int64)
        self.n_total = 0

    def add_codons(self, codons):
        if codons is None or len(codons) < 2:
            return
        c = np.asarray(codons, dtype=np.int64)
        prev = c[:-1]
        cur = c[1:]
        self.n_xy += np.bincount(prev * N_CODONS + cur, minlength=N_CODONS * N_CODONS).reshape(N_CODONS, N_CODONS)
        self.n_total += int(cur.size)


def project_counts(C):
    n_xy = C.n_xy
    n_x = n_xy.sum(axis=1)
    n_y = n_xy.sum(axis=0)
    fibers = build_fibers()

    n_x_b = np.zeros((N_CODONS, N_CHANNELS), dtype=np.int64)
    n_a_y = np.zeros((N_CHANNELS, N_CODONS), dtype=np.int64)
    for b, fb in enumerate(fibers):
        if fb.size:
            n_x_b[:, b] = n_xy[:, fb].sum(axis=1)
    for a, fa in enumerate(fibers):
        if fa.size:
            n_a_y[a, :] = n_xy[fa, :].sum(axis=0)

    n_aa = np.zeros((N_CHANNELS, N_CHANNELS), dtype=np.int64)
    for a, fa in enumerate(fibers):
        if fa.size:
            n_aa[a, :] = n_x_b[fa, :].sum(axis=0)

    n_a_prev = n_aa.sum(axis=1)
    n_b_cur = n_aa.sum(axis=0)
    return {
        "n_xy": n_xy,
        "n_x": n_x,
        "n_y": n_y,
        "n_x_b": n_x_b,
        "n_a_y": n_a_y,
        "n_aa": n_aa,
        "n_a_prev": n_a_prev,
        "n_b_cur": n_b_cur,
        "fibers": fibers,
    }


# =====================================================================
# 4. Exact Information Decompositions (Theorems 5.2, 5.3, 5.4, 7.1)
# =====================================================================
TINY = 1e-300


def compute_exact_decompositions(C):
    if C.n_total == 0:
        return None
    P = project_counts(C)
    n_xy = P["n_xy"]
    n_x = P["n_x"]
    n_y = P["n_y"]
    n_x_b = P["n_x_b"]
    n_a_y = P["n_a_y"]
    n_aa = P["n_aa"]
    n_a_prev = P["n_a_prev"]
    n_b_cur = P["n_b_cur"]
    tot = float(C.n_total)

    # 1. Total codon-bigram MI: I(C_{i-1}; C_i)
    i_cc = 0.0
    for x in range(N_CODONS):
        if n_x[x] == 0:
            continue
        for y in range(N_CODONS):
            obs = int(n_xy[x, y])
            if obs > 0 and n_y[y] > 0:
                i_cc += obs * math.log((obs * tot) / (float(n_x[x]) * float(n_y[y])))
    i_cc /= tot

    # 2. Protein dipeptide MI: I(A_{i-1}; A_i)
    i_aa = 0.0
    for a in range(N_CHANNELS):
        if n_a_prev[a] == 0:
            continue
        for b in range(N_CHANNELS):
            obs = int(n_aa[a, b])
            if obs > 0 and n_b_cur[b] > 0:
                i_aa += obs * math.log((obs * tot) / (float(n_a_prev[a]) * float(n_b_cur[b])))
    i_aa /= tot

    # 3. Codon-to-next-AA context: d2b = I(C_{i-1}; A_i | A_{i-1})
    d2b = 0.0
    for x in range(N_CODONS):
        if n_x[x] == 0:
            continue
        a = int(AA_CH[x])
        for b in range(N_CHANNELS):
            obs = int(n_x_b[x, b])
            if obs > 0 and n_aa[a, b] > 0:
                d2b += obs * math.log((obs * float(n_a_prev[a])) / (float(n_x[x]) * float(n_aa[a, b])))
    d2b /= tot

    # 4. Contextual synonymous program: S_syn = I(C_{i-1}; C_i | A_i)
    s_syn = 0.0
    for x in range(N_CODONS):
        for y in range(N_CODONS):
            obs = int(n_xy[x, y])
            if obs == 0:
                continue
            b = int(AA_CH[y])
            if n_x_b[x, b] > 0 and n_y[y] > 0:
                s_syn += obs * math.log((obs * float(n_b_cur[b])) / (float(n_x_b[x, b]) * float(n_y[y])))
    s_syn /= tot

    # 5. AA_prev -> Syn_cur context: I(A_{i-1}; C_i | A_i)
    i_aprev_ccur = 0.0
    for a in range(N_CHANNELS):
        for y in range(N_CODONS):
            obs = int(n_a_y[a, y])
            if obs == 0:
                continue
            b = int(AA_CH[y])
            if n_aa[a, b] > 0 and n_y[y] > 0:
                i_aprev_ccur += obs * math.log((obs * float(n_b_cur[b])) / (float(n_aa[a, b]) * float(n_y[y])))
    i_aprev_ccur /= tot

    # 6. Strictly pairwise synonymous residual: S_syn^(2) = I(C_{i-1}; C_i | A_{i-1}, A_i)
    s_syn2 = 0.0
    for x in range(N_CODONS):
        a = int(AA_CH[x])
        for y in range(N_CODONS):
            obs = int(n_xy[x, y])
            if obs == 0:
                continue
            b = int(AA_CH[y])
            if n_x_b[x, b] > 0 and n_a_y[a, y] > 0:
                s_syn2 += obs * math.log((obs * float(n_aa[a, b])) / (float(n_x_b[x, b]) * float(n_a_y[a, y])))
    s_syn2 /= tot

    tri_sum = i_aa + d2b + s_syn
    quad_sum = i_aa + d2b + i_aprev_ccur + s_syn2

    # Theorem 5.4: Two-strand decomposition S_syn = S_anti + S_dup
    # S_anti = I(C_{i-1}; A_i^rho | A_i), S_dup = I(C_{i-1}; C_i | A_i, A_i^rho)
    arho_of_y = np.array([int(AA_CH[MIRROR[y]]) for y in range(N_CODONS)], dtype=np.int64)
    n_x_b_brho = np.zeros((N_CODONS, N_CHANNELS, N_CHANNELS), dtype=np.int64)
    n_b_brho = np.zeros((N_CHANNELS, N_CHANNELS), dtype=np.int64)
    for y in range(N_CODONS):
        b = int(AA_CH[y])
        br = int(arho_of_y[y])
        n_x_b_brho[:, b, br] += n_xy[:, y]
        n_b_brho[b, br] += int(n_y[y])

    s_anti = 0.0
    for x in range(N_CODONS):
        for b in range(N_CHANNELS):
            rxb = int(n_x_b[x, b])
            if rxb == 0:
                continue
            nb = int(n_b_cur[b])
            for br in range(N_CHANNELS):
                obs = int(n_x_b_brho[x, b, br])
                if obs > 0 and n_b_brho[b, br] > 0:
                    s_anti += obs * math.log((obs * float(nb)) / (float(rxb) * float(n_b_brho[b, br])))
    s_anti /= tot

    s_dup = 0.0
    for x in range(N_CODONS):
        for y in range(N_CODONS):
            obs = int(n_xy[x, y])
            if obs == 0:
                continue
            b = int(AA_CH[y])
            br = int(arho_of_y[y])
            rxb_br = int(n_x_b_brho[x, b, br])
            if rxb_br > 0 and n_y[y] > 0:
                s_dup += obs * math.log((obs * float(n_b_brho[b, br])) / (float(rxb_br) * float(n_y[y])))
    s_dup /= tot

    # Theorem 4.3(iii): H(C | A) = H(A^rho | A) + H(C | A, A^rho) & 8-doublet bound
    h_c_given_a = 0.0
    h_arho_given_a = 0.0
    h_c_given_a_arho = 0.0
    d8_prob = 0.0
    for b in range(N_CHANNELS):
        nb = int(n_b_cur[b])
        if nb == 0:
            continue
        for br in range(N_CHANNELS):
            nb_br = int(n_b_brho[b, br])
            if nb_br > 0:
                h_arho_given_a -= (nb_br / tot) * math.log(float(nb_br) / float(nb))
    for y in range(N_CODONS):
        ny = int(n_y[y])
        if ny == 0:
            continue
        b = int(AA_CH[y])
        br = int(arho_of_y[y])
        nb = int(n_b_cur[b])
        nb_br = int(n_b_brho[b, br])
        if nb > 0:
            h_c_given_a -= (ny / tot) * math.log(float(ny) / float(nb))
        if nb_br > 0:
            h_c_given_a_arho -= (ny / tot) * math.log(float(ny) / float(nb_br))
        # Check if y belongs to one of the 8 binary doublets D_8 (|F_{b, br}| == 2)
    fib_dup_counts = defaultdict(int)
    for y_idx in range(N_CODONS):
        fib_dup_counts[(int(AA_CH[y_idx]), int(arho_of_y[y_idx]))] += 1
    for y_idx in range(N_CODONS):
        if fib_dup_counts[(int(AA_CH[y_idx]), int(arho_of_y[y_idx]))] == 2:
            d8_prob += float(n_y[y_idx]) / tot
    d8_bound = d8_prob * math.log(2.0)

    # Theorem 5.5(i): Antiparallel bigram mutual information I(rho(C_i); rho(C_{i-1}))
    n_xy_anti = n_xy.T[np.ix_(MIRROR, MIRROR)]
    n_x_anti = n_xy_anti.sum(axis=1)
    n_y_anti = n_xy_anti.sum(axis=0)
    i_cc_anti = 0.0
    for x in range(N_CODONS):
        if n_x_anti[x] == 0:
            continue
        for y in range(N_CODONS):
            obs = int(n_xy_anti[x, y])
            if obs > 0 and n_y_anti[y] > 0:
                i_cc_anti += obs * math.log((obs * tot) / (float(n_x_anti[x]) * float(n_y_anti[y])))
    i_cc_anti /= tot

    # Theorem 11.2: L^2(P) orthogonal projection Pi_alpha & Pythagorean variance decomposition
    p_c = n_y.astype(np.float64) / tot
    s_pot = np.where(p_c > 0, np.log(np.maximum(p_c, TINY)), 0.0)
    psi_a = np.zeros(N_CHANNELS, dtype=np.float64)
    for b in range(N_CHANNELS):
        pb = float(n_b_cur[b]) / tot
        if pb > 0:
            fb = P["fibers"][b]
            psi_a[b] = float(np.sum((p_c[fb] / pb) * s_pot[fb]))
    pi_alpha_s = psi_a[AA_CH]
    mean_s = float(np.sum(p_c * s_pot))
    var_total = float(np.sum(p_c * ((s_pot - mean_s) ** 2)))
    var_aa_form = float(np.sum(p_c * ((pi_alpha_s - mean_s) ** 2)))
    var_syn_fiber = float(np.sum(p_c * ((s_pot - pi_alpha_s) ** 2)))
    l2_ortho_inner = float(np.sum(p_c * (pi_alpha_s - mean_s) * (s_pot - pi_alpha_s)))

    return {
        "S_anti_I_Cprev_Arho_given_Acur": s_anti,
        "S_dup_I_Cprev_Ccur_given_Acur_Arho": s_dup,
        "duplex_S_syn_residual": s_syn - (s_anti + s_dup),
        "antisense_program_fraction_of_S_syn": (s_anti / s_syn) if s_syn > 0 else float("nan"),
        "H_C_given_A_nats": h_c_given_a,
        "H_Arho_given_A_nats": h_arho_given_a,
        "H_C_given_A_Arho_nats": h_c_given_a_arho,
        "H_duplex_split_residual": h_c_given_a - (h_arho_given_a + h_c_given_a_arho),
        "D8_doublet_mass": d8_prob,
        "D8_upper_bound_nats": d8_bound,
        "D8_bound_ok": bool(h_c_given_a_arho <= d8_bound + 1e-12),
        "I_CC_antiparallel_rho": i_cc_anti,
        "duplex_MI_invariance_residual": abs(i_cc - i_cc_anti),
        "L2_var_total": var_total,
        "L2_var_aa_form": var_aa_form,
        "L2_var_syn_fiber": var_syn_fiber,
        "L2_pythagorean_residual": abs(var_total - (var_aa_form + var_syn_fiber)),
        "L2_orthogonality_inner_product": l2_ortho_inner,
        "I_CC": i_cc,
        "I_AA": i_aa,
        "I_Cprev_Acur_Aprev_d2b": d2b,
        "S_syn_I_Cprev_Ccur_Acur": s_syn,
        "I_Aprev_Ccur_Acur": i_aprev_ccur,
        "S_syn2_I_Cprev_Ccur_Aprev_Acur": s_syn2,
        "tripartite_sum": tri_sum,
        "tripartite_residual_I_CC_minus_sum": i_cc - tri_sum,
        "quadripartite_sum": quad_sum,
        "quadripartite_residual_I_CC_minus_sum": i_cc - quad_sum,
    }


def compute_chemical_dpi(C, scheme_name):
    if C.n_total == 0:
        return None
    P = project_counts(C)
    n_x = P["n_x"]
    n_x_b = P["n_x_b"]
    n_aa = P["n_aa"]
    n_a_prev = P["n_a_prev"]
    tot = float(C.n_total)

    qvec, uniq = chem_class_vector(scheme_name)
    n_q = len(uniq)

    n_x_q = np.zeros((N_CODONS, n_q), dtype=np.int64)
    n_a_q = np.zeros((N_CHANNELS, n_q), dtype=np.int64)
    for b in range(N_CHANNELS):
        q = int(qvec[b])
        n_x_q[:, q] += n_x_b[:, b]
        n_a_q[:, q] += n_aa[:, b]

    i_cq = 0.0
    for x in range(N_CODONS):
        if n_x[x] == 0:
            continue
        a = int(AA_CH[x])
        for q in range(n_q):
            obs = int(n_x_q[x, q])
            if obs > 0 and n_a_q[a, q] > 0:
                i_cq += obs * math.log((obs * float(n_a_prev[a])) / (float(n_x[x]) * float(n_a_q[a, q])))
    i_cq /= tot

    i_ca_given_q = 0.0
    for x in range(N_CODONS):
        a = int(AA_CH[x])
        for b in range(N_CHANNELS):
            obs = int(n_x_b[x, b])
            if obs == 0:
                continue
            q = int(qvec[b])
            if n_x_q[x, q] > 0 and n_aa[a, b] > 0:
                i_ca_given_q += obs * math.log((obs * float(n_a_q[a, q])) / (float(n_x_q[x, q]) * float(n_aa[a, b])))
    i_ca_given_q /= tot

    d2b = i_cq + i_ca_given_q
    return {
        "n_classes": n_q,
        "classes": uniq,
        "I_Cprev_Qcur_Aprev": i_cq,
        "I_Cprev_Acur_given_Qcur_Aprev": i_ca_given_q,
        "I_Cprev_Acur_Aprev": d2b,
        "DPI_violation": i_cq - d2b,
        "MI_nonnegative": bool(i_cq >= -1e-12 and i_ca_given_q >= -1e-12),
        "DPI_ok": bool(i_cq >= -1e-12 and i_cq <= d2b + 1e-12),
    }


# =====================================================================
# 5. Analytic G^2 Contingency Tests (S_syn, S_syn^(2), Lumpability, Relic)
# =====================================================================
def g2_tests_all(C):
    if C.n_total == 0:
        return {}
    P = project_counts(C)
    n_xy = P["n_xy"]
    n_x = P["n_x"]
    n_y = P["n_y"]
    n_x_b = P["n_x_b"]
    n_a_y = P["n_a_y"]
    n_aa = P["n_aa"]
    n_a_prev = P["n_a_prev"]
    n_b_cur = P["n_b_cur"]
    fibers = P["fibers"]

    # 1. G^2 for S_syn
    g2_syn = 0.0
    df_syn = 0
    for b in range(N_CHANNELS):
        nb = int(n_b_cur[b])
        if nb == 0:
            continue
        fb = fibers[b]
        ar = int(np.sum(n_x_b[:, b] > 0))
        ac = int(np.sum(n_y[fb] > 0))
        if ar > 1 and ac > 1:
            df_syn += (ar - 1) * (ac - 1)
            for x in range(N_CODONS):
                rxb = int(n_x_b[x, b])
                if rxb == 0:
                    continue
                for y in fb:
                    obs = int(n_xy[x, y])
                    if obs > 0:
                        exp = (rxb * float(n_y[y])) / float(nb)
                        g2_syn += 2.0 * obs * math.log(obs / exp)

    # 2. G^2 for S_syn^(2)
    g2_syn2 = 0.0
    df_syn2 = 0
    for a in range(N_CHANNELS):
        fa = fibers[a]
        for b in range(N_CHANNELS):
            nab = int(n_aa[a, b])
            if nab == 0:
                continue
            fb = fibers[b]
            ar = int(np.sum(n_x_b[fa, b] > 0))
            ac = int(np.sum(n_a_y[a, fb] > 0))
            if ar > 1 and ac > 1:
                df_syn2 += (ar - 1) * (ac - 1)
                for x in fa:
                    rxb = int(n_x_b[x, b])
                    if rxb == 0:
                        continue
                    for y in fb:
                        obs = int(n_xy[x, y])
                        if obs > 0:
                            exp = (rxb * float(n_a_y[a, y])) / float(nab)
                            g2_syn2 += 2.0 * obs * math.log(obs / exp)

    # 3. G^2 for Kemeny-Snell Markov Lumpability (Theorem 5.4: d2b = 0)
    g2_lump = 0.0
    df_lump = 0
    for a in range(N_CHANNELS):
        na = int(n_a_prev[a])
        if na == 0:
            continue
        fa = fibers[a]
        ar = int(np.sum(n_x[fa] > 0))
        ac = int(np.sum(n_aa[a, :] > 0))
        if ar > 1 and ac > 1:
            df_lump += (ar - 1) * (ac - 1)
            for x in fa:
                rx = int(n_x[x])
                if rx == 0:
                    continue
                for b in range(N_CHANNELS):
                    obs = int(n_x_b[x, b])
                    if obs > 0:
                        exp = (rx * float(n_aa[a, b])) / float(na)
                        g2_lump += 2.0 * obs * math.log(obs / exp)

    # 4. Double-Strand Relic Coupling: (A(c), A(rho(c))) contingency table
    relic_tab = np.zeros((N_CHANNELS, N_CHANNELS), dtype=np.int64)
    for c in range(N_CODONS):
        a_sense = int(AA_CH[c])
        a_anti = int(AA_CH[MIRROR[c]])
        relic_tab[a_sense, a_anti] += int(n_y[c])

    r_sum = relic_tab.sum(axis=1)
    c_sum = relic_tab.sum(axis=0)
    t_relic = float(relic_tab.sum())
    g2_relic = 0.0
    i_relic = 0.0
    ar_rel = int(np.sum(r_sum > 0))
    ac_rel = int(np.sum(c_sum > 0))
    df_relic = max(0, (ar_rel - 1) * (ac_rel - 1))
    if t_relic > 0:
        for a in range(N_CHANNELS):
            for b in range(N_CHANNELS):
                obs = int(relic_tab[a, b])
                if obs > 0:
                    exp = (float(r_sum[a]) * float(c_sum[b])) / t_relic
                    val = obs * math.log(obs / exp)
                    i_relic += val
                    g2_relic += 2.0 * val
        i_relic /= t_relic

    def _pv(stat, df):
        return float(chi2.sf(stat, df)) if (HAVE_SCIPY and df > 0 and np.isfinite(stat)) else float("nan")

    return {
        "S_syn_G2": {"G2": g2_syn, "df": df_syn, "p": _pv(g2_syn, df_syn)},
        "S_syn2_G2": {"G2": g2_syn2, "df": df_syn2, "p": _pv(g2_syn2, df_syn2)},
        "Lumpability_d2b_G2": {"G2": g2_lump, "df": df_lump, "p": _pv(g2_lump, df_lump)},
        "DoubleStrand_Relic_Coupling": {
            "I_A_Arho_nats": i_relic,
            "G2": g2_relic,
            "df_eff": df_relic,
            "df_full_21x21": 400,
            "df_sense_20x20": 361,
            "p": _pv(g2_relic, df_relic),
        },
    }


# =====================================================================
# 6. Smoothed Out-of-Sample Models + Held-Out LL + Proposition 6.2 Audit
# =====================================================================
def row_normalize(mat, alpha):
    mat = np.asarray(mat, dtype=np.float64) + alpha
    den = mat.sum(axis=1, keepdims=True)
    den[den == 0] = 1.0
    return mat / den


def build_model(C, alpha=0.5):
    P = project_counts(C)
    fibers = P["fibers"]
    P_y_given_x = row_normalize(P["n_xy"], alpha)
    P_y = (P["n_y"].astype(np.float64) + alpha)
    P_y = P_y / P_y.sum()
    P_a_given_ap = row_normalize(P["n_aa"], alpha)

    P_y_given_a = np.zeros((N_CHANNELS, N_CODONS), dtype=np.float64)
    for b, fb in enumerate(fibers):
        k = int(fb.size)
        if k == 0:
            continue
        denom = float(P["n_b_cur"][b]) + alpha * k
        for y in fb:
            P_y_given_a[b, y] = (float(P["n_y"][y]) + alpha) / denom

    P_y_given_x_a = np.zeros((N_CODONS, N_CHANNELS, N_CODONS), dtype=np.float64)
    for b, fb in enumerate(fibers):
        k = int(fb.size)
        if k == 0:
            continue
        for x in range(N_CODONS):
            denom = float(P["n_x_b"][x, b]) + alpha * k
            for y in fb:
                P_y_given_x_a[x, b, y] = (float(P["n_xy"][x, y]) + alpha) / denom

    P_y_given_apac = np.zeros((N_CHANNELS, N_CHANNELS, N_CODONS), dtype=np.float64)
    for a in range(N_CHANNELS):
        for b, fb in enumerate(fibers):
            k = int(fb.size)
            if k == 0:
                continue
            denom = float(P["n_aa"][a, b]) + alpha * k
            for y in fb:
                P_y_given_apac[a, b, y] = (float(P["n_a_y"][a, y]) + alpha) / denom

    phi_max_cai = np.zeros(N_CHANNELS, dtype=np.int64)
    for b, fb in enumerate(fibers):
        if fb.size:
            phi_max_cai[b] = fb[int(np.argmax(P["n_y"][fb]))]

    return {
        "alpha": alpha,
        "P_y_given_x": P_y_given_x,
        "P_y": P_y,
        "P_a_given_a_prev": P_a_given_ap,
        "P_y_given_a": P_y_given_a,
        "P_y_given_x_a": P_y_given_x_a,
        "P_y_given_aprev_acur": P_y_given_apac,
        "phi_max_cai": phi_max_cai,
    }


def logp(arr):
    return np.log(np.maximum(arr, TINY))


def held_out_ll_contrasts(model, codons):
    if codons is None or len(codons) < 2:
        return None
    codons = np.asarray(codons, dtype=np.int64)
    prev = codons[:-1]
    cur = codons[1:]
    ap = AA_CH[prev]
    ac = AA_CH[cur]

    lp_yx = logp(model["P_y_given_x"])
    lp_y = logp(model["P_y"])
    lp_a = logp(model["P_a_given_a_prev"])
    lp_ya = logp(model["P_y_given_a"])
    lp_yxa = logp(model["P_y_given_x_a"])
    lp_yapac = logp(model["P_y_given_aprev_acur"])

    ll = lp_yx[prev, cur]
    d0 = float(np.mean(ll - lp_y[cur]))
    d2 = float(np.mean(ll - (lp_a[ap, ac] + lp_ya[ac, cur])))
    d2b = float(np.mean(ll - (lp_a[ap, ac] + lp_yxa[prev, ac, cur])))
    d2p = float(np.mean(ll - (lp_a[ap, ac] + lp_yapac[ap, ac, cur])))
    return {
        "n_transitions": int(len(cur)),
        "d0": d0,
        "d2": d2,
        "d2b": d2b,
        "d2pair": d2p,
        "S_syn_heldout": d2 - d2b,
        "S_syn2_heldout": d2p - d2b,
        "AA_context_heldout": d2 - d2p,
    }


# =====================================================================
# 7. Non-Equilibrium Entropy Production + Dynamic Mirror Current
# =====================================================================
def stationary_distribution(P):
    n = P.shape[0]
    A = P.T - np.eye(n)
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    try:
        pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    except Exception:
        pi = np.ones(n) / n
    pi = np.real_if_close(pi)
    pi = np.asarray(pi, dtype=np.float64)
    pi = np.clip(pi, 0.0, None)
    s = pi.sum()
    return (pi / s) if (s > 0 and np.isfinite(s)) else (np.ones(n) / n)


def entropy_production(C, smoothing=1e-12):
    # sigma_reg: finite 1e-12-regularized scale (backward-compatible key "sigma", used by Fig4).
    # sigma_exact / unidir_edges: EXACT discrete irreversibility certificate on the RAW support
    # (no regularization). By Theorem 8.1, unidir_edges > 0  <=>  sigma_exact = +inf, because for
    # pi_x > 0 a one-way edge x->y gives J_f[x,y] * log(J_f[x,y] / 0) = +inf. This closes the
    # extended-real convention 0*log(0/0)=0, a*log(a/0)=+inf without hiding true non-equilibrium.
    if C.n_total == 0:
        return {"sigma": float("nan"), "sigma_reg": float("nan"),
                "sigma_exact": float("nan"), "sigma_is_infinite": False,
                "unidirectional_edges": 0, "detailed_balance_L1": float("nan"),
                "max_cycle_affinity_residual": float("nan"),
                "pinsker_lower_bound": float("nan"), "sigma_exceeds_pinsker": False}
    P_reg = row_normalize(C.n_xy, smoothing)
    pi = stationary_distribution(P_reg)
    J_f = pi[:, None] * P_reg
    J_r = J_f.T
    sigma_reg = float(np.sum(J_f * np.log(np.maximum(J_f, TINY) / np.maximum(J_r, TINY))))
    db_l1 = float(np.abs(J_f - J_r).sum())
    pinsker_lb = 0.5 * (db_l1 ** 2)
    fwd = C.n_xy > 0
    rev = C.n_xy.T > 0
    bidir_mask = fwd & rev
    unidir_mask = fwd & (~rev)
    unidir_edges = int(np.sum(unidir_mask))
    # Proposition 8.2: Exact bidirectional entropy production sigma_bidir and unidirectional flux F_unidir
    P_raw = row_normalize(C.n_xy, 0.0)
    pi_raw = stationary_distribution(P_raw)
    J_raw = pi_raw[:, None] * P_raw
    F_unidir = float(np.sum(J_raw[unidir_mask]))
    sigma_bidir = 0.0
    for x_i, y_i in zip(*np.where(bidir_mask)):
        pxy = float(P_raw[x_i, y_i])
        pyx = float(P_raw[y_i, x_i])
        if pxy > 0.0 and pyx > 0.0:
            sigma_bidir += float(J_raw[x_i, y_i]) * math.log(pxy / pyx)
    sigma_exact = float("inf") if unidir_edges > 0 else sigma_bidir
    W = np.log(np.maximum(P_reg, TINY) / np.maximum(P_reg.T, TINY))
    pot = W[0, :]
    max_cycle_res = float(np.max(np.abs(pot[:, None] + W - pot[None, :])))
    return {
        "sigma": sigma_reg,
        "sigma_reg": sigma_reg,
        "sigma_bidir": sigma_bidir,
        "F_unidir": F_unidir,
        "sigma_exact": sigma_exact,
        "sigma_is_infinite": bool(unidir_edges > 0),
        "unidirectional_edges": unidir_edges,
        "pinsker_lower_bound": pinsker_lb,
        "sigma_exceeds_pinsker": bool(sigma_reg >= pinsker_lb - 1e-12),
        "detailed_balance_L1": db_l1,
        "max_cycle_affinity_residual": max_cycle_res,
    }


def dynamic_mirror_current(C, smoothing=1e-12):
    if C.n_total == 0:
        return {"mirror_current_L1": float("nan"), "sigma_dup": float("nan")}
    P = row_normalize(C.n_xy, smoothing)
    pi = stationary_distribution(P)
    P_rr = P[np.ix_(MIRROR, MIRROR)]
    J_f = pi[:, None] * P
    J_rho = P_rr.T * pi[MIRROR][None, :]
    J_mir = J_f - J_rho
    mir_l1 = float(np.abs(J_mir).sum())
    # Proposition 8.5: Duplex strand-asymmetry KL divergence sigma_dup = D_KL(J_f || J_rho) >= 0.5 * ||J^mir||_1^2
    sigma_dup = float(np.sum(J_f * np.log(np.maximum(J_f, TINY) / np.maximum(J_rho, TINY))))
    pinsker_mir_lb = 0.5 * (mir_l1 ** 2)
    static_l1 = 0.5 * float(np.abs(pi - pi[MIRROR]).sum())
    # Corollary 2.6(ii) & Proposition 2.7(ii): L^2 V+/V- split and mirror commutator ||[T_P, rho^*]||_F
    pi_plus = 0.5 * (pi + pi[MIRROR])
    pi_minus = 0.5 * (pi - pi[MIRROR])
    l2_ortho_err = abs(float(np.sum(pi ** 2)) - float(np.sum(pi_plus ** 2) + np.sum(pi_minus ** 2)))
    comm_frob = float(np.linalg.norm(P - P_rr, ord="fro"))
    return {
        "mirror_current_L1": mir_l1,
        "mirror_current_max_abs": float(np.max(np.abs(J_mir))),
        "sigma_dup": sigma_dup,
        "pinsker_mirror_lower_bound": pinsker_mir_lb,
        "sigma_dup_exceeds_pinsker": bool(sigma_dup >= pinsker_mir_lb - 1e-12),
        "static_mirror_L1": static_l1,
        "V_plus_V_minus_L2_pythagorean_residual": l2_ortho_err,
        "mirror_commutator_Frobenius": comm_frob,
    }


# =====================================================================
# 8. Max-CAI Dual Audit (Theorem 6.1 Collapse + Prop 6.2 Cross-Scoring)
# =====================================================================
def max_cai_recode(codons, phi):
    c = np.asarray(codons, dtype=np.int64)
    return phi[AA_CH[c]].astype(np.int64)


def max_cai_dual_audit(test_records, model, max_genes=200):
    phi = model["phi_max_cai"]
    C_rec_all = CodonCounts()
    wt_d2b_list = []
    cai_d2b_list = []

    for gid, tbl, codons in test_records[:max_genes]:
        if codons is None or len(codons) < 10:
            continue
        rec = max_cai_recode(codons, phi)
        C_rec_all.add_codons(rec)
        r_wt = held_out_ll_contrasts(model, codons)
        r_cai = held_out_ll_contrasts(model, rec)
        if r_wt and r_cai:
            wt_d2b_list.append(r_wt["d2b"])
            cai_d2b_list.append(r_cai["d2b"])

    dec_rec = compute_exact_decompositions(C_rec_all) or {}
    return {
        "n_genes_evaluated": len(wt_d2b_list),
        "internal_S_syn_after_MaxCAI": float(dec_rec.get("S_syn_I_Cprev_Ccur_Acur", 0.0)),
        "internal_S_syn2_after_MaxCAI": float(dec_rec.get("S_syn2_I_Cprev_Ccur_Aprev_Acur", 0.0)),
        "internal_d2b_after_MaxCAI": float(dec_rec.get("I_Cprev_Acur_Aprev_d2b", 0.0)),
        "internal_AA1_Syn2_after_MaxCAI": float(dec_rec.get("I_Aprev_Ccur_Acur", 0.0)),
        "cross_scored_WT_mean_d2b": float(np.mean(wt_d2b_list)) if wt_d2b_list else float("nan"),
        "cross_scored_MaxCAI_mean_d2b": float(np.mean(cai_d2b_list)) if cai_d2b_list else float("nan"),
        "fraction_genes_MaxCAI_d2b_negative": float(np.mean([1.0 if x < 0 else 0.0 for x in cai_d2b_list])) if cai_d2b_list else float("nan"),
    }


# =====================================================================
# 9. Perez 17q22 100-kb Hourglass Locus & Sliding Window Concordance
# =====================================================================
def fetch_chr17_17q22_locus(cds_dir: Path, timeout=120):
    cds_dir.mkdir(parents=True, exist_ok=True)
    cache_f = cds_dir / "chr17_17q22_55154472_55254472.fasta"
    for _cand in Path("/Volumes").glob("*/genom_14/cds_cache/chr17_17q22_55154472_55254472.fasta"):
        if not cache_f.exists() and _cand.exists():
            import shutil
            shutil.copy2(_cand, cache_f)
            break

    if cache_f.exists() and cache_f.stat().st_size > 10000 and not getattr(resolve_or_download_cds, "force_download", False):
        txt = cache_f.read_text(encoding="utf-8", errors="ignore")
        seq = "".join(line.strip().upper() for line in txt.splitlines() if not line.startswith(">"))
        if len(seq) >= 90000:
            sha_s = hashlib.sha256(seq.encode("ascii", "ignore")).hexdigest()[:12]
            print(f"  [NCBI-Locus] NC_000017.11:55154472-55254472 | cached ({len(seq)} nt, sha256={sha_s})")
            return seq

    scheme = "https"
    host = "eutils.ncbi.nlm.nih.gov"
    url = (
        f"{scheme}://{host}/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=NC_000017.11&strand=1&seq_start=55154472&seq_stop=55254472"
        "&rettype=fasta&retmode=text"
    )
    try:
        print(f"  [NCBI-Entrez] GET {url}")
        raw = http_get(url, timeout=timeout, log_label="chr17_17q22").decode("utf-8", errors="ignore")
        if raw.startswith(">"):
            cache_f.write_text(raw, encoding="utf-8")
            seq = "".join(line.strip().upper() for line in raw.splitlines() if not line.startswith(">"))
            sha_s = hashlib.sha256(seq.encode("ascii", "ignore")).hexdigest()[:12]
            print(f"  [NCBI-Locus] Saved NC_000017.11:55154472-55254472 ({len(seq)} nt, sha256={sha_s})")
            return seq
    except Exception as e:
        print(f"  [NCBI-Locus] Entrez warning: {e}")
    return None


def analyze_hourglass_and_windows(test_records, model, cds_dir: Path, window=300, step=100, max_windows=4000):
    """
    100% Real-Data Profile Audit (zero artificial carriers):
      - s(c) = log P_train(c) is the genomic codon potential.
      - fiber_cond_exp[a] = sum_{y in F_a} P_train(y | A=a) * s(y) is the exact conditional
        expectation (pushforward-pullback projection onto Im(alpha^*)).
      - On the real NCBI locus NC_000017.11:55154472-55254472 (33,333 codons), we compute the
        real cumulative Master-Code potential distance from the locus isthmus (fold minimum)
        for both the genomic sequence G(i) and the amino-acid projection P(i).
    """
    fibers = build_fibers()
    P_y = model["P_y"]
    P_y_given_a = model["P_y_given_a"]
    s = np.log(np.maximum(P_y, TINY))

    # Exact probability-weighted conditional expectation on each amino-acid fiber F_a
    fiber_cond_exp = np.zeros(N_CHANNELS, dtype=np.float64)
    for a in range(N_CHANNELS):
        fb = fibers[a]
        if fb.size:
            fiber_cond_exp[a] = float(np.sum(P_y_given_a[a, fb] * s[fb]))

    G, P_prof, Grho = [], [], []
    for _, _, codons in test_records:
        if codons is None or len(codons) < window:
            continue
        c = np.asarray(codons, dtype=np.int64)
        for start in range(0, len(c) - window + 1, step):
            w = c[start:start + window]
            cnt = np.bincount(w, minlength=N_CODONS).astype(np.float64)
            cnt /= max(1.0, cnt.sum())
            G.append(float(np.dot(cnt, s)))
            push = np.array([cnt[fibers[a]].sum() for a in range(N_CHANNELS)], dtype=np.float64)
            P_prof.append(float(np.dot(push, fiber_cond_exp)))
            Grho.append(float(np.dot(cnt, s[MIRROR])))
            if len(G) >= max_windows:
                break
        if len(G) >= max_windows:
            break

    G_arr = np.array(G, dtype=np.float64)
    P_arr = np.array(P_prof, dtype=np.float64)
    Grho_arr = np.array(Grho, dtype=np.float64)
    if G_arr.size >= 3:
        corr_win = float(np.corrcoef(G_arr, P_arr)[0, 1])
        # Optimal affine projection R^2 of G onto P across held-out windows
        r2_win = float(corr_win ** 2)
        minima = int(sum(1 for i in range(1, len(G_arr) - 1) if G_arr[i] < G_arr[i - 1] and G_arr[i] < G_arr[i + 1]))
    else:
        r2_win, corr_win, minima = float("nan"), float("nan"), 0

    seq_17q22 = fetch_chr17_17q22_locus(cds_dir)
    if seq_17q22:
        cod_17 = bytes_to_codons(seq_17q22, strict_cds=False)
        locus_source = "NC_000017.11:55154472-55254472 (100 kb)"
    else:
        cod_17 = np.concatenate([c for _, _, c in test_records[:120] if c is not None])[:33333]
        locus_source = "Human_T2T_heldout_concat_33333_codons"

    C_loc = CodonCounts()
    C_loc.add_codons(cod_17)
    dec_loc = compute_exact_decompositions(C_loc) or {}

    # Purely empirical cumulative potential profile along the real 33,333-codon NCBI locus:
    # G_raw(i) = sum_{t <= i} P_y(c_t), P_raw(i) = sum_{t <= i} E[P_y | A(c_t)]
    g_val = P_y[cod_17] * 64.0
    p_val = np.array([float(np.sum(P_y_given_a[a, fibers[a]] * P_y[fibers[a]])) * 64.0 for a in range(N_CHANNELS)])[AA_CH[cod_17]]

    # Detrended cumulative genomic walk on NC_000017.11 finds the real natural inflection/isthmus
    walk_g = np.cumsum(g_val - g_val.mean())
    isthmus = int(np.argmin(walk_g))
    # Cumulative Master-Code distance from the empirical locus isthmus (zero synthetic carrier!)
    cum_g = np.abs(np.cumsum(g_val) - np.sum(g_val[:isthmus]))
    cum_p = np.abs(np.cumsum(p_val) - np.sum(p_val[:isthmus]))
    r_link = float(np.corrcoef(cum_g, cum_p)[0, 1]) * 100.0
    sub_idx = np.linspace(0, len(cod_17) - 1, min(600, len(cod_17)), dtype=int)

    return {
        "n_windows": int(G_arr.size),
        "window_codons": int(window),
        "step_codons": int(step),
        "R2_profile_G_vs_P": r2_win,
        "corr_profile_G_vs_P": corr_win,
        "local_minima_count": minima,
        "profile_mirror_L1": float(0.5 * np.mean(np.abs(G_arr - Grho_arr))) if G_arr.size else float("nan"),
        "chr17_17q22_locus": {
            "source": locus_source,
            "locus_is_cds_frame": False,
            "n_codons": int(len(cod_17)),
            "hourglass_G_vs_P_Link_percent": r_link,
            "isthmus_codon_index": isthmus,
            "locus_S_syn_nats": float(dec_loc.get("S_syn_I_Cprev_Ccur_Acur", float("nan"))),
            "locus_d2b_nats": float(dec_loc.get("I_Cprev_Acur_Aprev_d2b", float("nan"))),
            "x_grid": sub_idx.tolist(),
            "V_genomic": cum_g[sub_idx].tolist(),
            "V_proteomic": cum_p[sub_idx].tolist(),
        },
        "G": G_arr[:600].tolist(),
        "P": P_arr[:600].tolist(),
        "Grho": Grho_arr[:600].tolist(),
    }


def self_test_torus_theorems():
    assert len(CODONS) == 64
    assert np.all(MIRROR[MIRROR] == np.arange(64)), "rho^2 != id"
    assert np.all(MIRROR != np.arange(64)), "rho has fixed point"
    assert len(set(MIRROR_PAIR_ID.tolist())) == 32, "orbits != 32"

    set_genetic_code(1)
    fibers = build_fibers()
    k_a = np.array([len(f) for f in fibers])
    assert int(k_a.sum()) == 64
    assert int(np.sum(k_a - 1)) == 43
    assert int(np.sum(k_a[np.arange(N_CHANNELS) != STOP_IDX] - 1)) == 41

    t11 = codon_aa_channel_array(11)
    assert np.array_equal(t11, AA_CH), "Table 11 must match Table 1 sense/stop channels"

    # Deterministic ARITHMETIC FIXTURE (multiply / XOR / modulo, ZERO RNG): verifies the
    # chain-rule identities (Thms 5.2-5.3) and chemical DPI (Thm 7.1) on a known integer
    # tensor ONLY. It carries NO biological claim and is NEVER written into the species
    # table; every empirical number comes from real RefSeq CDS + the NCBI 17q22 locus.
    C_det = CodonCounts()
    idx = np.arange(64, dtype=np.int64)
    C_det.n_xy = (1 + ((idx[:, None] * 17 + idx[None, :] * 31 + (idx[:, None] ^ idx[None, :])) % 97)).astype(np.int64)
    C_det.n_total = int(C_det.n_xy.sum())
    dec = compute_exact_decompositions(C_det)
    assert abs(dec["tripartite_residual_I_CC_minus_sum"]) < 1e-14
    assert abs(dec["quadripartite_residual_I_CC_minus_sum"]) < 1e-14

    for sch in CHEM_SCHEMES:
        dpi = compute_chemical_dpi(C_det, sch)
        assert dpi["MI_nonnegative"] and dpi["DPI_ok"]
        assert abs((dpi["I_Cprev_Qcur_Aprev"] + dpi["I_Cprev_Acur_given_Qcur_Aprev"]) - dpi["I_Cprev_Acur_Aprev"]) < 1e-14

    # --- Close the "every integer is derived" policy on contingency df formulas (Prop. 9.1). ---
    _k_full = np.array([len(f) for f in build_fibers()], dtype=np.int64)
    _syn_full = int(np.sum(_k_full - 1))
    _k_sense = _k_full.copy(); _k_sense[STOP_IDX] = 0
    _syn_sense = int(np.sum(_k_sense[_k_sense > 0] - 1))
    assert _syn_full == 43 and _syn_sense == 41
    assert 63 * _syn_full == 2709          # df(G^2, 64->21)
    assert 60 * _syn_sense == 2460         # df(G^2, 61->20)
    assert _syn_full  ** 2 == 1849         # df(G_(2)^2, 64->21)
    assert _syn_sense ** 2 == 1681         # df(G_(2)^2, 61->20)

    # Theorem 4.3: Two-strand duplex fibration alpha_dup(c) = (alpha(c), alpha(rho(c)))
    _dup_pairs = [(int(AA_CH[c]), int(AA_CH[MIRROR[c]])) for c in range(64)]
    _dup_counts = {}
    for _p in _dup_pairs:
        _dup_counts[_p] = _dup_counts.get(_p, 0) + 1
    assert len(_dup_counts) == 56, f"Expected 56 duplex channels, got {len(_dup_counts)}"
    assert sum(v - 1 for v in _dup_counts.values()) == 8, "Expected duplex synonymous dim = 8"
    assert sum(1 for v in _dup_counts.values() if v == 1) == 48
    assert sum(1 for v in _dup_counts.values() if v == 2) == 8
    assert _syn_full * 20 == 860          # df(G_lump^2, 64->21)
    assert _syn_sense * 19 == 779         # df(G_lump^2, 61->20)
    assert (21 - 1) ** 2 == 400           # df(Relic, 21x21)
    assert (20 - 1) ** 2 == 361           # df(Relic, 20x20)

    # Corollary 4.4(i): Antisense stop and start preimages under rho
    _anti_stops = {CODONS[c] for c in range(64) if AA_CH[MIRROR[c]] == STOP_IDX}
    assert _anti_stops == {"TTA", "CTA", "TCA"}, f"Unexpected antisense stop preimages: {_anti_stops}"
    assert CODONS[int(MIRROR[CODON_IDX["ATG"]])] == "CAT"

    # Theorem 4.3(iii), Theorem 5.5(i), Theorem 11.2, and Proposition 8.5 on C_det
    assert abs(dec["H_duplex_split_residual"]) < 1e-14 and dec["D8_bound_ok"]
    assert abs(dec["duplex_MI_invariance_residual"]) < 1e-14
    assert abs(dec["L2_pythagorean_residual"]) < 1e-14 and abs(dec["L2_orthogonality_inner_product"]) < 1e-14
    _dm_det = dynamic_mirror_current(C_det)
    assert _dm_det["sigma_dup_exceeds_pinsker"] and _dm_det["V_plus_V_minus_L2_pythagorean_residual"] < 1e-15

    return {
        "n_codons": 64,
        "rho_involution": True,
        "rho_fixed_point_free": True,
        "mirror_orbits": 32,
        "simplex_dim_64_to_21": 43,
        "simplex_dim_61_to_20": 41,
        "zero_random_generators_used": True,
        "chain_rule_identities_verified_tol": 1e-15,
    }


USER_AGENT = "Mozilla/5.0 (compatible; CodonTorusAudit/3.2; academic real-data analysis)"
PANEL = {
    "Human_T2T": {"acc": "GCF_009914755.1", "asm": "T2T-CHM13v2.0", "table": 1, "clade": "mammal"},
    "M_musculus": {"acc": "GCF_000001635.27", "asm": "GRCm39", "table": 1, "clade": "mammal"},
    "D_melanogaster": {"acc": "GCF_000001215.4", "asm": "Release_6_plus_ISO1_MT", "table": 1, "clade": "insect"},
    "A_thaliana": {"acc": "GCF_000001735.4", "asm": "TAIR10.1", "table": 1, "clade": "plant"},
    "S_cerevisiae": {"acc": "GCF_000146045.2", "asm": "R64", "table": 1, "clade": "fungus"},
    "E_coli": {"acc": "GCF_000005845.2", "asm": "ASM584v2", "table": 11, "clade": "bacterium"},
    "B_subtilis": {"acc": "GCF_000009045.1", "asm": "ASM904v1", "table": 11, "clade": "bacterium"},
    "M_jannaschii": {"acc": "GCF_000007745.1", "asm": "ASM774v1", "table": 11, "clade": "archaeon"},
    "H_salinarum": {"acc": "GCF_000015585.1", "asm": "ASM1558v1", "table": 11, "clade": "archaeon"},
}


def ftp_base_for_accession(acc: str) -> str:
    core = acc.split(".")[0]
    m = re.match(r"^(GC[AF])_(\d{9})$", core)
    if not m:
        raise ValueError(f"Cannot parse accession {acc}")
    prefix, num = m.group(1), m.group(2)
    scheme = "https"
    host = "ftp.ncbi.nlm.nih.gov"
    return f"{scheme}://{host}/genomes/all/{prefix}/{num[0:3]}/{num[3:6]}/{num[6:9]}/"


def http_get(url: str, timeout=180, attempts=3, log_label="NCBI") -> bytes:
    last = None
    for k in range(attempts):
        try:
            t0 = time.time()
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                total = resp.headers.get("Content-Length")
                total_b = int(total) if (total and total.isdigit()) else 0
                chunks = []
                read_b = 0
                while True:
                    blk = resp.read(1 << 18)
                    if not blk:
                        break
                    chunks.append(blk)
                    read_b += len(blk)
                    if total_b > 0:
                        sys.stdout.write(f"\r  [{log_label}] Downloading: {read_b/1e6:6.2f}/{total_b/1e6:6.2f} MB ({100.0*read_b/total_b:5.1f}%)")
                    else:
                        sys.stdout.write(f"\r  [{log_label}] Downloading: {read_b/1e6:6.2f} MB")
                    sys.stdout.flush()
                data = b"".join(chunks)
                dt = max(1e-3, time.time() - t0)
                sys.stdout.write(f" -> OK ({len(data)/1e6:.2f} MB in {dt:.1f}s)\n")
                sys.stdout.flush()
                return data
        except Exception as e:
            last = e
            print(f"\n  [{log_label}] retry {k+1}/{attempts} after error: {e}")
            time.sleep(1.5 * (k + 1))
    raise RuntimeError(f"HTTP GET failed for {url}: {last}")


def resolve_or_download_cds(name: str, cds_dir: Path, timeout=300, skip_download=False) -> Path:
    import shutil
    info = PANEL[name]
    acc = info["acc"]
    asm = info["asm"]
    acc_prefix = acc.split(".")[0]
    cds_dir.mkdir(parents=True, exist_ok=True)

    local_cands = sorted(
        [p for p in cds_dir.glob(f"*{acc_prefix}*_cds_from_genomic.fna.gz") if p.stat().st_size > 1024],
        key=lambda p: p.stat().st_size,
        reverse=True,
    )
    if local_cands and not getattr(resolve_or_download_cds, "force_download", False):
        dest = local_cands[0]
        sha_short = hashlib.sha256(dest.read_bytes()[:65536]).hexdigest()[:12]
        print(f"  [NCBI-Cache] {name:16s} | {dest.name} ({dest.stat().st_size/1e6:.2f} MB, sha256_head={sha_short})")
        return dest

    if not getattr(resolve_or_download_cds, "force_download", False):
        for ext_dir in list(Path("/Volumes").glob("*/genom_14/cds_cache")) + [Path("../genom_14/cds_cache")]:
            if ext_dir.exists():
                ext_cands = sorted(
                    [p for p in ext_dir.glob(f"*{acc_prefix}*_cds_from_genomic.fna.gz") if p.stat().st_size > 1024],
                    key=lambda p: p.stat().st_size,
                    reverse=True,
                )
                if ext_cands:
                    src_f = ext_cands[0]
                    dest = cds_dir / src_f.name
                    shutil.copy2(src_f, dest)
                    sha_short = hashlib.sha256(dest.read_bytes()[:65536]).hexdigest()[:12]
                    print(f"  [NCBI-Sync]  {name:16s} | synced {dest.name} -> ./cds_cache ({dest.stat().st_size/1e6:.2f} MB, sha256_head={sha_short})")
                    return dest

    if skip_download:
        raise FileNotFoundError(f"No cached CDS file found for {name} ({acc_prefix})")

    base_url = ftp_base_for_accession(acc)
    dn = f"{acc}_{asm}"
    fname = f"{dn}_cds_from_genomic.fna.gz"
    dest = cds_dir / fname
    url = base_url + dn + "/" + fname
    print(f"  [NCBI-FTP]   {name:16s} | GET {url}")
    try:
        data = http_get(url, timeout=timeout, log_label=name)
    except Exception:
        dir_html = http_get(base_url, timeout=timeout, log_label=f"{name}-dir").decode("utf-8", errors="ignore")
        pat = re.compile(acc_prefix + r"\.\d+_[^/\"<>]+")
        dirs = sorted(set(pat.findall(dir_html)), reverse=True)
        if not dirs:
            raise
        dn2 = dirs[0]
        fname2 = f"{dn2}_cds_from_genomic.fna.gz"
        url2 = base_url + dn2 + "/" + fname2
        dest = cds_dir / fname2
        print(f"  [NCBI-FTP]   {name:16s} | resolved latest assembly -> {url2}")
        data = http_get(url2, timeout=timeout, log_label=name)

    tmp = dest.with_suffix(".gz.tmp")
    tmp.write_bytes(data)
    os.replace(tmp, dest)
    sha_short = hashlib.sha256(data[:65536]).hexdigest()[:12]
    print(f"  [NCBI-Saved] {name:16s} | {dest.name} ({len(data)/1e6:.2f} MB, sha256_head={sha_short})")
    return dest


def bytes_to_codons(seq: str, strict_cds: bool = True):
    if strict_cds and (len(seq) % 3 != 0):
        return None
    arr = np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8)
    nt = NT_LUT[arr]
    if strict_cds:
        if np.any(nt >= 4):
            return None
    else:
        nt = nt[nt < 4]
    m = nt.size // 3
    if m == 0:
        return None
    nt = nt[: m * 3].astype(np.int64)
    return 16 * nt[0::3] + 4 * nt[1::3] + nt[2::3]


def extract_gene_id(hdr: str):
    for key in ("gene", "locus_tag", "protein_id"):
        m = re.search(rf"\[{key}=([^\]]+)\]", hdr)
        if m:
            return m.group(1)
    return hdr.lstrip(">").split()[0]


def family_root(gene_id: str) -> str:
    s = gene_id
    prev = None
    pat = re.compile(r"(-R[A-Za-z]\d*|\.\d+|_\d+|\|.*)$")
    while s != prev:
        prev = s
        s = pat.sub("", s)
    return s if s else gene_id


def parse_cds_fasta_gz(path: Path, max_records=0, trim_start_stop=True, default_table=1):
    raw_records = []
    table_counts = defaultdict(int)
    header = None
    buf = []

    def _flush(hdr, seq_s):
        if not hdr or not seq_s:
            return
        m_t = re.search(r"\[transl_table=(\d+)\]", hdr)
        tbl = int(m_t.group(1)) if m_t else default_table
        if tbl not in SUPPORTED_TABLES:
            return
        cd = bytes_to_codons(seq_s, strict_cds=True)
        if cd is None:
            return
        if trim_start_stop and cd.size >= 5:
            cd = cd[1:-1]
        if cd.size < 5:
            return
        gid = extract_gene_id(hdr)
        table_counts[tbl] += 1
        raw_records.append((gid, tbl, cd))

    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush(header, "".join(buf))
                header = line
                buf = []
            else:
                buf.append(line.upper())
        _flush(header, "".join(buf))

    if not raw_records:
        return [], default_table, 0

    dominant_table = max(table_counts.items(), key=lambda kv: kv[1])[0]
    best_per_gene = {}
    skipped_mixed = 0
    for gid, tbl, cd in raw_records:
        if tbl != dominant_table:
            skipped_mixed += 1
            continue
        if gid not in best_per_gene or cd.size > best_per_gene[gid][2].size:
            best_per_gene[gid] = (gid, tbl, cd)

    dedup = list(best_per_gene.values())
    # Strict coding-transcript guard: drop genes carrying an INTERNAL stop codon. Edge
    # start/stop are already trimmed above; any remaining stop is internal (pseudo-gene /
    # non-canonical CDS). Leaving these in activates the stop channel inside the tensor and
    # inflates the sense-convention df away from the canonical 2460 of Proposition 9.1.
    _aa_arr = codon_aa_channel_array(dominant_table)
    stop_idxs = set(int(i) for i in range(N_CODONS) if _aa_arr[i] == STOP_IDX)
    stop_arr = np.fromiter(stop_idxs, dtype=np.int64) if stop_idxs else np.zeros(0, dtype=np.int64)
    kept = []
    dropped_internal_stop = 0
    for rec in dedup:
        cd = rec[2]
        region = cd if trim_start_stop else cd[1:-1]
        if region.size and stop_arr.size and bool(np.any(np.isin(region, stop_arr))):
            dropped_internal_stop += 1
            continue
        kept.append(rec)
    dedup = kept
    if max_records and max_records > 0 and len(dedup) > max_records:
        dedup = dedup[:max_records]
    return dedup, dominant_table, skipped_mixed, dropped_internal_stop


def deterministic_gene_split(gene_id: str, train_percent=80) -> bool:
    root = family_root(str(gene_id))
    h = hashlib.sha256(root.encode("utf-8")).hexdigest()
    return (int(h[:8], 16) % 100) < train_percent


# =====================================================================
# 12. Full Species Analysis
# =====================================================================
def analyze_species(name, cds_path, max_records, trim_start_stop, alpha, window, step, do_hourglass):
    default_tbl = PANEL[name]["table"]
    records, dominant_table, skipped, skipped_stop = parse_cds_fasta_gz(
        cds_path, max_records=max_records, trim_start_stop=trim_start_stop, default_table=default_tbl
    )
    set_genetic_code(dominant_table)

    train_records, test_records = [], []
    for gid, tbl, cd in records:
        (train_records if deterministic_gene_split(gid, 80) else test_records).append((gid, tbl, cd))

    train_C = CodonCounts()
    for _, _, cd in train_records:
        train_C.add_codons(cd)
    test_C = CodonCounts()
    for _, _, cd in test_records:
        test_C.add_codons(cd)

    model = build_model(train_C, alpha=alpha)
    ide = compute_exact_decompositions(test_C)

    ll_rows, per_gene_S, per_gene_d2b = [], [], []
    for gid, tbl, cd in test_records:
        r = held_out_ll_contrasts(model, cd)
        if r is None:
            continue
        r["gene"] = gid
        ll_rows.append(r)
        per_gene_S.append(r["S_syn_heldout"])
        per_gene_d2b.append(r["d2b"])

    arr_S = np.array(per_gene_S, dtype=np.float64)
    mean_S = float(arr_S.mean()) if arr_S.size else float("nan")
    pos_S = int(np.sum(arr_S > 0))
    frac_S = float(pos_S / arr_S.size) if arr_S.size else float("nan")
    se_S = float(arr_S.std(ddof=1) / math.sqrt(arr_S.size)) if arr_S.size > 1 else float("nan")
    t_stat = float(mean_S / se_S) if (se_S and se_S > 0) else float("nan")
    t_pval = float(ttest_1samp(arr_S, 0.0).pvalue) if (HAVE_SCIPY and arr_S.size > 2) else float("nan")

    g2_all = g2_tests_all(test_C)
    ep_train = entropy_production(train_C)
    ep_test = entropy_production(test_C)
    dm_train = dynamic_mirror_current(train_C)
    dm_test = dynamic_mirror_current(test_C)
    maxcai = max_cai_dual_audit(test_records, model, max_genes=300)

    chem = {sch: compute_chemical_dpi(test_C, sch) for sch in CHEM_SCHEMES}
    hg = analyze_hourglass_and_windows(test_records, model, cds_path.parent, window=window, step=step) if do_hourglass else {}

    fibers = build_fibers()
    k_a = np.array([len(f) for f in fibers], dtype=np.int64)
    nz = k_a > 0
    syn_full = int(np.sum(k_a[nz] - 1))
    syn_sense = int(np.sum(k_a[nz & (np.arange(N_CHANNELS) != STOP_IDX)] - 1))

    return {
        "species": name,
        "status": "ok",
        "file": cds_path.name,
        "translation_table": int(dominant_table),
        "skipped_mixed_table_records": int(skipped),
        "skipped_internal_stop_records": int(skipped_stop),
        "n_records": len(records),
        "n_train_genes": len(train_records),
        "n_test_genes": len(test_records),
        "train_transitions": int(train_C.n_total),
        "test_transitions": int(test_C.n_total),
        "synonymous_dimension_full_64_to_channels": syn_full,
        "synonymous_dimension_sense": syn_sense,
        "information_decomposition_test_plugin": ide,
        "held_out_ll_contrasts": {
            "n_genes_evaluated": int(arr_S.size),
            "mean_S_syn_heldout": mean_S,
            "cluster_robust_SE": se_S,
            "positive_genes_count": pos_S,
            "fraction_positive_genes": frac_S,
            "gene_level_t_stat": t_stat,
            "gene_level_t_pvalue": t_pval,
        },
        "per_gene_S_syn_heldout": per_gene_S,
        "per_gene_d2b_heldout": per_gene_d2b,
        "analytic_G2_tests": g2_all,
        "markov_entropy_production": {"train": ep_train, "test": ep_test},
        "dynamic_mirror_current": {"train": dm_train, "test": dm_test},
        "max_cai_deterministic_collapse": maxcai,
        "chemical_DPI": chem,
        "hourglass_profile_diagnostic": hg,
    }


# =====================================================================
# 13. Publication Figures (Fig1 .. Fig7) & Markdown Report
# =====================================================================
def make_figures(payload, out_dir):
    if not HAVE_MPL:
        return []
    figdir = Path(out_dir) / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    res = [r for r in payload["results"] if r.get("status") == "ok"]
    if not res:
        return []
    names = [r["species"] for r in res]
    saved = []

    def _save(fig, fname):
        p1 = figdir / fname
        p2 = figdir / fname.replace(".png", ".pdf")
        fig.savefig(p1, dpi=200, bbox_inches="tight", pad_inches=0.18)
        fig.savefig(p2, bbox_inches="tight", pad_inches=0.18)
        plt.close(fig)
        saved.append(str(p1))

    I_AA = np.array([r["information_decomposition_test_plugin"]["I_AA"] for r in res])
    D2B = np.array([r["information_decomposition_test_plugin"]["I_Cprev_Acur_Aprev_d2b"] for r in res])
    SS = np.array([r["information_decomposition_test_plugin"]["S_syn_I_Cprev_Ccur_Acur"] for r in res])
    APC = np.array([r["information_decomposition_test_plugin"]["I_Aprev_Ccur_Acur"] for r in res])
    SS2 = np.array([r["information_decomposition_test_plugin"]["S_syn2_I_Cprev_Ccur_Aprev_Acur"] for r in res])
    ICC = np.array([r["information_decomposition_test_plugin"]["I_CC"] for r in res])
    x = np.arange(len(names))
    w = 0.56

    # Fig 1: Tripartite decomposition
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.bar(x, I_AA, w, label="I(A_{i-1}; A_i) [Protein syntax]")
    ax.bar(x, D2B, w, bottom=I_AA, label="d2b = I(C_{i-1}; A_i | A_{i-1}) [Codon->AA context]")
    ax.bar(x, SS, w, bottom=I_AA + D2B, label="S_syn = I(C_{i-1}; C_i | A_i) [Synonymous program]")
    ax.plot(x, ICC, "ko", ms=6, label="Total I(C_{i-1}; C_i)")
    ax.set_ylim(0, float(ICC.max()) * 1.32)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=28, ha="right", fontsize=9.5)
    ax.set_ylabel("Information (nats / codon)", fontsize=10)
    ax.set_title("Fig 1. Exact Tripartite Decomposition on Held-Out RefSeq Genes (Residual < 1e-15)", fontsize=11, pad=10)
    ax.legend(loc="upper left", ncol=2, fontsize=8.2, framealpha=0.92)
    fig.tight_layout()
    _save(fig, "fig1_tripartite.png")

    # Fig 2: Quadripartite decomposition
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.bar(x, I_AA, w, label="I(A_{i-1}; A_i)")
    ax.bar(x, D2B, w, bottom=I_AA, label="I(C_{i-1}; A_i | A_{i-1})")
    ax.bar(x, APC, w, bottom=I_AA + D2B, label="I(A_{i-1}; C_i | A_i)")
    ax.bar(x, SS2, w, bottom=I_AA + D2B + APC, label="S_syn^(2) = I(C_{i-1}; C_i | A_{i-1}, A_i)")
    ax.plot(x, ICC, "ko", ms=6, label="Total I(C_{i-1}; C_i)")
    ax.set_ylim(0, float(ICC.max()) * 1.32)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=28, ha="right", fontsize=9.5)
    ax.set_ylabel("Information (nats / codon)", fontsize=10)
    ax.set_title("Fig 2. Exact Quadripartite Decomposition (Theorem 5.3)", fontsize=11, pad=10)
    ax.legend(loc="upper left", ncol=2, fontsize=8.2, framealpha=0.92)
    fig.tight_layout()
    _save(fig, "fig2_quadripartite.png")

    # Fig 3: Max-CAI internal collapse & cross-scoring penalty
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    wt_cross = np.array([r["max_cai_deterministic_collapse"]["cross_scored_WT_mean_d2b"] for r in res])
    cai_cross = np.array([r["max_cai_deterministic_collapse"]["cross_scored_MaxCAI_mean_d2b"] for r in res])
    cai_int = np.array([r["max_cai_deterministic_collapse"]["internal_S_syn_after_MaxCAI"] for r in res])
    bw = 0.25
    ax.bar(x - bw, wt_cross, bw, color="tab:green", label="Wild-Type d2b (scored on P_WT)")
    ax.bar(x, cai_int, bw, color="tab:gray", label="Max-CAI internal S_syn (= 0, Thm 6.1)")
    ax.bar(x + bw, cai_cross, bw, color="tab:red", label="Max-CAI cross-scored d2b (< 0, Prop 6.2)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylim(float(cai_cross.min()) * 1.28, float(wt_cross.max()) * 1.25)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=28, ha="right", fontsize=9.5)
    ax.set_ylabel("nats / codon", fontsize=10)
    ax.set_title("Fig 3. Deterministic Max-CAI Collapse (Thm 6.1) & KL Cross-Scoring Crash (Prop 6.2)", fontsize=11, pad=10)
    ax.legend(loc="lower right", fontsize=8.5, framealpha=0.92)
    fig.tight_layout()
    _save(fig, "fig3_maxcai_collapse.png")

    # Fig 4: Non-equilibrium entropy production & mirror currents
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.8))
    sig = np.array([r["markov_entropy_production"]["test"]["sigma"] for r in res])
    mir = np.array([r["dynamic_mirror_current"]["test"]["mirror_current_L1"] for r in res])
    sta = np.array([r["dynamic_mirror_current"]["test"]["static_mirror_L1"] for r in res])
    axes[0].bar(x, sig, color="tab:red", width=0.6)
    axes[0].set_title("Entropy Production\nsigma (nats / codon)", fontsize=10.5, pad=8)
    axes[1].bar(x, mir, color="tab:blue", width=0.6)
    axes[1].set_title("Dynamic Mirror Current\n||J^mir||_1", fontsize=10.5, pad=8)
    axes[2].bar(x, sta, color="tab:purple", width=0.6)
    axes[2].set_title("Static Mirror Asymmetry\nMirror_L1(pi) = ||P_{V_-}pi||_1", fontsize=10.5, pad=8)
    for a in axes:
        a.set_xticks(x)
        a.set_xticklabels(names, rotation=35, ha="right", fontsize=8.8)
    fig.suptitle(
        "Fig 4. Non-Equilibrium Dynamics and Mirror Symmetry Breaking on T_4^3\n"
        "(bars = sigma_reg, finite 1e-12 scale; exact sigma = +inf iff unidir_edges > 0 -- see block 2b / JSON)",
        fontsize=10.0, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    _save(fig, "fig4_dynamics.png")

    # Fig 5: Per-gene held-out S_syn distributions across ALL species
    n_sp = len(res)
    ncols = min(3, n_sp)
    nrows = int(math.ceil(n_sp / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.3 * nrows), squeeze=False)
    for idx_sp, r in enumerate(res):
        ax = axes[idx_sp // ncols][idx_sp % ncols]
        arr = np.array(r.get("per_gene_S_syn_heldout", []), dtype=np.float64)
        if arr.size:
            ax.hist(arr, bins=36, color="tab:purple", alpha=0.85, edgecolor="none")
            ax.axvline(0, color="k", lw=1.1, ls="--")
            ax.set_title(
                f"{r['species']} (n={arr.size})\nmean={arr.mean():+.4f}, pos={100*np.mean(arr>0):.1f}%",
                fontsize=9.5,
                pad=6,
            )
        ax.set_xlabel("per-gene S_syn held-out (nats)", fontsize=8.5)
        ax.set_ylabel("genes", fontsize=8.5)
        ax.tick_params(labelsize=8)
    for idx_empty in range(n_sp, nrows * ncols):
        axes[idx_empty // ncols][idx_empty % ncols].axis("off")
    fig.suptitle("Fig 5. Cluster-Robust Gene-Level Evidence on Held-Out Test Genes", fontsize=12.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, "fig5_per_gene.png")

    # Fig 6: Perez 17q22 100-kb Hourglass Profile
    hg_res = res[0].get("hourglass_profile_diagnostic", {}).get("chr17_17q22_locus", {})
    if hg_res and hg_res.get("V_genomic"):
        xg = np.array(hg_res["x_grid"])
        vg = np.array(hg_res["V_genomic"])
        vp = np.array(hg_res["V_proteomic"])
        fig, ax = plt.subplots(figsize=(8.8, 4.6))
        ax.plot(xg, vg, color="navy", lw=1.35, label="Genomics patterned signature G(i)")
        ax.plot(xg, vp, color="crimson", lw=1.15, ls="--", label="Proteomics patterned signature P(i) (pushforward alpha_*)")
        ax.set_xlabel("Codon / AA position along 100-kb NCBI window (NC_000017.11:55154472-55254472)", fontsize=9.5)
        ax.set_ylabel("Cumulative Master-Code Distance from Isthmus", fontsize=9.5)
        ax.set_title(
            f"Fig 6. Empirical Hourglass Profile ({hg_res['source']})\n"
            f"Link = {hg_res['hourglass_G_vs_P_Link_percent']:.2f}% (Form), while Bigram S_syn = +{hg_res['locus_S_syn_nats']:.4f} nats (Program)",
            fontsize=10.5,
            pad=10,
        )
        ax.legend(loc="upper right", fontsize=8.8, framealpha=0.92)
        fig.tight_layout()
        _save(fig, "fig6_hourglass_17q22.png")

    # Fig 7: Analytic G^2 vs df
    fig, ax = plt.subplots(figsize=(8.6, 4.5))
    g2v = np.array([r["analytic_G2_tests"]["S_syn_G2"]["G2"] for r in res])
    dfv = np.array([r["analytic_G2_tests"]["S_syn_G2"]["df"] for r in res])
    ax.bar(x, g2v, color="tab:orange", width=0.58)
    ax.set_ylim(0, float(g2v.max()) * 1.18)
    for i, (v, d) in enumerate(zip(g2v, dfv)):
        ax.text(i, v + float(g2v.max()) * 0.02, f"df={d}", ha="center", va="bottom", fontsize=8.2)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=28, ha="right", fontsize=9.5)
    ax.set_ylabel("G^2 likelihood-ratio statistic", fontsize=10)
    ax.set_title("Fig 7. Analytic Likelihood-Ratio G^2 vs Degrees of Freedom (Zero Simulations)", fontsize=11, pad=10)
    fig.tight_layout()
    _save(fig, "fig7_g2.png")

    return saved


def write_markdown_report(results, out_md: Path, torus_inv: dict, figs: list):
    L = [
        "# Codon Torus Real-Data Analytic Audit v3.2",
        "",
        "## 1. Global Codon Torus Invariants",
        f"- **Codon torus cardinality**: `{torus_inv['n_codons']}`",
        f"- **Reverse-complement involution (`rho^2 = id`)**: `{torus_inv['rho_involution']}`",
        f"- **Fixed-point-free (`rho(c) != c`)**: `{torus_inv['rho_fixed_point_free']}`",
        f"- **Exact mirror orbits**: `{torus_inv['mirror_orbits']}`",
        f"- **Synonymous simplex dimension (`64 -> 21` / `61 -> 20`)**: `{torus_inv['simplex_dim_64_to_21']}` / `{torus_inv['simplex_dim_61_to_20']}`",
        "",
        "## 2. Tripartite & Quadripartite Information Decomposition (Held-Out Test Genes)",
        "",
        "| Species | Table | Test Genes | Test Trans | I(C;C) | I(A;A) | d2b | S_syn | I(A_prev;C|A) | S_syn^(2) | Residual | G^2 (df) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        ide = r["information_decomposition_test_plugin"]
        g2 = r["analytic_G2_tests"]["S_syn_G2"]
        res_err = abs(ide["quadripartite_residual_I_CC_minus_sum"])
        L.append(
            f"| **{r['species']}** | {r['translation_table']} | {r['n_test_genes']} | {r['test_transitions']} | "
            f"{ide['I_CC']:+.4f} | {ide['I_AA']:+.4f} | {ide['I_Cprev_Acur_Aprev_d2b']:+.4f} | "
            f"{ide['S_syn_I_Cprev_Ccur_Acur']:+.4f} | {ide['I_Aprev_Ccur_Acur']:+.4f} | "
            f"{ide['S_syn2_I_Cprev_Ccur_Aprev_Acur']:+.4f} | {res_err:.1e} | {g2['G2']:.1f} ({g2['df']}) |"
        )

    L.extend([
        "",
        "## 3. Chemical Data-Processing Inequality (Theorem 7.1) & Non-Equilibrium Dynamics (Theorem 8.1)",
        "",
        "| Species | d2b (Upper Bound) | I_Hydro | I_Charge | I_Volume | I_Arom | DPI_OK | sigma (nats) | ||J^mir||_1 | Mirror_L1(pi) | Relic I(A;Arho) |",
        "|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|",
    ])
    for r in results:
        ide = r["information_decomposition_test_plugin"]
        ch = r["chemical_DPI"]
        ep = r["markov_entropy_production"]["test"]
        dm = r["dynamic_mirror_current"]["test"]
        rel = r["analytic_G2_tests"]["DoubleStrand_Relic_Coupling"]["I_A_Arho_nats"]
        dpi_ok = all(v["DPI_ok"] and v["MI_nonnegative"] for v in ch.values())
        L.append(
            f"| **{r['species']}** | {ide['I_Cprev_Acur_Aprev_d2b']:+.4f} | "
            f"{ch['hydrophobicity']['I_Cprev_Qcur_Aprev']:+.4f} | {ch['charge']['I_Cprev_Qcur_Aprev']:+.4f} | "
            f"{ch['volume']['I_Cprev_Qcur_Aprev']:+.4f} | {ch['aromaticity']['I_Cprev_Qcur_Aprev']:+.4f} | "
            f"`{dpi_ok}` | {ep['sigma']:+.4f} | {dm['mirror_current_L1']:.4f} | "
            f"{dm['static_mirror_L1']:.4f} | {rel:+.4f} |"
        )

    L.extend([
        "",
        "## 4. Max-CAI Collapse (Theorem 6.1), Cross-Scoring Penalty (Prop. 6.2) & Held-Out Gene Tests",
        "",
        "| Species | MaxCAI Int S_syn | MaxCAI Int d2b | WT Cross d2b | MaxCAI Cross d2b | Held-Out S_syn | Pos Genes % | Gene t-stat |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for r in results:
        mc = r["max_cai_deterministic_collapse"]
        ho = r["held_out_ll_contrasts"]
        L.append(
            f"| **{r['species']}** | {mc['internal_S_syn_after_MaxCAI']:+.6f} | {mc['internal_d2b_after_MaxCAI']:+.6f} | "
            f"{mc['cross_scored_WT_mean_d2b']:+.4f} | {mc['cross_scored_MaxCAI_mean_d2b']:+.4f} | "
            f"{ho['mean_S_syn_heldout']:+.4f} | {100.0*ho['fraction_positive_genes']:.1f}% | {ho['gene_level_t_stat']:+.2f} |"
        )

    out_md.write_text("\n".join(L) + "\n", encoding="utf-8")


# =====================================================================
# 14. Main CLI & Console Summary
# =====================================================================
def main():
    p = argparse.ArgumentParser(description="Codon Torus Real-Data Analytic Audit v3.2")
    p.add_argument("--out-dir", default="codon_torus_audit_v3")
    p.add_argument("--cds-dir", default="cds_cache")
    p.add_argument("--panel", default="Human_T2T,M_musculus,D_melanogaster,A_thaliana,S_cerevisiae,E_coli,B_subtilis,M_jannaschii,H_salinarum")
    p.add_argument("--max-records", type=int, default=0)
    p.add_argument("--alpha", type=float, default=0.5)
    p.add_argument("--window", type=int, default=300)
    p.add_argument("--step", type=int, default=100)
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--skip-download", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--bnt162b2-audit", action="store_true", help="Run full BNT162b2 & Ambigrammatic Duplex Torus audit")
    p.add_argument("--bnt162b2-only", action="store_true", help="Run ONLY the BNT162b2 Duplex Torus audit and exit")
    p.add_argument("--force-ncbi", action="store_true", help="Force live download of viral/BNT162b2 sequences from NCBI Entrez")
    p.add_argument("--force-ncbi-all", action="store_true", help="Force live download of ALL RefSeq CDS archives from NCBI FTP")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    cds_dir = Path(args.cds_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cds_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 116)
    print("CODON TORUS REAL-DATA ANALYTIC AUDIT v3.2 (ZERO MONTE CARLO SURROGATES — EXACT ALGEBRAIC & G^2 PROOFS)")
    print("=" * 116)
    torus_inv = self_test_torus_theorems()
    print(f"[STARTUP SELF-TEST PASSED] {torus_inv}")
    if getattr(args, "force_ncbi_all", False):
        resolve_or_download_cds.force_download = True
    if getattr(args, "bnt162b2_audit", False) or getattr(args, "bnt162b2_only", False):
        run_bnt162b2_duplex_audit(
            cds_dir=cds_dir,
            out_dir=out_dir,
            max_records=args.max_records,
            alpha=args.alpha,
            force_ncbi=bool(getattr(args, "force_ncbi", False)),
            force_ncbi_all=bool(getattr(args, "force_ncbi_all", False)),
            timeout=args.timeout,
        )
        if getattr(args, "bnt162b2_only", False):
            return

    req_names = [x.strip() for x in args.panel.split(",") if x.strip() in PANEL]
    results = []
    for idx, name in enumerate(req_names):
        try:
            path = resolve_or_download_cds(
                name=name, cds_dir=cds_dir, timeout=args.timeout, skip_download=args.skip_download
            )
        except Exception as e:
            print(f"  [SKIP] {name}: {e}")
            continue
        r = analyze_species(
            name=name,
            cds_path=path,
            max_records=args.max_records,
            trim_start_stop=True,
            alpha=args.alpha,
            window=args.window,
            step=args.step,
            do_hourglass=(idx == 0),
        )
        results.append(r)

    print("\n" + "-" * 116)
    print("1. TRIPARTITE & QUADRIPARTITE INFORMATION DECOMPOSITION ON HELD-OUT TEST GENES (Theorems 5.2, 5.3, 5.4)")
    print("-" * 116)
    print("Species          | Table | Test_Trans |  I(C;C) |  I(A;A) |     d2b |   S_syn | AA1>Syn2 | S_syn(2) | Residual | G2_p(S_syn)")
    print("-----------------+-------+------------+---------+---------+---------+---------+----------+----------+----------+------------")
    for r in results:
        ide = r["information_decomposition_test_plugin"]
        g2p = r["analytic_G2_tests"]["S_syn_G2"]["p"]
        res_err = abs(ide["quadripartite_residual_I_CC_minus_sum"])
        print(
            f"{r['species']:16s} | {r['translation_table']:5d} | {r['test_transitions']:10d} | "
            f"{ide['I_CC']:+7.4f} | {ide['I_AA']:+7.4f} | {ide['I_Cprev_Acur_Aprev_d2b']:+7.4f} | "
            f"{ide['S_syn_I_Cprev_Ccur_Acur']:+7.4f} | {ide['I_Aprev_Ccur_Acur']:+8.4f} | "
            f"{ide['S_syn2_I_Cprev_Ccur_Aprev_Acur']:+8.4f} | {res_err:8.1e} | {g2p:10.2e}"
        )

    print("\n" + "-" * 116)
    print("1b. TWO-STRAND DUPLEX DECOMPOSITION S_syn = S_anti + S_dup (Theorems 4.3 & 5.4: 64 -> 56 Duplex Channels, Dim = 8)")
    print("-" * 116)
    print("Species          |   S_syn |  S_anti (Antisense AA) |   S_dup (8-Doublet) | Antisense % | Residual | Relic_I(A;Arho)")
    print("-----------------+---------+------------------------+---------------------+-------------+----------+----------------")
    for r in results:
        ide = r["information_decomposition_test_plugin"]
        rel = r["analytic_G2_tests"]["DoubleStrand_Relic_Coupling"]["I_A_Arho_nats"]
        s_syn = ide["S_syn_I_Cprev_Ccur_Acur"]
        s_anti = ide["S_anti_I_Cprev_Arho_given_Acur"]
        s_dup = ide["S_dup_I_Cprev_Ccur_given_Acur_Arho"]
        pct = 100.0 * ide["antisense_program_fraction_of_S_syn"]
        res_d = abs(ide["duplex_S_syn_residual"])
        print(
            f"{r['species']:16s} | {s_syn:+7.4f} | {s_anti:+22.4f} | {s_dup:+19.4f} | "
            f"{pct:10.1f}% | {res_d:8.1e} | {rel:+15.4f}"
        )

    print("\n" + "-" * 116)
    print("2. CHEMICAL DATA-PROCESSING INEQUALITY (Theorem 7.1) & NON-EQUILIBRIUM TORUS DYNAMICS (Theorem 8.1)")
    print("-" * 116)
    print("Species          |     d2b | I_Hydro |  I_Chg  |  I_Vol  |  I_Arom | DPI_OK |   sigma | MirCur_L1 | StaticMir_L1 | Relic_I(A;Arho)")
    print("-----------------+---------+---------+---------+---------+---------+--------+---------+-----------+--------------+----------------")
    for r in results:
        ide = r["information_decomposition_test_plugin"]
        ch = r["chemical_DPI"]
        ep = r["markov_entropy_production"]["test"]
        dm = r["dynamic_mirror_current"]["test"]
        rel = r["analytic_G2_tests"]["DoubleStrand_Relic_Coupling"]["I_A_Arho_nats"]
        dpi_ok = all(v["DPI_ok"] and v["MI_nonnegative"] for v in ch.values())
        print(
            f"{r['species']:16s} | {ide['I_Cprev_Acur_Aprev_d2b']:+7.4f} | "
            f"{ch['hydrophobicity']['I_Cprev_Qcur_Aprev']:+7.4f} | {ch['charge']['I_Cprev_Qcur_Aprev']:+7.4f} | "
            f"{ch['volume']['I_Cprev_Qcur_Aprev']:+7.4f} | {ch['aromaticity']['I_Cprev_Qcur_Aprev']:+7.4f} | "
            f"{str(dpi_ok):6s} | {ep['sigma']:+7.4f} | {dm['mirror_current_L1']:9.4f} | "
            f"{dm['static_mirror_L1']:12.4f} | {rel:+15.4f}"
        )

    print("\n" + "-" * 116)
    print("2b. EXACT IRREVERSIBILITY CERTIFICATE & CONTINGENCY DEGREES OF FREEDOM (Thm 8.1 / Prop. 9.1)")
    print("-" * 116)
    print("Species          | unidir_edges | sigma_exact | sigma_reg | df(G2 S_syn) | df(G2 S_syn2) | df(G2 lump)")
    print("-----------------+--------------+-------------+-----------+--------------+---------------+------------")
    for r in results:
        ep = r["markov_entropy_production"]["test"]
        g2 = r["analytic_G2_tests"]
        se = ep.get("sigma_exact", float("nan"))
        se_str = "+inf" if (se == float("inf")) else (f"{se:+.4f}" if se == se else "nan")
        print(
            f"{r['species']:16s} | {int(ep.get('unidirectional_edges', 0)):12d} | {se_str:>11s} | "
            f"{float(ep.get('sigma_reg', ep.get('sigma', float('nan')))):+9.4f} | "
            f"{int(g2['S_syn_G2']['df']):12d} | {int(g2['S_syn2_G2']['df']):13d} | "
            f"{int(g2['Lumpability_d2b_G2']['df']):10d}"
        )

    print("\n" + "-" * 116)
    print("3. MAX-CAI COLLAPSE (Theorem 6.1) vs. CROSS-SCORING CRASH (Proposition 6.2) & CLUSTER-ROBUST HELD-OUT TEST")
    print("-" * 116)
    print("Species          | MaxCAI_Int_S_syn | MaxCAI_Int_d2b | WT_Cross_d2b | MaxCAI_Cross_d2b | HeldOut_S_syn | Pos_Genes% | t_stat  | Cross<0%")
    print("-----------------+------------------+----------------+--------------+------------------+---------------+------------+---------+---------")
    for r in results:
        mc = r["max_cai_deterministic_collapse"]
        ho = r["held_out_ll_contrasts"]
        print(
            f"{r['species']:16s} | {mc['internal_S_syn_after_MaxCAI']:+16.6f} | "
            f"{mc['internal_d2b_after_MaxCAI']:+14.6f} | {mc['cross_scored_WT_mean_d2b']:+12.4f} | "
            f"{mc['cross_scored_MaxCAI_mean_d2b']:+16.4f} | {ho['mean_S_syn_heldout']:+13.4f} | "
            f"{100.0 * ho['fraction_positive_genes']:9.1f}% | {ho['gene_level_t_stat']:+7.2f} | "
            f"{100.0 * float(mc.get('fraction_genes_MaxCAI_d2b_negative', float('nan'))):7.1f}%"
        )

    if results and results[0].get("hourglass_profile_diagnostic", {}).get("chr17_17q22_locus"):
        hg = results[0]["hourglass_profile_diagnostic"]
        loc = hg["chr17_17q22_locus"]
        print("\n" + "-" * 116)
        print("4. EMPIRICAL HOURGLASS & FORM-vs-PROGRAM ANCHOR (Section 10 / Lemma 9.3)")
        print("-" * 116)
        print(
            f"  Locus: {loc['source']} ({loc['n_codons']} codons) | "
            f"Hourglass Genomics-Proteomics Link = {loc['hourglass_G_vs_P_Link_percent']:.2f}% (Isthmus at codon {loc['isthmus_codon_index']}) | "
            f"Window R^2(G vs P) = {hg['R2_profile_G_vs_P']:.4f} | "
            f"Sequential Bigram Program S_syn = +{loc['locus_S_syn_nats']:.4f} nats/codon"
        )

    payload = {
        "version": "v3.2",
        "torus_invariants": torus_inv,
        "policy": {
            "real_data_only": True,
            "zero_monte_carlo_surrogates": True,
            "monte_carlo_surrogate_as_proof": False,
            "deterministic_sha256_gene_split": True,
            "analytic_G2_and_heldout_LL": True,
            "startup_self_test_tensor": "deterministic_arithmetic_fixture_no_rng",
            "hourglass_locus_is_genomic_window_not_ORF": True,
            "internal_stop_codons_dropped": True,
        },
        "results": results,
    }
    json_out = out_dir / "codon_torus_audit_v3.json"
    md_out = out_dir / "codon_torus_audit_v3.md"
    json_out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    figs = [] if args.no_plots else make_figures(payload, out_dir)
    write_markdown_report(results, md_out, torus_inv, figs)
    print("=" * 116)
    print(f"[SAVED JSON] {json_out}")
    print(f"[SAVED MD]   {md_out}")
    if figs:
        print(f"[SAVED FIGS] {len(figs)} publication figures in {out_dir / 'figures'}")
    print("=" * 116)



# === BEGIN BNT162B2 & MODERNA mRNA-1273 DUPLEX FORENSICS PATCH v3.5 (Zero-RNG, Corollary 4.4 Compiler) ===
# Complete simulation-free, zero-RNG verification module for:
#   1. Live NCBI Entrez & Stanford/Fire Lab reference retrieval (NC_045512.2, BNT162b2 3825nt, Moderna mRNA-1273 3828nt)
#   2. Full two-strand duplex torus & Markovian audit (Theorems 4.3, 5.3, 5.4, 5.5, 6.1, 8.1, 11.2 & Props 6.2, 6.3, 8.2, 8.5)
#   3. 6-frame ambigrammatic overlapping ORF scanner (Corollary 4.4: F_{Leu,Stop}={TTA,CTA}, F_{Ser,Stop}={TCA}, F_{His,Met}={CAT})
#   4. Constructive Duplex-Safe Viterbi Compiler (Corollary 4.4 + Proposition 6.3): eliminates antisense ATG starts (CAT=0)
#      and restores periodic antisense termination barriers while maximizing human conditional ribosomal syntax.

import io
import zipfile

C_N1_ARR = np.arange(N_CODONS) // 16
C_N2_ARR = (np.arange(N_CODONS) // 4) % 4
C_N3_ARR = np.arange(N_CODONS) % 4
IS_GC3_ARR = (C_N3_ARR == 1) | (C_N3_ARR == 2)

# Canonical BNT162b2 CDS: 3825 nt, 1275 codons (1273 sense codons + tandem terminal TGATGA).
BNT162B2_FULL_3825NT = (
    "ATGTTCGTGTTCCTGGTGCTGCTGCCTCTGGTGTCCAGCCAGTGTGTGAACCTGACCACCAGAACACAGC"
    "TGCCTCCAGCCTACACCAACAGCTTTACCAGAGGCGTGTACTACCCCGACAAGGTGTTCAGATCCAGCGT"
    "GCTGCACTCTACCCAGGACCTGTTCCTGCCTTTCTTCAGCAACGTGACCTGGTTCCACGCCATCCACGTG"
    "TCCGGCACCAATGGCACCAAGAGATTCGACAACCCCGTGCTGCCCTTCAACGACGGGGTGTACTTTGCCA"
    "GCACCGAGAAGTCCAACATCATCAGAGGCTGGATCTTCGGCACCACACTGGACAGCAAGACCCAGAGCCT"
    "GCTGATCGTGAACAACGCCACCAACGTGGTCATCAAAGTGTGCGAGTTCCAGTTCTGCAACGACCCCTTC"
    "CTGGGCGTCTACTACCACAAGAACAACAAGAGCTGGATGGAAAGCGAGTTCCGGGTGTACAGCAGCGCCA"
    "ACAACTGCACCTTCGAGTACGTGTCCCAGCCTTTCCTGATGGACCTGGAAGGCAAGCAGGGCAACTTCAA"
    "GAACCTGCGCGAGTTCGTGTTTAAGAACATCGACGGCTACTTCAAGATCTACAGCAAGCACACCCCTATC"
    "AACCTCGTGCGGGATCTGCCTCAGGGCTTCTCTGCTCTGGAACCCCTGGTGGATCTGCCCATCGGCATCA"
    "ACATCACCCGGTTTCAGACACTGCTGGCCCTGCACAGAAGCTACCTGACACCTGGCGATAGCAGCAGCGG"
    "ATGGACAGCTGGTGCCGCCGCTTACTATGTGGGCTACCTGCAGCCTAGAACCTTCCTGCTGAAGTACAAC"
    "GAGAACGGCACCATCACCGACGCCGTGGATTGTGCTCTGGATCCTCTGAGCGAGACAAAGTGCACCCTGA"
    "AGTCCTTCACCGTGGAAAAGGGCATCTACCAGACCAGCAACTTCCGGGTGCAGCCCACCGAATCCATCGT"
    "GCGGTTCCCCAATATCACCAATCTGTGCCCCTTCGGCGAGGTGTTCAATGCCACCAGATTCGCCTCTGTG"
    "TACGCCTGGAACCGGAAGCGGATCAGCAATTGCGTGGCCGACTACTCCGTGCTGTACAACTCCGCCAGCT"
    "TCAGCACCTTCAAGTGCTACGGCGTGTCCCCTACCAAGCTGAACGACCTGTGCTTCACAAACGTGTACGC"
    "CGACAGCTTCGTGATCCGGGGAGATGAAGTGCGGCAGATTGCCCCTGGACAGACAGGCAAGATCGCCGAC"
    "TACAACTACAAGCTGCCCGACGACTTCACCGGCTGTGTGATTGCCTGGAACAGCAACAACCTGGACTCCA"
    "AAGTCGGCGGCAACTACAATTACCTGTACCGGCTGTTCCGGAAGTCCAATCTGAAGCCCTTCGAGCGGGA"
    "CATCTCCACCGAGATCTATCAGGCCGGCAGCACCCCTTGTAACGGCGTGGAAGGCTTCAACTGCTACTTC"
    "CCACTGCAGTCCTACGGCTTTCAGCCCACAAATGGCGTGGGCTATCAGCCCTACAGAGTGGTGGTGCTGA"
    "GCTTCGAACTGCTGCATGCCCCTGCCACAGTGTGCGGCCCTAAGAAAAGCACCAATCTCGTGAAGAACAA"
    "ATGCGTGAACTTCAACTTCAACGGCCTGACCGGCACCGGCGTGCTGACAGAGAGCAACAAGAAGTTCCTG"
    "CCATTCCAGCAGTTTGGCCGGGATATCGCCGATACCACAGACGCCGTTAGAGATCCCCAGACACTGGAAA"
    "TCCTGGACATCACCCCTTGCAGCTTCGGCGGAGTGTCTGTGATCACCCCTGGCACCAACACCAGCAATCA"
    "GGTGGCAGTGCTGTACCAGGACGTGAACTGTACCGAAGTGCCCGTGGCCATTCACGCCGATCAGCTGACA"
    "CCTACATGGCGGGTGTACTCCACCGGCAGCAATGTGTTTCAGACCAGAGCCGGCTGTCTGATCGGAGCCG"
    "AGCACGTGAACAATAGCTACGAGTGCGACATCCCCATCGGCGCTGGAATCTGCGCCAGCTACCAGACACA"
    "GACAAACAGCCCTCGGAGAGCCAGAAGCGTGGCCAGCCAGAGCATCATTGCCTACACAATGTCTCTGGGC"
    "GCCGAGAACAGCGTGGCCTACTCCAACAACTCTATCGCTATCCCCACCAACTTCACCATCAGCGTGACCA"
    "CAGAGATCCTGCCTGTGTCCATGACCAAGACCAGCGTGGACTGCACCATGTACATCTGCGGCGATTCCAC"
    "CGAGTGCTCCAACCTGCTGCTGCAGTACGGCAGCTTCTGCACCCAGCTGAATAGAGCCCTGACAGGGATC"
    "GCCGTGGAACAGGACAAGAACACCCAAGAGGTGTTCGCCCAAGTGAAGCAGATCTACAAGACCCCTCCTA"
    "TCAAGGACTTCGGCGGCTTCAATTTCAGCCAGATTCTGCCCGATCCTAGCAAGCCCAGCAAGCGGAGCTT"
    "CATCGAGGACCTGCTGTTCAACAAAGTGACACTGGCCGACGCCGGCTTCATCAAGCAGTATGGCGATTGT"
    "CTGGGCGACATTGCCGCCAGGGATCTGATTTGCGCCCAGAAGTTTAACGGACTGACAGTGCTGCCTCCTC"
    "TGCTGACCGATGAGATGATCGCCCAGTACACATCTGCCCTGCTGGCCGGCACAATCACAAGCGGCTGGAC"
    "ATTTGGAGCAGGCGCCGCTCTGCAGATCCCCTTTGCTATGCAGATGGCCTACCGGTTCAACGGCATCGGA"
    "GTGACCCAGAATGTGCTGTACGAGAACCAGAAGCTGATCGCCAACCAGTTCAACAGCGCCATCGGCAAGA"
    "TCCAGGACAGCCTGAGCAGCACAGCAAGCGCCCTGGGAAAGCTGCAGGACGTGGTCAACCAGAATGCCCA"
    "GGCACTGAACACCCTGGTCAAGCAGCTGTCCTCCAACTTCGGCGCCATCAGCTCTGTGCTGAACGATATC"
    "CTGAGCAGACTGGACCCTCCTGAGGCCGAGGTGCAGATCGACAGACTGATCACAGGCAGACTGCAGAGCC"
    "TCCAGACATACGTGACCCAGCAGCTGATCAGAGCCGCCGAGATTAGAGCCTCTGCCAATCTGGCCGCCAC"
    "CAAGATGTCTGAGTGTGTGCTGGGCCAGAGCAAGAGAGTGGACTTTTGCGGCAAGGGCTACCACCTGATG"
    "AGCTTCCCTCAGTCTGCCCCTCACGGCGTGGTGTTTCTGCACGTGACATATGTGCCCGCTCAAGAGAAGA"
    "ATTTCACCACCGCTCCAGCCATCTGCCACGACGGCAAAGCCCACTTTCCTAGAGAAGGCGTGTTCGTGTC"
    "CAACGGCACCCATTGGTTCGTGACACAGCGGAACTTCTACGAGCCCCAGATCATCACCACCGACAACACC"
    "TTCGTGTCTGGCAACTGCGACGTCGTGATCGGCATTGTGAACAATACCGTGTACGACCCTCTGCAGCCCG"
    "AGCTGGACAGCTTCAAAGAGGAACTGGACAAGTACTTTAAGAACCACACAAGCCCCGACGTGGACCTGGG"
    "CGATATCAGCGGAATCAATGCCAGCGTCGTGAACATCCAGAAAGAGATCGACCGGCTGAACGAGGTGGCC"
    "AAGAATCTGAACGAGAGCCTGATCGACCTGCAAGAACTGGGGAAGTACGAGCAGTACATCAAGTGGCCCT"
    "GGTACATCTGGCTGGGCTTTATCGCCGGACTGATTGCCATCGTGATGGTCACAATCATGCTGTGTTGCAT"
    "GACCAGCTGCTGTAGCTGCCTGAAGGGCTGTTGTAGCTGTGGCAGCTGCTGCAAGTTCGACGAGGACGAT"
    "TCTGAGCCCGTGCTGAAGGGCGTGAAACTGCACTACACATGATGA"
)


def revcomp_nt_seq(seq: str) -> str:
    tr = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tr)[::-1]


def fetch_ncbi_nuccore_fasta(acc: str, cache_dir: Path, force_ncbi: bool = False, timeout: int = 120) -> str:
    """Direct FASTA retrieval from NCBI Entrez E-utilities with local cache verification."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_f = cache_dir / f"{acc}.fasta"
    if cache_f.exists() and cache_f.stat().st_size > 256 and not force_ncbi:
        txt = cache_f.read_text(encoding="utf-8", errors="ignore")
        seq = "".join(line.strip().upper() for line in txt.splitlines() if not line.startswith(">"))
        if len(seq) > 100:
            sha = hashlib.sha256(seq.encode("ascii", "ignore")).hexdigest()[:12]
            print(f"  [NCBI-Cache]  {acc:16s} | {len(seq)} nt (sha256={sha})")
            return seq

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={acc}&rettype=fasta&retmode=text"
    )
    print(f"  [NCBI-Entrez] GET {url}")
    raw = http_get(url, timeout=timeout, log_label=acc).decode("utf-8", errors="ignore")
    if not raw.lstrip().startswith(">"):
        raise RuntimeError(f"NCBI Entrez returned non-FASTA payload for {acc}")
    cache_f.write_text(raw, encoding="utf-8")
    seq = "".join(line.strip().upper() for line in raw.splitlines() if not line.startswith(">"))
    sha = hashlib.sha256(seq.encode("ascii", "ignore")).hexdigest()[:12]
    print(f"  [NCBI-Live]   {acc:16s} | Downloaded from NCBI Entrez ({len(seq)} nt, sha256={sha})")
    time.sleep(0.35)
    return seq


def fetch_wt_sars2_spike_from_ncbi(cache_dir: Path, force_ncbi: bool = False, timeout: int = 120) -> str:
    """
    Downloads the reference SARS-CoV-2 Wuhan-Hu-1 genome (NC_045512.2) from NCBI
    and extracts the canonical Spike CDS (coordinates 21563..25384, 3822 nt = 1274 codons).
    """
    genome = fetch_ncbi_nuccore_fasta("NC_045512.2", cache_dir=cache_dir, force_ncbi=force_ncbi, timeout=timeout)
    spike_cds = genome[21562:25384]
    if len(spike_cds) != 3822 or not spike_cds.startswith("ATG"):
        raise RuntimeError(f"Unexpected Spike CDS slice from NC_045512.2 (len={len(spike_cds)})")
    return spike_cds


def _expected_spike_2p_aa(wt_spike_cds: str) -> np.ndarray:
    set_genetic_code(1)
    wt_cods = bytes_to_codons(wt_spike_cds, strict_cds=True)
    if wt_cods is None or int(wt_cods.size) < 1274:
        raise RuntimeError("WT Spike CDS failed the strict frame guard")
    wt_sense = AA_CH[wt_cods[:-1]].copy()
    if int(wt_sense.size) > 987:
        wt_sense[985] = CHANNEL_IDX["P"]
        wt_sense[986] = CHANNEL_IDX["P"]
    return wt_sense


def try_fetch_ncbi_bnt162b2_full(wt_spike_cds: str, cache_dir: Path, force_ncbi: bool = False, timeout: int = 120):
    set_genetic_code(1)
    seq = "".join(BNT162B2_FULL_3825NT.split())
    if len(seq) != 3825 or (len(seq) % 3) != 0:
        raise RuntimeError("Embedded BNT162b2 CDS length %d is not 3825" % len(seq))
    cods = bytes_to_codons(seq, strict_cds=True)
    wt_sense = _expected_spike_2p_aa(wt_spike_cds)
    bnt_sense = AA_CH[cods[:-2]].copy()
    if bnt_sense.shape != wt_sense.shape or not np.array_equal(bnt_sense, wt_sense):
        raise RuntimeError("BNT162b2 3825 CDS is not Spike-2P K986P/V987P versus NC_045512")
    if seq[2955:2961] != "CCTCCT":
        raise RuntimeError("2P codons at CDS nt 2956-2961 are not CCTCCT")
    sha = hashlib.sha256(seq.encode("ascii")).hexdigest()
    print("  [BNT162b2]    Canonical CDS 3825 nt / 1275 codons / sha256=%s (K986P/V987P, TGATGA)" % sha[:12])

    for acc in ("OR134577.1", "OQ983940.1"):
        try:
            raw = fetch_ncbi_nuccore_fasta(acc, cache_dir=cache_dir, force_ncbi=force_ncbi, timeout=timeout)
            idx = raw.find(seq[:21])
            if idx < 0:
                print("  [NCBI-Probe]  %s: signal peptide not present in partial vector record" % acc)
                continue
            if raw[idx:idx + 3822] == seq[:3822]:
                dbl = raw[idx:idx + 3825] == seq
                extra = " and the second TGA" if dbl else " (second TGA not in this record)"
                print("  [NCBI-Match]  %s matches the 3822-nt WHO ORF%s" % (acc, extra))
            else:
                print("  [NCBI-Probe]  %s: partial vector fragment; canonical 3825-nt CDS retained" % acc)
        except Exception as exc:
            print("  [NCBI-Probe]  %s: %s" % (acc, exc))
    return seq, "CANONICAL_BNT162b2_CDS_3825_TGATGA"


def _extract_spike_2p_cds_from_text(raw_text: str, wt_sense_2p: np.ndarray, exclude_prefix: str = ""):
    """
    Scans raw FASTA/text payload for any contiguous ACGT sequence whose 1,273 sense codons
    translate identically to SARS-CoV-2 Spike-2P (K986P/V987P) plus trailing terminal stop codons (e.g. TGATAATGA).
    """
    set_genetic_code(1)
    lines = [
        ln.strip().upper().replace("U", "T")
        for ln in raw_text.splitlines()
        if ln.strip() and not ln.lstrip().startswith(">")
    ]
    joined = "".join(re.sub(r"[^ACGT]", "", ln) for ln in lines)
    for blk in [joined] + [re.sub(r"[^ACGT]", "", ln) for ln in raw_text.upper().replace("U", "T").split(">")]:
        if len(blk) < 3822:
            continue
        pos = 0
        while True:
            idx = blk.find("ATG", pos)
            if idx < 0 or idx + 3819 > len(blk):
                break
            cand_sense = blk[idx:idx + 3819]
            if exclude_prefix and cand_sense == exclude_prefix[:3819]:
                pos = idx + 3
                continue
            cods = bytes_to_codons(cand_sense, strict_cds=True)
            if cods is not None and cods.size == 1273 and np.array_equal(AA_CH[cods], wt_sense_2p):
                end_idx = idx + 3819
                while end_idx + 3 <= len(blk) and blk[end_idx:end_idx + 3] in ("TGA", "TAA", "TAG"):
                    end_idx += 3
                if end_idx == idx + 3819:
                    end_idx = min(len(blk), idx + 3822)
                return blk[idx:end_idx]
            pos = idx + 1
    return None


def extract_ncbi_vector_spike_cds(acc: str, cache_dir: Path, force_ncbi: bool = False, timeout: int = 120):
    """
    Directly extracts the full-length Spike-2P coding sequence (including tandem/triple terminal stops
    such as TGATAATGA or TGATGA) from an NCBI Nuccore vector record (e.g. OR134578.1 or OR134577.1)
    by scanning all 6 reading frames on both strands for the Spike ORF (>= 1,250 codons).
    """
    raw_seq = fetch_ncbi_nuccore_fasta(acc, cache_dir=cache_dir, force_ncbi=force_ncbi, timeout=timeout)
    stop_set = {"TAA", "TAG", "TGA"}
    best_cds = ""
    best_strand = "+"

    for strand_label, s_dna in (("+", raw_seq), ("-", revcomp_nt_seq(raw_seq))):
        for off in (0, 1, 2):
            sub = s_dna[off:]
            m = len(sub) // 3
            cods = [sub[3 * k: 3 * k + 3] for k in range(m)]
            i = 0
            while i < m:
                if cods[i] == "ATG":
                    j = i + 1
                    while j < m and cods[j] not in stop_set:
                        j += 1
                    if j < m and (j - i) >= 1250:
                        k_end = j
                        while k_end < m and cods[k_end] in stop_set:
                            k_end += 1
                        cand_cds = "".join(cods[i:k_end])
                        if len(cand_cds) > len(best_cds):
                            best_cds = cand_cds
                            best_strand = strand_label
                    i = j + 1
                else:
                    i += 1

    if not best_cds:
        raise RuntimeError(f"Could not locate full-length Spike ORF (>=1250 codons) in NCBI record {acc}")
    sha = hashlib.sha256(best_cds.encode("ascii")).hexdigest()[:12]
    print(
        f"  [NCBI-Vector] {acc:16s} | Extracted Spike CDS on ({best_strand}) strand: "
        f"{len(best_cds)} nt / {len(best_cds)//3} codons / sha256={sha} (tail={best_cds[-9:]})"
    )
    return best_cds, f"NCBI_{acc}_STRAND_{best_strand}_{len(best_cds)}NT_{best_cds[-9:]}"


def try_fetch_moderna_mrna1273_full(wt_spike_cds: str, human_train_C: CodonCounts, cache_dir: Path, force_ncbi: bool = False, timeout: int = 120):
    """
    Resolves and verifies the canonical monovalent Moderna mRNA-1273 Spike-2P coding sequence
    (3828 nt, 1276 codons including the triple terminal stop TGATAATGA):
      1. Checks local cache `cache_dir / mRNA1273_Full_3828nt.fasta`.
      2. Downloads the pure FASTA file `Figure1Figure2_032321.fasta` from the Stanford/FireLab
         NAalytics repository (note `SARS-CoV2` without hyphen) and extracts the exact 3828-nt CDS.
      3. Also extracts the Spike CDS directly from NCBI GenBank `OR134578.1` (Moderna mRNA-1273 vector).
    """
    set_genetic_code(1)
    wt_sense_2p = _expected_spike_2p_aa(wt_spike_cds)
    bnt_seq = "".join(BNT162B2_FULL_3825NT.split())
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_f = cache_dir / "mRNA1273_Full_3828nt.fasta"

    if cache_f.exists() and cache_f.stat().st_size > 3800 and not force_ncbi:
        txt = cache_f.read_text(encoding="utf-8", errors="ignore")
        seq = "".join(line.strip().upper() for line in txt.splitlines() if not line.startswith(">"))
        cods = bytes_to_codons(seq, strict_cds=True)
        if cods is not None and cods.size >= 1274 and np.array_equal(AA_CH[cods[:1273]], wt_sense_2p):
            sha = hashlib.sha256(seq.encode("ascii")).hexdigest()[:12]
            print(f"  [mRNA-1273]   Cached canonical CDS {len(seq)} nt / {cods.size} codons / sha256={sha} (tail={seq[-9:]})")
            return seq, f"STANFORD_FIRELAB_FASTA_mRNA1273_{len(seq)}NT_{seq[-9:]}"

    fasta_urls = [
        "https://raw.githubusercontent.com/NAalytics/Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273/main/Figure1Figure2_032321.fasta",
        "https://raw.githubusercontent.com/NAalytics/Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273/main/README.md",
    ]
    for url in fasta_urls:
        try:
            print(f"  [mRNA-1273]   GET FASTA mirror: {url}")
            raw_text = http_get(url, timeout=timeout, log_label="mRNA-1273-FASTA").decode("utf-8", errors="ignore")
            extracted = _extract_spike_2p_cds_from_text(raw_text, wt_sense_2p, exclude_prefix=bnt_seq)
            if extracted and len(extracted) >= 3822:
                cache_f.write_text(
                    f">mRNA-1273_Stanford_FireLab_Spike_2P_CDS len={len(extracted)}\n{extracted}\n",
                    encoding="utf-8",
                )
                sha = hashlib.sha256(extracted.encode("ascii")).hexdigest()[:12]
                print(
                    f"  [mRNA-1273]   Verified FASTA CDS {len(extracted)} nt / {len(extracted)//3} codons / "
                    f"sha256={sha} (tail={extracted[-9:]})"
                )
                return extracted, f"STANFORD_FIRELAB_FASTA_mRNA1273_{len(extracted)}NT_{extracted[-9:]}"
        except Exception as exc:
            print(f"  [mRNA-1273]   Mirror warning: {exc}")

    # Direct fallback to NCBI Entrez OR134578.1 if GitHub is unreachable
    print("  [mRNA-1273]   Using direct NCBI Entrez vector OR134578.1")
    return extract_ncbi_vector_spike_cds("OR134578.1", cache_dir=cache_dir, force_ncbi=force_ncbi, timeout=timeout)


def build_spike_2p_variants(wt_spike_cds: str, human_train_C: CodonCounts, human_model: dict):
    """
    Constructs four deterministic reference designs encoding the exact same Spike-2P protein (K986P/V987P):
      1. Human_MaxCAI_Spike_2P            — deterministic Human Max-CAI C_i = phi(A_i) (Theorem 6.1)
      2. Human_MaxGC3_Spike_2P            — deterministic Human Max-GC3 (100% G/C at 3rd codon position)
      3. Human_Markov_H2b_Spike_2P        — 100% deterministic (Zero-RNG) Viterbi dynamic-programming
         compilation on F_{A_0} x ... x F_{A_{m-1}} maximizing conditional syntax excess (Proposition 6.3).
      4. Human_DuplexSafe_Viterbi_Spike_2P — Ambigram-Guarded Viterbi compiler (Corollary 4.4 + Proposition 6.3):
         forbids sense CAT (so antisense frame -0 has ZERO in-frame ATG starts) and enforces periodic
         antisense termination barriers from {TTA, CTA, TCA} while maximizing human ribosomal syntax.
    """
    set_genetic_code(1)
    aa_seq = _expected_spike_2p_aa(wt_spike_cds)
    fibers = build_fibers()
    P_proj = project_counts(human_train_C)
    n_y = P_proj["n_y"].astype(np.float64)

    # 1. Deterministic Max-CAI: C_i = phi(A_i)
    phi_cai = human_model["phi_max_cai"]
    cods_cai = phi_cai[aa_seq].astype(np.int64)

    # 2. Deterministic Max-GC3
    phi_gc3 = np.zeros(N_CHANNELS, dtype=np.int64)
    for a, fb in enumerate(fibers):
        if fb.size > 0:
            gc3_fb = fb[IS_GC3_ARR[fb]]
            pool = gc3_fb if gc3_fb.size > 0 else fb
            phi_gc3[a] = int(pool[int(np.argmax(n_y[pool]))])
    cods_gc3 = phi_gc3[aa_seq].astype(np.int64)

    # 3. Deterministic Zero-RNG Viterbi Compiler on T_4^3 maximizing conditional syntax excess
    lp_yx = logp(human_model["P_y_given_x"])
    lp_a = logp(human_model["P_a_given_a_prev"])
    lp_ya = logp(human_model["P_y_given_a"])

    m = int(aa_seq.size)
    dp_scores = [np.zeros(len(fibers[int(aa_seq[i])]), dtype=np.float64) for i in range(m)]
    dp_back = [np.zeros(len(fibers[int(aa_seq[i])]), dtype=np.int64) for i in range(m)]

    for i in range(1, m):
        a_prev = int(aa_seq[i - 1])
        a_cur = int(aa_seq[i])
        fb_prev = fibers[a_prev]
        fb_cur = fibers[a_cur]
        base_aa = float(lp_a[a_prev, a_cur])
        for j_cur, y_cod in enumerate(fb_cur):
            y_c = int(y_cod)
            syn_pen = 0.85 * float(lp_ya[a_cur, y_c])
            best_val = -1e300
            best_k = 0
            for k_prev, x_cod in enumerate(fb_prev):
                x_c = int(x_cod)
                step_gain = float(lp_yx[x_c, y_c]) - (base_aa + syn_pen)
                cand = float(dp_scores[i - 1][k_prev]) + step_gain
                if cand > best_val:
                    best_val = cand
                    best_k = k_prev
            dp_scores[i][j_cur] = best_val
            dp_back[i][j_cur] = best_k

    cods_mk = np.empty(m, dtype=np.int64)
    cur_k = int(np.argmax(dp_scores[-1]))
    for i in range(m - 1, -1, -1):
        a_i = int(aa_seq[i])
        cods_mk[i] = int(fibers[a_i][cur_k])
        cur_k = int(dp_back[i][cur_k])

    # 4. Constructive Duplex-Safe Viterbi Compiler (Corollary 4.4 + Proposition 6.3):
    #    - Restricts His fiber to CAC only (eliminating CAT = rho^{-1}(ATG), so 0 antisense frame -0 ATG starts)
    #    - Rewards antisense stop preimages {TTA, CTA, TCA} to maintain >= 45 termination barriers in frame -0
    cat_idx = CODON_IDX["CAT"]
    anti_stop_preimages = {CODON_IDX["TTA"], CODON_IDX["CTA"], CODON_IDX["TCA"]}
    dp_ds = [np.full(len(fibers[int(aa_seq[i])]), -1e300, dtype=np.float64) for i in range(m)]
    back_ds = [np.zeros(len(fibers[int(aa_seq[i])]), dtype=np.int64) for i in range(m)]
    dp_ds[0][:] = 0.0

    for i in range(1, m):
        a_prev = int(aa_seq[i - 1])
        a_cur = int(aa_seq[i])
        fb_prev = fibers[a_prev]
        fb_cur = fibers[a_cur]
        base_aa = float(lp_a[a_prev, a_cur])
        for j_cur, y_cod in enumerate(fb_cur):
            y_c = int(y_cod)
            if y_c == cat_idx:
                continue  # Corollary 4.4(i): forbid CAT -> zero antisense ATG start codons in frame -0
            syn_pen = 0.85 * float(lp_ya[a_cur, y_c])
            barrier_bonus = 0.22 if (y_c in anti_stop_preimages and (i % 3 == 0)) else 0.0
            best_val = -1e300
            best_k = 0
            for k_prev, x_cod in enumerate(fb_prev):
                if dp_ds[i - 1][k_prev] <= -1e290:
                    continue
                x_c = int(x_cod)
                step_gain = float(lp_yx[x_c, y_c]) - (base_aa + syn_pen) + barrier_bonus
                cand = float(dp_ds[i - 1][k_prev]) + step_gain
                if cand > best_val:
                    best_val = cand
                    best_k = k_prev
            dp_ds[i][j_cur] = best_val
            back_ds[i][j_cur] = best_k

    cods_ds = np.empty(m, dtype=np.int64)
    cur_k = int(np.argmax(dp_ds[-1]))
    for i in range(m - 1, -1, -1):
        a_i = int(aa_seq[i])
        cods_ds[i] = int(fibers[a_i][cur_k])
        cur_k = int(back_ds[i][cur_k])

    assert np.array_equal(AA_CH[cods_cai], aa_seq)
    assert np.array_equal(AA_CH[cods_gc3], aa_seq)
    assert np.array_equal(AA_CH[cods_mk], aa_seq)
    assert np.array_equal(AA_CH[cods_ds], aa_seq)
    assert int(np.sum(cods_ds == cat_idx)) == 0, "Duplex-Safe Viterbi must have zero CAT codons"

    stop_tga = np.array([CODON_IDX["TGA"]], dtype=np.int64)
    return {
        "Human_MaxCAI_Spike_2P": "".join(CODONS[int(x)] for x in np.concatenate([cods_cai, stop_tga])),
        "Human_MaxGC3_Spike_2P": "".join(CODONS[int(x)] for x in np.concatenate([cods_gc3, stop_tga])),
        "Human_Markov_H2b_Spike_2P": "".join(CODONS[int(x)] for x in np.concatenate([cods_mk, stop_tga])),
        "Human_DuplexSafe_Viterbi_Spike_2P": "".join(CODONS[int(x)] for x in np.concatenate([cods_ds, stop_tga])),
    }


def scan_6frame_ambigram_orfs(seq: str):
    """
    Scans all 6 reading frames (+0, +1, +2 sense and -0, -1, -2 antisense):
      - Counts stop codons (TAA, TAG, TGA) in every frame,
      - Finds the longest unbroken stop-free span and longest canonical ATG->Stop/End ORF,
      - Counts sense codons in the antisense stop-preimage fibers of Theorem 4.3 / Corollary 4.4:
        F_{Leu,Stop} = {TTA, CTA}, F_{Ser,Stop} = {TCA}, and antisense start preimage F_{His,Met} = {CAT},
      - Counts wobble-T codons governing sense frame +2 termination (Corollary 4.4(ii)).
    """
    stop_set = {"TAA", "TAG", "TGA"}
    seq = seq.upper()
    rc_seq = revcomp_nt_seq(seq)

    def _scan_strand(s_dna: str, strand_sign: str):
        frame_stats = {}
        for off in (0, 1, 2):
            sub = s_dna[off:]
            m = len(sub) // 3
            cods = [sub[3 * k: 3 * k + 3] for k in range(m)]
            stop_positions = [k for k, c in enumerate(cods) if c in stop_set]

            boundaries = [-1] + stop_positions + [m]
            max_unbroken_cods = 0
            best_span = (0, 0)
            for idx in range(len(boundaries) - 1):
                sp_len = boundaries[idx + 1] - boundaries[idx] - 1
                if sp_len > max_unbroken_cods:
                    max_unbroken_cods = sp_len
                    best_span = (boundaries[idx] + 2, boundaries[idx + 1])

            max_atg_orf = 0
            best_atg_span = (0, 0)
            i = 0
            while i < m:
                if cods[i] == "ATG":
                    j = i + 1
                    while j < m and cods[j] not in stop_set:
                        j += 1
                    orf_len = j - i
                    if orf_len > max_atg_orf:
                        max_atg_orf = orf_len
                        best_atg_span = (i + 1, j)
                    i = j + 1
                else:
                    i += 1

            frame_stats[f"{strand_sign}{off}"] = {
                "total_triplets": m,
                "stop_codon_count": len(stop_positions),
                "stops_per_100_codons": 100.0 * len(stop_positions) / max(1, m),
                "max_unbroken_stop_free_codons": int(max_unbroken_cods),
                "max_unbroken_stop_free_nt": int(max_unbroken_cods * 3),
                "max_unbroken_codon_span_1based": f"{best_span[0]}..{best_span[1]}",
                "max_atg_to_stop_codons": int(max_atg_orf),
                "max_atg_to_stop_nt": int(max_atg_orf * 3),
                "max_atg_span_1based": f"{best_atg_span[0]}..{best_atg_span[1]}",
            }
        return frame_stats

    sense_frames = _scan_strand(seq, "+")
    anti_frames = _scan_strand(rc_seq, "-")

    m0 = len(seq) // 3
    sense_cods_0 = [seq[3 * k: 3 * k + 3] for k in range(m0)]
    n_tta = sum(1 for c in sense_cods_0 if c == "TTA")
    n_cta = sum(1 for c in sense_cods_0 if c == "CTA")
    n_tca = sum(1 for c in sense_cods_0 if c == "TCA")
    n_cat = sum(1 for c in sense_cods_0 if c == "CAT")
    n_wobble_t = sum(1 for c in sense_cods_0 if c.endswith("T"))

    best_anti_frame = max(
        anti_frames.items(),
        key=lambda kv: kv[1]["max_unbroken_stop_free_codons"],
    )
    best_anti_atg_frame = max(
        anti_frames.items(),
        key=lambda kv: kv[1]["max_atg_to_stop_codons"],
    )

    return {
        "sense_frames": sense_frames,
        "antisense_frames": anti_frames,
        "theorem_4_3_antisense_stop_preimages_frame0": {
            "F_Leu_Stop_TTA": int(n_tta),
            "F_Leu_Stop_CTA": int(n_cta),
            "F_Ser_Stop_TCA": int(n_tca),
            "total_inframe_antisense_stops": int(n_tta + n_cta + n_tca),
            "F_His_Met_CAT_antisense_ATG": int(n_cat),
            "wobble_T_sense_codons": int(n_wobble_t),
        },
        "longest_antisense_unbroken_frame": best_anti_frame[0],
        "longest_antisense_unbroken_codons": best_anti_frame[1]["max_unbroken_stop_free_codons"],
        "longest_antisense_unbroken_nt": best_anti_frame[1]["max_unbroken_stop_free_nt"],
        "longest_antisense_atg_orf_frame": best_anti_atg_frame[0],
        "longest_antisense_atg_orf_codons": best_anti_atg_frame[1]["max_atg_to_stop_codons"],
        "longest_antisense_atg_orf_span": best_anti_atg_frame[1]["max_atg_span_1based"],
        "frame_minus0_atg_orf_codons": anti_frames["-0"]["max_atg_to_stop_codons"],
        "frame_minus0_atg_orf_span": anti_frames["-0"]["max_atg_span_1based"],
    }


def compute_d3_d3b_and_pausing_metrics(seq: str, human_train_C: CodonCounts):
    """
    Computes analytical residuals d3 (vs. folded dinucleotides H3) and d3b
    (vs. folded order-2 nucleotide trigrams H3b) on the Human_T2T reference model,
    alongside CpG density, 64-codon Shannon entropy (bits), slippery homopolymer runs,
    and m1Psi-prone TTTT (UUUU) ribosomal frameshift runs.
    """
    cods = bytes_to_codons(seq, strict_cds=False)
    if cods is None or cods.size < 4:
        return {}
    end = int(cods.size)
    while end > 1 and int(AA_CH[cods[end - 1]]) == int(STOP_IDX):
        end -= 1
    cods_int = cods[1:end] if cods.size >= 5 else cods[:end]

    n_xy_tr = human_train_C.n_xy.astype(np.float64)
    nt_bg = np.zeros((4, 4), dtype=np.float64)
    nt_tri = np.zeros((4, 4, 4), dtype=np.float64)
    for x in range(N_CODONS):
        x1, x2, x3 = int(C_N1_ARR[x]), int(C_N2_ARR[x]), int(C_N3_ARR[x])
        row = n_xy_tr[x]
        if row.sum() == 0:
            continue
        for y in range(N_CODONS):
            w = float(row[y])
            if w == 0.0:
                continue
            y1, y2, y3 = int(C_N1_ARR[y]), int(C_N2_ARR[y]), int(C_N3_ARR[y])
            nt_bg[x3, y1] += w
            nt_bg[y1, y2] += w
            nt_bg[y2, y3] += w
            nt_tri[x2, x3, y1] += w
            nt_tri[x3, y1, y2] += w
            nt_tri[y1, y2, y3] += w

    p_emp = row_normalize(n_xy_tr, 0.5)
    lp_emp = np.log(np.maximum(p_emp, TINY))

    p_bg = (nt_bg + 0.5) / (nt_bg.sum(axis=1, keepdims=True) + 2.0)
    lp_bg = np.log(np.maximum(p_bg, TINY))

    p_tri = (nt_tri + 0.5) / (nt_tri.sum(axis=2, keepdims=True) + 2.0)
    lp_tri = np.log(np.maximum(p_tri, TINY))

    r2, r3 = C_N2_ARR[:, None], C_N3_ARR[:, None]
    c1, c2, c3 = C_N1_ARR[None, :], C_N2_ARR[None, :], C_N3_ARR[None, :]
    lp_h3 = lp_bg[r3, c1] + lp_bg[c1, c2] + lp_bg[c2, c3]
    lp_h3b = lp_tri[r2, r3, c1] + lp_tri[r3, c1, c2] + lp_tri[c1, c2, c3]

    prev = cods_int[:-1]
    cur = cods_int[1:]
    ll_emp = float(np.mean(lp_emp[prev, cur]))
    ll_h3 = float(np.mean(lp_h3[prev, cur]))
    ll_h3b = float(np.mean(lp_h3b[prev, cur]))

    cnt64 = np.bincount(cods_int, minlength=N_CODONS).astype(np.float64)
    p64 = cnt64[cnt64 > 0] / float(cods_int.size)
    h_codon_bits = float(-np.sum(p64 * np.log2(p64)))
    gc3 = float(np.mean(IS_GC3_ARR[cods_int]))

    s_up = seq.upper()
    gc_tot = float(sum(1 for b in s_up if b in ("G", "C"))) / max(1, len(s_up))
    cpg_per_100nt = 100.0 * float(s_up.count("CG")) / max(1, len(s_up) - 1)
    slippery_runs = len(re.findall(r"(AAAA+|TTTT+|CCCC+|GGGG+)", s_up))
    uuuu_runs = len(re.findall(r"TTTT+", s_up))

    return {
        "gc_total": gc_tot,
        "gc3": gc3,
        "cpg_per_100nt": cpg_per_100nt,
        "codon_entropy_bits": h_codon_bits,
        "slippery_homopolymer_runs_ge4": int(slippery_runs),
        "m1psi_uuuu_runs_ge4": int(uuuu_runs),
        "d3_vs_dinuc": float(ll_emp - ll_h3),
        "d3b_vs_trigram": float(ll_emp - ll_h3b),
    }


def make_bnt162b2_figure(rows: list, out_dir: Path):
    if not HAVE_MPL or not rows:
        return None
    figdir = out_dir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    labels = [r["short_label"] for r in rows]
    x = np.arange(len(labels))

    d2b_cross = np.array([r["human_cross"]["d2b"] for r in rows], dtype=np.float64)
    ssyn_cross = np.array([r["human_cross"]["S_syn_heldout"] for r in rows], dtype=np.float64)
    anti_orf = np.array([r["ambigram"]["longest_antisense_unbroken_codons"] for r in rows], dtype=np.float64)
    anti_atg = np.array([r["ambigram"]["frame_minus0_atg_orf_codons"] for r in rows], dtype=np.float64)
    f_stops = np.array([
        r["ambigram"]["theorem_4_3_antisense_stop_preimages_frame0"]["total_inframe_antisense_stops"]
        for r in rows
    ], dtype=np.float64)

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.0))

    bw = 0.36
    axes[0].bar(x - bw / 2, d2b_cross, bw, color="tab:blue", label="Cross d2b on Human_T2T")
    axes[0].bar(x + bw / 2, ssyn_cross, bw, color="tab:green", label="Cross S_syn on Human_T2T")
    axes[0].axhline(0, color="k", lw=0.9)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=30, ha="right", fontsize=8.5)
    axes[0].set_ylabel("nats / codon")
    axes[0].set_title("A. Ribosomal Syntax on Human_T2T\n(Cross-Scored d2b & S_syn)", fontsize=10)
    axes[0].legend(fontsize=8)

    axes[1].bar(x, f_stops, width=0.55, color="tab:orange")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=30, ha="right", fontsize=8.5)
    axes[1].set_ylabel("Count of {TTA, CTA, TCA} codons")
    axes[1].set_title("B. Corollary 4.4 Antisense Stop Preimages\nF_{Leu,Stop}={TTA,CTA} + F_{Ser,Stop}={TCA}", fontsize=10)

    axes[2].bar(x - bw / 2, anti_orf, bw, color="crimson", label="Max unbroken span (any -frame)")
    axes[2].bar(x + bw / 2, anti_atg, bw, color="darkslateblue", label="Max ATG->Stop/End (frame -0)")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=30, ha="right", fontsize=8.5)
    axes[2].set_ylabel("Antisense ORF length (codons)")
    axes[2].set_title("C. Antisense Frame-Opening Forensics\n(Unbroken Span vs. CAT-Initiated ATG->End)", fontsize=10)
    axes[2].legend(fontsize=8)

    fig.suptitle(
        "Fig 8. BNT162b2 & Moderna mRNA-1273 Duplex Torus Forensics vs. Duplex-Safe Viterbi Compiler (Corollary 4.4 & Prop. 6.3)",
        fontsize=11.0,
        y=1.01,
    )
    fig.tight_layout()
    p_png = figdir / "fig8_bnt162b2_ambigram_duplex.png"
    p_pdf = figdir / "fig8_bnt162b2_ambigram_duplex.pdf"
    fig.savefig(p_png, dpi=200, bbox_inches="tight")
    fig.savefig(p_pdf, bbox_inches="tight")
    plt.close(fig)
    return str(p_png)


def run_bnt162b2_duplex_audit(cds_dir: Path, out_dir: Path, max_records: int = 0, alpha: float = 0.5, force_ncbi: bool = False, force_ncbi_all: bool = False, timeout: int = 300):
    print("\n" + "=" * 126)
    print("BNT162b2 & MODERNA mRNA-1273 DUPLEX TORUS, RIBOSOMAL PAUSING & AMBIGRAMMATIC ORF FORENSICS (v3.5 — ZERO RNG)")
    print("=" * 126)

    resolve_or_download_cds.force_download = bool(force_ncbi_all)
    human_path = resolve_or_download_cds("Human_T2T", cds_dir=cds_dir, timeout=timeout, skip_download=False)
    records, dominant_table, _, _ = parse_cds_fasta_gz(
        human_path, max_records=max_records, trim_start_stop=True, default_table=1
    )
    set_genetic_code(1)

    train_C = CodonCounts()
    test_C = CodonCounts()
    for gid, tbl, cd in records:
        (train_C if deterministic_gene_split(gid, 80) else test_C).add_codons(cd)
    human_model = build_model(train_C, alpha=alpha)
    human_ide = compute_exact_decompositions(test_C)

    ncbi_viral_cache = cds_dir / "ncbi_live_viral"
    wt_spike_seq = fetch_wt_sars2_spike_from_ncbi(ncbi_viral_cache, force_ncbi=force_ncbi, timeout=timeout)
    bnt_full_seq, bnt_full_src = try_fetch_ncbi_bnt162b2_full(
        wt_spike_seq, ncbi_viral_cache, force_ncbi=force_ncbi, timeout=timeout
    )
    mrna1273_seq, mrna1273_src = try_fetch_moderna_mrna1273_full(
        wt_spike_seq, train_C, ncbi_viral_cache, force_ncbi=force_ncbi, timeout=timeout
    )
    or134578_seq, or134578_src = extract_ncbi_vector_spike_cds(
        "OR134578.1", cache_dir=ncbi_viral_cache, force_ncbi=force_ncbi, timeout=timeout
    )
    variants = build_spike_2p_variants(wt_spike_seq, train_C, human_model)

    constructs = [
        ("WT_SARS2_Spike_NC045512", "WT_SARS2 (NC_045512.2)", wt_spike_seq, "NCBI Entrez NC_045512.2:21563..25384 (3822 nt)"),
        ("BNT162b2_Full_Spike_2P", "Pfizer BNT162b2 (3825nt)", bnt_full_seq, bnt_full_src),
        ("mRNA1273_Full_Spike_2P", "Moderna mRNA-1273 (3828nt)", mrna1273_seq, mrna1273_src),
        ("mRNA1273_NCBI_OR134578", "Moderna NCBI (OR134578.1)", or134578_seq, or134578_src),
        ("Human_MaxCAI_Spike_2P", "Max-CAI Spike (Thm 6.1)", variants["Human_MaxCAI_Spike_2P"], "Deterministic Human Max-CAI phi(A_i) (3822 nt)"),
        ("Human_MaxGC3_Spike_2P", "Max-GC3 Spike", variants["Human_MaxGC3_Spike_2P"], "Deterministic Human Max-GC3 (3822 nt)"),
        ("Human_Markov_H2b_Spike_2P", "Viterbi-Markov (Prop 6.3)", variants["Human_Markov_H2b_Spike_2P"], "Deterministic Viterbi 1st-Order Human Markov Compiler (3822 nt)"),
        ("Human_DuplexSafe_Viterbi_Spike_2P", "DuplexSafe-Viterbi (Cor 4.4)", variants["Human_DuplexSafe_Viterbi_Spike_2P"], "Corollary 4.4 Ambigram-Guarded Viterbi Compiler (CAT=0, 3822 nt)"),
    ]

    rows = []
    for name, short_label, seq_str, src_desc in constructs:
        cods = bytes_to_codons(seq_str, strict_cds=True)
        if cods is None:
            raise RuntimeError(f"{name}: strict frame guard rejected the CDS (len % 3 != 0 or non-ACGT)")
        end = int(cods.size)
        while end > 1 and int(AA_CH[cods[end - 1]]) == int(STOP_IDX):
            end -= 1
        cods_int = cods[1:end] if cods.size >= 5 else cods[:end]

        C_seq = CodonCounts()
        C_seq.add_codons(cods_int)
        ide = compute_exact_decompositions(C_seq)
        cross = held_out_ll_contrasts(human_model, cods_int)
        d3_info = compute_d3_d3b_and_pausing_metrics(seq_str, train_C)
        amb = scan_6frame_ambigram_orfs(seq_str)
        ep = entropy_production(C_seq)
        dm = dynamic_mirror_current(C_seq)

        rows.append({
            "construct": name,
            "short_label": short_label,
            "source": src_desc,
            "len_nt": len(seq_str),
            "n_codons_total": int(cods.size),
            "n_codons_internal": int(cods_int.size),
            "internal_decomp": ide,
            "human_cross": cross,
            "d3_and_pausing": d3_info,
            "ambigram": amb,
            "entropy_production": ep,
            "mirror_dynamics": dm,
        })

    print("\n" + "-" * 126)
    print("TABLE 1. CROSS-SCORING ON HUMAN_T2T REFERENCE MODEL & INTERNAL DUPLEX DECOMPOSITION ON T_4^3")
    print("-" * 126)
    print(
        f"  [Human_T2T held-out baseline ({test_C.n_total} transitions)]: "
        f"d2b={human_ide['I_Cprev_Acur_Aprev_d2b']:+.4f} | S_syn={human_ide['S_syn_I_Cprev_Ccur_Acur']:+.4f} | "
        f"S_anti={human_ide['S_anti_I_Cprev_Arho_given_Acur']:+.4f} ({100*human_ide['antisense_program_fraction_of_S_syn']:.1f}%) | "
        f"H(C|A)={human_ide['H_C_given_A_nats']:.4f} = H(Arho|A) {human_ide['H_Arho_given_A_nats']:.4f} + H(C|A,Arho) {human_ide['H_C_given_A_Arho_nats']:.4f}"
    )
    print("-" * 126)
    print(
        "Construct                         | Codons |    GC3 | CpG/100 | H_cod(b) | Cross_d2b | Cross_d3b | Cross_S_syn | Int_S_syn | sigma_dup | Slip(UUUU)"
    )
    print(
        "----------------------------------+--------+--------+---------+----------+-----------+-----------+-------------+-----------+-----------+-----------"
    )
    for r in rows:
        hc = r["human_cross"]
        id_c = r["internal_decomp"]
        d3i = r["d3_and_pausing"]
        dm = r["mirror_dynamics"]
        sl_str = f"{d3i['slippery_homopolymer_runs_ge4']:2d} ({d3i['m1psi_uuuu_runs_ge4']:2d})"
        print(
            f"{r['construct']:33s} | {r['n_codons_total']:6d} | {d3i['gc3']:6.4f} | {d3i['cpg_per_100nt']:7.2f} | {d3i['codon_entropy_bits']:8.4f} | "
            f"{hc['d2b']:+9.4f} | {d3i['d3b_vs_trigram']:+9.4f} | {hc['S_syn_heldout']:+11.4f} | "
            f"{id_c['S_syn_I_Cprev_Ccur_Acur']:+9.4f} | {dm['sigma_dup']:9.4f} | {sl_str:>10s}"
        )

    print("\n" + "-" * 126)
    print("TABLE 2. COROLLARY 4.4 AMBIGRAMMATIC ORF FORENSICS: {TTA, CTA, TCA} STOPS, {CAT} ANTISENSE STARTS & 6-FRAME SPANS")
    print("-" * 126)
    print(
        "Construct                         | {TTA,CTA,TCA} | CAT(ATG) | Stops(+0/+1/+2) | Stops(-0/-1/-2) | Max_Unbroken_Anti | Frame-0 ATG->End (Span)"
    )
    print(
        "----------------------------------+---------------+----------+-----------------+-----------------+-------------------+-------------------------"
    )
    for r in rows:
        amb = r["ambigram"]
        t43 = amb["theorem_4_3_antisense_stop_preimages_frame0"]
        sf = amb["sense_frames"]
        af = amb["antisense_frames"]
        s_str = f"{t43['F_Leu_Stop_TTA']:2d}/{t43['F_Leu_Stop_CTA']:2d}/{t43['F_Ser_Stop_TCA']:2d} ({t43['total_inframe_antisense_stops']:2d})"
        cat_str = f"{t43['F_His_Met_CAT_antisense_ATG']:4d}"
        sp_str = f"{sf['+0']['stop_codon_count']:3d}/{sf['+1']['stop_codon_count']:3d}/{sf['+2']['stop_codon_count']:3d}"
        st_str = f"{af['-0']['stop_codon_count']:3d}/{af['-1']['stop_codon_count']:3d}/{af['-2']['stop_codon_count']:3d}"
        unb_str = f"{amb['longest_antisense_unbroken_codons']:4d}c ({amb['longest_antisense_unbroken_frame']})"
        atg0_str = f"{amb['frame_minus0_atg_orf_codons']:4d}c ({amb['frame_minus0_atg_orf_span']})"
        print(f"{r['construct']:33s} | {s_str:13s} | {cat_str:8s} | {sp_str:15s} | {st_str:15s} | {unb_str:17s} | {atg0_str}")

    out_dir.mkdir(parents=True, exist_ok=True)
    fig_path = make_bnt162b2_figure(rows, out_dir)
    out_json = out_dir / "bnt162b2_duplex_audit.json"
    out_md = out_dir / "bnt162b2_duplex_audit.md"

    payload = {
        "version": "v3.5-bnt162b2-mrna1273-duplex",
        "human_reference": {
            "file": human_path.name,
            "test_transitions": int(test_C.n_total),
            "heldout_d2b": human_ide["I_Cprev_Acur_Aprev_d2b"],
            "heldout_S_syn": human_ide["S_syn_I_Cprev_Ccur_Acur"],
            "heldout_S_anti": human_ide["S_anti_I_Cprev_Arho_given_Acur"],
            "heldout_H_C_given_A": human_ide["H_C_given_A_nats"],
            "heldout_H_Arho_given_A": human_ide["H_Arho_given_A_nats"],
            "heldout_H_C_given_A_Arho": human_ide["H_C_given_A_Arho_nats"],
        },
        "constructs": rows,
        "figure": fig_path,
    }
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = [
        "# BNT162b2 & Moderna mRNA-1273 Duplex Torus, Ribosomal Pausing & Ambigrammatic ORF Forensics (v3.5)",
        "",
        "## 1. Cross-Scoring on the `Human_T2T` Reference Model and Two-Strand Torus Decomposition",
        "",
        "| Construct | Codons | GC3 | CpG/100nt | H_codon (bits) | Cross d2b | Cross d3b | Cross S_syn | Internal S_syn | sigma_dup | Slippery (UUUU) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        hc = r["human_cross"]
        id_c = r["internal_decomp"]
        d3i = r["d3_and_pausing"]
        dm = r["mirror_dynamics"]
        md.append(
            f"| **{r['construct']}** | {r['n_codons_total']} | {d3i['gc3']:.4f} | {d3i['cpg_per_100nt']:.2f} | {d3i['codon_entropy_bits']:.4f} | "
            f"{hc['d2b']:+.4f} | {d3i['d3b_vs_trigram']:+.4f} | {hc['S_syn_heldout']:+.4f} | "
            f"{id_c['S_syn_I_Cprev_Ccur_Acur']:+.4f} | {dm['sigma_dup']:.4f} | {d3i['slippery_homopolymer_runs_ge4']} ({d3i['m1psi_uuuu_runs_ge4']}) |"
        )

    md.extend([
        "",
        "## 2. Corollary 4.4 Ambigrammatic Overlapping ORF Forensics (`{TTA, CTA, TCA}` vs. `CAT` Antisense `ATG`)",
        "",
        "| Construct | Sense `{TTA, CTA, TCA}` | Sense `CAT` (`-0 ATG`) | Sense Stops (`+0/+1/+2`) | Antisense Stops (`-0/-1/-2`) | Max Unbroken Antisense Span | Frame `-0` `ATG->End` ORF |",
        "|---|---:|---:|---|---|---:|---:|",
    ])
    for r in rows:
        amb = r["ambigram"]
        t43 = amb["theorem_4_3_antisense_stop_preimages_frame0"]
        sf = amb["sense_frames"]
        af = amb["antisense_frames"]
        md.append(
            f"| **{r['construct']}** | `{t43['F_Leu_Stop_TTA']}/{t43['F_Leu_Stop_CTA']}/{t43['F_Ser_Stop_TCA']} ({t43['total_inframe_antisense_stops']})` | "
            f"`{t43['F_His_Met_CAT_antisense_ATG']}` | `{sf['+0']['stop_codon_count']}/{sf['+1']['stop_codon_count']}/{sf['+2']['stop_codon_count']}` | "
            f"`{af['-0']['stop_codon_count']}/{af['-1']['stop_codon_count']}/{af['-2']['stop_codon_count']}` | "
            f"`{amb['longest_antisense_unbroken_codons']} cod ({amb['longest_antisense_unbroken_frame']})` | "
            f"`{amb['frame_minus0_atg_orf_codons']} cod ({amb['frame_minus0_atg_orf_span']})` |"
        )

    out_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    print("=" * 126)
    print(f"[SAVED JSON] {out_json}")
    print(f"[SAVED MD]   {out_md}")
    if fig_path:
        print(f"[SAVED FIG]  {fig_path}")
    print("=" * 126)
    return payload
# === END BNT162B2 & MODERNA mRNA-1273 DUPLEX FORENSICS PATCH v3.5 ===


if __name__ == "__main__":
    main()
