#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
CODON MEMORY SPECTRUM  v3.0  —  codon_memory_spectrum_v3.py
================================================================================
Unified, autonomous, zero-RNG analytic engine & publication figure suite for
sequential memory on the double-stranded codon torus T_4^3 across 100 RefSeq
genomes (merging engine v2.4 + replot_memory_v24.py into 1 standalone script).

--------------------------------------------------------------------------------
CHANGELOG v2.4 -> v3.0
--------------------------------------------------------------------------------
[UNIFIED] Single autonomous script: merges replot_memory_v24.py directly into
  the engine with --replot-only JSON, --dpi, --plot-format, and --plot-mode
  (auto / lines / clade / both), producing publication-grade 300-DPI Figs 1-7.
  Supports both v3.0 and legacy v2.4 JSON payloads seamlessly.

[BIO-100] Expanded curated cellular biosphere catalog from 54 to 100 RefSeq
  species (--biosphere-100) spanning 12 taxonomic clades and FIVE distinct NCBI
  translation tables (Tables 1, 4, 6, 11, 12):
    * Table 4 (TGA->W, q_s=62, 3721 df): M. mycoides, M. pneumoniae, M. genitalium
    * Table 6 (TAA/TAG->Q, q_s=63, 3844 df, 2 duplex quartets): Tetrahymena, Paramecium
    * Table 12 (CTG->S, b_1=25 invariant): C. albicans, D. hansenii, S. stipitis
  Includes automatic NCBI E-utilities (esearch -> esummary -> FtpPath_RefSeq)
  fallback resolution + manifest_v30.json cache + --resume checkpointing.

[THM-AB]  EXACT ZERO-RNG PERIODICITY & C2H2 ZINC-FINGER STRATIFICATION OF
  DIRECTED KLEIN FLUXES (resolves Paper III Policy 10.1).
  Translates each test CDS via its genetic code and classifies genes with zero
  external annotations into C2H2 Kruppel zinc-finger arrays (C.{2,4}C.{12}H.{3,5}H
  repeated >=2x, ~28-aa period), collagen G-X-Y repeats, and non-ZNF genes.
  Computes exact log-detrended DFT power spectra over k=5..60, exact Delta_k=14
  dip depths at k in {14,28,42,56} vs flanks, and compares full vs non-ZNF vs
  ZNF-only cohorts to test causally whether periodic dips live in ZNF arrays.

[THM-AD]  ISOCHORE-FREE SCREEN DECOMPOSITION & EXACT SCREEN-STRATIFIED KL CHAIN
  RULE AGAINST E_het^(k).
  Answers whether the long-range k=30 plateau in D_3'wobble and Delta_x1 is
  driven by between-gene GC3 isochores or true within-gene memory:
    (1) Per-layer Miller-Madow bias correction (Lemma 7.1) minus exact E_het^(k)
        screen MI gives isochore-stripped intra-gene excess l_intra^MM(k).
    (2) Exact forward & reverse 3-term screen KL chain rules (all terms >= 0,
        exact to <1e-14, exact antiparallel duplex duality under R):
          D_KL(P || P_het) = K_bridge + K_lead|SO + K_wob|XO
                           = K_bridge + K_wob|SO  + K_lead|SY.
        Reduces identically to Theorem T when P_het -> P_ind = p_X x p_Y!
        K_lead|SY(P || P_het) >= 0 is the exact isochore-free analogue of Delta_x1,
        and K_nonbridge = D_KL(P || P_het) - K_bridge >= 0 is exact isochore-free
        memory invisible to any 16x16 dinucleotide screen.

[THM-AC]  EXACT ZERO-RNG BIOSPHERE PHYLOGENY & MULTI-CODE FINGERPRINT MATRIX.
  Computes S x S pairwise geometric distances across all analyzed taxa in four
  orthogonal channels: (i) scalar GC3, (ii) static form (W_++ + W_+-),
  (iii) directed program (W_-+ + W_--), and (iv) 93.75% non-bridge screen
  complement, plus exact Leave-One-Out 1-NN clade purity and exact cross-species
  held-out d2b_transfer / d3b_transfer matrices.

[FIX-52] Corrected stale FIX-40 docstring note (Tables 3 & 5 have q_s=62, not 61).

--------------------------------------------------------------------------------
LEGACY CHANGELOG (v1.0 -> v2.4)
--------------------------------------------------------------------------------
CODON MEMORY SPECTRUM  v2.1  —  codon_memory_spectrum_v1.py
================================================================================
Zero-RNG analytic engine for sequential memory in real coding DNA on T_4^3.
--------------------------------------------------------------------------------
CHANGELOG v2.0 -> v2.1
--------------------------------------------------------------------------------
[FIX-40] EXACT SENSE-CONVENTION DEGREES OF FREEDOM OF THE 16x16 SCREEN.
v2.0 carried df_delta_x1_sense = 16*3*60 = 2880.  That number mixes two
incompatible conventions: a FULL-torus left margin (4 leading bases in
every one of the 16 suffixes => 48) with a SENSE right margin (60).
Deleting the three stop triplets TAA=(3,0,0), TAG=(3,0,2), TGA=(3,2,0)
removes codons on BOTH sides of the pair table:
  suffix AA (s=0): TAA gone -> R_0 = 3
  suffix AG (s=2): TAG gone -> R_2 = 3
  suffix GA (s=8): TGA gone -> R_8 = 3
  prefix TA (o=12): TAA,TAG gone -> C_12 = 2
  prefix TG (o=14): TGA gone     -> C_14 = 3
hence sum_s (R_s-1) = sum_o (C_o-1) = q_s - 16 = 61 - 16 = 45  (NOT 48),
and the exact sense partition of the SAME 3600 df as Theorem G.2 is
      bridge 225 + lead 675 + wobble 675 + d13 2025 = 3600 = 60^2,
      Delta_x1 = 675 + 2025 = 2700,   I_suffix16 = 225 + 675 = 900,
      (G.2 cross-check: 361 + 779 + 2040 + 420 = 3600).
The formula is now computed GENERICALLY from the translation table:
  Table 1/11/12: q_s=61 -> 3600 ; Table 2: q_s=60 -> 3481 ;
  Table 3/4/5: q_s=62 -> 3721 ; Table 6: q_s=63 -> 3844.
Kernel parameter accounting in the sense convention is also corrected:
      675 + 675 + 2025 + 240 (bridge A) + 45 (wobble R) = 3660 = 61*60,
because R(y3|o) has C_o-1 free parameters per box, i.e. 45 and not 48.
NO empirical estimator changes: this is pure df bookkeeping.

[THM-T.5] FIVE-TERM HELD-OUT DECOMPOSITION AND EXACT PHASE CANCELLATION.
Substituting Thm T.2 (Delta_x1 = D_5'lead + D_1<->3) and Thm T.4
(D_suffix = D_3'wobble + D_bridge_kl + D_wobble_kl) into FIX-39 gives
  d3b(k) = D_5'lead + D_3'wobble + D_1<->3 + D_bridge_kl + D_wobble_kl - D_gen.
The three MI layers and D_gen do NOT depend on the trigram null, hence
  d3b - d3b_phase = (D_bridge_kl - D_bridge_kl_phase)
                  + (D_wobble_kl - D_wobble_kl_phase)   EXACTLY.
Moreover FIX-37 gives C_ph[2] = codon marginal, so
  R_wobble_phase(y3|o) = P_train(Y3 | Y1,Y2) identically,
and D_wobble_kl_phase reduces to a pure train->test marginal shift
(zero on train==test up to Dirichlet smoothing; verified with alpha=1e-9).

[THM-U]  DUPLEX MEANING OF THE SCREEN AND THE 8+8 BOX PARTITION.
O(c)=(c1,c2) is the SENSE codon box, c3 the SENSE wobble; the antisense
codon rho(c)=(k(c3),k(c2),k(c1)) has box determined by S(c)=(c2,c3) and
wobble k(c1).  Therefore
  I_bridge      = I(antisense box_t ; sense box_{t+k}),
  D_3'wobble    = memory of the SENSE wobble at fixed sense box,
  D_5'lead      = memory of the ANTISENSE wobble at fixed antisense box,
  D_1<->3       = pure long-range wobble-wobble coupling of the two strands.
In the standard code the 16 boxes split 8+8 into unsplit 4-boxes
  B4 = {CT,GT,TC,CC,AC,GC,CG,GG}   (wobble choice is 100% synonymous)
and split boxes {TT,AT,TA,CA,AA,GA,TG,AG}, with S4 = rho~(B4) =
{AG,AC,GA,GG,GT,GC,CG,CC}.  Both layers split additively and the parts
map into each other COMPONENTWISE under the antiparallel involution:
  (D_5'lead^S4)^anti = D_3'wobble^B4 ,  (D_5'lead^Ssplit)^anti = D_3'wobble^Bsplit.
Masks are derived from the actual translation table (Table 12 splits CT).

--------------------------------------------------------------------------------
RETAINED (v1.0 -> v2.0)
--------------------------------------------------------------------------------
FIX-37  exact linear projection of the nucleotide trigram tensor:
        C_ph[2]=marginal, C_ph[0]/C_ph[1] = projections of B_1,
        n + 2(n-1) = 3n-2 exactly; H3b at EVERY lag is a deterministic
        functional of (B_1^train, marginal^train).
FIX-38  frame cohorts built ONCE for all lags (Theorem O).
FIX-39  exact held-out three-term identity
        d3b = Delta_x1_test + D_suffix_test - D_gen, all terms >= 0.
THM-T   two-sided 16x16 Markov screen, four-term MI split
        225+720+720+2304 = 3969 = 63^2 (full torus), duplex duality,
        three-term KL split of D_suffix, rank-16 bottleneck,
        exact semigroup Q_1^k == Q_k (homogeneous and 3-phase).
A/A.2/A.3/A.4  bilinear analytic_h4; J_g=ceil(n/L)-1; gamma; R_junc;
        cohort bias; Psi_LR = I + D_mix_LR - D_mismatch; Jeffreys identity
        and the FALSE J>=2I counterexample guard.
B       held-out memory ladder + coverage diagnostics.
C/C.2/C.3 exact Markov lag prediction; hypergeometric between-gene null;
        packed window-matched null on the identical pair set.
D       exact zero-RNG H5 (mass invariant) / H6 / H7 on a COMMON cohort.
E/E.2   Delta_H4 = gamma[Psi + d_pos - d_isochore - d_endpoint - d_finite].
F/F.2   sorted Perron spectrum; exact chi2 proxy; biorthogonal
        reconstruction; kappa_j >= 1 non-normality weights.
G/G.2   lag-resolved tripartite + duplex split; exact df partition
        361+779+2040+420 = 3600 = 60^2; effective df; gated enrichment.
H/H.3   tower of exact layer nulls (2F BLAS calls); B << E.
I/K     one-sided four-term KL chain rule; D8 localisation.
J       canonical rho-invariant CpG / Leu quartet split of D8.
L       exact time-reversal symmetry of the F=1 null; JS >= A^2/2.
M/N/O   Klein quadrants (signed-safe); symmetric duplex tower; frame lock
        with THREE-period de Bruijn B(4,6) null (4^6 == 1 mod 3).
P/Q/Q.2/R/S  (alpha,M) transversality; V_4 equivariance and both commutator
        invariants; order-2 nt propagator T^(3k-1), M^(k-1)T0T1.
FIX-1..36 full 64x64 kernel; None-safe formatting; bilinear form; Psi sign
        guard; locus_tag collapse guard; cache auto-discovery; Candida
        Table 12; dual kernel mode; H5 cohort; mass invariant; df_eff;
        exact tower diagonal; enrichment gate; start/stop accounting;
        frame alignment; lag-list evaluation.
--------------------------------------------------------------------------------
POLICY
--------------------------------------------------------------------------------
* real data only (NCBI RefSeq CDS); ZERO RNG in every proof path
* strict CDS frame guard; internal-stop pseudogenes dropped
* one-gene-one-vote longest isoform; deterministic SHA-256 split with fallback
* all ladder orders scored on IDENTICAL target positions
* H5/H6/H7 and every tower null on the IDENTICAL gene cohort
* cluster-robust gene-level t inference (no bootstrap RNG)
* --mc-check runs a seeded shuffle ONLY to validate Theorem A
Author: Victor Logvinovich (ORCID 0009-0007-2040-6955)
Jean-Claude Perez (ORCID 000-0001-6446-2042)
================================================================================
"""
import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
import re
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

try:
    from scipy.stats import ttest_1samp
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

VERSION = "3.0"
TINY = 1e-300
SELFTEST_LONG_N = 200000
MASS_TOL = 1e-12
UA = "codon-memory-spectrum/3.0 (academic real-data analysis)"


# =====================================================================
# 0. Safe formatting
# =====================================================================
def fnum(x, nd=3, sign=False):
    if x is None:
        return "n/a"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "n/a"
    if not math.isfinite(v):
        return "n/a"
    return f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"


def _finite(x):
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


# =====================================================================
# 1. Genetic code, mirror involution, duplex fibers, CpG quartets,
#    screen boxes and exact screen degrees of freedom
# =====================================================================
BASES = "ACGT"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODON_IDX = {c: i for i, c in enumerate(CODONS)}
NC = 64
AA_STRING = "ACDEFGHIKLMNPQRSTVWY*"
AA_IDX = {ch: i for i, ch in enumerate(AA_STRING)}
NAA = 21
STOP = AA_IDX["*"]

NT_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(BASES):
    NT_LUT[ord(_b)] = _i

TABLE_1 = {
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
    1: {}, 11: {},
    2: {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},
    3: {"CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "ATA": "M", "TGA": "W"},
    4: {"TGA": "W"},
    5: {"ATA": "M", "TGA": "W", "AGA": "S", "AGG": "S"},
    6: {"TAA": "Q", "TAG": "Q"},
    12: {"CTG": "S"},
}
SUPPORTED_TABLES = set(TABLE_OVERRIDES)


def aa_array(table_id):
    if int(table_id) not in SUPPORTED_TABLES:
        raise ValueError(f"unsupported transl_table {table_id}")
    t = dict(TABLE_1)
    t.update(TABLE_OVERRIDES[int(table_id)])
    return np.array([AA_IDX[t[c]] for c in CODONS], dtype=np.int64)


def fibers_of(aa):
    return [np.flatnonzero(aa == a).astype(np.int64) for a in range(NAA)]


AA_TABLE1 = aa_array(1)

_COMP = str.maketrans("ACGT", "TGCA")
MIRROR = np.array([CODON_IDX[c.translate(_COMP)[::-1]] for c in CODONS], dtype=np.int64)
MIRROR_PAIR = np.minimum(np.arange(NC), MIRROR)
_uniq_pairs = np.unique(MIRROR_PAIR)
MIRROR32 = np.zeros(NC, dtype=np.int64)
for _i, _p in enumerate(_uniq_pairs):
    MIRROR32[MIRROR_PAIR == _p] = _i

WOBBLE4 = (np.arange(NC) % 4).astype(np.int64)
GC3BIN = ((WOBBLE4 == 1) | (WOBBLE4 == 2)).astype(np.int64)
IS_GC3 = GC3BIN.astype(bool)
HAS_CPG = np.array([("CG" in c) for c in CODONS], dtype=bool)

C_N1 = np.arange(NC) // 16
C_N2 = (np.arange(NC) // 4) % 4
C_N3 = np.arange(NC) % 4
SUF16 = 4 * C_N2 + C_N3          # [THM-S/T] suffix dinucleotide of the LEFT codon
OUT16 = 4 * C_N1 + C_N2          # [THM-T]   prefix dinucleotide of the RIGHT codon

# rho~ : sense suffix box s=(c2,c3)  ->  antisense prefix box (k(c3),k(c2))
RHO_BOX = np.array([4 * (3 - (s % 4)) + (3 - (s // 4)) for s in range(16)], dtype=np.int64)

_SUF_IDX = [np.flatnonzero(SUF16 == s) for s in range(16)]
_OUT_IDX = [np.flatnonzero(OUT16 == o) for o in range(16)]
_BOX_CACHE = {}


def screen_df_bounds(sense_mask=None):
    """[FIX-40] EXACT full-support df of the four screen layers.

    sense_mask=None  -> full 64-codon torus:  225 + 720 + 720 + 2304 = 3969 = 63^2
    sense_mask=(aa!=STOP), Table 1 -> 225 + 675 + 675 + 2025 = 3600 = 60^2.

    With R_s = #codons in suffix class s and C_o = #codons in prefix class o,
        a = sum_s (R_s - 1),  b = sum_o (C_o - 1),
        bridge = (n_s-1)(n_o-1), lead = a(n_o-1), wobble = (n_s-1)b, d13 = ab,
    and the total is always (q-1)^2 with q = #codons in the support.
    """
    if sense_mask is None:
        sm = np.ones(NC, dtype=bool)
    else:
        sm = np.asarray(sense_mask, dtype=bool)
    Rs = np.array([int(np.count_nonzero(sm[ix])) for ix in _SUF_IDX], dtype=np.int64)
    Co = np.array([int(np.count_nonzero(sm[ix])) for ix in _OUT_IDX], dtype=np.int64)
    nsa = int(np.count_nonzero(Rs > 0))
    noa = int(np.count_nonzero(Co > 0))
    a = int(np.sum(np.maximum(Rs - 1, 0)))
    b = int(np.sum(np.maximum(Co - 1, 0)))
    ls = max(nsa - 1, 0)
    rs = max(noa - 1, 0)
    d = {"bridge": ls * rs, "lead": a * rs, "wobble": ls * b, "d13": a * b}
    d["sum"] = int(d["bridge"] + d["lead"] + d["wobble"] + d["d13"])
    q = int(np.count_nonzero(sm))
    d["q_codons"] = q
    d["sum_left_df"] = a
    d["sum_right_df"] = b
    d["n_suffix_active"] = nsa
    d["n_prefix_active"] = noa
    d["delta_x1"] = int(d["lead"] + d["d13"])
    d["suffix16"] = int(d["bridge"] + d["wobble"])
    d["identity_ok"] = bool(d["sum"] == (q - 1) ** 2)
    d["kernel_A_df"] = int(nsa * rs)
    d["kernel_R_df"] = int(b)
    d["kernel_total_df"] = int(d["lead"] + d["wobble"] + d["d13"]
                               + d["kernel_A_df"] + d["kernel_R_df"])
    d["kernel_identity_ok"] = bool(d["kernel_total_df"] == q * (q - 1))
    return d


def box_partition(aa):
    """[THM-U] canonical 8+8 box partition of the standard code (table-aware).

    unsplit_prefix[o] : alpha(y) is independent of y3 inside prefix box o
                        -> the SENSE wobble choice is 100% synonymous.
    unsplit_suffix[s] : alpha(rho(x)) is independent of x1 inside suffix box s
                        -> the ANTISENSE wobble choice is 100% synonymous.
    Always unsplit_suffix[s] == unsplit_prefix[RHO_BOX[s]].
    """
    aa = np.asarray(aa, dtype=np.int64)
    key = aa.tobytes()
    hit = _BOX_CACHE.get(key)
    if hit is not None:
        return hit
    ar = aa[MIRROR]
    unp = np.array([np.unique(aa[ix]).size == 1 for ix in _OUT_IDX], dtype=bool)
    uns = np.array([np.unique(ar[ix]).size == 1 for ix in _SUF_IDX], dtype=bool)
    out = {
        "unsplit_prefix": unp,
        "unsplit_suffix": uns,
        "rho_box": RHO_BOX,
        "rho_consistent": bool(np.array_equal(uns, unp[RHO_BOX])),
        "prefix_codon_mask": unp[OUT16],
        "suffix_codon_mask": uns[SUF16],
        "B4_boxes": [BASES[o // 4] + BASES[o % 4] for o in range(16) if unp[o]],
        "Bsplit_boxes": [BASES[o // 4] + BASES[o % 4] for o in range(16) if not unp[o]],
        "S4_boxes": [BASES[s // 4] + BASES[s % 4] for s in range(16) if uns[s]],
        "Ssplit_boxes": [BASES[s // 4] + BASES[s % 4] for s in range(16) if not uns[s]],
        "n_unsplit_prefix": int(unp.sum()),
        "n_split_prefix": int((~unp).sum()),
        "n_unsplit_suffix": int(uns.sum()),
        "n_split_suffix": int((~uns).sum()),
        "n_codons_B4": int(np.count_nonzero(unp[OUT16])),
        "n_codons_Bsplit": int(np.count_nonzero(~unp[OUT16])),
    }
    _BOX_CACHE[key] = out
    return out


def duplex_layers(aa):
    ar = aa[MIRROR]
    key = aa.astype(np.int64) * NAA + ar.astype(np.int64)
    uk = np.unique(key)
    lay = np.searchsorted(uk, key).astype(np.int64)
    F = int(uk.size)
    sizes = np.array([int(np.sum(lay == a)) for a in range(F)], dtype=np.int64)
    per_codon = sizes[lay]
    return lay, F, per_codon, (per_codon > 1)


def cpg_quartet_partition(aa):
    """[THM-J] canonical rho-invariant split D8 = D4_CpG sqcup D4_Leu."""
    lay, F, sizes, d8 = duplex_layers(aa)
    cpg_mask = np.zeros(NC, dtype=bool)
    leu_mask = np.zeros(NC, dtype=bool)
    n_cpg, n_leu, n_mixed = 0, 0, 0
    doublets = []
    for a in range(F):
        fa = np.flatnonzero(lay == a)
        if fa.size < 2:
            continue
        ncpg = int(HAS_CPG[fa].sum())
        doublets.append({"codons": [CODONS[int(c)] for c in fa], "n_cpg": ncpg,
                         "channels": [AA_STRING[int(aa[fa[0]])],
                                      AA_STRING[int(aa[MIRROR[fa[0]]])]]})
        if ncpg == 1:
            cpg_mask[fa] = True
            n_cpg += 1
        elif ncpg == 0:
            leu_mask[fa] = True
            n_leu += 1
        else:
            n_mixed += 1
    return {"cpg_mask": cpg_mask, "leu_mask": leu_mask,
            "n_cpg_doublets": n_cpg, "n_leu_doublets": n_leu,
            "n_mixed_doublets": n_mixed,
            "cpg_quartet_rho_invariant": bool(np.array_equal(cpg_mask[MIRROR], cpg_mask)),
            "leu_quartet_rho_invariant": bool(np.array_equal(leu_mask[MIRROR], leu_mask)),
            "doublets": doublets}


ALPHABETS = {
    "codon64":  (64, None,      4),
    "aa20":     (21, "AA",      6),
    "mirror32": (32, MIRROR32,  4),
    "wobble4":  (4,  WOBBLE4,  12),
    "gc3bin2":  (2,  GC3BIN,   20),
}


def encode(codons, alphabet, aa_ch):
    if alphabet == "codon64":
        return codons.astype(np.int64)
    if alphabet == "aa20":
        return aa_ch[codons].astype(np.int64)
    return ALPHABETS[alphabet][1][codons].astype(np.int64)


# =====================================================================
# 2. Parsing + robust family split
# =====================================================================
def extract_gene_id(hdr):
    for key in ("gene", "locus_tag", "protein_id"):
        m = re.search(rf"\[{key}=([^\]]+)\]", hdr)
        if m:
            return m.group(1)
    return hdr.lstrip(">").split()[0]


_ROOT_PAT = re.compile(r"(-R[A-Za-z]\d*|\.\d+|_[0-9]{1,2}|\|.*)$")


def family_root(gid):
    s, prev = str(gid), None
    while s != prev:
        prev = s
        s = _ROOT_PAT.sub("", s)
    return s or str(gid)


def _hash_train(key, pct):
    h = hashlib.sha256(str(key).encode("utf-8")).hexdigest()
    return (int(h[:8], 16) % 100) < pct


def assign_split(gids, pct=80, min_frac=0.10, collapse_ratio=0.5, min_roots=20):
    gids = [str(g) for g in gids]
    n = len(gids)
    roots = [family_root(g) for g in gids]
    n_roots = len(set(roots))
    mode, keys = "family_root", roots
    if n_roots < max(min_roots, collapse_ratio * n):
        mode, keys = "gid_fallback_root_collapse", gids
    mask = np.array([_hash_train(k, pct) for k in keys], dtype=bool)
    n_tr, n_te = int(mask.sum()), int((~mask).sum())
    if n and min(n_tr, n_te) < min_frac * n:
        mode, keys = "gid_fallback_split_imbalance", gids
        mask = np.array([_hash_train(k, pct) for k in keys], dtype=bool)
        n_tr, n_te = int(mask.sum()), int((~mask).sum())
    return mask, {"mode": mode, "n_genes": n, "n_family_roots": n_roots,
                  "n_train": n_tr, "n_test": n_te,
                  "train_fraction": (n_tr / n) if n else None}


def parse_cds(path, max_genes=0, trim=True, default_table=1, min_codons=8, gene_order="hash"):
    """[FIX-28/29] record-own transl_table; kept counters summed AFTER the slice."""
    raw, tables = [], {}
    stats = dict(n_records=0, rej_nonacgt=0, rej_frame=0, rej_short=0,
                 rej_mixed_table=0, rej_internal_stop=0, kept=0,
                 n_codons_kept=0, max_gene_codons=0, median_gene_codons=0,
                 n_atg_start_raw=0, n_stop_end_raw=0,
                 n_atg_start_kept=0, n_stop_end_kept=0)
    hdr, buf = None, []
    aa_cache = {t: aa_array(t) for t in SUPPORTED_TABLES}

    def flush(h, parts):
        if not h or not parts:
            return
        s = "".join(parts)
        stats["n_records"] += 1
        arr = np.frombuffer(s.encode("ascii", "ignore"), dtype=np.uint8)
        nt = NT_LUT[arr]
        if np.any(nt >= 4):
            stats["rej_nonacgt"] += 1
            return
        if nt.size % 3 != 0:
            stats["rej_frame"] += 1
            return
        cd = (16 * nt[0::3] + 4 * nt[1::3] + nt[2::3]).astype(np.int64)
        mt = re.search(r"\[transl_table=(\d+)\]", h)
        tb = int(mt.group(1)) if mt else default_table
        if tb not in SUPPORTED_TABLES:
            tb = default_table
        aa_rec = aa_cache[tb]
        atg = bool(cd.size >= 2 and int(cd[0]) == CODON_IDX["ATG"])
        stp = bool(cd.size >= 2 and int(aa_rec[cd[-1]]) == STOP)
        stats["n_atg_start_raw"] += int(atg)
        stats["n_stop_end_raw"] += int(stp)
        if trim and cd.size >= 5:
            cd = cd[1:-1]
        if cd.size < min_codons:
            stats["rej_short"] += 1
            return
        tables[tb] = tables.get(tb, 0) + 1
        raw.append((extract_gene_id(h), tb, cd, atg, stp))

    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush(hdr, buf)
                hdr, buf = line, []
            else:
                buf.append(line.upper())
        flush(hdr, buf)

    if not raw:
        return [], default_table, stats
    dom = max(tables.items(), key=lambda kv: kv[1])[0]
    stop_ids = np.flatnonzero(aa_cache[dom] == STOP)
    best = {}
    for gid, tb, cd, atg, stp in raw:
        if tb != dom:
            stats["rej_mixed_table"] += 1
            continue
        if gid not in best or cd.size > best[gid][2].size:
            best[gid] = (gid, tb, cd, atg, stp)
    recs4 = []
    for gid, tb, cd, atg, stp in best.values():
        region = cd if trim else (cd[1:-1] if cd.size >= 2 else cd)
        if stop_ids.size and np.isin(region, stop_ids).any():
            stats["rej_internal_stop"] += 1
            continue
        recs4.append((gid, cd, atg, stp))
    if max_genes and max_genes > 0 and len(recs4) > max_genes:
        if gene_order == "hash":
            recs4.sort(key=lambda r: (hashlib.sha256(r[0].encode("utf-8")).hexdigest(), r[0]))
            recs4 = recs4[:max_genes]
        else:
            recs4.sort(key=lambda r: r[0])
            recs4 = recs4[:max_genes]
    recs4.sort(key=lambda r: r[0])
    stats["n_atg_start_kept"] = int(sum(int(r[2]) for r in recs4))
    stats["n_stop_end_kept"] = int(sum(int(r[3]) for r in recs4))
    recs = [(r[0], r[1]) for r in recs4]
    stats["kept"] = len(recs)
    if recs:
        ln = np.array([cd.size for _, cd in recs])
        stats["n_codons_kept"] = int(ln.sum())
        stats["max_gene_codons"] = int(ln.max())
        stats["median_gene_codons"] = int(np.median(ln))
    return recs, dom, stats


# =====================================================================
# 3. Curated 100-Species Biosphere Catalog, Cache Discovery & NCBI Resolver
# =====================================================================
# Format: species_key -> (known_gcf_or_None, table_id, clade, organism_name)
# First 54 cellular species match v2.4 exactly (100% cache compatible).
# Next 46 species expand to 100 cellular taxa across 12 clades & 5 NCBI tables (1,4,6,11,12).
BIOSPHERE_100 = {
    # --- Core 15 (v2.4 PANEL) ---
    "Human_T2T":             ("GCF_009914755", 1,  "mammal",     "Homo sapiens"),
    "M_musculus":            ("GCF_000001635", 1,  "mammal",     "Mus musculus"),
    "D_melanogaster":        ("GCF_000001215", 1,  "insect",     "Drosophila melanogaster"),
    "A_thaliana":            ("GCF_000001735", 1,  "plant",      "Arabidopsis thaliana"),
    "S_cerevisiae":          ("GCF_000146045", 1,  "fungus",     "Saccharomyces cerevisiae"),
    "E_coli":                ("GCF_000005845", 11, "bacterium",  "Escherichia coli"),
    "B_subtilis":            ("GCF_000009045", 11, "bacterium",  "Bacillus subtilis"),
    "M_jannaschii":          ("GCF_000007745", 11, "archaeon",   "Methanocaldococcus jannaschii"),
    "H_salinarum":           ("GCF_000015585", 11, "archaeon",   "Halobacterium salinarum"),
    "C_elegans":             ("GCF_000002985", 1,  "nematode",   "Caenorhabditis elegans"),
    "Zebrafish":             ("GCF_000002035", 1,  "fish",       "Danio rerio"),
    "Chicken":               ("GCF_016699485", 1,  "bird",       "Gallus gallus"),
    "Thermus":               ("GCF_000008125", 11, "bacterium",  "Thermus thermophilus"),
    "Dictyostelium":         ("GCF_000004695", 1,  "protist",    "Dictyostelium discoideum"),
    "P_aeruginosa":          ("GCF_000006765", 11, "bacterium",  "Pseudomonas aeruginosa"),
    # --- Remaining 39 cellular species from v2.4 EXTENDED (total = 54) ---
    "Chimp":                 ("GCF_002880755", 1,  "mammal",     "Pan troglodytes"),
    "Rattus_norvegicus":     ("GCF_015227675", 1,  "mammal",     "Rattus norvegicus"),
    "Macaca_mulatta":        ("GCF_003339765", 1,  "mammal",     "Macaca mulatta"),
    "Gorilla_gorilla":       ("GCF_000151905", 1,  "mammal",     "Gorilla gorilla"),
    "Bos_taurus":            ("GCF_002263795", 1,  "mammal",     "Bos taurus"),
    "Sus_scrofa":            ("GCF_000003025", 1,  "mammal",     "Sus scrofa"),
    "Canis_familiaris":      ("GCF_000002285", 1,  "mammal",     "Canis lupus familiaris"),
    "Felis_catus":           ("GCF_000181335", 1,  "mammal",     "Felis catus"),
    "Equus_caballus":        ("GCF_002863925", 1,  "mammal",     "Equus caballus"),
    "Oryctolagus_cunic":     ("GCF_000003625", 1,  "mammal",     "Oryctolagus cuniculus"),
    "Xenopus_tropicalis":    ("GCF_000004195", 1,  "amphibian",  "Xenopus tropicalis"),
    "Takifugu_rubripes":     ("GCF_901000725", 1,  "fish",       "Takifugu rubripes"),
    "Anopheles_gambiae":     ("GCF_000005575", 1,  "insect",     "Anopheles gambiae"),
    "Apis_mellifera":        ("GCF_003254395", 1,  "insect",     "Apis mellifera"),
    "Bombyx_mori":           ("GCF_014905235", 1,  "insect",     "Bombyx mori"),
    "C_briggsae":            ("GCF_000004555", 1,  "nematode",   "Caenorhabditis briggsae"),
    "Oryza_sativa":          ("GCF_001433935", 1,  "plant",      "Oryza sativa"),
    "Zea_mays":              ("GCF_000005005", 1,  "plant",      "Zea mays"),
    "Solanum_lycopers":      ("GCF_000188115", 1,  "plant",      "Solanum lycopersicum"),
    "Glycine_max":           ("GCF_000004515", 1,  "plant",      "Glycine max"),
    "Populus_trichocarpa":   ("GCF_000002775", 1,  "plant",      "Populus trichocarpa"),
    "S_pombe":               ("GCF_000002945", 1,  "fungus",     "Schizosaccharomyces pombe"),
    "Candida_albicans":      ("GCF_000182965", 12, "fungus",     "Candida albicans"),
    "Aspergillus_nidulans":  ("GCF_000149205", 1,  "fungus",     "Aspergillus nidulans"),
    "Neurospora_crassa":     ("GCF_000182925", 1,  "fungus",     "Neurospora crassa"),
    "Plasmodium_falcip":     ("GCF_000002765", 1,  "protist",    "Plasmodium falciparum"),
    "Trypanosoma_brucei":    ("GCF_000002445", 1,  "protist",    "Trypanosoma brucei"),
    "Leishmania_major":      ("GCF_000002725", 1,  "protist",    "Leishmania major"),
    "Toxoplasma_gondii":     ("GCF_000006565", 1,  "protist",    "Toxoplasma gondii"),
    "Chlamydomonas_rein":    ("GCF_000002595", 1,  "plant",      "Chlamydomonas reinhardtii"),
    "Staphylococcus_aureus": ("GCF_000013425", 11, "bacterium",  "Staphylococcus aureus"),
    "Streptococcus_pyog":    ("GCF_000006785", 11, "bacterium",  "Streptococcus pyogenes"),
    "M_tuberculosis":        ("GCF_000195955", 11, "bacterium",  "Mycobacterium tuberculosis"),
    "Salmonella_enterica":   ("GCF_000006945", 11, "bacterium",  "Salmonella enterica"),
    "Lactococcus_lactis":    ("GCF_000006865", 11, "bacterium",  "Lactococcus lactis"),
    "Helicobacter_pylori":   ("GCF_000008525", 11, "bacterium",  "Helicobacter pylori"),
    "Bacteroides_fragilis":  ("GCF_000025985", 11, "bacterium",  "Bacteroides fragilis"),
    "C_difficile":           ("GCF_000009205", 11, "bacterium",  "Clostridioides difficile"),
    "Synechocystis_PCC6803": ("GCF_000009725", 11, "bacterium",  "Synechocystis sp. PCC 6803"),
    # --- 46 New Biosphere-100 Cellular Species (auto-resolved via NCBI RefSeq E-utilities) ---
    # Mammals (+8 -> 20 total; adds marsupial + monotreme!)
    "Pan_paniscus":          (None, 1,  "mammal",     "Pan paniscus"),
    "Pongo_abelii":          (None, 1,  "mammal",     "Pongo abelii"),
    "Cavia_porcellus":       (None, 1,  "mammal",     "Cavia porcellus"),
    "Mesocricetus_auratus":  (None, 1,  "mammal",     "Mesocricetus auratus"),
    "Ovis_aries":            (None, 1,  "mammal",     "Ovis aries"),
    "Capra_hircus":          (None, 1,  "mammal",     "Capra hircus"),
    "Monodelphis_domestica": (None, 1,  "mammal",     "Monodelphis domestica"),
    "Ornithorhynchus_anat":  (None, 1,  "mammal",     "Ornithorhynchus anatinus"),
    # Birds, Reptiles, Fishes (+8 -> 12 non-mammal vertebrates)
    "Taeniopygia_guttata":   (None, 1,  "bird",       "Taeniopygia guttata"),
    "Anas_platyrhynchos":    (None, 1,  "bird",       "Anas platyrhynchos"),
    "Anolis_carolinensis":   (None, 1,  "reptile",    "Anolis carolinensis"),
    "Chrysemys_picta":       (None, 1,  "reptile",    "Chrysemys picta bellii"),
    "Salmo_salar":           (None, 1,  "fish",       "Salmo salar"),
    "Oryzias_latipes":       (None, 1,  "fish",       "Oryzias latipes"),
    "Latimeria_chalumnae":   (None, 1,  "fish",       "Latimeria chalumnae"),
    "Callorhinchus_milii":   (None, 1,  "fish",       "Callorhinchus milii"),
    # Invertebrates & Chordate Bridge (+6 -> 12 total)
    "Aedes_aegypti":         (None, 1,  "insect",     "Aedes aegypti"),
    "Tribolium_castaneum":   (None, 1,  "insect",     "Tribolium castaneum"),
    "Daphnia_pulex":         (None, 1,  "invert",     "Daphnia pulex"),
    "Octopus_bimaculoides":  (None, 1,  "invert",     "Octopus bimaculoides"),
    "Strongylocentrotus":    (None, 1,  "invert",     "Strongylocentrotus purpuratus"),
    "Ciona_intestinalis":    (None, 1,  "invert",     "Ciona intestinalis"),
    # Plants & Bryophytes (+5 -> 12 total)
    "Brachypodium_distach":  (None, 1,  "plant",      "Brachypodium distachyon"),
    "Sorghum_bicolor":       (None, 1,  "plant",      "Sorghum bicolor"),
    "Vitis_vinifera":        (None, 1,  "plant",      "Vitis vinifera"),
    "Brassica_rapa":         (None, 1,  "plant",      "Brassica rapa"),
    "Physcomitrium_patens":  (None, 1,  "plant",      "Physcomitrium patens"),
    # Fungi (+5 -> 10 total; includes 2 MORE Table-12 yeasts!)
    "Debaryomyces_hansenii": (None, 12, "fungus",     "Debaryomyces hansenii"),
    "Scheffersomyces_stip":  (None, 12, "fungus",     "Scheffersomyces stipitis"),
    "Yarrowia_lipolytica":   (None, 1,  "fungus",     "Yarrowia lipolytica"),
    "Cryptococcus_neoform":  (None, 1,  "fungus",     "Cryptococcus neoformans"),
    "Ustilago_maydis":       (None, 1,  "fungus",     "Ustilago maydis"),
    # Protists & Ciliates (+4 -> 9 total; includes 2 Table-6 ciliates with duplex QUARTETS!)
    "Tetrahymena_thermoph":  (None, 6,  "protist",    "Tetrahymena thermophila"),
    "Paramecium_tetraur":    (None, 6,  "protist",    "Paramecium tetraurelia"),
    "Giardia_lamblia":       (None, 1,  "protist",    "Giardia lamblia"),
    "Thalassiosira_pseudo":  (None, 1,  "protist",    "Thalassiosira pseudonana"),
    # Bacteria & Table-4 Mycoplasmas (+6 -> 19 total; includes 3 Table-4 genomes!)
    "Mycoplasma_mycoides":   (None, 4,  "bacterium",  "Mycoplasma mycoides"),
    "Mycoplasma_pneumoniae": (None, 4,  "bacterium",  "Mycoplasma pneumoniae"),
    "Mycoplasma_genitalium": (None, 4,  "bacterium",  "Mycoplasma genitalium"),
    "Streptomyces_coelic":   (None, 11, "bacterium",  "Streptomyces coelicolor"),
    "Deinococcus_radiodur":  (None, 11, "bacterium",  "Deinococcus radiodurans"),
    "Vibrio_cholerae":       (None, 11, "bacterium",  "Vibrio cholerae"),
    # Archaea (+4 -> 6 total; adds 100C/122C extremophiles & Crenarchaeota)
    "Pyrococcus_furiosus":   (None, 11, "archaeon",   "Pyrococcus furiosus"),
    "Methanopyrus_kandleri": (None, 11, "archaeon",   "Methanopyrus kandleri"),
    "Saccharolobus_solfat":  (None, 11, "archaeon",   "Saccharolobus solfataricus"),
    "Methanosarcina_acet":   (None, 11, "archaeon",   "Methanosarcina acetivorans"),
}

PANEL = {nm: (acc, tbl) for nm, (acc, tbl, _cl, _org) in list(BIOSPHERE_100.items())[:15]}
EXTENDED = {acc: (nm, tbl) for nm, (acc, tbl, _cl, _org) in list(BIOSPHERE_100.items())[15:] if acc}
EXTENDED.update({
    "GCF_000008895": ("Human_mito", 2),
    "GCF_009858895": ("SARS_CoV_2", 1),
    "GCF_000002025": ("SARS_CoV_2_alt", 1),
})
ACC2NAME = {acc: nm for nm, (acc, _t, _c, _o) in BIOSPHERE_100.items() if acc}
ACC2TABLE = {acc: t for _nm, (acc, t, _c, _o) in BIOSPHERE_100.items() if acc}
ACC2CLADE = {acc: cl for _nm, (acc, _t, cl, _o) in BIOSPHERE_100.items() if acc}
NAME2CLADE = {nm: cl for nm, (_a, _t, cl, _o) in BIOSPHERE_100.items()}
NAME2TABLE = {nm: t for nm, (_a, t, _c, _o) in BIOSPHERE_100.items()}
for _acc, (_nm, _t) in EXTENDED.items():
    ACC2NAME.setdefault(_acc, _nm)
    ACC2TABLE.setdefault(_acc, _t)
    NAME2TABLE.setdefault(_nm, _t)


def load_manifest(cds_dir):
    mf = Path(cds_dir) / "manifest_v30.json"
    if mf.exists():
        try:
            data = json.loads(mf.read_text(encoding="utf-8"))
            for acc, rec in data.items():
                nm = rec.get("species")
                if nm:
                    ACC2NAME[acc] = nm
                    if rec.get("table") is not None:
                        ACC2TABLE[acc] = int(rec["table"])
                        NAME2TABLE.setdefault(nm, int(rec["table"]))
                    if rec.get("clade"):
                        ACC2CLADE[acc] = rec["clade"]
                        NAME2CLADE.setdefault(nm, rec["clade"])
            return data
        except Exception:
            return {}
    return {}


def save_manifest_entry(cds_dir, acc, species, table, clade, filename):
    mf = Path(cds_dir) / "manifest_v30.json"
    data = {}
    if mf.exists():
        try:
            data = json.loads(mf.read_text(encoding="utf-8"))
        except Exception:
            data = {}
    data[acc] = {"species": species, "table": int(table), "clade": clade, "file": str(filename)}
    ACC2NAME[acc] = species
    ACC2TABLE[acc] = int(table)
    ACC2CLADE[acc] = clade
    NAME2CLADE[species] = clade
    NAME2TABLE[species] = int(table)
    mf.parent.mkdir(parents=True, exist_ok=True)
    tmp = mf.with_name(mf.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2), encoding="utf-8")
    os.replace(tmp, mf)


def discover_cache(cds_dir):
    load_manifest(cds_dir)
    roots = [Path(cds_dir), Path("../genom_14/cds_cache"), Path("./cds_cache")]
    try:
        roots += sorted(Path("/Volumes").glob("*/genom_14/cds_cache"))
        roots += sorted(Path("/Volumes").glob("*/GitHub/cds_cache"))
    except Exception:
        pass
    out, seen = {}, set()
    for r in roots:
        try:
            rr = r.resolve()
        except Exception:
            continue
        if rr in seen or not rr.is_dir():
            continue
        seen.add(rr)
        load_manifest(rr)
        for p in sorted(rr.glob("*_cds_from_genomic.fna.gz")):
            if p.stat().st_size < 1024:
                continue
            m = re.search(r"(GCF_\d+)", p.name)
            if m:
                acc = m.group(1)
                if acc not in out or p.stat().st_size > out[acc].stat().st_size:
                    out[acc] = p
    return out


def ftp_base_for_accession(acc):
    tail = acc.split("_", 1)[1] if "_" in acc else acc
    parts = [tail[i:i + 3] for i in range(0, len(tail), 3)]
    return "https" + "://" + "ftp.ncbi.nlm.nih.gov" + "/genomes/all/GCF/" + "/".join(parts) + "/"


def _download_cds_url(url, dest, timeout=180):
    if dest.exists() and dest.stat().st_size > 1024:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    if not data:
        return None
    tmp = dest.with_suffix(".part")
    tmp.write_bytes(data)
    with gzip.open(tmp, "rb") as g:
        if not g.read(1).startswith(b">"):
            raise ValueError("not FASTA")
    os.replace(tmp, dest)
    return dest


def auto_download(acc, dest_dir, timeout=180):
    dest_dir.mkdir(parents=True, exist_ok=True)
    base = ftp_base_for_accession(acc)
    try:
        req = urllib.request.Request(base, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            html = r.read().decode("utf-8", "ignore")
    except Exception as e:
        print(f"    [auto-download] {acc}: listing failed ({e})")
        return None
    dirs = sorted({l.rstrip("/").split("/")[-1] for l in re.findall(r'href="([^"]+)"', html)
                   if l.rstrip("/").split("/")[-1].startswith(acc)}, reverse=True)
    for d in dirs:
        fn = f"{d}_cds_from_genomic.fna.gz"
        dest = dest_dir / fn
        try:
            got = _download_cds_url(base + d + "/" + fn, dest, timeout=timeout)
            if got is not None:
                print(f"    [auto-download] {acc}: {fn} ({got.stat().st_size/1e6:.1f} MB)")
                return got
        except Exception as e:
            print(f"    [auto-download] {acc}/{d}: {e}")
    return None


def auto_download_by_organism(species_key, organism_name, table_id, clade, dest_dir, timeout=180):
    """Resolve official latest RefSeq CDS archive via NCBI E-utilities esearch -> esummary."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    eutils = "https" + "://" + "eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    alt_orgs = [organism_name]
    if organism_name == "Mycoplasma pneumoniae":
        alt_orgs.append("Mycoplasmoides pneumoniae")
    queries = []
    for org_q in alt_orgs:
        queries += [
            f'"{org_q}"[Organism] AND "latest refseq"[filter] AND "reference genome"[filter]',
            f'"{org_q}"[Organism] AND "latest refseq"[filter] AND "representative genome"[filter]',
            f'"{org_q}"[Organism] AND "latest refseq"[filter]',
        ]
    ids = []
    for q in queries:
        try:
            url = eutils + "esearch.fcgi?db=assembly&retmode=json&retmax=5&term=" + urllib.parse.quote(q)
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=min(timeout, 60)) as r:
                js = json.loads(r.read().decode("utf-8", "ignore"))
            ids = js.get("esearchresult", {}).get("idlist", [])
            time.sleep(0.35)
            if ids:
                break
        except Exception as e:
            print(f"    [eutils-search] {species_key}: {e}")
            time.sleep(0.5)
    if not ids:
        print(f"    [eutils] {species_key} ({organism_name}): no RefSeq assembly found")
        return None
    try:
        url = eutils + "esummary.fcgi?db=assembly&retmode=json&id=" + ",".join(ids[:5])
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=min(timeout, 60)) as r:
            sjs = json.loads(r.read().decode("utf-8", "ignore"))
        time.sleep(0.35)
    except Exception as e:
        print(f"    [eutils-summary] {species_key}: {e}")
        return None
    res_map = sjs.get("result", {})
    for aid in ids[:5]:
        rec = res_map.get(str(aid), {})
        ftp = rec.get("ftppath_refseq", "")
        if not ftp:
            continue
        https_dir = ftp.replace("ftp://", "https://").rstrip("/")
        dname = https_dir.split("/")[-1]
        m = re.search(r"(GCF_\d+)", dname)
        acc = m.group(1) if m else dname
        fn = f"{dname}_cds_from_genomic.fna.gz"
        dest = dest_dir / fn
        try:
            got = _download_cds_url(https_dir + "/" + fn, dest, timeout=timeout)
            if got is not None:
                save_manifest_entry(dest_dir, acc, species_key, table_id, clade, fn)
                print(f"    [eutils-download] {species_key} -> {fn} ({got.stat().st_size/1e6:.1f} MB)")
                return got
        except Exception as e:
            print(f"    [eutils-download] {species_key}/{dname}: {e}")
    return None


# =====================================================================
# 4. Context models + dual kernel
# =====================================================================
def ctx_at(s, positions, order, q):
    codes = np.zeros(positions.size, dtype=np.int64)
    for j in range(order):
        codes += s[positions - order + j].astype(np.int64) * (q ** (order - 1 - j))
    return codes


class OrderModel:
    __slots__ = ("order", "q", "uniq", "logP", "p_marg", "logp_marg",
                 "n_ctx", "n_train_sym")

    def __init__(self, order, q):
        self.order, self.q = order, q
        self.uniq = np.empty(0, dtype=np.int64)
        self.logP = np.empty((0, q), dtype=np.float64)
        self.p_marg = np.ones(q) / q
        self.logp_marg = np.log(self.p_marg)
        self.n_ctx = 0
        self.n_train_sym = 0


def _lookup_rows(models, order, codes, q):
    out = np.empty((codes.size, q), dtype=np.float64)
    todo = np.arange(codes.size)
    for o in range(order, -1, -1):
        if todo.size == 0:
            break
        if o == 0:
            out[todo] = models[0].p_marg[None, :]
            break
        m = models[o]
        if m.uniq.size == 0:
            continue
        c = codes[todo] % (q ** o)
        loc = np.searchsorted(m.uniq, c)
        ok = m.uniq[np.minimum(loc, m.uniq.size - 1)] == c
        if ok.any():
            out[todo[ok]] = np.exp(m.logP[loc[ok]])
        todo = todo[~ok]
    return out


def build_models(seqs, q, max_order, kappa=8.0, alpha=0.5,
                 min_ctx_count=3, max_ctx=200000):
    models = [OrderModel(p, q) for p in range(max_order + 1)]
    cnt = np.zeros(q, dtype=np.float64)
    for s in seqs:
        cnt += np.bincount(s, minlength=q)
    p0 = (cnt + alpha) / (cnt.sum() + alpha * q)
    models[0].p_marg = p0
    models[0].logp_marg = np.log(np.maximum(p0, TINY))
    models[0].n_train_sym = int(cnt.sum())
    for p in range(1, max_order + 1):
        cl, tl = [], []
        for s in seqs:
            if s.size <= p:
                continue
            cl.append(ctx_at(s, np.arange(p, s.size), p, q))
            tl.append(s[p:])
        if not cl:
            continue
        ctx = np.concatenate(cl)
        tgt = np.concatenate(tl)
        del cl, tl
        uniq, inv = np.unique(ctx, return_inverse=True)
        tot = np.bincount(inv, minlength=uniq.size)
        keep = np.flatnonzero(tot >= min_ctx_count)
        if keep.size > max_ctx:
            keep = np.sort(keep[np.argsort(-tot[keep])[:max_ctx]])
        if keep.size == 0:
            del ctx, tgt, uniq, inv, tot
            continue
        remap = np.full(uniq.size, -1, dtype=np.int64)
        remap[keep] = np.arange(keep.size)
        inv_sel = remap[inv]
        sel = inv_sel >= 0
        counts = np.bincount(inv_sel[sel] * q + tgt[sel],
                             minlength=keep.size * q).reshape(keep.size, q).astype(np.float64)
        u = uniq[keep]
        models[p].n_train_sym = int(sel.sum())
        del ctx, tgt, uniq, inv, tot, inv_sel, sel, remap
        Ppar = _lookup_rows(models, p - 1, u % (q ** (p - 1)), q)
        counts += kappa * Ppar
        counts /= counts.sum(axis=1, keepdims=True)
        del Ppar
        np.log(np.maximum(counts, TINY), out=counts)
        models[p].uniq = u
        models[p].logP = counts
        models[p].n_ctx = int(u.size)
    return models


def full_first_order_logp(models, q):
    P = np.tile(models[0].p_marg.astype(np.float64), (q, 1))
    m = models[1]
    if m.uniq.size:
        P[m.uniq] = np.exp(m.logP)
    P = np.maximum(P, TINY)
    P /= P.sum(axis=1, keepdims=True)
    return np.log(P)


def bigram_counts(seqs, q):
    B = np.zeros((q, q), dtype=np.float64)
    for s in seqs:
        if s.size < 2:
            continue
        B += np.bincount(s[:-1] * q + s[1:], minlength=q * q).reshape(q, q)
    return B


def laplace_logp(B, alpha=0.1):
    q = B.shape[0]
    P = (B + alpha) / (B.sum(axis=1, keepdims=True) + q * alpha)
    return np.log(np.maximum(P, TINY))


def build_kernel(train_seqs, q, mode="backoff", kappa=8.0, laplace_alpha=0.1):
    if mode == "laplace":
        return laplace_logp(bigram_counts(train_seqs, q), laplace_alpha)
    m = build_models(train_seqs, q, 1, kappa=kappa, min_ctx_count=1)
    return full_first_order_logp(m, q)


def seq_logprob(models, s, order, eval_start, q, return_order=False):
    N = s.size - eval_start
    if N <= 0:
        return (np.empty(0), np.empty(0, dtype=np.int64)) if return_order else np.empty(0)
    out = np.empty(N, dtype=np.float64)
    used = np.full(N, -1, dtype=np.int64)
    todo = np.arange(N)
    tgt = s[eval_start:]
    for o in range(order, -1, -1):
        if todo.size == 0:
            break
        if o == 0:
            out[todo] = models[0].logp_marg[tgt[todo]]
            used[todo] = 0
            todo = np.empty(0, dtype=np.int64)
            break
        m = models[o]
        if m.uniq.size == 0:
            continue
        codes = ctx_at(s, eval_start + todo, o, q)
        loc = np.searchsorted(m.uniq, codes)
        ok = m.uniq[np.minimum(loc, m.uniq.size - 1)] == codes
        if ok.any():
            hit = todo[ok]
            out[hit] = m.logP[loc[ok], tgt[hit]]
            used[hit] = o
        todo = todo[~ok]
    if todo.size:
        out[todo] = models[0].logp_marg[tgt[todo]]
        used[todo] = 0
    return (out, used) if return_order else out


# =====================================================================
# 5. THEOREM A / A.2 / A.3 / A.4 / E / E.2
# =====================================================================
def junction_count(n, L):
    if L < 1 or n // L < 2:
        return 0
    return (n + L - 1) // L - 1


def _gene_boundary_terms(s, logP, L):
    n = s.size
    m = n // L
    if m < 2:
        return None
    q = logP.shape[0]
    tail = n - m * L
    ends = s[L - 1:m * L:L]
    starts = s[0:m * L:L]
    D = float(np.sum(logP[s[L - 1:(m - 1) * L:L], s[L:m * L:L]]))
    if tail > 0:
        D += float(logP[s[m * L - 1], s[m * L]])
    ce = np.bincount(ends, minlength=q).astype(np.float64)
    cs = np.bincount(starts, minlength=q).astype(np.float64)
    N = (float(ce @ logP @ cs) - float(np.sum(logP[ends, starts]))) / m
    if tail > 0:
        N += float(ce @ logP[:, s[m * L]]) / m
    J = (m - 1) + (1 if tail > 0 else 0)
    Ncomp = J * float((ce / m) @ logP @ (cs / m))
    pg = np.bincount(s, minlength=q).astype(np.float64) / float(n)
    Niso = J * float(pg @ logP @ pg)
    return J, D, N, Ncomp, Niso


def analytic_h4(seqs, logP, L, collect_terms=False):
    if L < 1:
        return None
    real_sum, shuf_sum, ntrans, jtot = 0.0, 0.0, 0, 0
    Dsum = Nsum = NcompSum = NisoSum = 0.0
    n_genes_used = 0
    for s in seqs:
        n = s.size
        if n < 2:
            continue
        ll = float(np.sum(logP[s[:-1], s[1:]]))
        real_sum += ll
        ntrans += n - 1
        t = _gene_boundary_terms(s, logP, L)
        if t is None:
            shuf_sum += ll
            continue
        J, D, N, Nc, Ni = t
        n_genes_used += 1
        jtot += J
        Dsum += D
        Nsum += N
        NcompSum += Nc
        NisoSum += Ni
        shuf_sum += ll - D + N
    if ntrans == 0:
        return None
    out = {"L": int(L), "ll_real": real_sum / ntrans,
           "ll_shuf_expected": shuf_sum / ntrans,
           "delta_H4_analytic": (real_sum - shuf_sum) / ntrans,
           "n_transitions": int(ntrans), "n_genes_shuffled": int(n_genes_used),
           "junctions_total": int(jtot), "gamma_L": jtot / ntrans}
    if collect_terms:
        out["_Dsum"], out["_Nsum"] = Dsum, Nsum
        out["_NcompSum"], out["_NisoSum"] = NcompSum, NisoSum
    return out


def theorem_E(seqs, logP, L, psi_dict, h4rec):
    if psi_dict is None or h4rec is None or h4rec["junctions_total"] == 0:
        return None
    q = logP.shape[0]
    B = bigram_counts(seqs, q)
    N = B.sum()
    pL, pR = B.sum(axis=1) / N, B.sum(axis=0) / N
    mean_real = float(np.sum((B / N) * logP))
    mean_indep = float(np.sum((pL[:, None] * pR[None, :]) * logP))
    J = h4rec["junctions_total"]
    Dbar = h4rec["_Dsum"] / J
    Nbar = h4rec["_Nsum"] / J
    Ncomp = h4rec["_NcompSum"] / J
    Niso = h4rec["_NisoSum"] / J
    d_pos = Dbar - mean_real
    d_within = Ncomp - mean_indep
    d_isochore = Niso - mean_indep
    d_endpoint = Ncomp - Niso
    d_finite = Nbar - Ncomp
    psi = psi_dict["Psi_LR"]
    gamma = h4rec["gamma_L"]
    recon = gamma * (psi + d_pos - d_within - d_finite)
    return {"L": int(L), "gamma_L": gamma, "Psi_LR": psi,
            "Dbar": Dbar, "Nbar": Nbar, "Ncomp": Ncomp, "Niso": Niso,
            "d_pos": d_pos, "d_within": d_within,
            "d_isochore": d_isochore, "d_endpoint": d_endpoint, "d_finite": d_finite,
            "within_split_abs_err": abs(d_within - (d_isochore + d_endpoint)),
            "delta_reconstructed": recon, "delta_direct": h4rec["delta_H4_analytic"],
            "identity_abs_err": abs(recon - h4rec["delta_H4_analytic"]),
            "R_junc_from_terms": (1.0 + (d_pos - d_within - d_finite) / psi) if psi != 0 else None}


def cohort_decomposition(seqs, logP, L):
    all_ll, all_n, keep_ll, keep_n, sh_ll, sh_n = 0.0, 0, 0.0, 0, 0.0, 0
    for s in seqs:
        n = s.size
        if n < 2:
            continue
        ll = float(np.sum(logP[s[:-1], s[1:]]))
        all_ll += ll
        all_n += n - 1
        t = _gene_boundary_terms(s, logP, L)
        if t is None:
            continue
        J, D, N, _, _ = t
        keep_ll += ll
        keep_n += n - 1
        sh_ll += ll - D + N
        sh_n += n - 1
    if all_n == 0 or keep_n == 0 or sh_n == 0:
        return None
    lr_all, lr_keep, ls_keep = all_ll / all_n, keep_ll / keep_n, sh_ll / sh_n
    return {"L": int(L), "frac_genes_kept": keep_n / all_n,
            "B_cohort": lr_all - lr_keep, "delta_matched": lr_keep - ls_keep,
            "delta_v14_8": lr_all - ls_keep}


def _mi_from_counts(B, mm=False):
    N = float(B.sum())
    if N <= 0:
        return float("nan"), 0, 0, 0, 0
    rx, ry = B.sum(axis=1), B.sum(axis=0)
    den = np.maximum(rx[:, None] * ry[None, :], TINY)
    mask = B > 0
    ratio = np.where(mask, B * N / den, 1.0)
    I = float(np.sum(np.where(mask, (B / N) * np.log(np.maximum(ratio, TINY)), 0.0)))
    mXY = int(np.count_nonzero(B))
    mX, mY = int(np.count_nonzero(rx)), int(np.count_nonzero(ry))
    if mm:
        I -= (mXY - mX - mY + 1) / (2.0 * N)
    return I, mXY, mX, mY, int(N)


def _mi_terms(obs, den, N):
    """Elementwise plug-in MI contributions over an arbitrary projected table."""
    m = obs > 0
    r = np.where(m, obs / np.maximum(den, TINY), 1.0)
    return np.where(m, (obs / N) * np.log(np.maximum(r, TINY)), 0.0)


def _mi_tab(obs, den, N):
    """Plug-in MI term sum over an arbitrary projected table."""
    return float(np.sum(_mi_terms(obs, den, N)))


def _kl_joint(P, Q):
    m = P > 0
    return float(np.sum(P[m] * np.log(P[m] / np.maximum(Q[m], TINY))))


def js_divergence(P, Q):
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    P = P / P.sum(); Q = Q / Q.sum()
    M = 0.5 * (P + Q)
    return 0.5 * _kl_joint(P, M) + 0.5 * _kl_joint(Q, M)


def psi_constant(seqs, logP):
    q = logP.shape[0]
    B = bigram_counts(seqs, q)
    marg = np.zeros(q)
    for s in seqs:
        if s.size >= 2:
            marg += np.bincount(s, minlength=q)
    N = B.sum()
    if N == 0:
        return None
    Phat = B / N
    pL, pR = Phat.sum(axis=1), Phat.sum(axis=0)
    pi = marg / marg.sum()
    Qk = np.exp(logP)
    I, mXY, mX, mY, _ = _mi_from_counts(B)
    Dmix_LR = Dmis = 0.0
    for x in np.flatnonzero(pL > 0):
        Dmix_LR += pL[x] * _kl_joint(pR, Qk[x])
        Dmis += pL[x] * _kl_joint(Phat[x] / pL[x], Qk[x])
    Psi_LR = float(np.sum((Phat - pL[:, None] * pR[None, :]) * logP))
    Psi_pi = float(np.sum((Phat - pi[:, None] * pi[None, :]) * logP))
    delta_edge = float(np.sum((pL[:, None] * pR[None, :] - pi[:, None] * pi[None, :]) * logP))
    return {"Psi_LR": Psi_LR, "Psi_pi": Psi_pi, "Psi": Psi_LR,
            "I_bigram_plugin": I, "I_bigram_MM": I - (mXY - mX - mY + 1) / (2.0 * N),
            "D_mix_LR": Dmix_LR, "D_mismatch": Dmis, "delta_edge": delta_edge,
            "identity_abs_err": abs(Psi_LR - (I + Dmix_LR - Dmis)),
            "edge_abs_err": abs(Psi_pi - (Psi_LR + delta_edge)),
            "n_transitions": int(N)}


def jeffreys_of_counts(B):
    N = B.sum()
    P = B / N
    pL, pR = P.sum(1), P.sum(0)
    prod = pL[:, None] * pR[None, :]
    I = _kl_joint(P, prod)
    Jrev = _kl_joint(prod, P)
    return I, Jrev, I + Jrev


def mc_h4_check(seqs, logP, L, reps, seed):
    rng = np.random.default_rng(seed)
    real, sh, ntr = 0.0, 0.0, 0
    for s in seqs:
        if s.size < 2:
            continue
        real += float(np.sum(logP[s[:-1], s[1:]]))
        ntr += s.size - 1
    if ntr == 0 or reps <= 0:
        return None
    for _ in range(reps):
        acc = 0.0
        for s in seqs:
            n = s.size
            if n < 2:
                continue
            m = n // L
            if m < 2:
                acc += float(np.sum(logP[s[:-1], s[1:]]))
                continue
            out = s[:m * L].reshape(m, L)[rng.permutation(m)].reshape(-1)
            if n > m * L:
                out = np.concatenate([out, s[m * L:]])
            acc += float(np.sum(logP[out[:-1], out[1:]]))
        sh += acc
    return {"L": int(L), "delta_H4_montecarlo": (real - sh / reps) / ntr,
            "reps": int(reps), "seed": int(seed)}


# =====================================================================
# 6. THEOREM M  Klein four-group quadrants (signed-safe)
# =====================================================================
def klein_traces(q=NC, mirror=MIRROR):
    x = np.arange(q)
    tr_id = q * q
    tr_T = int(np.sum(x[:, None] == x[None, :]))
    tr_R = int(np.sum(mirror[x][:, None] == x[None, :]))
    tr_S = int(np.sum((mirror[x][:, None] == x[None, :]) &
                      (mirror[x][None, :] == x[None, :])))
    return {"id": tr_id, "T": tr_T, "R": tr_R, "S": tr_S}


def klein_dims(q=NC, mirror=MIRROR):
    t = klein_traces(q, mirror)
    d = {}
    for eT in (+1, -1):
        for eR in (+1, -1):
            key = ("+" if eT > 0 else "-") + ("+" if eR > 0 else "-")
            d["W_" + key] = (t["id"] + eT * t["T"] + eR * t["R"] + eT * eR * t["S"]) // 4
    d["sum"] = sum(v for k, v in d.items() if k != "sum")
    return d, t


def klein_quadrants_of_matrix(M):
    """V_4 quadrant norms of ANY real 64x64 matrix (SIGNED, zero-mass allowed)."""
    M = np.asarray(M, dtype=np.float64)
    TP = M.T
    RP = M.T[np.ix_(MIRROR, MIRROR)]
    SP = M[np.ix_(MIRROR, MIRROR)]
    Q = {"pp": 0.25 * (M + TP + RP + SP),
         "pm": 0.25 * (M + TP - RP - SP),
         "mp": 0.25 * (M - TP + RP - SP),
         "mm": 0.25 * (M - TP - RP + SP)}
    n2 = {k: float(np.sum(v * v)) for k, v in Q.items()}
    n_tot = float(np.sum(M * M))
    recon = Q["pp"] + Q["pm"] + Q["mp"] + Q["mm"]
    qt = 0.25 * float(np.sum((M - TP) ** 2))
    qm = 0.25 * float(np.sum((M - RP) ** 2))
    return {"norm2_total": n_tot,
            "norm2_W_pp": n2["pp"],
            "norm2_W_pm_static_strand": n2["pm"],
            "norm2_W_mp_sym_time_flux": n2["mp"],
            "norm2_W_mm_duplex_time_flux": n2["mm"],
            "directed_program_norm2": n2["mp"] + n2["mm"],
            "reconstruction_max_abs_err": float(np.abs(recon - M).max()),
            "pythagorean_residual": abs(n_tot - sum(n2.values())),
            "time_flux_split_residual": abs(qt - (n2["mp"] + n2["mm"])),
            "mirror_flux_split_residual": abs(qm - (n2["pm"] + n2["mm"])),
            "dynamical_share_of_mirror_flux_F2":
                (n2["mm"] / (n2["pm"] + n2["mm"])) if (n2["pm"] + n2["mm"]) > 0 else None}


def theorem_M_klein_quadrants(B):
    N = float(np.sum(B))
    if N <= 0:
        return None
    out = klein_quadrants_of_matrix(np.asarray(B, float) / N)
    dims, tr = klein_dims()
    out["dims_exact"] = dims
    out["traces"] = tr
    return out


# =====================================================================
# 7. THEOREM C / C.2 / C.3 / F / F.2 / L
# =====================================================================
def lag_counts(seqs, q, k):
    B = np.zeros((q, q), dtype=np.float64)
    for s in seqs:
        if s.size <= k:
            continue
        B += np.bincount(s[:-k] * q + s[k:], minlength=q * q).reshape(q, q)
    return B


def het_expected_counts(seqs, q, k):
    G = len(seqs)
    C = np.zeros((G, q), dtype=np.float64)
    w = np.zeros(G, dtype=np.float64)
    for i, s in enumerate(seqs):
        n = s.size
        C[i] = np.bincount(s, minlength=q)
        if n > k and n > 1:
            w[i] = (n - k) / (n * (n - 1.0))
    return C.T @ (w[:, None] * C) - np.diag(C.T @ w)


def theorem_L_diagnostics(B, E, label="het"):
    N = float(B.sum())
    if N <= 0:
        return None
    Esym_err = float(np.abs(E - E.T).max())
    A_obs = float(0.5 * np.abs(B - B.T).sum() / N)
    P = B / N
    js = js_divergence(P, P.T)
    return {"null_label": label,
            "null_symmetry_max_abs_err": Esym_err,
            "null_is_symmetric": bool(Esym_err < 1e-9),
            "null_A_time": float(0.5 * np.abs(E - E.T).sum() / max(E.sum(), TINY)),
            "null_JS_time": js_divergence(E, E.T) if E.sum() > 0 else None,
            "A_time_obs": A_obs, "JS_time_obs": js,
            "JS_ge_half_A2": bool(js >= 0.5 * A_obs ** 2 - 1e-12),
            "JS_minus_half_A2": js - 0.5 * A_obs ** 2}


def window_cache(seqs, q, W, max_windows=150000):
    full, tails = [], []
    for s in seqs:
        n = s.size
        nf = n // W
        if nf:
            full.append(s[:nf * W].reshape(nf, W))
        if n > nf * W:
            t = s[nf * W:]
            if t.size >= 2:
                tails.append(t)
    S_full = np.concatenate(full, axis=0) if full else np.empty((0, W), dtype=np.int64)
    total = int(S_full.shape[0]) + len(tails)
    stride = 1
    if S_full.shape[0] > max_windows:
        stride = int(math.ceil(S_full.shape[0] / max_windows))
        S_full = S_full[::stride]
    if S_full.shape[0] == 0 and not tails:
        return None
    C_full = np.zeros((S_full.shape[0], q), dtype=np.float64)
    for j in range(q):
        C_full[:, j] = np.sum(S_full == j, axis=1)
    C_tail = np.array([np.bincount(t, minlength=q) for t in tails], dtype=np.float64) \
        if tails else np.empty((0, q))
    L_tail = np.array([t.size for t in tails], dtype=np.int64) if tails else np.empty(0, np.int64)
    return {"S_full": S_full, "C_full": C_full, "tails": tails,
            "C_tail": C_tail, "L_tail": L_tail, "W": int(W),
            "n_windows_total": total, "n_full_used": int(S_full.shape[0]),
            "n_tails": len(tails), "stride": stride}


def window_pairs_cached(wc, q, k):
    W, S = wc["W"], wc["S_full"]
    obs = np.zeros((q, q), dtype=np.float64)
    if S.shape[0] and W > k:
        idx = (S[:, :-k] * q + S[:, k:]).ravel()
        obs += np.bincount(idx, minlength=q * q).reshape(q, q)
    parts_C, parts_w = [], []
    if S.shape[0] and W > k:
        parts_C.append(wc["C_full"])
        parts_w.append(np.full(S.shape[0], (W - k) / (W * (W - 1.0))))
    for i, t in enumerate(wc["tails"]):
        L = int(wc["L_tail"][i])
        if L <= k:
            continue
        obs += np.bincount(t[:-k] * q + t[k:], minlength=q * q).reshape(q, q)
        parts_C.append(wc["C_tail"][i:i + 1])
        parts_w.append(np.array([(L - k) / (L * (L - 1.0))]))
    if not parts_C:
        return obs, np.zeros((q, q))
    C = np.concatenate(parts_C, axis=0)
    w = np.concatenate(parts_w)
    return obs, C.T @ (w[:, None] * C) - np.diag(C.T @ w)


def kl_of_joint(B, E):
    N1, N2 = B.sum(), E.sum()
    if N1 <= 0 or N2 <= 0:
        return float("nan")
    return _kl_joint(B / N1, np.maximum(E / N2, TINY))


def lagged_mi_profile(seqs, q, kmax):
    out = []
    for k in range(1, kmax + 1):
        B = lag_counts(seqs, q, k)
        Iobs, mXY, mX, mY, N = _mi_from_counts(B, mm=True)
        Iraw, *_ = _mi_from_counts(B, mm=False)
        if N == 0:
            out.append({"lag": k, "N": 0})
            continue
        E = het_expected_counts(seqs, q, k)
        Ihet, *_ = _mi_from_counts(E, mm=False)
        kl = kl_of_joint(B, E)
        out.append({"lag": k, "N": N, "support_cells": mXY,
                    "I_plugin": Iraw, "I_MM": Iobs, "MM_correction": Iraw - Iobs,
                    "I_het": Ihet,
                    "het_fraction_of_obs": (Ihet / Iobs) if Iobs > 0 else None,
                    "delta_I_intra": Iobs - Ihet,
                    "intra_fraction_bounded": ((Iobs - Ihet) / Iobs) if Iobs > 0 else None,
                    "KL_obs_vs_het": kl,
                    "eta_KL_over_Iobs": (kl / Iobs) if Iobs > 0 else None,
                    "theorem_L": theorem_L_diagnostics(B, E, "het_gene"),
                    "theorem_M": theorem_M_klein_quadrants(B)})
    return out


def stationary(P):
    n = P.shape[0]
    A = P.T - np.eye(n)
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    pi, *_ = np.linalg.lstsq(A, b, rcond=None)
    pi = np.real_if_close(np.asarray(pi, dtype=np.float64))
    pi = np.clip(pi, 0.0, None)
    s = pi.sum()
    return pi / s if s > 0 else np.ones(n) / n


def markov_lag_prediction(P, kmax, do_spectral=True):
    q = P.shape[0]
    pi = stationary(P)
    ev = np.linalg.eigvals(P)
    mods = np.sort(np.abs(ev))[::-1]
    lam2 = float(mods[1]) if mods.size > 1 else 0.0
    xi = float(-1.0 / math.log(lam2)) if 0 < lam2 < 1 else float("inf")
    J1 = pi[:, None] * P
    rev_defect = float(np.abs(J1 - J1.T).sum())
    Padj = (P.T * pi[None, :]) / np.maximum(pi[:, None], TINY)
    comm = P @ Padj - Padj @ P
    nrm = float(np.linalg.norm(P, ord="fro")) ** 2
    normal_defect = float(np.linalg.norm(comm, ord="fro") / max(nrm, TINY))
    ones_pi = np.ones((q, 1)) * pi[None, :]
    wgt = np.where(pi[None, :] > 0, pi[:, None] / np.maximum(pi[None, :], TINY), 0.0)
    spec = None
    if do_spectral:
        try:
            lam, V = np.linalg.eig(P)
            U = np.linalg.inv(V)
            ip = int(np.argmax(np.abs(lam)))
            keep = [j for j in range(q) if j != ip]
            Vk, Uk, lk = V[:, keep], U[keep, :], lam[keep]
            Gv = (Vk.conj().T * pi[None, :]) @ Vk
            Gu = (Uk.conj() / np.maximum(pi[None, :], TINY)) @ Uk.T
            kappa = np.real(np.diag(Gv) * np.diag(Gu))
            spec = {"lam": lk, "GvGu": Gv * Gu, "kappa": kappa,
                    "cond_V": float(np.linalg.cond(V)),
                    "kappa_min": float(kappa.min()) if kappa.size else None,
                    "kappa_mean": float(kappa.mean()) if kappa.size else None,
                    "kappa_max": float(kappa.max()) if kappa.size else None}
        except Exception:
            spec = None
    E1 = P - ones_pi
    E1 -= ones_pi @ E1
    E1 -= (E1 @ np.ones((q, 1))) * pi[None, :]
    Pk = np.eye(q)
    Ek = np.eye(q) - ones_pi
    out = []
    for k in range(1, kmax + 1):
        Pk = Pk @ P
        Ek = Ek @ E1
        Ek -= ones_pi @ Ek
        Ek -= (Ek @ np.ones((q, 1))) * pi[None, :]
        Iex, *_ = _mi_from_counts(pi[:, None] * Pk, mm=False)
        chi2 = 0.5 * float(np.sum(wgt * Ek * Ek))
        row = {"lag": k, "I_markov": Iex, "I_chi2_proxy": chi2,
               "chi2_ratio": (Iex / chi2) if chi2 > 0 else None,
               "I_reversible_spectral_sum": 0.5 * float(np.sum(mods[1:] ** (2.0 * k)))}
        if spec is not None:
            Lk = spec["lam"] ** k
            Mspec = np.outer(np.conj(Lk), Lk) * spec["GvGu"]
            row["chi2_biorthogonal_full"] = 0.5 * float(np.real(Mspec.sum()))
            row["chi2_diagonal_kappa"] = 0.5 * float(
                np.sum((np.abs(spec["lam"]) ** (2.0 * k)) * spec["kappa"]))
            row["biorthogonal_rel_err"] = abs(row["chi2_biorthogonal_full"] - chi2) / max(chi2, TINY)
        out.append(row)
    meta = {"lambda2_modulus": lam2, "xi_codons": xi,
            "reversibility_defect_L1": rev_defect,
            "is_reversible": bool(rev_defect < 1e-10),
            "normality_defect_relF": normal_defect,
            "is_normal_in_L2pi": bool(normal_defect < 1e-9),
            "spectral_sum_valid": bool(rev_defect < 1e-10),
            "n_complex_eigenvalues": int(np.sum(np.abs(ev.imag) > 1e-12))}
    if spec is not None:
        meta.update({k: spec[k] for k in ("cond_V", "kappa_min", "kappa_mean", "kappa_max")})
    return out, meta


# =====================================================================
# 8. THEOREM G / G.2 / J / P
# =====================================================================
def duplex_df_table(aa, support=None):
    lay, F, sizes, _ = duplex_layers(aa)
    if support is None:
        cod, chans = np.arange(NC), np.arange(NAA)
    else:
        cod = np.flatnonzero(support)
        chans = np.unique(aa[cod])
    q, nA = int(cod.size), int(chans.size)
    ka = np.array([int(np.sum(aa[cod] == a)) for a in chans])
    syn = int(np.sum(ka - 1))
    joint = int(len(np.unique(lay[cod])))
    df_AA = (nA - 1) ** 2
    df_d2b = syn * (nA - 1)
    df_anti = (q - 1) * (joint - nA)
    df_dup = (q - 1) * (q - joint)
    return {"n_codons": q, "n_AA_channels": nA, "n_joint_channels": joint,
            "syn_dim_single_strand": syn, "syn_dim_duplex": q - joint,
            "absorbed_by_antisense": syn - (q - joint),
            "df_AA": df_AA, "df_d2b": df_d2b, "df_anti": df_anti, "df_dup": df_dup,
            "df_syn": df_anti + df_dup, "df_total": (q - 1) ** 2,
            "df_sum_check": df_AA + df_d2b + df_anti + df_dup,
            "df_anti_share_of_syn": (df_anti / (df_anti + df_dup)) if (df_anti + df_dup) else None}


def theorem_P_transversality(aa):
    key = aa.astype(np.int64) * 32 + MIRROR32
    uniq, cnts = np.unique(key, return_counts=True)
    self_pairs, diag = [], np.zeros(NAA, dtype=np.int64)
    for c in range(NC):
        r = int(MIRROR[c])
        if int(aa[c]) == int(aa[r]):
            diag[int(aa[c])] += 1
            if c < r:
                self_pairs.append([CODONS[c], CODONS[r], AA_STRING[int(aa[c])]])
    return {"n_joint_AA_mirror_classes": int(uniq.size),
            "join_is_identity": bool(uniq.size == NC),
            "join_synonymous_dim": int(NC - uniq.size),
            "max_join_class_size": int(cnts.max()),
            "n_self_complementary_aa_mirror_pairs": len(self_pairs),
            "self_complementary_pairs": self_pairs,
            "relic_diagonal_total": int(diag.sum()),
            "note": ("100% mirror capture is an INJECTIVITY fact (32 labels vs 21 "
                     "for alpha_rho), not evidence of a privileged mirror channel")}


def _s_syn_of_counts(B, aa):
    N = float(B.sum())
    if N <= 0:
        return float("nan")
    q = B.shape[0]
    xs = np.arange(q)
    ny = B.sum(axis=0)
    nxb = np.bincount((xs[:, None] * NAA + aa[None, :]).ravel(),
                      weights=B.ravel(), minlength=q * NAA).reshape(q, NAA)
    nb = np.bincount(aa, weights=ny, minlength=NAA)
    bY = aa[None, :]
    den = np.maximum(nxb[:, aa] * ny[None, :], TINY) / np.maximum(nb[bY], TINY)
    return _mi_tab(B, den, N)


def lag_duplex_decomposition(seqs, q, k, aa, d8mask, df=None, cpg=None):
    B = lag_counts(seqs, q, k)
    N = float(B.sum())
    if N <= 0:
        return None
    ar = aa[MIRROR]
    nx, ny = B.sum(axis=1), B.sum(axis=0)
    xs = np.arange(q)
    nab = np.bincount((aa[:, None] * NAA + aa[None, :]).ravel(),
                      weights=B.ravel(), minlength=NAA * NAA).reshape(NAA, NAA)
    nxb = np.bincount((xs[:, None] * NAA + aa[None, :]).ravel(),
                      weights=B.ravel(), minlength=q * NAA).reshape(q, NAA)
    na, nb = nab.sum(axis=1), nab.sum(axis=0)
    nbr = np.bincount(aa * NAA + ar, weights=ny,
                      minlength=NAA * NAA).reshape(NAA, NAA)
    nxbr = np.bincount((xs[:, None] * NAA * NAA + aa[None, :] * NAA + ar[None, :]).ravel(),
                       weights=B.ravel(), minlength=q * NAA * NAA).reshape(q, NAA, NAA)
    Icc = _mi_tab(B, np.maximum(nx[:, None] * ny[None, :], TINY) / N, N)
    Iaa = _mi_tab(nab, np.maximum(na[:, None] * nb[None, :], TINY) / N, N)
    d2b = _mi_tab(nxb, np.maximum(nx[:, None] * nab[aa, :], TINY)
                  / np.maximum(na[aa[:, None]], TINY), N)
    bY, rY = aa[None, :], ar[None, :]
    Ssyn = _mi_tab(B, np.maximum(nxb[:, aa] * ny[None, :], TINY)
                   / np.maximum(nb[bY], TINY), N)
    t_anti = np.where(nxbr > 0, (nxbr / N) * np.log(np.maximum(
        nxbr / np.maximum(np.maximum(nxb[:, :, None] * nbr[None, :, :], TINY)
                          / np.maximum(nb[None, :, None], TINY), TINY), TINY)), 0.0)
    Santi = float(t_anti.sum())
    nxbr_xy = nxbr[:, aa, ar]
    den_dup = np.maximum(nxbr_xy * ny[None, :], TINY) / np.maximum(nbr[bY, rY], TINY)
    m_dup = B > 0
    t_dup = np.where(m_dup, (B / N) * np.log(np.maximum(
        np.where(m_dup, B / np.maximum(den_dup, TINY), 1.0), TINY)), 0.0)
    Sdup = float(t_dup.sum())
    Banti = B.T[np.ix_(MIRROR, MIRROR)]
    Ianti, *_ = _mi_from_counts(Banti, mm=False)
    _layD, _FD, _perD, _mulD = duplex_layers(aa)
    _logkD = np.log(np.maximum(_perD.astype(np.float64), 1.0))
    d8mass = float(np.sum(ny[d8mask]) / N)
    d8bound = float(np.sum(ny * _logkD) / N)
    res = {"lag": k, "N": int(N), "I_CC": Icc, "I_AA": Iaa, "d2b": d2b,
           "S_syn": Ssyn, "S_anti": Santi, "S_dup": Sdup,
           "tripartite_residual": Icc - (Iaa + d2b + Ssyn),
           "duplex_residual": Ssyn - (Santi + Sdup),
           "antisense_fraction_raw": (Santi / Ssyn) if Ssyn > 0 else None,
           "I_CC_antiparallel": Ianti,
           "duplex_invariance_abs_err": abs(Icc - Ianti),
           "D8_mass": d8mass, "D8_bound_nats": d8bound,
           "D8_bound_log2_legacy": d8mass * math.log(2.0),
           "D8_has_quartets": bool(np.any(_perD > 2)),
           "D8_bound_ok": bool(Sdup <= d8bound + 1e-12)}
    # ---- [THM-P] mirror-orbit capture ----
    NM = 32
    lay_m = MIRROR32
    nbm = np.zeros((NAA, NM), dtype=np.float64)
    np.add.at(nbm, (aa, lay_m), ny)
    nxbm = np.bincount((xs[:, None] * NAA * NM + aa[None, :] * NM
                        + lay_m[None, :]).ravel(), weights=B.ravel(),
                       minlength=q * NAA * NM).reshape(q, NAA, NM)
    S_mirror = _mi_tab(nxbm, np.maximum(nxb[:, :, None] * nbm[None, :, :], TINY)
                       / np.maximum(nb[None, :, None], TINY), N)
    nxbm_xy = nxbm[:, aa, lay_m]
    S_res_m = _mi_tab(B, np.maximum(nxbm_xy * ny[None, :], TINY)
                      / np.maximum(nbm[aa[None, :], lay_m[None, :]], TINY), N)
    res["S_mirror_orbit"] = S_mirror
    res["mirror_capture_residual"] = Ssyn - S_mirror
    res["mirror_cond_residual"] = S_res_m
    res["mirror_capture_fraction"] = (S_mirror / Ssyn) if Ssyn > 0 else None
    if cpg is not None:
        cm, lm = cpg["cpg_mask"], cpg["leu_mask"]
        Sdup_cpg = float(t_dup[:, cm].sum())
        Sdup_leu = float(t_dup[:, lm].sum())
        mass_cpg = float(ny[cm].sum() / N)
        mass_leu = float(ny[lm].sum() / N)
        H_arho_a = 0.0
        for b in range(NAA):
            if nb[b] <= 0:
                continue
            for r in range(NAA):
                if nbr[b, r] > 0:
                    H_arho_a -= (nbr[b, r] / N) * math.log(nbr[b, r] / nb[b])
        H_c_a_arho = 0.0
        for y in range(q):
            if ny[y] > 0:
                d = nbr[aa[y], ar[y]]
                if d > 0:
                    H_c_a_arho -= (ny[y] / N) * math.log(ny[y] / d)
        u_anti = (Santi / H_arho_a) if H_arho_a > 1e-15 else None
        u_dup = (Sdup / H_c_a_arho) if H_c_a_arho > 1e-15 else None
        res.update({
            "S_dup_cpg": Sdup_cpg, "S_dup_leu": Sdup_leu,
            "S_dup_split_abs_err": abs(Sdup - (Sdup_cpg + Sdup_leu)),
            "D8_cpg_mass": mass_cpg, "D8_leu_mass": mass_leu,
            "cpg_bound_ok": bool(Sdup_cpg <= float(np.sum(ny[cm] * _logkD[cm]) / N) + 1e-12),
            "leu_bound_ok": bool(Sdup_leu <= float(np.sum(ny[lm] * _logkD[lm]) / N) + 1e-12),
            "cpg_share_of_S_dup": (Sdup_cpg / Sdup) if Sdup > 0 else None,
            "H_Arho_given_A": H_arho_a, "H_C_given_A_Arho": H_c_a_arho,
            "capacity_util_anti": u_anti, "capacity_util_dup": u_dup,
            "R_cap": (u_anti / u_dup) if (u_anti is not None and u_dup not in (None, 0)) else None})
    if df:
        da_full, dd_full = df["df_anti"], df["df_dup"]
        da_eff = 0
        for b in range(NAA):
            R = int(np.count_nonzero(nxb[:, b] > 0))
            Cc = int(np.count_nonzero(nbr[b, :] > 0))
            if R > 1 and Cc > 1:
                da_eff += (R - 1) * (Cc - 1)
        dd_eff = 0
        for b in range(NAA):
            for r in range(NAA):
                if nbr[b, r] <= 0:
                    continue
                R = int(np.count_nonzero(nxbr[:, b, r] > 0))
                ys = np.flatnonzero((aa == b) & (ar == r))
                Cc = int(np.count_nonzero(ny[ys] > 0))
                if R > 1 and Cc > 1:
                    dd_eff += (R - 1) * (Cc - 1)
        out = {"df_anti_full": da_full, "df_dup_full": dd_full,
               "df_anti_eff": da_eff, "df_dup_eff": dd_eff}
        out["enrich_per_df_raw"] = ((Santi / da_eff) / (Sdup / dd_eff)) \
            if (da_eff and dd_eff and Sdup > 0) else None
        for tag, da, dd in (("full", da_full, dd_full), ("eff", da_eff, dd_eff)):
            ba, bd = da / (2.0 * N), dd / (2.0 * N)
            ac, dc = Santi - ba, Sdup - bd
            sc = ac + dc
            out[f"MM_bias_anti_{tag}"] = ba
            out[f"MM_bias_dup_{tag}"] = bd
            out[f"S_anti_MM_{tag}"] = ac
            out[f"S_dup_MM_{tag}"] = dc
            out[f"anti_fraction_MM_{tag}"] = (ac / sc) if sc > 0 else None
            gated = bool(Santi > 2.0 * ba and Sdup > 2.0 * bd and ac > 0 and dc > 0)
            out[f"enrich_gate_ok_{tag}"] = gated
            out[f"enrich_per_df_{tag}"] = ((ac / da) / (dc / dd)) \
                if (gated and da and dd) else None
        out["df_anti_share_full"] = da_full / (da_full + dd_full) if (da_full + dd_full) else None
        out["df_anti_share_eff"] = da_eff / (da_eff + dd_eff) if (da_eff + dd_eff) else None
        res.update(out)
    return res


# =====================================================================
# 9. THEOREM H / H.3 / K  tower of exact layer nulls
# =====================================================================
def layer_tower_expected(seqs, q, k, layer, F, max_genes=0):
    S = seqs if max_genes <= 0 else seqs[:max_genes]
    G = len(S)
    fibers = [np.flatnonzero(layer == a) for a in range(F)]
    B = lag_counts(S, q, k)
    meta = {"n_genes": G, "F": F, "lag": k, "n_genes_input": len(seqs),
            "cohort_truncated": bool(G != len(seqs))}
    if G == 0 or B.sum() <= 0:
        return np.zeros((q, q)), B, meta
    C = np.zeros((G, q), dtype=np.float64)
    for g, s in enumerate(S):
        C[g] = np.bincount(s, minlength=q)
    Nlay = np.zeros((G, F), dtype=np.float64)
    for a in range(F):
        if fibers[a].size:
            Nlay[:, a] = C[:, fibers[a]].sum(axis=1)
    M = np.zeros((G, F, F), dtype=np.float64)
    for g, s in enumerate(S):
        if s.size <= k:
            continue
        Ls = layer[s]
        M[g] = np.bincount(Ls[:-k] * F + Ls[k:], minlength=F * F).reshape(F, F)
    V = np.zeros((G, q), dtype=np.float64)
    Wm = np.zeros((G, q), dtype=np.float64)
    Wd = np.zeros((G, q), dtype=np.float64)
    diag_exact = np.zeros(q, dtype=np.float64)
    for a in range(F):
        fa = fibers[a]
        if fa.size == 0:
            continue
        Na = Nlay[:, a]
        ok1 = Na >= 1
        ok2 = Na >= 2
        inv1 = np.zeros(G); inv1[ok1] = 1.0 / Na[ok1]
        inv2 = np.zeros(G); inv2[ok2] = 1.0 / (Na[ok2] - 1.0)
        V[:, fa] = C[:, fa] * inv1[:, None]
        Wm[:, fa] = C[:, fa] * inv1[:, None]
        Wd[:, fa] = C[:, fa] * inv2[:, None]
        if fa.size == 1:
            diag_exact[int(fa[0])] = float(M[ok2, a, a].sum())
        else:
            wd = np.zeros(G)
            wd[ok2] = M[ok2, a, a] / (Na[ok2] * (Na[ok2] - 1.0))
            cx = C[:, fa]
            diag_exact[fa] = np.sum(wd[:, None] * cx * (cx - 1.0), axis=0)
    E = np.zeros((q, q), dtype=np.float64)
    for a in range(F):
        fa = fibers[a]
        if fa.size == 0:
            continue
        VaT = V[:, fa].T
        m_off = M[:, a, :].copy()
        m_off[:, a] = 0.0
        E[fa, :] = VaT @ (m_off[:, layer] * Wm)
        if fa.size >= 2:
            E[np.ix_(fa, fa)] = VaT @ (M[:, a, a][:, None] * Wd[:, fa])
    np.fill_diagonal(E, diag_exact)
    meta["min_E_value"] = float(E.min())
    meta["no_clipping_required"] = bool(E.min() >= 0.0)
    return E, B, meta


def _layer_tower_reference(seqs, q, k, layer, F):
    """Slow F^2 reference implementation, used ONLY as a regression fixture."""
    G = len(seqs)
    fibers = [np.flatnonzero(layer == a) for a in range(F)]
    B = lag_counts(seqs, q, k)
    C = np.zeros((G, q), dtype=np.float64)
    for g, s in enumerate(seqs):
        C[g] = np.bincount(s, minlength=q)
    Nlay = np.zeros((G, F), dtype=np.float64)
    for a in range(F):
        if fibers[a].size:
            Nlay[:, a] = C[:, fibers[a]].sum(axis=1)
    M = np.zeros((G, F, F), dtype=np.float64)
    for g, s in enumerate(seqs):
        if s.size <= k:
            continue
        Ls = layer[s]
        M[g] = np.bincount(Ls[:-k] * F + Ls[k:], minlength=F * F).reshape(F, F)
    V = np.zeros((G, q)); Wm = np.zeros((G, q)); Wd = np.zeros((G, q))
    dg = np.zeros(q)
    for a in range(F):
        fa = fibers[a]
        if fa.size == 0:
            continue
        Na = Nlay[:, a]
        o1, o2 = Na >= 1, Na >= 2
        i1 = np.zeros(G); i1[o1] = 1.0 / Na[o1]
        i2 = np.zeros(G); i2[o2] = 1.0 / (Na[o2] - 1.0)
        V[:, fa] = C[:, fa] * i1[:, None]
        Wm[:, fa] = C[:, fa] * i1[:, None]
        Wd[:, fa] = C[:, fa] * i2[:, None]
        if fa.size == 1:
            dg[int(fa[0])] = float(M[o2, a, a].sum())
        else:
            wd = np.zeros(G); wd[o2] = M[o2, a, a] / (Na[o2] * (Na[o2] - 1.0))
            cx = C[:, fa]
            dg[fa] = np.sum(wd[:, None] * cx * (cx - 1.0), axis=0)
    E = np.zeros((q, q))
    for a in range(F):
        fa = fibers[a]
        if fa.size == 0:
            continue
        for b in range(F):
            fb = fibers[b]
            if fb.size == 0:
                continue
            w = M[:, a, b]
            if not np.any(w):
                continue
            Wb = Wd[:, fb] if b == a else Wm[:, fb]
            E[np.ix_(fa, fb)] = V[:, fa].T @ (w[:, None] * Wb)
    np.fill_diagonal(E, dg)
    return E, B


def tower_projection_check(E, B, layer, F):
    idx = (layer[:, None] * F + layer[None, :]).ravel()
    El = np.bincount(idx, weights=E.ravel(), minlength=F * F).reshape(F, F)
    Bl = np.bincount(idx, weights=B.ravel(), minlength=F * F).reshape(F, F)
    return float(np.abs(El - Bl).max()), float(abs(E.sum() - B.sum()))


def lemma_H3_check(E, B):
    m = B > 0
    if not m.any():
        return None
    Ev, Bv = E[m], B[m]
    return {"n_support_cells": int(m.sum()),
            "min_E_on_support": float(Ev.min()),
            "min_E_over_B_on_support": float((Ev / Bv).min()),
            "absolute_continuity_holds": bool(Ev.min() > 0.0),
            "n_zero_E_on_support": int(np.count_nonzero(Ev <= 0.0))}


def theorem_K_localisation(E, B, layer):
    sizes = np.bincount(layer, minlength=int(layer.max()) + 1)
    single = sizes[layer] == 1
    ss = single[:, None] & single[None, :]
    N1, N2 = B.sum(), E.sum()
    if N1 <= 0 or N2 <= 0:
        return None
    P, Q = B / N1, np.maximum(E / N2, TINY)
    contrib = np.where(P > 0, P * np.log(np.maximum(P, TINY) / Q), 0.0)
    tot = float(contrib.sum())
    ss_kl = float(contrib[ss].sum())
    return {"n_singleton_codons": int(single.sum()),
            "n_singleton_pair_cells": int(ss.sum()),
            "singleton_cell_max_abs_err": float(np.abs(E[ss] - B[ss]).max()) if ss.any() else 0.0,
            "KL_total": tot, "KL_on_singleton_pairs": ss_kl,
            "KL_localised_on_D8": bool(abs(ss_kl) < 1e-12),
            "D8_KL_fraction": (1.0 - ss_kl / tot) if tot > 0 else None}


# =====================================================================
# 10. THEOREM Q / Q.2  V_4 equivariance, both commutator invariants
# =====================================================================
def _layer_avg_matrix(layer, F):
    sizes = np.bincount(layer, minlength=F).astype(np.float64)
    D = (layer[:, None] == layer[None, :]).astype(np.float64)
    D /= np.maximum(sizes[layer][:, None], TINY)
    return D


def klein_commutator_audit(aa, probe=None):
    """[THM-Q/Q.2] exact 1-index and 2-index commutator invariants.
    S_2 = Tr(D R D R) = sum_(a,b) m_ab^2/(k_a k_b)
    ||[D,R]||_F^2   = 2C   - 2 S_2      (operator on R^64)
    ||[Pi,R]||_HS^2 = 2C^2 - 2 S_2^2    (superoperator on R^{64x64})
    Table 1 AA: C=21, S_2=569/72 -> 943/36 and 1962383/2592.  All COMPUTED."""
    Rr = np.zeros((NC, NC))
    Rr[np.arange(NC), MIRROR] = 1.0
    layD, FD, _, _ = duplex_layers(aa)
    if probe is None:
        ix = np.arange(NC)
        probe = (1 + ((ix[:, None] * 17 + ix[None, :] * 31
                       + (ix[:, None] ^ ix[None, :])) % 97)).astype(float)
        probe /= probe.sum()
    TP = probe.T
    RP = probe.T[np.ix_(MIRROR, MIRROR)]
    SP = probe[np.ix_(MIRROR, MIRROR)]
    quads = {"pp": 0.25 * (probe + TP + RP + SP),
             "pm": 0.25 * (probe + TP - RP - SP),
             "mp": 0.25 * (probe - TP + RP - SP),
             "mm": 0.25 * (probe - TP - RP + SP)}
    out = {}
    for name, lay, F in (("AA", aa, NAA), ("dup", layD, FD), ("mirror", MIRROR32, 32)):
        D = _layer_avg_matrix(lay, F)
        cR = float(np.linalg.norm(D @ Rr - Rr @ D, ord="fro"))
        cT = float(np.linalg.norm((D @ probe @ D).T - D @ probe.T @ D, ord="fro"))
        sizes = np.bincount(lay, minlength=F).astype(np.float64)
        kx, krx = sizes[lay], sizes[lay[MIRROR]]
        same = (lay[:, None] == lay[None, :]) & \
               (lay[MIRROR][:, None] == lay[MIRROR][None, :])
        mx = same.sum(axis=1).astype(np.float64)
        C = float(np.count_nonzero(sizes))
        S2 = float(np.sum(mx / np.maximum(kx * krx, TINY)))
        pred = 2.0 * C - 2.0 * S2
        S1 = 0.0
        for a in range(F):
            fa = np.flatnonzero(lay == a)
            if fa.size == 0:
                continue
            for b in range(F):
                fb = np.flatnonzero(lay == b)
                if fb.size == 0:
                    continue
                mab = int(np.sum(lay[MIRROR[fa]] == b))
                if mab:
                    S1 += mab / (sizes[a] * sizes[b])
        hs2 = 2.0 * C * C - 2.0 * S2 * S2
        py = 0.0
        for qm in quads.values():
            a_ = D @ qm @ D
            b_ = qm - a_
            py = max(py, abs(float(np.sum(qm * qm))
                             - float(np.sum(a_ * a_)) - float(np.sum(b_ * b_))))
        out[name] = {"n_layers": int(C),
                     "commutator_R_fro": cR,
                     "commutator_T_fro": cT,
                     "S2_trace_DRDR": S2,
                     "S1_linear_weight_sum": float(S1),
                     "commutator_identity_pred_sq": pred,
                     "commutator_identity_abs_err": abs(cR * cR - pred),
                     "superoperator_commutator_HS_sq": hs2,
                     "superoperator_commutator_HS": math.sqrt(max(hs2, 0.0)),
                     "commutes_with_V4": bool(cR < 1e-12),
                     "quadrant_pythagoras_max_err": py}
    tau = np.full(FD, -1, dtype=np.int64)
    ok_tau = True
    for c in range(NC):
        u, v = int(layD[c]), int(layD[MIRROR[c]])
        if tau[u] >= 0 and tau[u] != v:
            ok_tau = False
        tau[u] = v
    nfix = int(np.sum(tau == np.arange(FD)))
    tr = {"id": FD * FD, "T": FD, "R": FD, "S": nfix}
    full, _ = klein_dims()
    dch, dker = {}, {}
    for eT in (+1, -1):
        for eR in (+1, -1):
            k = "W_" + ("+" if eT > 0 else "-") + ("+" if eR > 0 else "-")
            dch[k] = (tr["id"] + eT * tr["T"] + eR * tr["R"] + eT * eR * tr["S"]) // 4
            dker[k] = full[k] - dch[k]
    out["duplex_channel_traces"] = tr
    out["duplex_rho_well_defined"] = bool(ok_tau)
    out["n_duplex_channels"] = FD
    out["dims_full"] = {k: full[k] for k in dch}
    out["dims_duplex_image"] = dch
    out["dims_D8_kernel"] = dker
    out["dims_image_sum"] = int(sum(dch.values()))
    out["dims_kernel_sum"] = int(sum(dker.values()))
    return out


# =====================================================================
# 11. THEOREM I  one-sided four-term KL chain rule
# =====================================================================
def kl_four_way(B, E, aa):
    N1, N2 = float(B.sum()), float(E.sum())
    if N1 <= 0 or N2 <= 0:
        return None
    P = B / N1
    Q = np.maximum(E / N2, TINY)
    Q = Q / Q.sum()
    ar = aa[MIRROR]
    xs = np.arange(len(aa))

    def proj(M, r, c, nr, nc):
        return np.bincount((r[:, None] * nc + c[None, :]).ravel(),
                           weights=M.ravel(), minlength=nr * nc).reshape(nr, nc)

    Pab, Qab = proj(P, aa, aa, NAA, NAA), proj(Q, aa, aa, NAA, NAA)
    Pxb, Qxb = proj(P, xs, aa, NC, NAA), proj(Q, xs, aa, NC, NAA)
    key3 = (xs[:, None] * NAA * NAA + aa[None, :] * NAA + ar[None, :]).ravel()
    Pxbr = np.bincount(key3, weights=P.ravel(), minlength=NC * NAA * NAA).reshape(NC, NAA, NAA)
    Qxbr = np.bincount(key3, weights=Q.ravel(), minlength=NC * NAA * NAA).reshape(NC, NAA, NAA)
    D_AA = _kl_joint(Pab, np.maximum(Qab, TINY))
    den2 = np.where(Pab[aa, :] > 0, Qxb * Pab[aa, :] / np.maximum(Qab[aa, :], TINY), TINY)
    m2 = Pxb > 0
    D_d2b = float(np.sum(Pxb[m2] * np.log(Pxb[m2] / np.maximum(den2[m2], TINY))))
    Pxb3 = Pxb[:, :, None]
    Qxb3 = np.maximum(Qxb[:, :, None], TINY)
    den3 = np.where(Pxb3 > 0, Qxbr * Pxb3 / Qxb3, TINY)
    m3 = Pxbr > 0
    D_anti = float(np.sum(Pxbr[m3] * np.log(Pxbr[m3] / np.maximum(den3[m3], TINY))))
    Pxbr_xy = Pxbr[:, aa, ar]
    Qxbr_xy = np.maximum(Qxbr[:, aa, ar], TINY)
    den4 = np.where(Pxbr_xy > 0, Q * Pxbr_xy / Qxbr_xy, TINY)
    m4 = P > 0
    D_dup = float(np.sum(P[m4] * np.log(P[m4] / np.maximum(den4[m4], TINY))))
    total = _kl_joint(P, Q)
    ssum = D_AA + D_d2b + D_anti + D_dup
    return {"KL_total": total, "D_AA": D_AA, "D_d2b": D_d2b,
            "D_anti": D_anti, "D_dup": D_dup, "sum": ssum,
            "residual": total - ssum,
            "all_nonnegative": bool(min(D_AA, D_d2b, D_anti, D_dup) > -1e-12),
            "anti_share_clean": (D_anti / (D_anti + D_dup)) if (D_anti + D_dup) > 0 else None}


# =====================================================================
# 12. THEOREM N  symmetric two-sided duplex tower
# =====================================================================
def kl_symmetric_duplex_tower(B, E, aa):
    N1, N2 = float(B.sum()), float(E.sum())
    if N1 <= 0 or N2 <= 0:
        return None
    P = B / N1
    Q = np.maximum(E / N2, TINY)
    Q = Q / Q.sum()
    layD, FD, _, _ = duplex_layers(aa)

    def proj2d(M, lay, F):
        idx = (lay[:, None] * F + lay[None, :]).ravel()
        return np.bincount(idx, weights=M.ravel(), minlength=F * F).reshape(F, F)

    P_aa, Q_aa = proj2d(P, aa, NAA), proj2d(Q, aa, NAA)
    P_dup, Q_dup = proj2d(P, layD, FD), proj2d(Q, layD, FD)
    kl_tot = _kl_joint(P, Q)
    kl_dup = _kl_joint(P_dup, np.maximum(Q_dup, TINY))
    kl_aa = _kl_joint(P_aa, np.maximum(Q_aa, TINY))
    D_AA_sym = kl_aa
    D_anti_sym = kl_dup - kl_aa
    D_D8_sym = kl_tot - kl_dup
    ssum = D_AA_sym + D_anti_sym + D_D8_sym
    syn_clean = D_anti_sym + D_D8_sym
    return {"KL_total": kl_tot, "KL_pi_dup": kl_dup, "KL_pi_AA": kl_aa,
            "D_AA_sym": D_AA_sym, "D_anti_sym": D_anti_sym, "D_D8_sym": D_D8_sym,
            "sum": ssum, "residual": kl_tot - ssum,
            "all_nonnegative": bool(min(D_AA_sym, D_anti_sym, D_D8_sym) > -1e-12),
            "pi_AA_max_abs_diff": float(np.abs(P_aa - Q_aa).max()),
            "pi_dup_max_abs_diff": float(np.abs(P_dup - Q_dup).max()),
            "anti_share_sym_clean": (D_anti_sym / syn_clean) if syn_clean > 1e-15 else None}


# =====================================================================
# 13. THEOREM O  reading-frame phase locking  [FIX-33 + FIX-38]
# =====================================================================
def build_frame_cohorts(seqs, aa):
    """[FIX-38] build the three ALIGNED frame cohorts ONCE for all lags.
    [FIX-33] every frame reads exactly m_common = n_g - 1 triplets."""
    seqs_f = [[], [], []]
    trip = [0, 0, 0]
    stops = [0, 0, 0]
    stop_mask = (aa == STOP)
    for s in seqs:
        ng = int(s.size)
        if ng < 3:
            continue
        nt = np.empty(ng * 3, dtype=np.int64)
        nt[0::3] = s // 16
        nt[1::3] = (s // 4) % 4
        nt[2::3] = s % 4
        m_common = ng - 1
        span = 3 * m_common
        for phi in (0, 1, 2):
            sub = nt[phi:phi + span]
            if sub.size < span:
                continue
            cd = (16 * sub[0::3] + 4 * sub[1::3] + sub[2::3]).astype(np.int64)
            seqs_f[phi].append(cd)
            trip[phi] += int(m_common)
            stops[phi] += int(stop_mask[cd].sum())
    return {"seqs_f": seqs_f, "triplet_counts": trip, "stop_counts": stops,
            "n_genes_in": len(seqs),
            "n_genes_used": [len(x) for x in seqs_f]}


def theorem_O_from_cohorts(fc, k, aa, return_counts=False):
    """[THM-O] on CACHED frame cohorts; identical maths to the per-call version."""
    k = int(k)
    B_frames = [lag_counts(sq, NC, k) for sq in fc["seqs_f"]]
    if any(Bf.sum() <= 0 for Bf in B_frames):
        return None
    P_frames = [Bf / Bf.sum() for Bf in B_frames]
    P_bar = (P_frames[0] + P_frames[1] + P_frames[2]) / 3.0
    omega_phase = sum(_kl_joint(Pf, P_bar) for Pf in P_frames) / 3.0
    pX = [Pf.sum(axis=1) for Pf in P_frames]
    pX_bar = (pX[0] + pX[1] + pX[2]) / 3.0
    omega_marg = sum(_kl_joint(px, pX_bar) for px in pX) / 3.0
    omega_cond = 0.0
    for phi in (0, 1, 2):
        Pf, pxf = P_frames[phi], pX[phi]
        for x in np.flatnonzero(pxf > 0):
            omega_cond += pxf[x] * _kl_joint(Pf[x] / pxf[x],
                                             P_bar[x] / max(pX_bar[x], TINY))
    omega_cond /= 3.0
    frames = {}
    for phi in (0, 1, 2):
        I_cc, *_ = _mi_from_counts(B_frames[phi], mm=True)
        frames[f"frame_plus{phi}"] = {
            "I_CC_MM": I_cc, "S_syn": _s_syn_of_counts(B_frames[phi], aa),
            "n_triplets": fc["triplet_counts"][phi],
            "n_pairs": int(B_frames[phi].sum()),
            "stops_per_100_triplets": 100.0 * fc["stop_counts"][phi]
            / max(1, fc["triplet_counts"][phi])}
    s0 = frames["frame_plus0"]["S_syn"]
    N_tot = float(np.mean([Bf.sum() for Bf in B_frames]))
    q_eff = int(max(np.count_nonzero(Bf.sum(axis=1)) for Bf in B_frames))
    m_eff = int(max(np.count_nonzero(Bf) for Bf in B_frames))
    df_marg = 2 * (q_eff - 1)
    df_cond = 2 * q_eff * (q_eff - 1)
    df_full = 2 * (m_eff - 1)
    Bu = B_frames[0] + B_frames[1] + B_frames[2]
    q_u = int(np.count_nonzero(Bu.sum(axis=1) > 0))
    m_u = int(np.count_nonzero(Bu > 0))
    df_marg_eff = 2 * max(q_u - 1, 0)
    df_cond_eff = 2 * (m_u - q_u)
    df_full_eff = 2 * max(m_u - 1, 0)
    df_add_err = abs(df_full_eff - (df_marg_eff + df_cond_eff))
    den = 6.0 * max(N_tot, 1.0)
    npair = [int(Bf.sum()) for Bf in B_frames]
    res = {"lag": k,
           "Omega_phase_JS_nats": float(omega_phase),
           "Omega_marg_JS_nats": float(omega_marg),
           "Omega_cond_JS_nats": float(omega_cond),
           "chain_rule_abs_err": abs(omega_phase - (omega_marg + omega_cond)),
           "cond_nonnegative": bool(omega_cond >= -1e-12),
           "frames": frames, "frame_pair_counts": npair,
           "frame_pair_counts_equal": bool(npair[0] == npair[1] == npair[2]),
           "frame0_max_S_syn": bool(s0 > frames["frame_plus1"]["S_syn"]
                                    and s0 > frames["frame_plus2"]["S_syn"]),
           "frame0_stops_zero_by_construction": True,
           "N_tot_mean": N_tot, "q_eff": q_eff, "m_eff_cells": m_eff,
           "df_phase_marg": df_marg, "df_phase_cond": df_cond, "df_phase_full": df_full,
           "df_phase_marg_eff": df_marg_eff, "df_phase_cond_eff": df_cond_eff,
           "df_phase_full_eff": df_full_eff, "df_phase_additivity_abs_err": df_add_err,
           "q_union": q_u, "m_union": m_u,
           "df_phase_full_support_additivity_ok": bool(df_full == df_marg + df_cond),
           "bias_marg": df_marg / den, "bias_cond": df_cond / den,
           "bias_phase": df_full / den,
           "G2_phase_cond": 6.0 * N_tot * float(omega_cond),
           "omega_cond_above_bias": bool(omega_cond > df_cond / den),
           "bias_caveat": ("chi-square df is a BIAS SCALE only; frames are dependent "
                           "overlapping readings, not independent samples"),
           "note": ("frame-0 internal stops are 0 by parser construction; only "
                    "frames +1/+2 are informative. frame0_max_S_syn is an empirical "
                    "prediction, not a theorem. Omega_cond>0 alone does not exclude "
                    "phase-dependent dinucleotide inertia -- see Theorems R, S, T.")}
    if return_counts:
        res["_B_frames"] = B_frames
    return res


def theorem_O_multi_lag(seqs, lags, aa):
    """[FIX-38] ONE frame build, many lags."""
    fc = build_frame_cohorts(seqs, aa)
    out = []
    for k in sorted({int(x) for x in lags}):
        o = theorem_O_from_cohorts(fc, k, aa)
        if o is not None:
            out.append(o)
    return out


def theorem_O_frame_phase_audit(seqs, k, aa, return_counts=False):
    """Backward-compatible single-lag wrapper over the cached-cohort engine."""
    return theorem_O_from_cohorts(build_frame_cohorts(seqs, aa), int(k), aa,
                                  return_counts=return_counts)


def de_bruijn(k, n):
    """Deterministic de Bruijn sequence B(k,n) (FKMM Lyndon-word construction)."""
    a = [0] * (k * n)
    seq = []

    def db(t, p):
        if t > n:
            if n % p == 0:
                seq.extend(a[1:p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)

    db(1, 1)
    return seq


def debruijn_phase_null(n_periods=3, order=6):
    db = de_bruijn(4, order)
    nt = np.array(db * n_periods + db[:order], dtype=np.int64)
    m = nt.size // 3
    cd = (16 * nt[0:3 * m:3] + 4 * nt[1:3 * m:3] + nt[2:3 * m:3]).astype(np.int64)
    return cd, {"period_len": 4 ** order, "n_periods": n_periods,
                "period_len_mod3": (4 ** order) % 3, "n_nt": int(nt.size),
                "n_codons": int(m)}


# =====================================================================
# 14. THEOREM R / S / T / T.5 / U  order-2 nt propagator, screens, splits
# =====================================================================
def _nt_order2_counts(seqs):
    """Naive reference implementation (per-gene nucleotide unpacking)."""
    C_all = np.zeros((4, 4, 4), dtype=np.float64)
    C_ph = np.zeros((3, 4, 4, 4), dtype=np.float64)
    for s in seqs:
        nt = np.empty(s.size * 3, dtype=np.int64)
        nt[0::3] = s // 16
        nt[1::3] = (s // 4) % 4
        nt[2::3] = s % 4
        if nt.size < 3:
            continue
        a, b, c = nt[:-2], nt[1:-1], nt[2:]
        np.add.at(C_all, (a, b, c), 1.0)
        np.add.at(C_ph, ((np.arange(a.size) + 2) % 3, a, b, c), 1.0)
    return C_all, C_ph


def _nt_order2_counts_projected(seqs):
    """[FIX-37] EXACT linear projection: no nucleotide unpacking at all.
    C_ph[2] = codon marginal                (n_g triples)
    C_ph[0] = B_1 projected over x1,(y2,y3) (n_g-1)
    C_ph[1] = B_1 projected over (x1,x2),y3 (n_g-1)
    Total n + 2(n-1) = 3n-2 == nt.size-2, so boundaries match EXACTLY."""
    marg = np.zeros(NC, dtype=np.float64)
    B1 = np.zeros((NC, NC), dtype=np.float64)
    for s in seqs:
        marg += np.bincount(s, minlength=NC)
        if s.size >= 2:
            B1 += np.bincount(s[:-1] * NC + s[1:], minlength=NC * NC).reshape(NC, NC)
    C_ph = np.zeros((3, 4, 4, 4), dtype=np.float64)
    C_ph[2] = marg.reshape(4, 4, 4)
    M0 = np.bincount((SUF16[:, None] * 4 + C_N1[None, :]).ravel(),
                     weights=B1.ravel(), minlength=16 * 4).reshape(16, 4)
    C_ph[0] = M0.reshape(4, 4, 4)
    M1 = np.bincount((C_N3[:, None] * 16 + OUT16[None, :]).ravel(),
                     weights=B1.ravel(), minlength=4 * 16).reshape(4, 16)
    C_ph[1] = M1.reshape(4, 4, 4)
    return C_ph.sum(axis=0), C_ph


def _cond3(C, alpha):
    P = np.asarray(C, dtype=np.float64) + alpha
    return P / np.maximum(P.sum(axis=-1, keepdims=True), TINY)


def _to_T(P):
    T = np.zeros((16, 16), dtype=np.float64)
    for u in range(4):
        for v in range(4):
            for w in range(4):
                T[4 * u + v, 4 * v + w] = P[u, v, w]
    return T


def build_h3b_lag_kernels(train_seqs, lags, alpha=0.1, projected=True):
    """[THM-R/S/T] EXACT semigroup Q_1^k == Q_k, computed only for `lags`.
    Q_k       = (T_nt^(3k-1))[SUF16(x), OUT16(y)] * P_nt(y3|y1,y2)
    Q_k^phase = (M^(k-1) T0 T1)[SUF16(x), OUT16(y)] * P2(y3|y1,y2)
    [THM-T] A_bridge = T^(3k-1) is the 16x16 screen; R_wobble = P_nt(y3|o).
    [THM-T.5] R_wobble_phase == P_train(Y3|Y1,Y2) exactly (FIX-37)."""
    lags = sorted({int(k) for k in lags if int(k) >= 1})
    if projected:
        C_all, C_ph = _nt_order2_counts_projected(train_seqs)
    else:
        C_all, C_ph = _nt_order2_counts(train_seqs)
    P_all = _cond3(C_all, alpha)
    P_ph = [_cond3(C_ph[i], alpha) for i in range(3)]
    T_all = _to_T(P_all)
    T0, T1, T2 = (_to_T(P_ph[i]) for i in range(3))
    M_codon = T0 @ T1 @ T2
    last = P_all[C_N1, C_N2, C_N3]
    last2 = P_ph[2][C_N1, C_N2, C_N3]
    R_wobble = np.zeros((16, 4), dtype=np.float64)
    R_wobble2 = np.zeros((16, 4), dtype=np.float64)
    for o in range(16):
        R_wobble[o] = P_all[o // 4, o % 4, :]
        R_wobble2[o] = P_ph[2][o // 4, o % 4, :]
    A1 = T_all @ T_all
    A1p = T0 @ T1
    Q1 = A1[np.ix_(SUF16, OUT16)] * last[None, :]
    Q1p = A1p[np.ix_(SUF16, OUT16)] * last2[None, :]
    out = {}
    for k in lags:
        Ak = np.linalg.matrix_power(T_all, 3 * k - 1)
        Q = Ak[np.ix_(SUF16, OUT16)] * last[None, :]
        Ap = np.linalg.matrix_power(M_codon, k - 1) @ T0 @ T1
        Qp = Ap[np.ix_(SUF16, OUT16)] * last2[None, :]
        sg = float(np.abs(np.linalg.matrix_power(Q1, k) - Q).max())
        sgp = float(np.abs(np.linalg.matrix_power(Q1p, k) - Qp).max())
        rowvar = 0.0
        for s_ in range(16):
            idx = np.flatnonzero(SUF16 == s_)
            if idx.size > 1:
                rowvar = max(rowvar, float(np.abs(Q[idx] - Q[idx[0]]).max()))
        out[k] = {"Q": Q, "Q_phase": Qp,
                  "A_bridge": Ak, "A_bridge_phase": Ap,
                  "row_sum_err": float(np.abs(Q.sum(axis=1) - 1.0).max()),
                  "row_sum_err_phase": float(np.abs(Qp.sum(axis=1) - 1.0).max()),
                  "bridge_row_sum_err": float(np.abs(Ak.sum(axis=1) - 1.0).max()),
                  "bridge_row_sum_err_phase": float(np.abs(Ap.sum(axis=1) - 1.0).max()),
                  "semigroup_closure_err": sg,
                  "semigroup_closure_err_phase": sgp,
                  "suffix_row_max_dev": rowvar,
                  "rank_Q": int(np.linalg.matrix_rank(Q, tol=1e-10)),
                  "rank_Q_phase": int(np.linalg.matrix_rank(Qp, tol=1e-10))}
    ev1 = np.linalg.eigvals(Q1)
    meta = {"T_all": T_all, "T0": T0, "T1": T1, "T2": T2,
            "P_all": P_all, "P_ph": P_ph, "M_codon": M_codon,
            "R_wobble": R_wobble, "R_wobble_phase": R_wobble2,
            "Q1": Q1, "Q1_phase": Q1p, "lags": lags,
            "n_zero_eigenvalues_Q1": int(np.sum(np.abs(ev1) < 1e-10)),
            "rank_Q1": int(np.linalg.matrix_rank(Q1, tol=1e-10)),
            "model_free_params_homogeneous": 16 * 3,
            "model_free_params_phase": 3 * 16 * 3,
            "unlinked_screen_free_params": 16 * 15 + 16 * 3,
            "full_kernel_free_params": NC * (NC - 1),
            "param_note": ("240+48=288 is the UNLINKED screen family {A[s,o]*R(y3|o)} "
                           "on the FULL torus; H3b is a 48-dim (homogeneous) / 144-dim "
                           "(phase) submanifold because A_k and R share the SAME P_nt. "
                           "In the SENSE convention R has sum_o (C_o-1) = q_s-16 free "
                           "parameters (45 for Table 1), so 675+675+2025+240+45 = 3660 "
                           "= 61*60.")}
    return out, meta


def two_sided_screen_split(B, aa=None):
    N = float(B.sum())
    if N <= 0:
        return None
    if aa is None:
        aa = AA_TABLE1
    bp = box_partition(aa)
    b4c = bp["prefix_codon_mask"]
    s4c = bp["suffix_codon_mask"]
    xs = np.arange(NC)
    nx, ny = B.sum(axis=1), B.sum(axis=0)
    Bso = np.bincount((SUF16[:, None] * 16 + OUT16[None, :]).ravel(),
                      weights=B.ravel(), minlength=256).reshape(16, 16)
    Bxo = np.bincount((xs[:, None] * 16 + OUT16[None, :]).ravel(),
                      weights=B.ravel(), minlength=NC * 16).reshape(NC, 16)
    Bsy = np.bincount((SUF16[:, None] * NC + xs[None, :]).ravel(),
                      weights=B.ravel(), minlength=16 * NC).reshape(16, NC)
    ns, no = Bso.sum(axis=1), Bso.sum(axis=0)
    I_cc = _mi_tab(B, np.maximum(nx[:, None] * ny[None, :], TINY) / N, N)
    I_br = _mi_tab(Bso, np.maximum(ns[:, None] * no[None, :], TINY) / N, N)
    den_lead = np.maximum(nx[:, None] * Bso[SUF16, :], TINY) / np.maximum(ns[SUF16][:, None], TINY)
    t_lead = _mi_terms(Bxo, den_lead, N)
    D_lead = float(t_lead.sum())
    den_wob = np.maximum(Bso[:, OUT16] * ny[None, :], TINY) / np.maximum(no[OUT16][None, :], TINY)
    t_wob = _mi_terms(Bsy, den_wob, N)
    D_wob = float(t_wob.sum())
    D_13 = _mi_tab(B, np.maximum(Bxo[:, OUT16] * Bsy[SUF16, :], TINY) / np.maximum(Bso[np.ix_(SUF16, OUT16)], TINY), N)
    D_lead_S4 = float(t_lead[s4c, :].sum())
    D_lead_Ssplit = float(t_lead[~s4c, :].sum())
    D_wob_B4 = float(t_wob[:, b4c].sum())
    D_wob_Bsplit = float(t_wob[:, ~b4c].sum())
    df_br = max(0, (int(np.count_nonzero(ns > 0)) - 1)) * max(0, (int(np.count_nonzero(no > 0)) - 1))
    df_lead = 0
    for s_ in range(16):
        if ns[s_] <= 0:
            continue
        R = int(np.count_nonzero(Bxo[SUF16 == s_, :].sum(axis=1) > 0))
        Cc = int(np.count_nonzero(Bso[s_, :] > 0))
        if R > 1 and Cc > 1:
            df_lead += (R - 1) * (Cc - 1)
    df_wob = 0
    for o_ in range(16):
        if no[o_] <= 0:
            continue
        R = int(np.count_nonzero(Bso[:, o_] > 0))
        Cc = int(np.count_nonzero(ny[OUT16 == o_] > 0))
        if R > 1 and Cc > 1:
            df_wob += (R - 1) * (Cc - 1)
    df_13 = 0
    for s_ in range(16):
        for o_ in range(16):
            if Bso[s_, o_] <= 0:
                continue
            R = int(np.count_nonzero(Bxo[(SUF16 == s_), o_] > 0))
            Cc = int(np.count_nonzero(Bsy[s_, (OUT16 == o_)] > 0))
            if R > 1 and Cc > 1:
                df_13 += (R - 1) * (Cc - 1)
    df_full = screen_df_bounds(None)
    df_sense = screen_df_bounds(aa != STOP)
    return {"N": N, "I_CC": I_cc,
            "I_bridge": I_br, "D_5prime_lead": D_lead,
            "D_3prime_wobble": D_wob, "D_1_3": D_13,
            "four_term_residual": I_cc - (I_br + D_lead + D_wob + D_13),
            "Delta_x1": D_lead + D_13, "I_suffix16": I_br + D_wob,
            "all_nonnegative": bool(min(I_br, D_lead, D_wob, D_13) > -1e-12),
            "D_5prime_lead_S4": D_lead_S4,
            "D_5prime_lead_Ssplit": D_lead_Ssplit,
            "D_3prime_wobble_B4": D_wob_B4,
            "D_3prime_wobble_Bsplit": D_wob_Bsplit,
            "box_split_lead_residual": D_lead - (D_lead_S4 + D_lead_Ssplit),
            "box_split_wobble_residual": D_wob - (D_wob_B4 + D_wob_Bsplit),
            "df_full": df_full,
            "df_sense": df_sense,
            "df_eff": {"bridge": df_br, "lead": df_lead, "wobble": df_wob,
                       "d13": df_13, "sum": df_br + df_lead + df_wob + df_13},
            "H3b_visible_df_fraction_full": 225.0 / 3969.0,
            "H3b_visible_df_fraction_sense": 225.0 / 3600.0,
            "Bso": Bso, "Bxo": Bxo, "Bsy": Bsy,
            "ns": ns, "no": no, "nx": nx, "ny": ny}


def screen_antiparallel_dual(B, aa=None):
    """[THM-T.3 / THM-U] under (X,Y)->(rho(Y),rho(X)):
    lead and wobble SWAP (also their box parts and their effective df),
    while bridge and D_1<->3 are invariant."""
    Ba = B.T[np.ix_(MIRROR, MIRROR)]
    f = two_sided_screen_split(B, aa)
    a = two_sided_screen_split(Ba, aa)
    if f is None or a is None:
        return None
    keys = ("I_bridge", "D_5prime_lead", "D_3prime_wobble", "D_1_3",
            "D_5prime_lead_S4", "D_5prime_lead_Ssplit",
            "D_3prime_wobble_B4", "D_3prime_wobble_Bsplit")
    fe, ae = f["df_eff"], a["df_eff"]
    return {"bridge_abs_err": abs(f["I_bridge"] - a["I_bridge"]),
            "d13_abs_err": abs(f["D_1_3"] - a["D_1_3"]),
            "lead_vs_wobble_abs_err": abs(f["D_5prime_lead"] - a["D_3prime_wobble"]),
            "wobble_vs_lead_abs_err": abs(f["D_3prime_wobble"] - a["D_5prime_lead"]),
            "I_CC_abs_err": abs(f["I_CC"] - a["I_CC"]),
            "U_leadS4_vs_wobB4_abs_err": abs(f["D_5prime_lead_S4"] - a["D_3prime_wobble_B4"]),
            "U_leadSsplit_vs_wobBsplit_abs_err":
                abs(f["D_5prime_lead_Ssplit"] - a["D_3prime_wobble_Bsplit"]),
            "U_wobB4_vs_leadS4_abs_err": abs(f["D_3prime_wobble_B4"] - a["D_5prime_lead_S4"]),
            "U_wobBsplit_vs_leadSsplit_abs_err":
                abs(f["D_3prime_wobble_Bsplit"] - a["D_5prime_lead_Ssplit"]),
            "df_bridge_abs_err": int(abs(fe["bridge"] - ae["bridge"])),
            "df_d13_abs_err": int(abs(fe["d13"] - ae["d13"])),
            "df_lead_vs_wobble_abs_err": int(abs(fe["lead"] - ae["wobble"])),
            "df_wobble_vs_lead_abs_err": int(abs(fe["wobble"] - ae["lead"])),
            "sense": {k: f[k] for k in keys},
            "antisense": {k: a[k] for k in keys}}


def suffix16_split(B, aa=None):
    """[THM-S] two-term split (retained; Thm T refines it into four terms).
    [FIX-40] sense-convention df now EXACT: (q_s-16)*(q_s-1) = 45*60 = 2700."""
    sp = two_sided_screen_split(B, aa)
    if sp is None:
        return None
    dff, dfs = sp["df_full"], sp["df_sense"]
    return {"I_CC": sp["I_CC"], "I_suffix16": sp["I_suffix16"],
            "delta_x1": sp["Delta_x1"],
            "delta_x1_nonnegative": bool(sp["Delta_x1"] >= -1e-12),
            "df_delta_x1_full_torus": dff["delta_x1"],
            "df_delta_x1_sense": dfs["delta_x1"],
            "df_suffix16_full_torus": dff["suffix16"],
            "df_suffix16_sense": dfs["suffix16"],
            "df_sense_partition": dfs,
            "df_delta_x1_eff": sp["df_eff"]["lead"] + sp["df_eff"]["d13"],
            "n_x": sp["nx"], "n_s": sp["ns"], "B_s": sp["Bsy"], "N": sp["N"]}


def h3b_lag_profile(train_seqs, test_seqs, lags, alpha=0.1, aa=None):
    """[THM-R/S/T/T.5/U + FIX-37/39/40] only the REQUESTED lags are evaluated.
    Reports per lag:
      * held-out d3b / d3b_phase (test counts under the train kernel);
      * [FIX-39] EXACT three-term identity
            d3b = Delta_x1_test + D_suffix_test - D_gen;
      * [THM-T] four-term two-sided screen split of I_test + antiparallel dual;
      * [THM-T.4] three-term split D_suffix = D_3'wobble + D_bridge + D_wobble;
      * [THM-T.5] EXACT five-term identity
            d3b = D_5'lead + D_3'wobble + D_1<->3 + D_bridge + D_wobble - D_gen
        and the EXACT phase cancellation
            d3b - d3b_phase = (D_bridge - D_bridge_ph) + (D_wobble - D_wobble_ph);
      * [THM-U] additive 8+8 box decomposition of both wobble layers;
      * the v1.9 in-sample model-joint identity, retained as a cross-check.
    """
    lags = sorted({int(k) for k in lags if int(k) >= 1})
    if aa is None:
        aa = AA_TABLE1
    ker, meta = build_h3b_lag_kernels(train_seqs, lags, alpha)
    Rw = meta["R_wobble"]
    Rw2 = meta["R_wobble_phase"]
    rep = np.array([np.flatnonzero(SUF16 == s_)[0] for s_ in range(16)])
    rows = []
    for k in lags:
        Btr = lag_counts(train_seqs, NC, k)
        Bte = lag_counts(test_seqs, NC, k)
        Ntr, Nte = float(Btr.sum()), float(Bte.sum())
        if Nte <= 0 or Ntr <= 0:
            continue
        Ptr = (Btr + alpha) / (Btr.sum(axis=1, keepdims=True) + NC * alpha)
        lp = np.log(np.maximum(Ptr, TINY))
        Q, Qp = ker[k]["Q"], ker[k]["Q_phase"]
        Ak, Akp = ker[k]["A_bridge"], ker[k]["A_bridge_phase"]
        lq, lqp = np.log(np.maximum(Q, TINY)), np.log(np.maximum(Qp, TINY))
        d3b = float(np.sum(Bte * (lp - lq)) / Nte)
        d3bp = float(np.sum(Bte * (lp - lqp)) / Nte)
        # ---- [THM-T] four-term screen split on TEST + antiparallel dual ----
        sp_te = two_sided_screen_split(Bte, aa)
        dual = screen_antiparallel_dual(Bte, aa)
        # ---- [FIX-39] EXACT held-out three-term identity ----
        nx_te, ns_te, no_te = sp_te["nx"], sp_te["ns"], sp_te["no"]
        Pte = np.zeros((NC, NC))
        okx = nx_te > 0
        Pte[okx] = Bte[okx] / nx_te[okx, None]
        D_gen = 0.0
        for x in np.flatnonzero(okx):
            D_gen += (nx_te[x] / Nte) * _kl_joint(Pte[x], Ptr[x])
        Psuf = np.zeros((16, NC))
        oks = ns_te > 0
        Psuf[oks] = sp_te["Bsy"][oks] / ns_te[oks, None]
        Qrow = Q[rep]
        D_suf = 0.0
        for s_ in np.flatnonzero(oks):
            D_suf += (ns_te[s_] / Nte) * _kl_joint(Psuf[s_], Qrow[s_])
        ho_res = abs(d3b - (sp_te["Delta_x1"] + D_suf - D_gen))
        # ---- [THM-T.4] three-term split of D_suffix ----
        Pbr = np.zeros((16, 16))
        Pbr[oks] = sp_te["Bso"][oks] / ns_te[oks, None]
        D_bridge = 0.0
        for s_ in np.flatnonzero(oks):
            D_bridge += (ns_te[s_] / Nte) * _kl_joint(Pbr[s_], Ak[s_])
        oko = no_te > 0
        Pw = np.zeros((16, 4))
        for o_ in np.flatnonzero(oko):
            v = sp_te["ny"][OUT16 == o_]
            Pw[o_] = v / max(v.sum(), TINY)
        D_wob = 0.0
        for o_ in np.flatnonzero(oko):
            D_wob += (no_te[o_] / Nte) * _kl_joint(Pw[o_], Rw[o_])
        t4_res = abs(D_suf - (sp_te["D_3prime_wobble"] + D_bridge + D_wob))
        # ---- [THM-T.5] phase-kernel KL terms, five-term identity, cancellation ----
        D_bridge_p = 0.0
        for s_ in np.flatnonzero(oks):
            D_bridge_p += (ns_te[s_] / Nte) * _kl_joint(Pbr[s_], Akp[s_])
        D_wob_p = 0.0
        for o_ in np.flatnonzero(oko):
            D_wob_p += (no_te[o_] / Nte) * _kl_joint(Pw[o_], Rw2[o_])
        Rw2k = wobble_kernel_from_positions(train_seqs, k, alpha)
        D_wob_p_al = 0.0
        for o_ in np.flatnonzero(oko):
            D_wob_p_al += (no_te[o_] / Nte) * _kl_joint(Pw[o_], Rw2k[o_])
        D_wob_edge = D_wob_p - D_wob_p_al
        D_wob_shift = 0.0
        for o_ in np.flatnonzero(oko):
            D_wob_shift += (no_te[o_] / Nte) * _kl_joint(Rw2k[o_], Rw2[o_])
        D_suf_p = sp_te["D_3prime_wobble"] + D_bridge_p + D_wob_p
        nonbridge = sp_te["D_5prime_lead"] + sp_te["D_3prime_wobble"] + sp_te["D_1_3"]
        t5_res = abs(d3b - (nonbridge + D_bridge + D_wob - D_gen))
        t5p_res = abs(d3bp - (nonbridge + D_bridge_p + D_wob_p - D_gen))
        t5d_res = abs((d3b - d3bp) - ((D_bridge - D_bridge_p) + (D_wob - D_wob_p)))
        ho_p_res = abs(d3bp - (sp_te["Delta_x1"] + D_suf_p - D_gen))
        V23 = thm_V_screen_split(Bte, aa)
        W23 = thm_W_wobble_split(Pw, Rw, no_te, Nte, aa)
        W23p = thm_W_wobble_split(Pw, Rw2, no_te, Nte, aa)
        Y23 = thm_Y_diamond(Bte, aa, sense=True)
        Z23 = thm_Z_tower(Bte, aa, sense=True)
        AA24 = thm_AA_matrix(Bte, aa, sense=True)
        AAd24 = thm_AA_d2b(Bte, aa, sense=True)
        # ---- retained v1.9 in-sample model-joint cross-check ----
        nx_tr = Btr.sum(axis=1)
        Bmodel = nx_tr[:, None] * Ptr
        sp_md = two_sided_screen_split(Bmodel, aa)
        ns_md = sp_md["ns"]
        Pm = np.zeros((16, NC))
        okm = ns_md > 0
        Pm[okm] = sp_md["Bsy"][okm] / ns_md[okm, None]
        kl_md = 0.0
        for s_ in np.flatnonzero(okm):
            kl_md += (ns_md[s_] / Ntr) * _kl_joint(Pm[s_], Qrow[s_])
        d3b_plug = float(np.sum(Bmodel * (lp - lq)) / Ntr)
        rows.append({"lag": k, "N_test": int(Nte), "N_train": int(Ntr),
                     "d3b": d3b, "d3b_phase": d3bp,
                     "row_sum_err": ker[k]["row_sum_err"],
                     "row_sum_err_phase": ker[k]["row_sum_err_phase"],
                     "bridge_row_sum_err": ker[k]["bridge_row_sum_err"],
                     "bridge_row_sum_err_phase": ker[k]["bridge_row_sum_err_phase"],
                     "semigroup_closure_err_k": ker[k]["semigroup_closure_err"],
                     "semigroup_closure_err_phase_k": ker[k]["semigroup_closure_err_phase"],
                     "rank_Q": ker[k]["rank_Q"], "rank_Q_phase": ker[k]["rank_Q_phase"],
                     "suffix_row_max_dev": ker[k]["suffix_row_max_dev"],
                     "mm_bias_scale": (NC * NC - 1) / (2.0 * Nte),
                     # THM-T
                     "I_CC_test": sp_te["I_CC"],
                     "I_bridge": sp_te["I_bridge"],
                     "D_5prime_lead": sp_te["D_5prime_lead"],
                     "D_3prime_wobble": sp_te["D_3prime_wobble"],
                     "D_1_3": sp_te["D_1_3"],
                     "thmT_four_term_abs_err": abs(sp_te["four_term_residual"]),
                     "Delta_x1_test": sp_te["Delta_x1"],
                     "I_suffix16_test": sp_te["I_suffix16"],
                     "I_nonbridge_test": nonbridge,
                     "df_eff": sp_te["df_eff"], "df_full": sp_te["df_full"],
                     "df_sense": sp_te["df_sense"],
                     "antiparallel_dual": dual,
                     # THM-U
                     "D_5prime_lead_S4": sp_te["D_5prime_lead_S4"],
                     "D_5prime_lead_Ssplit": sp_te["D_5prime_lead_Ssplit"],
                     "D_3prime_wobble_B4": sp_te["D_3prime_wobble_B4"],
                     "D_3prime_wobble_Bsplit": sp_te["D_3prime_wobble_Bsplit"],
                     "thmU_box_split_lead_abs_err": abs(sp_te["box_split_lead_residual"]),
                     "thmU_box_split_wobble_abs_err": abs(sp_te["box_split_wobble_residual"]),
                     # FIX-39 exact held-out identity
                     "D_suffix_test": D_suf, "D_gen_train_test": D_gen,
                     "thmT_heldout_identity_abs_err": ho_res,
                     # THM-T.4
                     "D_bridge_kl": D_bridge, "D_wobble_kl": D_wob,
                     "thmT_suffix_split_abs_err": t4_res,
                     # THM-T.5
                     "D_suffix_test_phase": D_suf_p,
                     "D_bridge_kl_phase": D_bridge_p, "D_wobble_kl_phase": D_wob_p,
                     "D_wobble_kl_phase_aligned": D_wob_p_al,
                     "D_wobble_phase_5prime_edge": D_wob_edge,
                     "thmT5_five_term_abs_err": t5_res,
                     "thmT5_phase_five_term_abs_err": t5p_res,
                     "thmT5_phase_diff_abs_err": t5d_res,
                     "thmT_heldout_identity_phase_abs_err": ho_p_res,
                     "D_wobble_kernel_shift_nonneg": D_wob_shift,
                     "thmV": V23, "thmW": W23, "thmW_phase": W23p,
                     "thmY": Y23, "thmZ": Z23,
                     "thmAA": AA24, "thmAA_d2b": AAd24,
                     # retained in-sample
                     "d3b_plugin_insample": d3b_plug,
                     "delta_x1_insample": sp_md["Delta_x1"],
                     "kl_suffix_insample": kl_md,
                     "thmS_identity_abs_err": abs(d3b_plug - (sp_md["Delta_x1"] + kl_md)),
                     "identity_note": ("held-out: d3b = Delta_x1_test + D_suffix_test "
                                       "- D_gen (EXACT); five-term: d3b = D_5'lead + "
                                       "D_3'wobble + D_1<->3 + D_bridge_kl + "
                                       "D_wobble_kl - D_gen (EXACT). in-sample: model "
                                       "joint n_x*P_train(y|x) (EXACT, independent check).")})
    return rows, meta


# =====================================================================
# 15. THEOREM B  memory ladder
# =====================================================================
def _tstat(d):
    if d.size < 2:
        return (float("nan"),) * 4
    se = float(d.std(ddof=1) / math.sqrt(d.size))
    mean = float(d.mean())
    t = mean / se if se > 0 else float("nan")
    pv = float(ttest_1samp(d, 0.0).pvalue) if (HAVE_SCIPY and d.size > 2) else float("nan")
    return mean, se, t, pv


def memory_ladder(train_seqs, test_seqs, q, max_order, kappa,
                  min_ctx_count, max_ctx, gene_gc3=None):
    models = build_models(train_seqs, q, max_order, kappa=kappa,
                          min_ctx_count=min_ctx_count, max_ctx=max_ctx)
    es = max_order
    idx = [i for i, s in enumerate(test_seqs) if s.size > es + 1]
    if not idx:
        return None
    usable = [test_seqs[i] for i in idx]
    gc3 = None if gene_gc3 is None else np.asarray([gene_gc3[i] for i in idx], float)
    ll = np.zeros((max_order + 1, len(usable)))
    cov = np.zeros_like(ll)
    nsym = np.zeros(len(usable), dtype=np.int64)
    for gi, s in enumerate(usable):
        n = s.size - es
        nsym[gi] = n
        for p in range(max_order + 1):
            lp, used = seq_logprob(models, s, p, es, q, return_order=True)
            ll[p, gi] = float(lp.sum()) / n
            cov[p, gi] = float(np.mean(used == p)) if p > 0 else 1.0
    ladder = []
    for p in range(max_order + 1):
        row = {"order": p,
               "ll_per_symbol_pooled": float(np.average(ll[p], weights=nsym)),
               "gene_mean_ll": float(ll[p].mean()),
               "n_contexts": int(models[p].n_ctx),
               "frac_symbols_at_full_order": float(np.average(cov[p], weights=nsym))}
        if p >= 1:
            m, se, t, pv = _tstat(ll[p] - ll[p - 1])
            row.update({"M_hat": m, "gene_SE": se, "t_stat": t, "p_value_two_sided": pv,
                        "ci95_low": m - 1.96 * se if _finite(se) else None,
                        "ci95_high": m + 1.96 * se if _finite(se) else None,
                        "frac_genes_positive": float(np.mean(ll[p] > ll[p - 1]))})
            mc, sc, tc, pc = _tstat(ll[p] - ll[1])
            row.update({"cum_gain_over_order1": mc, "cum_SE": sc, "cum_t": tc,
                        "cum_p_value": pc,
                        "cum_ci95_low": mc - 1.96 * sc if _finite(sc) else None,
                        "cum_ci95_high": mc + 1.96 * sc if _finite(sc) else None})
        ladder.append(row)
    lls = [r["ll_per_symbol_pooled"] for r in ladder]
    p_arg = int(np.argmax(lls))
    p_sig = 0
    for p in range(1, max_order + 1):
        lo = ladder[p].get("ci95_low")
        if _finite(lo) and lo > 0:
            p_sig = p
        else:
            break
    max_range = 0
    for p in range(1, max_order + 1):
        lo = ladder[p].get("cum_ci95_low")
        if _finite(lo) and lo > 0:
            max_range = p
    comp = None
    if gc3 is not None and p_arg >= 2 and gc3.size > 3:
        g = ll[p_arg] - ll[1]
        if float(np.std(g)) > 0 and float(np.std(gc3)) > 0:
            comp = {"pearson_r_gain_vs_gc3": float(np.corrcoef(g, gc3)[0, 1]),
                    "gene_gc3_mean": float(gc3.mean()),
                    "gene_gc3_sd": float(gc3.std(ddof=1))}
    return {"q": int(q), "max_order": int(max_order), "kappa": float(kappa),
            "n_test_genes": len(usable), "n_eval_symbols": int(nsym.sum()),
            "ladder": ladder, "p_star_argmax": p_arg, "p_star_significant": p_sig,
            "max_significant_range": max_range,
            "memory_beyond_first_order_nats": float(lls[p_arg] - lls[1]) if max_order >= 1 else 0.0,
            "coverage_at_p_star": ladder[p_arg]["frac_symbols_at_full_order"],
            "composition_association": comp}


# =====================================================================
# 16. THEOREM D  exact zero-RNG H5 / H6 / H7 on a COMMON cohort
# =====================================================================
def build_h2b_sampler_from_counts(B, fibers, alpha=0.1):
    bg = np.asarray(B, dtype=np.float64) + alpha
    Q = np.zeros((NC, NAA, NC), dtype=np.float64)
    for a, f in enumerate(fibers):
        if f.size == 0:
            continue
        Q[:, a, f] = bg[:, f] / np.maximum(bg[:, f].sum(axis=1, keepdims=True), TINY)
    return Q


def build_h2b_sampler(train_seqs, aa, fibers, alpha=0.1):
    return build_h2b_sampler_from_counts(bigram_counts(train_seqs, NC), fibers, alpha)


def h5_pair_tables(Q, logP, fibers):
    Qt, Rt = {}, {}
    for a, fa in enumerate(fibers):
        if fa.size == 0:
            continue
        for b, fb in enumerate(fibers):
            if fb.size == 0:
                continue
            M = Q[np.ix_(fa, [b], fb)][:, 0, :]
            Qt[(a, b)] = M
            Rt[(a, b)] = (M * logP[np.ix_(fa, fb)]).sum(axis=1)
    return Qt, Rt


def h5_analytic_gene(s, aa, fibers, Qt, Rt, mass_tol=MASS_TOL):
    n = s.size
    if n < 2:
        return None, 0.0
    a0 = int(aa[s[0]])
    f0 = fibers[a0]
    hit = np.flatnonzero(f0 == s[0])
    if hit.size == 0:
        return None, 0.0
    u = np.zeros(f0.size)
    u[int(hit[0])] = 1.0
    tot, ap, maxdev = 0.0, a0, 0.0
    for t in range(1, n):
        ac = int(aa[s[t]])
        key = (ap, ac)
        if key not in Qt:
            return None, maxdev
        tot += float(u @ Rt[key])
        u = u @ Qt[key]
        dev = abs(float(u.sum()) - 1.0)
        maxdev = max(maxdev, dev)
        if dev > mass_tol:
            return None, maxdev
        ap = ac
    return tot, maxdev


def layer_shuffle_analytic_gene(s, q, layer, F, logP, diag_logP):
    n = s.size
    if n < 2:
        return None
    lay = layer[s]
    N = np.bincount(lay, minlength=F).astype(np.float64)
    M = np.bincount(lay[:-1] * F + lay[1:], minlength=F * F).reshape(F, F).astype(np.float64)
    U = np.zeros((F, q), dtype=np.float64)
    np.add.at(U, (lay, s), 1.0)
    nz = N > 0
    U[nz] /= N[nz, None]
    W = U @ logP @ U.T
    dg = U @ diag_logP
    tot = 0.0
    for a in range(F):
        if not nz[a]:
            continue
        for b in range(F):
            if M[a, b] == 0.0:
                continue
            if a != b:
                tot += M[a, b] * W[a, b]
            elif N[a] >= 2:
                tot += M[a, b] * (N[a] / (N[a] - 1.0) * W[a, a] - dg[a] / (N[a] - 1.0))
    return tot


def h567_audit(test_seqs, logP, aa, fibers, Q, h5_max_genes=0):
    diag = np.diag(logP).copy()
    Qt, Rt = h5_pair_tables(Q, logP, fibers)
    cohort = test_seqs if h5_max_genes <= 0 else test_seqs[:h5_max_genes]
    per = {"H5_h2b_surrogate": {}, "H6_synonymous": {}, "H7_mirror": {}}
    devs, nviol = [], 0
    t0 = time.time()
    for gi, s in enumerate(cohort):
        if s.size < 2:
            continue
        nt = s.size - 1
        r = float(np.sum(logP[s[:-1], s[1:]]))
        e5, dev = h5_analytic_gene(s, aa, fibers, Qt, Rt)
        devs.append(dev)
        if e5 is None:
            nviol += 1
        else:
            per["H5_h2b_surrogate"][gi] = (r - e5) / nt
        e6 = layer_shuffle_analytic_gene(s, NC, aa, NAA, logP, diag)
        if e6 is not None:
            per["H6_synonymous"][gi] = (r - e6) / nt
        e7 = layer_shuffle_analytic_gene(s, NC, MIRROR32, 32, logP, diag)
        if e7 is not None:
            per["H7_mirror"][gi] = (r - e7) / nt
    common = sorted(set(per["H5_h2b_surrogate"]) &
                    set(per["H6_synonymous"]) &
                    set(per["H7_mirror"]))
    res = {}
    for name, dct in per.items():
        com = np.array([dct[g] for g in common], dtype=np.float64)
        if com.size < 3:
            continue
        m, se, t, pv = _tstat(com)
        res[name] = {"n_genes_own": len(dct), "n_genes_common": int(com.size),
                     "delta_per_transition_gene_mean": m,
                     "gene_SE": se, "t_stat": t, "p_value": pv,
                     "ci95_low": m - 1.96 * se if _finite(se) else None,
                     "ci95_high": m + 1.96 * se if _finite(se) else None,
                     "frac_genes_positive": float(np.mean(com > 0)),
                     "cohort": "COMMON_TO_H5_H6_H7",
                     "method": "EXACT_ANALYTIC_ZERO_RNG"}
    if res and "H5_h2b_surrogate" in res:
        res["H5_h2b_surrogate"].update({
            "max_markov_mass_deviation": float(max(devs)) if devs else 0.0,
            "n_mass_invariant_violations": int(nviol), "mass_tolerance": MASS_TOL})
    for v in res.values():
        v["seconds"] = time.time() - t0
    return res


# =====================================================================
# 17. Self-tests
# =====================================================================
# =====================================================================
# 16b. THEOREMS V / W / X  (v2.2)
# =====================================================================
KAPPA4 = np.array([3, 2, 1, 0], dtype=np.int64)


def wobble_kernel_from_positions(seqs, k, alpha=0.1):
    """[FIX-44] 16x4 wobble kernel from codons at positions t >= k."""
    cnt = np.zeros((16, 4), dtype=np.float64)
    for s in seqs:
        if s.size <= k:
            continue
        t = s[k:]
        cnt += np.bincount(OUT16[t] * 4 + WOBBLE4[t], minlength=64).reshape(16, 4)
    q = cnt + alpha
    return q / np.maximum(q.sum(axis=1, keepdims=True), TINY)


def duplex_doublet_partition(aa):
    """[THM-V] 8 binary doublets of alpha_dup, 2x2 Latin square."""
    aa = np.asarray(aa, dtype=np.int64)
    lay, F, sizes, d8 = duplex_layers(aa)
    bp = box_partition(aa)
    unp, uns = bp["unsplit_prefix"], bp["unsplit_suffix"]
    wob = np.zeros(NC, dtype=bool)
    lead = np.zeros(NC, dtype=bool)
    dbl, n_other = [], 0
    for a in range(F):
        fa = np.flatnonzero(lay == a)
        if fa.size != 2:
            continue
        c1, c2 = int(fa[0]), int(fa[1])
        sO = bool(OUT16[c1] == OUT16[c2])
        sS = bool(SUF16[c1] == SUF16[c2])
        kind = "wob" if sO else ("lead" if sS else "other")
        if kind == "other":
            n_other += 1
        nc = int(HAS_CPG[fa].sum())
        quart = "CpG" if nc == 1 else ("Leu" if nc == 0 else "mixed")
        if kind == "wob":
            wob[fa] = True
        elif kind == "lead":
            lead[fa] = True
        o, s = int(OUT16[c1]), int(SUF16[c1])
        dbl.append({"codons": [CODONS[c1], CODONS[c2]],
                    "channels": [AA_STRING[int(aa[c1])], AA_STRING[int(aa[MIRROR[c1]])]],
                    "kind": kind, "quartet": quart,
                    "box": (BASES[o // 4] + BASES[o % 4]) if sO else
                           ((BASES[s // 4] + BASES[s % 4]) if sS else None),
                    "in_unsplit": bool(unp[o]) if sO else (bool(uns[s]) if sS else None)})
    return {"doublets": dbl, "wob_mask": wob, "lead_mask": lead,
            "n_doublets": len(dbl), "n_other": n_other,
            "n_wob": int(wob.sum()), "n_lead": int(lead.sum()),
            "rho_exchanges": bool(np.array_equal(wob[MIRROR], lead))}


def wobble_injectivity(aa):
    """[THM-V] within-box injectivity of c3 -> A^rho and c1 -> A."""
    aa = np.asarray(aa, dtype=np.int64)
    ar = aa[MIRROR]
    injO = np.array([np.unique(ar[OUT16 == o]).size == 4 for o in range(16)], dtype=bool)
    injS = np.array([np.unique(aa[SUF16 == s]).size == 4 for s in range(16)], dtype=bool)
    return injO, injS


def _cmi_terms(Tab, marg_col, N, strata, nstrat, colmap=None, ncol=None):
    """p(g) * I(rows ; cols | stratum g), per stratum."""
    nrow = Tab.shape[0]
    out = np.zeros(nstrat, dtype=np.float64)
    for g in range(nstrat):
        ys = np.flatnonzero(strata == g)
        if ys.size == 0:
            continue
        ng = float(marg_col[ys].sum())
        if ng <= 0:
            continue
        if colmap is None:
            T = Tab[:, ys]
        else:
            T = np.zeros((nrow, ncol), dtype=np.float64)
            for j, y in enumerate(ys):
                T[:, int(colmap[y])] += Tab[:, y]
        rs, cs = T.sum(axis=1), T.sum(axis=0)
        out[g] = (ng / N) * _mi_tab(T, np.maximum(rs[:, None] * cs[None, :], TINY) / ng, ng)
    return out


def thm_V_screen_split(B, aa=None):
    """[THM-V] exact two-strand refinement of the THM-U box layers."""
    if aa is None:
        aa = AA_TABLE1
    aa = np.asarray(aa, dtype=np.int64)
    sp = two_sided_screen_split(B, aa)
    if sp is None:
        return None
    N = sp["N"]
    ar = aa[MIRROR]
    bp = box_partition(aa)
    b4, s4 = bp["unsplit_prefix"], bp["unsplit_suffix"]
    injO, injS = wobble_injectivity(aa)
    Bsy, Bxo = sp["Bsy"], sp["Bxo"]
    ny, nx = sp["ny"], sp["nx"]
    w1 = _cmi_terms(Bsy, ny, N, OUT16, 16)
    w2 = _cmi_terms(Bsy, ny, N, OUT16, 16, ar, NAA)
    w3 = _cmi_terms(Bsy, ny, N, OUT16 * NAA + ar, 16 * NAA).reshape(16, NAA).sum(axis=1)
    l1 = _cmi_terms(Bxo.T, nx, N, SUF16, 16)
    l2 = _cmi_terms(Bxo.T, nx, N, SUF16, 16, aa, NAA)
    l3 = _cmi_terms(Bxo.T, nx, N, SUF16 * NAA + aa, 16 * NAA).reshape(16, NAA).sum(axis=1)
    part = duplex_multi_partition(aa)
    KIND = {"wob": 0, "lead": 2, "mirror": 4, "mixed": 4}
    cls = np.full(NAA * NAA, -1, dtype=np.int64)
    for d in part["fibers"]:
        c = CODON_IDX[d["codons"][0]]
        key = int(aa[c]) * NAA + int(ar[c])
        cls[key] = KIND[d["kind"]] + (0 if d["quartet"] == "CpG" else 1)
    sd = _cmi_terms(B, ny, N, aa * NAA + ar, NAA * NAA)
    four = np.array([float(sd[cls == k].sum()) for k in range(6)])
    return {"D_3prime_wobble": sp["D_3prime_wobble"], "D_wob_total": float(w1.sum()),
            "D_5prime_lead": sp["D_5prime_lead"], "D_lead_total": float(l1.sum()),
            "chain_wob_max_abs_err": float(np.abs(w1 - (w2 + w3)).max()),
            "chain_lead_max_abs_err": float(np.abs(l1 - (l2 + l3)).max()),
            "wob_B4_total": float(w1[b4].sum()), "wob_B4_AA": float(w2[b4].sum()),
            "wob_B4_res": float(w3[b4].sum()),
            "wob_residual_on_injective_B4": float(w3[b4 & injO].sum()),
            "lead_S4_total": float(l1[s4].sum()), "lead_S4_AA": float(l2[s4].sum()),
            "lead_S4_res": float(l3[s4].sum()),
            "lead_residual_on_injective_S4": float(l3[s4 & injS].sum()),
            "S_dup_total": float(sd.sum()), "S_dup_four": four.tolist(),
            "S_dup_split_err": abs(float(four.sum()) - float(sd.sum())),
            "four_way_min": float(four.min())}


def thm_W_wobble_split(Pw, Rw, no, Nte, aa=None):
    """[THM-W] D_wobble = stop-leakage + intra-box AA-split + pure synonymous."""
    if aa is None:
        aa = AA_TABLE1
    aa = np.asarray(aa, dtype=np.int64)
    sense = (aa != STOP)
    b4 = box_partition(aa)["unsplit_prefix"]
    ts = ta = ty = 0.0
    per_box = []
    for o in range(16):
        if no[o] <= 0:
            continue
        w = float(no[o]) / Nte
        ys = np.flatnonzero(OUT16 == o)
        p = np.asarray(Pw[o], dtype=np.float64)
        r = np.asarray(Rw[o], dtype=np.float64)
        sm = sense[ys]
        Po = float(p[sm].sum())
        So = float(r[sm].sum())
        kl_bin = 0.0
        if Po > 0:
            kl_bin += Po * math.log(Po / max(So, TINY))
        if 1.0 - Po > 0:
            kl_bin += (1.0 - Po) * math.log((1.0 - Po) / max(1.0 - So, TINY))
        d_stop = w * kl_bin
        labs = aa[ys]
        d_aa = d_syn = 0.0
        if Po > 0:
            ps = np.where(sm, p / Po, 0.0)
            rs = np.where(sm, r / max(So, TINY), 0.0)
            for a in np.unique(labs[sm]):
                m = (labs == a) & sm
                pa, ra = float(ps[m].sum()), float(rs[m].sum())
                if pa > 0:
                    d_aa += Po * pa * math.log(pa / max(ra, TINY))
                    for j in np.flatnonzero(m):
                        if ps[j] > 0:
                            d_syn += Po * ps[j] * math.log((ps[j] / pa) / max(rs[j] / max(ra, TINY), TINY))
        if 1.0 - Po > 0:
            ust = ~sm
            pst = p[ust] / (1.0 - Po)
            rst = r[ust] / max(1.0 - So, TINY)
            for j in range(pst.size):
                if pst[j] > 0:
                    d_syn += (1.0 - Po) * float(pst[j]) * math.log(float(pst[j]) / max(float(rst[j]), TINY))
        ts += d_stop
        ta += w * d_aa
        ty += w * d_syn
        per_box.append({"box": BASES[o // 4] + BASES[o % 4], "S_sense": So,
                        "stop": d_stop, "aa_split": w * d_aa, "syn": w * d_syn,
                        "unsplit": bool(b4[o])})
    direct = 0.0
    for o in range(16):
        if no[o] > 0:
            direct += (float(no[o]) / Nte) * _kl_joint(np.asarray(Pw[o], float),
                                                       np.asarray(Rw[o], float))
    tot = ts + ta + ty
    return {"D_stop_leakage": ts, "D_AA_split": ta, "D_syn": ty, "sum": tot,
            "D_wobble_direct": direct, "residual": abs(tot - direct),
            "all_nonnegative": bool(min(ts, ta, ty) > -1e-12),
            "stop_on_B4": float(sum(b["stop"] for b in per_box if b["unsplit"])),
            "aa_on_B4": float(sum(b["aa_split"] for b in per_box if b["unsplit"])),
            "per_box": per_box}


def thm_X_screen_v4():
    """[THM-X] R D_S R = D_O, [Pi_screen,R]=0, Pi_C2 is full V_4."""
    R = np.zeros((NC, NC))
    R[np.arange(NC), MIRROR] = 1.0
    DS = _layer_avg_matrix(SUF16, 16)
    DO = _layer_avg_matrix(OUT16, 16)
    ix = np.arange(NC)
    probe = (1 + ((ix[:, None] * 17 + ix[None, :] * 31
                   + (ix[:, None] ^ ix[None, :])) % 97)).astype(float)
    probe /= probe.sum()
    Pi = lambda M: DS @ M @ DO
    Rop = lambda M: R @ M.T @ R
    Top = lambda M: M.T
    Sop = lambda M: M[np.ix_(MIRROR, MIRROR)]
    j = np.arange(256)
    img = RHO_BOX[j % 16] * 16 + RHO_BOX[j // 16]
    trR_img = int(np.count_nonzero(img == j))
    trR_full = int(np.sum(np.arange(NC)[None, :] == MIRROR[:, None]))
    dp, dm = (256 + trR_img) // 2, (256 - trR_img) // 2
    fp, fm = (NC * NC + trR_full) // 2, (NC * NC - trR_full) // 2
    DC2 = DS @ DO
    PiC = lambda M: DC2 @ M @ DC2
    t = np.arange(16)
    a, b = t // 4, t % 4
    trT = int(np.count_nonzero((b * 4 + a) == t))
    trRc = int(np.count_nonzero((KAPPA4[b] * 4 + KAPPA4[a]) == t))
    trS = int(np.count_nonzero((KAPPA4[a] * 4 + KAPPA4[b]) == t))
    return {"interleave_max_abs_err": max(float(np.abs(R @ DS @ R - DO).max()),
                                          float(np.abs(R @ DO @ R - DS).max())),
            "comm_DS_R_sq": float(np.linalg.norm(DS @ R - R @ DS, ord="fro") ** 2),
            "comm_DO_R_sq": float(np.linalg.norm(DO @ R - R @ DO, ord="fro") ** 2),
            "S2_trace_DS_DO": float(np.trace(DS @ DO)),
            "screen_R_commutator_max_abs_err": float(np.abs(Rop(Pi(probe)) - Pi(Rop(probe))).max()),
            "screen_T_commutator_max_abs_err": float(np.abs(Top(Pi(probe)) - Pi(Top(probe))).max()),
            "trace_R_on_image": trR_img, "dim_image_plus": dp, "dim_image_minus": dm,
            "dim_kernel_plus": fp - dp, "dim_kernel_minus": fm - dm,
            "C2_symmetry_err": float(np.abs(DC2 - DC2.T).max()),
            "C2_idempotence_err": float(np.abs(DC2 @ DC2 - DC2).max()),
            "C2_comm_R_err": float(np.abs(DC2 @ R - R @ DC2).max()),
            "C2_T_err": float(np.abs(Top(PiC(probe)) - PiC(Top(probe))).max()),
            "C2_R_err": float(np.abs(Rop(PiC(probe)) - PiC(Rop(probe))).max()),
            "C2_S_err": float(np.abs(Sop(PiC(probe)) - PiC(Sop(probe))).max()),
            "C2_traces": {"id": 16, "T": trT, "R": trRc, "S": trS},
            "C2_quadrant_dims": {"W_++": (16 + trT + trRc + trS) // 4,
                                 "W_+-": (16 + trT - trRc - trS) // 4,
                                 "W_-+": (16 - trT + trRc - trS) // 4,
                                 "W_--": (16 - trT - trRc + trS) // 4}}


def self_test_v22(ensG=None, aa1=None):
    ok = {}
    if aa1 is None:
        aa1 = AA_TABLE1
    X = thm_X_screen_v4()
    assert X["interleave_max_abs_err"] < 1e-14
    assert abs(X["comm_DS_R_sq"] - 24.0) < 1e-10
    assert abs(X["comm_DO_R_sq"] - 24.0) < 1e-10
    assert abs(X["S2_trace_DS_DO"] - 4.0) < 1e-12
    assert X["screen_R_commutator_max_abs_err"] < 1e-14
    assert X["screen_T_commutator_max_abs_err"] > 1e-6
    assert (X["dim_image_plus"], X["dim_image_minus"]) == (136, 120)
    assert (X["dim_kernel_plus"], X["dim_kernel_minus"]) == (1944, 1896)
    assert X["C2_symmetry_err"] < 1e-14 and X["C2_idempotence_err"] < 1e-14
    assert X["C2_comm_R_err"] < 1e-14
    assert max(X["C2_T_err"], X["C2_R_err"], X["C2_S_err"]) < 1e-14
    assert X["C2_quadrant_dims"] == {"W_++": 6, "W_+-": 4, "W_-+": 4, "W_--": 2}
    ok["thmX_screen_v4"] = X

    Pt = duplex_doublet_partition(aa1)
    assert Pt["n_doublets"] == 8 and Pt["n_other"] == 0
    assert Pt["n_wob"] == 8 and Pt["n_lead"] == 8 and Pt["rho_exchanges"]
    cells = {}
    for d in Pt["doublets"]:
        cells.setdefault((d["kind"], d["quartet"]), []).append((d["box"], d["in_unsplit"]))
    assert set(x[0] for x in cells[("wob", "CpG")]) == {"TC", "CC"}
    assert all(x[1] for x in cells[("wob", "CpG")])
    assert set(x[0] for x in cells[("wob", "Leu")]) == {"TA", "CA"}
    assert not any(x[1] for x in cells[("wob", "Leu")])
    assert set(x[0] for x in cells[("lead", "CpG")]) == {"GA", "GG"}
    assert all(x[1] for x in cells[("lead", "CpG")])
    assert set(x[0] for x in cells[("lead", "Leu")]) == {"TA", "TG"}
    assert not any(x[1] for x in cells[("lead", "Leu")])
    injO, injS = wobble_injectivity(aa1)
    assert np.array_equal(injO, injS[RHO_BOX])
    nm = lambda m: set(BASES[i // 4] + BASES[i % 4] for i in range(16) if m[i])
    assert nm(~injS) == {"GA", "GG", "TA", "TG"}
    assert nm(~injO) == {"TC", "CC", "TA", "CA"}
    ok["thmV_structure"] = {"n_wob": Pt["n_wob"], "n_lead": Pt["n_lead"],
                            "noninjective_prefix": sorted(nm(~injO)),
                            "noninjective_suffix": sorted(nm(~injS))}

    layD, FD, szD, d8m = duplex_layers(aa1)
    if ensG is not None:
        eV = eSd = 0.0
        for k in (1, 3, 7):
            Bk = lag_counts(ensG, 64, k)
            V = thm_V_screen_split(Bk, aa1)
            eV = max(eV, V["chain_wob_max_abs_err"], V["chain_lead_max_abs_err"],
                     abs(V["D_wob_total"] - V["D_3prime_wobble"]),
                     abs(V["D_lead_total"] - V["D_5prime_lead"]), V["S_dup_split_err"])
            assert V["wob_residual_on_injective_B4"] < 1e-14
            assert V["lead_residual_on_injective_S4"] < 1e-14
            assert V["four_way_min"] > -1e-12
            g = lag_duplex_decomposition(ensG, 64, k, aa1, d8m)
            eSd = max(eSd, abs(V["S_dup_total"] - g["S_dup"]))
        assert eV < 1e-13 and eSd < 1e-13, (eV, eSd)
        ok["thmV_numeric_max_abs_err"] = eV
        ok["thmV_Sdup_vs_thmG_max_abs_err"] = eSd

        ensG_sense = [s[aa1[s] != STOP] for s in ensG]
        _, mh = build_h3b_lag_kernels(ensG_sense, [1], alpha=0.1)
        _, mp = build_h3b_lag_kernels(ensG_sense, [1], alpha=1e-9)
        eW = 0.0
        for k in (1, 3):
            Bt = lag_counts(ensG_sense, 64, k)
            sp = two_sided_screen_split(Bt, aa1)
            no, Nte = sp["no"], float(Bt.sum())
            Pw = np.zeros((16, 4))
            for o in range(16):
                if no[o] > 0:
                    v = sp["ny"][OUT16 == o]
                    Pw[o] = v / max(v.sum(), TINY)
            W = thm_W_wobble_split(Pw, mh["R_wobble"], no, Nte, aa1)
            Wp = thm_W_wobble_split(Pw, mp["R_wobble_phase"], no, Nte, aa1)
            eW = max(eW, W["residual"], Wp["residual"])
            assert W["all_nonnegative"] and Wp["all_nonnegative"]
            assert W["stop_on_B4"] < 1e-14 and W["aa_on_B4"] < 1e-14
            assert W["D_stop_leakage"] > 1e-6
            assert Wp["D_stop_leakage"] < 1e-8
        assert eW < 1e-13, eW
        ok["thmW_three_way_max_abs_err"] = eW

        for o in theorem_O_multi_lag(ensG, [1, 5], aa1):
            assert o["df_phase_additivity_abs_err"] == 0
            assert o["df_phase_full_eff"] >= o["df_phase_full"]
        cdn, _ = debruijn_phase_null(3, 6)
        od = theorem_O_frame_phase_audit([cdn], 1, aa1)
        assert od["df_phase_additivity_abs_err"] == 0
        assert od["df_phase_full_eff"] == od["df_phase_full"] == 2 * 4095
        ok["fix42_phase_df_additivity"] = True
    ok["v22_all_passed"] = True
    return ok


# =====================================================================
# 16c. THEOREMS Y / Z  (v2.3)  —  exact bridge G.2 <-> T
# =====================================================================
def duplex_multi_partition(aa):
    """[FIX-48] generic multi-codon duplex fibers for ANY NCBI table.

    Table 6 has TWO QUARTETS (k=4), Table 5 a mirror doublet {TCT,AGA}.
    kinds: wob (same O), lead (same S), mirror (k=2, rho-paired), mixed (k>2).
    """
    aa = np.asarray(aa, dtype=np.int64)
    lay, F, per, multi = duplex_layers(aa)
    bp = box_partition(aa)
    unp, uns = bp["unsplit_prefix"], bp["unsplit_suffix"]
    out = []
    for a in range(F):
        fa = np.flatnonzero(lay == a)
        if fa.size < 2:
            continue
        sO = bool(np.unique(OUT16[fa]).size == 1)
        sS = bool(np.unique(SUF16[fa]).size == 1)
        kind = "wob" if sO else ("lead" if sS else ("mirror" if fa.size == 2 else "mixed"))
        nc = int(HAS_CPG[fa].sum())
        quart = "CpG" if nc == 1 else ("nonCpG" if nc == 0 else "mixedCpG")
        o, s = int(OUT16[fa[0]]), int(SUF16[fa[0]])
        out.append({"codons": [CODONS[int(c)] for c in fa], "size": int(fa.size),
                    "kind": kind, "quartet": quart,
                    "channels": [AA_STRING[int(aa[fa[0]])],
                                 AA_STRING[int(aa[MIRROR[fa[0]]])]],
                    "box": (BASES[o // 4] + BASES[o % 4]) if sO else
                           ((BASES[s // 4] + BASES[s % 4]) if sS else None),
                    "in_unsplit": bool(unp[o]) if sO else (bool(uns[s]) if sS else None)})
    sz = np.array([d["size"] for d in out], dtype=np.int64) if out else np.zeros(0, np.int64)
    ncod = int(sz.sum())
    return {"fibers": out, "n_multi": len(out),
            "n_doublets": int(np.sum(sz == 2)), "n_quartets": int(np.sum(sz == 4)),
            "n_codons_multi": ncod, "n_singletons": int(NC - ncod),
            "n_channels": int(F), "dim_dup": int(NC - F),
            "kinds": sorted({d["kind"] for d in out})}


def _bipartite_components(er, ec, nrow, ncol):
    par = list(range(nrow + ncol))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for r, c in zip(er, ec):
        a, b = find(int(r)), find(nrow + int(c))
        if a != b:
            par[a] = b
    return len({find(i) for i in range(nrow + ncol)})


def thm_Y_diamond(B, aa=None, sense=True):
    """[THM-Y] O–A bipartite forest, modular sigma-lattice diamond, and the
    EXACT balance  D_3'wobble + D_1<->3 - D_split_AA = S_syn - S_syn^hex_box.
    df: (q-1)(q-nO) - (q-1)(nE-nO) = (q-1)(q-nE) = (q-1)(q-nA) - (q-1)(nE-nA).
    """
    if aa is None:
        aa = AA_TABLE1
    aa = np.asarray(aa, dtype=np.int64)
    sm = (aa != STOP) if sense else np.ones(NC, dtype=bool)
    cod = np.flatnonzero(sm)
    ub, ia = np.unique(OUT16[cod], return_inverse=True)
    ua, ib = np.unique(aa[cod], return_inverse=True)
    nO, nA = int(ub.size), int(ua.size)
    ed = np.unique(ia * nA + ib)
    nE = int(ed.size)
    nC = _bipartite_components(ed // nA, ed % nA, nO, nA)
    q = int(cod.size)
    sp = two_sided_screen_split(B, aa)
    if sp is None:
        return None
    ny, N = sp["ny"], float(B.sum())
    lhs = sp["D_3prime_wobble"] + sp["D_1_3"]
    D_saa = float(_cmi_terms(B, ny, N, OUT16, 16, aa, NAA).sum())
    S_box = float(_cmi_terms(B, ny, N, OUT16 * NAA + aa, 16 * NAA).sum())
    S_syn = float(_cmi_terms(B, ny, N, aa, NAA).sum())
    S_hex = float(_cmi_terms(B, ny, N, aa, NAA, OUT16, 16).sum())
    df = {"I_given_O": (q - 1) * (q - nO), "split_AA": (q - 1) * (nE - nO),
          "syn_box": (q - 1) * (q - nE), "S_syn": (q - 1) * (q - nA),
          "hex_box": (q - 1) * (nE - nA)}
    return {"sense": bool(sense), "q": q, "n_O": nO, "n_A": nA, "n_OA": nE,
            "n_components": nC, "is_forest": bool(nC == nO + nA - nE),
            "modularity_residual": int((nO + nA) - (nE + nC)),
            "I_X_Y_given_O": lhs, "D_split_AA": D_saa, "S_syn_box": S_box,
            "S_syn": S_syn, "S_syn_hex_box": S_hex,
            "balance_abs_err": abs((lhs - D_saa) - (S_syn - S_hex)),
            "chain_abs_err": max(abs(lhs - (D_saa + S_box)),
                                 abs(S_syn - (S_hex + S_box))),
            "df": df,
            "df_balance_ok": bool(df["I_given_O"] - df["split_AA"] == df["syn_box"]
                                  == df["S_syn"] - df["hex_box"])}


def thm_Z_tower(B, aa=None, sense=True):
    """[THM-Z] tower sigma(O) < sigma(O,A) < sigma(O,A,A^rho) < sigma(Y),
    chain rule S_dup = S_dup^lead + S_dup^wob, and the CORRECTED diamond
        S_syn^hex_box + I(X;A^rho|A,O) = S_anti + S_dup^lead
    (180 + 2100 = 2040 + 240 = 2280 sense df; NOT 180 = 1 + 4).
    """
    if aa is None:
        aa = AA_TABLE1
    aa = np.asarray(aa, dtype=np.int64)
    ar = aa[MIRROR]
    sm = (aa != STOP) if sense else np.ones(NC, dtype=bool)
    cod = np.flatnonzero(sm)
    n1 = int(np.unique(OUT16[cod]).size)
    n2 = int(np.unique(OUT16[cod] * NAA + aa[cod]).size)
    n3 = int(np.unique((OUT16[cod] * NAA + aa[cod]) * NAA + ar[cod]).size)
    n4 = int(cod.size)
    ndup = int(np.unique(aa[cod] * NAA + ar[cod]).size)
    q = n4
    sp = two_sided_screen_split(B, aa)
    if sp is None:
        return None
    ny, N = sp["ny"], float(B.sum())
    dstr = aa * NAA + ar
    S_dup = float(_cmi_terms(B, ny, N, dstr, NAA * NAA).sum())
    S_lead = float(_cmi_terms(B, ny, N, dstr, NAA * NAA, OUT16, 16).sum())
    S_wob = float(_cmi_terms(B, ny, N, (OUT16 * NAA + aa) * NAA + ar,
                             16 * NAA * NAA).sum())
    S_anti = float(_cmi_terms(B, ny, N, aa, NAA, ar, NAA).sum())
    S_hex = float(_cmi_terms(B, ny, N, aa, NAA, OUT16, 16).sum())
    I_ar_AO = float(_cmi_terms(B, ny, N, aa * 16 + OUT16, NAA * 16, ar, NAA).sum())
    df = {"S_dup": (q - 1) * (q - ndup), "lead": (q - 1) * (n3 - ndup),
          "wob": (q - 1) * (q - n3), "hex_box": (q - 1) * (n2 - int(np.unique(aa[cod]).size)),
          "Arho_given_A_O": (q - 1) * (n3 - n2), "S_anti": (q - 1) * (ndup - int(np.unique(aa[cod]).size))}
    return {"sense": bool(sense), "tower": [n1, n2, n3, n4],
            "increments": [n2 - n1, n3 - n2, n4 - n3],
            "increments_sum_ok": bool((n2 - n1) + (n3 - n2) + (n4 - n3) == q - n1),
            "n_dup_channels": ndup, "S_dup": S_dup,
            "S_dup_lead": S_lead, "S_dup_wob": S_wob,
            "chain_abs_err": abs(S_dup - (S_lead + S_wob)),
            "S_anti": S_anti, "S_syn_hex_box": S_hex,
            "I_Arho_given_A_O": I_ar_AO,
            "diamond_abs_err": abs((S_hex + I_ar_AO) - (S_anti + S_lead)),
            "df": df,
            "df_ok": bool(df["lead"] + df["wob"] == df["S_dup"] and
                          df["hex_box"] + df["Arho_given_A_O"] == df["S_anti"] + df["lead"])}


def self_test_v23(ensG=None, aa1=None):
    ok = {}
    if aa1 is None:
        aa1 = AA_TABLE1
    # ---- [FIX-48] full 8-table duplex classification (hand-derived) ----
    EXP = {1: (56, 48, 8, 0, 8), 11: (56, 48, 8, 0, 8), 4: (56, 48, 8, 0, 8),
           2: (60, 56, 4, 0, 4), 3: (60, 56, 4, 0, 4), 5: (59, 54, 5, 0, 5),
           12: (58, 52, 6, 0, 6), 6: (54, 48, 4, 2, 10)}
    cls = {}
    for tid in sorted(SUPPORTED_TABLES):
        aat = aa_array(tid)
        pm = duplex_multi_partition(aat)
        got = (pm["n_channels"], pm["n_singletons"], pm["n_doublets"],
               pm["n_quartets"], pm["dim_dup"])
        assert got == EXP[tid], (tid, got, EXP[tid])
        lay, F, per, mul = duplex_layers(aat)
        assert int(np.sum(mul)) == pm["n_codons_multi"]
        assert int(F) == pm["n_channels"]
        cls[str(tid)] = {"F": got[0], "singletons": got[1], "doublets": got[2],
                         "quartets": got[3], "dim_dup": got[4], "kinds": pm["kinds"]}
    assert "mirror" in duplex_multi_partition(aa_array(5))["kinds"]
    assert "mixed" in duplex_multi_partition(aa_array(6))["kinds"]
    ok["fix48_duplex_classification"] = cls
    # ---- [THM-Y/Z] structure on ALL tables, sense and full ----
    yz = {}
    for tid in sorted(SUPPORTED_TABLES):
        aat = aa_array(tid)
        for sense in (True, False):
            cod = np.flatnonzero((aat != STOP) if sense else np.ones(NC, bool))
            ub, ia = np.unique(OUT16[cod], return_inverse=True)
            ua, ib = np.unique(aat[cod], return_inverse=True)
            nO, nA = int(ub.size), int(ua.size)
            ed = np.unique(ia * nA + ib)
            nE = int(ed.size)
            nC = _bipartite_components(ed // nA, ed % nA, nO, nA)
            assert nC == nO + nA - nE, (tid, sense, nO, nA, nE, nC)
            if tid == 1:
                assert (nO, nA, nE, nC) == ((16, 20, 23, 13) if sense else (16, 21, 25, 12))
                yz["sense" if sense else "full"] = [nO, nA, nE, nC]
    ok["thmY_forest_all_tables"] = True
    ok["thmY_table1_diamond"] = yz
    if ensG is None:
        ok["v23_all_passed"] = True
        return ok
    # ---- numeric Y / Z on the deterministic fixture ----
    eY = eZ = eD = 0.0
    for k in (1, 3, 7):
        Bk = lag_counts(ensG, 64, k)
        for sense in (True, False):
            Yv = thm_Y_diamond(Bk, aa1, sense)
            Zv = thm_Z_tower(Bk, aa1, sense)
            assert Yv["is_forest"] and Yv["modularity_residual"] == 0
            assert Yv["df_balance_ok"] and Zv["df_ok"] and Zv["increments_sum_ok"]
            eY = max(eY, Yv["balance_abs_err"], Yv["chain_abs_err"])
            eZ = max(eZ, Zv["chain_abs_err"])
            eD = max(eD, Zv["diamond_abs_err"])
            if sense:
                assert Yv["df"]["I_given_O"] == 2700 and Yv["df"]["split_AA"] == 420
                assert Yv["df"]["syn_box"] == 2280 and Yv["df"]["S_syn"] == 2460
                assert Yv["df"]["hex_box"] == 180
                assert Zv["tower"] == [16, 23, 58, 61]
                assert Zv["df"]["lead"] == 240 and Zv["df"]["wob"] == 180
                assert Zv["df"]["S_dup"] == 420
                assert Zv["df"]["hex_box"] + Zv["df"]["Arho_given_A_O"] == 2280
            else:
                assert Yv["df"]["I_given_O"] == 3024 and Yv["df"]["split_AA"] == 567
                assert Zv["tower"] == [16, 25, 60, 64]
    assert eY < 1e-13 and eZ < 1e-13 and eD < 1e-13, (eY, eZ, eD)
    ok["thmY_balance_max_abs_err"] = eY
    ok["thmZ_chain_max_abs_err"] = eZ
    ok["thmZ_diamond_max_abs_err"] = eD
    # ---- [FIX-48] Table 6 quartet bound regression ----
    aa6 = aa_array(6)
    _, _, _, m6 = duplex_layers(aa6)
    g6 = lag_duplex_decomposition(ensG, 64, 3, aa6, m6)
    assert g6["D8_has_quartets"] and g6["D8_bound_ok"]
    assert g6["D8_bound_nats"] > g6["D8_bound_log2_legacy"] - 1e-15
    g1 = lag_duplex_decomposition(ensG, 64, 3, aa1, duplex_layers(aa1)[3])
    assert abs(g1["D8_bound_nats"] - g1["D8_bound_log2_legacy"]) < 1e-15
    ok["fix48_table6_quartet_bound_ok"] = True
    # ---- V / W generic bins on Table 5 and 6 ----
    for tid in (5, 6, 12):
        aat = aa_array(tid)
        Bk = lag_counts(ensG, 64, 3)
        Vv = thm_V_screen_split(Bk, aat)
        assert Vv["S_dup_split_err"] < 1e-13, (tid, Vv["S_dup_split_err"])
        gg = lag_duplex_decomposition(ensG, 64, 3, aat, duplex_layers(aat)[3])
        assert abs(Vv["S_dup_total"] - gg["S_dup"]) < 1e-13
    ok["thmV_generic_tables_ok"] = True
    ok["v23_all_passed"] = True
    return ok


# =====================================================================
# 16d. THEOREM AA  (v2.4)  —  left tower, d2b split, self-dual 3x3 matrix
# =====================================================================
def _mi_cond2(B, rl, cl, cr, cc, ncr, ncc):
    """EXACT I( rl(X) ; cl(Y) | cr(X), cc(Y) ) by stratified chain rule."""
    N = float(B.sum())
    if N <= 0:
        return 0.0
    Sm = np.bincount((cr[:, None] * ncc + cc[None, :]).ravel(),
                     weights=B.ravel(), minlength=ncr * ncc).reshape(ncr, ncc)
    tot = 0.0
    for kr, kc in zip(*np.nonzero(Sm)):
        rows = np.flatnonzero(cr == kr)
        cols = np.flatnonzero(cc == kc)
        T = B[np.ix_(rows, cols)]
        n = float(Sm[kr, kc])
        ru, ri = np.unique(rl[rows], return_inverse=True)
        cu, ci = np.unique(cl[cols], return_inverse=True)
        agg = np.bincount((ri[:, None] * cu.size + ci[None, :]).ravel(), weights=T.ravel(),
                          minlength=ru.size * cu.size).reshape(ru.size, cu.size)
        rs, cs = agg.sum(axis=1), agg.sum(axis=0)
        tot += (n / N) * _mi_tab(agg, np.maximum(rs[:, None] * cs[None, :], TINY) / n, n)
    return float(tot)


def left_screen_structure(aa=None, sense=True):
    """[THM-AA.1] sigma(S,A) on the LEFT codon. NOT a forest (unlike Thm Y)."""
    if aa is None:
        aa = AA_TABLE1
    aa = np.asarray(aa, dtype=np.int64)
    sm = (aa != STOP) if sense else np.ones(NC, dtype=bool)
    cod = np.flatnonzero(sm)
    us, is_ = np.unique(SUF16[cod], return_inverse=True)
    ua, ia = np.unique(aa[cod], return_inverse=True)
    nS, nA = int(us.size), int(ua.size)
    ed = np.unique(is_ * nA + ia)
    nSA = int(ed.size)
    nC = _bipartite_components(ed // nA, ed % nA, nS, nA)
    q = int(cod.size)
    return {"sense": bool(sense), "q": q, "n_S": nS, "n_A": nA, "n_SA": nSA,
            "n_components": nC, "cycle_rank": int(nSA - (nS + nA) + nC),
            "is_forest": bool(nSA == nS + nA - nC),
            "modularity_residual": int((nS + nA) - (nSA + nC)),
            "n_lead_doublets": q - nSA, "tower": [nS, nSA, q],
            "increments": [nS - 1, nSA - nS, q - nSA]}


def _right_anti_structure(aa, sense=True):
    sm = (aa != STOP) if sense else np.ones(NC, dtype=bool)
    cod = np.flatnonzero(sm)
    ar = aa[MIRROR]
    nO = int(np.unique(OUT16[cod]).size)
    nOAr = int(np.unique((OUT16 * NAA + ar)[cod]).size)
    q = int(cod.size)
    return {"n_O": nO, "n_OAr": nOAr, "q": q,
            "tower": [nO, nOAr, q], "increments": [nO - 1, nOAr - nO, q - nOAr]}


def thm_AA_d2b(B, aa=None, sense=True):
    """[THM-AA.2] d2b = I(S;A'|A) + I(X;A'|A,S); 703+76=779 / 780+80=860."""
    if aa is None:
        aa = AA_TABLE1
    aa = np.asarray(aa, dtype=np.int64)
    idx = np.arange(NC)
    z = np.zeros(NC, dtype=np.int64)
    d2b = _mi_cond2(B, idx, aa, aa, z, NAA, 1)
    scr = _mi_cond2(B, SUF16, aa, aa, z, NAA, 1)
    lead = _mi_cond2(B, idx, aa, aa * 16 + SUF16, z, NAA * 16, 1)
    st = left_screen_structure(aa, sense)
    nSA, nA, q = st["n_SA"], st["n_A"], st["q"]
    df = {"screen": (nSA - nA) * (nA - 1), "lead": (q - nSA) * (nA - 1),
          "total": (q - nA) * (nA - 1)}
    return {"sense": bool(sense), "d2b": d2b, "d2b_screen": scr, "d2b_lead": lead,
            "residual": abs(d2b - (scr + lead)), "df": df,
            "df_ok": bool(df["screen"] + df["lead"] == df["total"]),
            "lead_share_df": (df["lead"] / df["total"]) if df["total"] else None}


def thm_AA_matrix(B, aa=None, sense=True, dual=True):
    """[THM-AA.3] 3x3 non-negative chain-rule matrix; M(anti) = M(sense)^T."""
    if aa is None:
        aa = AA_TABLE1
    aa = np.asarray(aa, dtype=np.int64)
    ar = aa[MIRROR]
    idx = np.arange(NC)
    z = np.zeros(NC, dtype=np.int64)
    Ls = [SUF16, SUF16 * NAA + aa, idx]
    Rs = [OUT16, OUT16 * NAA + ar, idx]
    Lc, Ln = [z, SUF16, SUF16 * NAA + aa], [1, 16, 16 * NAA]
    Rc, Rn = [z, OUT16, OUT16 * NAA + ar], [1, 16, 16 * NAA]

    def _mat(BB):
        M = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                M[i, j] = _mi_cond2(BB, Ls[i], Rs[j], Lc[i], Rc[j], Ln[i], Rn[j])
        return M

    M = _mat(B)
    stL = left_screen_structure(aa, sense)
    stR = _right_anti_structure(aa, sense)
    dL, dR = stL["increments"], stR["increments"]
    DF = np.outer(np.array(dL), np.array(dR))
    sp = two_sided_screen_split(B, aa)
    err_T = max(abs(M[0, 0] - sp["I_bridge"]),
                abs(M[1, 0] + M[2, 0] - sp["D_5prime_lead"]),
                abs(M[0, 1] + M[0, 2] - sp["D_3prime_wobble"]),
                abs(M[1, 1] + M[1, 2] + M[2, 1] + M[2, 2] - sp["D_1_3"]),
                abs(M.sum() - sp["I_CC"]),
                abs(M[0, :].sum() - sp["I_suffix16"]))
    derr = 0.0
    if dual:
        Ba = B.T[np.ix_(MIRROR, MIRROR)]
        derr = float(np.abs(_mat(Ba) - M.T).max())
    q = stL["q"]
    return {"sense": bool(sense), "M": M.tolist(), "df": DF.tolist(),
            "dL": dL, "dR": dR, "left_tower": stL["tower"],
            "right_tower": stR["tower"],
            "row_sums": M.sum(axis=1).tolist(), "col_sums": M.sum(axis=0).tolist(),
            "df_row_sums": DF.sum(axis=1).tolist(),
            "df_col_sums": DF.sum(axis=0).tolist(),
            "T_split_abs_err": float(err_T),
            "dual_transpose_abs_err": derr,
            "df_total_ok": bool(DF.sum() == (q - 1) ** 2),
            "M22_share_of_D13": float(M[1, 1] / max(M[1, 1] + M[1, 2] + M[2, 1] + M[2, 2], TINY)),
            "all_nonnegative": bool(M.min() > -1e-12)}


def self_test_v24(ensG=None, aa1=None):
    ok = {}
    if aa1 is None:
        aa1 = AA_TABLE1
    # ---- [FIX-51] cpg/leu masks cover ALL multi-fibers on every table ----
    cov = {}
    for tid in sorted(SUPPORTED_TABLES):
        aat = aa_array(tid)
        cp = cpg_quartet_partition(aat)
        _, _, _, mul = duplex_layers(aat)
        assert np.array_equal(cp["cpg_mask"] | cp["leu_mask"], mul), tid
        assert cp["cpg_quartet_rho_invariant"] and cp["leu_quartet_rho_invariant"], tid
        cov[str(tid)] = {"n_cpg": cp["n_cpg_doublets"], "n_nonCpG": cp["n_leu_doublets"],
                         "n_mixed": cp["n_mixed_doublets"]}
    assert cov["6"]["n_nonCpG"] == 2 and cov["6"]["n_cpg"] == 4
    ok["fix51_cpg_masks_all_tables"] = cov
    # ---- [THM-AA.1] left structure on all tables (NOT a forest) ----
    ls = {}
    for tid in sorted(SUPPORTED_TABLES):
        aat = aa_array(tid)
        for sense in (True, False):
            st = left_screen_structure(aat, sense)
            assert st["modularity_residual"] == -st["cycle_rank"]
            assert sum(st["increments"]) == st["q"] - 1
            if tid == 1:
                exp = (16, 20, 57, 4, 25) if sense else (16, 21, 60, 3, 26)
                got = (st["n_S"], st["n_A"], st["n_SA"], st["n_components"], st["cycle_rank"])
                assert got == exp, (sense, got, exp)
                assert not st["is_forest"]
                assert st["n_lead_doublets"] == 4
                ls["sense" if sense else "full"] = list(got)
    ok["thmAA_left_table1"] = ls
    if ensG is None:
        ok["v24_all_passed"] = True
        return ok
    # ---- numeric AA on the deterministic fixture ----
    eT = eD = e2 = 0.0
    for k in (1, 3, 7):
        Bk = lag_counts(ensG, 64, k)
        for sense in (True, False):
            A = thm_AA_matrix(Bk, aa1, sense)
            D = thm_AA_d2b(Bk, aa1, sense)
            assert A["all_nonnegative"] and A["df_total_ok"] and D["df_ok"]
            eT = max(eT, A["T_split_abs_err"])
            eD = max(eD, A["dual_transpose_abs_err"])
            e2 = max(e2, D["residual"])
            if sense:
                assert A["dL"] == [15, 41, 4] and A["dR"] == [15, 42, 3]
                assert A["df_row_sums"] == [900, 2460, 240]
                assert A["df_col_sums"] == [900, 2520, 180]
                assert A["df"][1][1] == 1722 and A["df"][2][2] == 12
                assert D["df"] == {"screen": 703, "lead": 76, "total": 779}
                assert A["left_tower"] == [16, 57, 61]
                assert A["right_tower"] == [16, 58, 61]
            else:
                assert A["dL"] == [15, 44, 4] and A["dR"] == [15, 44, 4]
                assert A["df_row_sums"] == [945, 2772, 252]
                assert A["df"][1][1] == 1936 and A["df"][2][2] == 16
                assert D["df"] == {"screen": 780, "lead": 80, "total": 860}
    assert eT < 1e-13 and eD < 1e-13 and e2 < 1e-13, (eT, eD, e2)
    ok["thmAA_T_split_max_abs_err"] = eT
    ok["thmAA_dual_transpose_max_abs_err"] = eD
    ok["thmAA_d2b_split_max_abs_err"] = e2
    # ---- Table 6 cpg split regression ----
    aa6 = aa_array(6)
    g6 = lag_duplex_decomposition(ensG, 64, 3, aa6, duplex_layers(aa6)[3],
                                  cpg=cpg_quartet_partition(aa6))
    assert g6["S_dup_split_abs_err"] < 1e-13, g6["S_dup_split_abs_err"]
    ok["fix51_table6_cpg_split_ok"] = True
    ok["v24_all_passed"] = True
    return ok



# =====================================================================
# 16e. THEOREMS AB, AD, AC  (v3.0) — Periodicity/ZNF, Isochore-Free
#      Screen KL Chain Rule, and Zero-RNG Biosphere Phylogeny Matrix
# =====================================================================
_ZNF_C2H2_RE = re.compile(r"C.{2,4}C.{12}H.{3,5}H")
_COLLAGEN_RE = re.compile(r"(?:G..){6,}")


def _translate_gene(s, aa):
    return "".join(AA_STRING[int(a)] for a in aa[s])


def _klein_directed_series(seqs, kmax):
    """Exact directed Klein energies ||W_-+||^2 and ||W_--||^2 over k=1..kmax."""
    mp = np.full(kmax, np.nan, dtype=np.float64)
    mm = np.full(kmax, np.nan, dtype=np.float64)
    pm = np.full(kmax, np.nan, dtype=np.float64)
    for k in range(1, kmax + 1):
        B = lag_counts(seqs, NC, k)
        if B.sum() <= 0:
            continue
        qM = theorem_M_klein_quadrants(B)
        pm[k - 1] = qM["norm2_W_pm_static_strand"]
        mp[k - 1] = qM["norm2_W_mp_sym_time_flux"]
        mm[k - 1] = qM["norm2_W_mm_duplex_time_flux"]
    return pm, mp, mm


def _dip14_and_dft(series, k_start=5):
    """Exact log10 dip depth at k in {14,28,42,56} vs +/-{1,2} flanks + DFT peak."""
    arr = np.asarray(series, dtype=np.float64)
    kmax = arr.size
    logv = np.where(np.isfinite(arr) & (arr > 0), np.log10(np.maximum(arr, 1e-30)), np.nan)
    dips = []
    for k0 in (14, 28, 42, 56):
        if k0 + 2 <= kmax and k0 - 2 >= 1:
            c = logv[k0 - 1]
            fl = [logv[j - 1] for j in (k0 - 2, k0 - 1, k0 + 1, k0 + 2) if np.isfinite(logv[j - 1])]
            if np.isfinite(c) and len(fl) == 4:
                dips.append(float(c - np.mean(fl)))
    dip14 = float(np.mean(dips)) if dips else None
    sub = logv[k_start - 1:]
    valid = np.isfinite(sub)
    peak_period, peak_snr, pwr14, pwr28 = None, None, None, None
    if int(valid.sum()) >= 20:
        x = np.arange(sub.size, dtype=np.float64)
        y = np.interp(x, x[valid], sub[valid])
        # Linear detrend (deterministic least squares)
        slope, intercept = np.polyfit(x, y, 1)
        yd = y - (slope * x + intercept)
        # Exact DFT evaluated on periods P in [3..30]
        n = yd.size
        periods = np.arange(3.0, min(31.0, n / 1.5), 0.5)
        t = np.arange(n, dtype=np.float64)
        pwrs = np.array([
            float(np.abs(np.sum(yd * np.exp(-2j * math.pi * t / P))) ** 2) / n
            for P in periods
        ])
        med = float(np.median(pwrs))
        ipk = int(np.argmax(pwrs))
        peak_period = float(periods[ipk])
        peak_snr = float(pwrs[ipk] / max(med, 1e-30))
        i14 = int(np.argmin(np.abs(periods - 14.0)))
        pwr14 = float(pwrs[i14] / max(med, 1e-30))
        if np.any(np.abs(periods - 28.0) < 0.6):
            i28 = int(np.argmin(np.abs(periods - 28.0)))
            pwr28 = float(pwrs[i28] / max(med, 1e-30))
    return {"dip14_log10": dip14, "dips_by_harmonic": dips,
            "dft_peak_period": peak_period, "dft_peak_snr": peak_snr,
            "dft_snr_at_14": pwr14, "dft_snr_at_28": pwr28}


def thm_AB_periodicity_znf_audit(test_seqs, aa, lag_profile=None, kmax=60):
    """[THM-AB] Exact zero-RNG periodicity & C2H2 Zinc-Finger stratification audit.
    Directly settles Paper III Policy 10.1 without external annotations."""
    znf_seqs, non_znf_seqs = [], []
    n_collagen = 0
    tot_codons = znf_codons = 0
    for s in test_seqs:
        prot = _translate_gene(s, aa)
        tot_codons += int(s.size)
        if _COLLAGEN_RE.search(prot):
            n_collagen += 1
        if len(_ZNF_C2H2_RE.findall(prot)) >= 2:
            znf_seqs.append(s)
            znf_codons += int(s.size)
        else:
            non_znf_seqs.append(s)

    if lag_profile is not None and len(lag_profile) >= kmax:
        mp_all = np.array([r["theorem_M"]["norm2_W_mp_sym_time_flux"] if r.get("theorem_M") else np.nan
                           for r in lag_profile[:kmax]], dtype=np.float64)
        mm_all = np.array([r["theorem_M"]["norm2_W_mm_duplex_time_flux"] if r.get("theorem_M") else np.nan
                           for r in lag_profile[:kmax]], dtype=np.float64)
    else:
        _, mp_all, mm_all = _klein_directed_series(test_seqs, kmax)

    st_mp_all = _dip14_and_dft(mp_all)
    st_mm_all = _dip14_and_dft(mm_all)

    strat = None
    verdict = "NO_SIGNIFICANT_DIP14"
    if len(znf_seqs) >= 3 and len(non_znf_seqs) >= 10:
        _, mp_non, mm_non = _klein_directed_series(non_znf_seqs, kmax)
        _, mp_znf, mm_znf = _klein_directed_series(znf_seqs, kmax)
        st_mp_non = _dip14_and_dft(mp_non)
        st_mm_non = _dip14_and_dft(mm_non)
        st_mp_znf = _dip14_and_dft(mp_znf)
        st_mm_znf = _dip14_and_dft(mm_znf)
        d_all = st_mp_all["dip14_log10"]
        d_non = st_mp_non["dip14_log10"]
        expl = (1.0 - d_non / d_all) if (d_all is not None and d_non is not None and abs(d_all) > 1e-6) else None
        if d_all is not None and d_all < -0.15:
            if d_non is not None and abs(d_non) < 0.35 * abs(d_all):
                verdict = "ZNF_C2H2_DRIVEN_DIP14"
            else:
                verdict = "NON_ZNF_PERIODIC_DIP14"
        elif st_mp_znf["dip14_log10"] is not None and st_mp_znf["dip14_log10"] < -0.25:
            verdict = "ZNF_SUBCOHORT_DIP14_PRESENT"
        strat = {
            "non_znf_mp": st_mp_non, "non_znf_mm": st_mm_non,
            "znf_only_mp": st_mp_znf, "znf_only_mm": st_mm_znf,
            "dip14_explained_fraction_mp": expl,
            "mp_all_series": mp_all.tolist(),
            "mp_non_znf_series": mp_non.tolist(),
            "mp_znf_only_series": mp_znf.tolist(),
        }
    elif st_mp_all["dip14_log10"] is not None and st_mp_all["dip14_log10"] < -0.15:
        verdict = "NON_ZNF_PERIODIC_DIP14"

    return {
        "kmax": int(kmax),
        "n_genes_total": len(test_seqs),
        "n_znf_c2h2_genes": len(znf_seqs),
        "znf_gene_fraction": len(znf_seqs) / max(1, len(test_seqs)),
        "znf_codon_fraction": znf_codons / max(1, tot_codons),
        "n_collagen_gxy_genes": int(n_collagen),
        "all_mp": st_mp_all, "all_mm": st_mm_all,
        "stratified": strat,
        "verdict": verdict,
    }


def screen_kl_chain_rule(B, E):
    """[THM-AD] EXACT forward & reverse 3-term screen-stratified KL chain rules:
        D_KL(P || Q) = K_bridge + K_lead|SO + K_wob|XO
                     = K_bridge + K_wob|SO  + K_lead|SY
    All terms >= 0. Reduces identically to Theorem T when Q = p_X x p_Y.
    Under antiparallel duplex inversion R:
        K_bridge invariant; K_lead|SO <-> K_wob|SO; K_wob|XO <-> K_lead|SY.
    """
    N1, N2 = float(np.sum(B)), float(np.sum(E))
    if N1 <= 0 or N2 <= 0:
        return None
    P = np.asarray(B, dtype=np.float64) / N1
    Q = np.maximum(np.asarray(E, dtype=np.float64) / N2, TINY)
    Q = Q / Q.sum()
    xs = np.arange(NC)

    def _projs(M):
        Mso = np.bincount((SUF16[:, None] * 16 + OUT16[None, :]).ravel(),
                          weights=M.ravel(), minlength=256).reshape(16, 16)
        Mxo = np.bincount((xs[:, None] * 16 + OUT16[None, :]).ravel(),
                          weights=M.ravel(), minlength=NC * 16).reshape(NC, 16)
        Msy = np.bincount((SUF16[:, None] * NC + xs[None, :]).ravel(),
                          weights=M.ravel(), minlength=16 * NC).reshape(16, NC)
        return Mso, Mxo, Msy

    Pso, Pxo, Psy = _projs(P)
    Qso, Qxo, Qsy = _projs(Q)
    kl_tot = _kl_joint(P, Q)
    K_br = _kl_joint(Pso, np.maximum(Qso, TINY))

    # Forward: SO -> X1 -> Y3
    m_xo = Pxo > 0
    den_l = np.where(m_xo, Qxo * Pso[SUF16, :] / np.maximum(Qso[SUF16, :], TINY), TINY)
    K_lead_SO = float(np.sum(Pxo[m_xo] * np.log(Pxo[m_xo] / np.maximum(den_l[m_xo], TINY))))
    m_xy = P > 0
    den_wx = np.where(m_xy, Q * Pxo[:, OUT16] / np.maximum(Qxo[:, OUT16], TINY), TINY)
    K_wob_XO = float(np.sum(P[m_xy] * np.log(P[m_xy] / np.maximum(den_wx[m_xy], TINY))))

    # Reverse: SO -> Y3 -> X1
    m_sy = Psy > 0
    den_w = np.where(m_sy, Qsy * Pso[:, OUT16] / np.maximum(Qso[:, OUT16], TINY), TINY)
    K_wob_SO = float(np.sum(Psy[m_sy] * np.log(Psy[m_sy] / np.maximum(den_w[m_sy], TINY))))
    den_ly = np.where(m_xy, Q * Psy[SUF16, :] / np.maximum(Qsy[SUF16, :], TINY), TINY)
    K_lead_SY = float(np.sum(P[m_xy] * np.log(P[m_xy] / np.maximum(den_ly[m_xy], TINY))))

    return {
        "KL_total": kl_tot,
        "K_bridge": K_br,
        "K_lead_given_SO": K_lead_SO,
        "K_wob_given_XO": K_wob_XO,
        "K_wob_given_SO": K_wob_SO,
        "K_lead_given_SY": K_lead_SY,
        "K_nonbridge": kl_tot - K_br,
        "nonbridge_share": ((kl_tot - K_br) / kl_tot) if kl_tot > 1e-15 else None,
        "fwd_residual": abs(kl_tot - (K_br + K_lead_SO + K_wob_XO)),
        "rev_residual": abs(kl_tot - (K_br + K_wob_SO + K_lead_SY)),
        "all_nonnegative": bool(min(K_br, K_lead_SO, K_wob_XO, K_wob_SO, K_lead_SY) > -1e-12),
    }


def thm_AD_isochore_screen_audit(test_seqs, lags, aa):
    """[THM-AD] Isochore-free screen decomposition against exact zero-RNG E_het^(k).
    Computes:
      * Lemma 7.1 Miller-Madow bias-corrected screen MI layers on B^(k),
      * Exact screen MI layers on the between-gene heterogeneity null E_het^(k),
      * Isochore-stripped within-gene excess l_intra^MM = l_MM(B) - l(E_het),
      * Exact 3-term forward/reverse screen-stratified KL(B^(k) || E_het^(k)).
    """
    rows = []
    for k in sorted({int(x) for x in lags if int(x) >= 1}):
        B = lag_counts(test_seqs, NC, k)
        N = float(B.sum())
        if N <= 0:
            continue
        E = het_expected_counts(test_seqs, NC, k)
        spB = two_sided_screen_split(B, aa)
        spE = two_sided_screen_split(E, aa)
        skl = screen_kl_chain_rule(B, E)
        df = spB["df_eff"]
        mm_br = spB["I_bridge"] - df["bridge"] / (2.0 * N)
        mm_ld = spB["D_5prime_lead"] - df["lead"] / (2.0 * N)
        mm_wb = spB["D_3prime_wobble"] - df["wobble"] / (2.0 * N)
        mm_13 = spB["D_1_3"] - df["d13"] / (2.0 * N)
        mm_dx1 = mm_ld + mm_13
        df_nonbr = int(df["lead"] + df["wobble"] + df["d13"])
        df_dx1 = int(df["lead"] + df["d13"])
        skl["df_nonbridge_eff"] = df_nonbr
        skl["df_lead_given_SY_eff"] = df_dx1
        skl["K_nonbridge_MM"] = skl["K_nonbridge"] - df_nonbr / (2.0 * N)
        skl["K_lead_given_SY_MM"] = skl["K_lead_given_SY"] - df_dx1 / (2.0 * N)
        rows.append({
            "lag": k, "N": int(N),
            # Bias-corrected observed screen layers (Lemma 7.1)
            "I_bridge_MM": mm_br, "D_5prime_lead_MM": mm_ld,
            "D_3prime_wobble_MM": mm_wb, "D_1_3_MM": mm_13, "Delta_x1_MM": mm_dx1,
            # Exact between-gene heterogeneity null screen layers (no MM subtraction!)
            "I_bridge_het": spE["I_bridge"], "D_5prime_lead_het": spE["D_5prime_lead"],
            "D_3prime_wobble_het": spE["D_3prime_wobble"], "D_1_3_het": spE["D_1_3"],
            "Delta_x1_het": spE["Delta_x1"],
            # Isochore-stripped intra-gene excess
            "I_bridge_intra_MM": mm_br - spE["I_bridge"],
            "D_5prime_lead_intra_MM": mm_ld - spE["D_5prime_lead"],
            "D_3prime_wobble_intra_MM": mm_wb - spE["D_3prime_wobble"],
            "D_1_3_intra_MM": mm_13 - spE["D_1_3"],
            "Delta_x1_intra_MM": mm_dx1 - spE["Delta_x1"],
            "wobble_isochore_share": (spE["D_3prime_wobble"] / mm_wb) if mm_wb > 1e-12 else None,
            "delta_x1_isochore_share": (spE["Delta_x1"] / mm_dx1) if mm_dx1 > 1e-12 else None,
            # Exact screen-stratified KL(B || E_het)
            "screen_kl": skl,
        })
    return rows


def thm_AC_biosphere_phylogeny(species_summaries, matrix_cache):
    """[THM-AC] Exact zero-RNG cross-species Klein-decomposed phylogeny & multi-code map.
    Compares Leave-One-Out 1-NN clade purity across 4 orthogonal channels:
      (1) scalar GC3, (2) static V4 form (W_++ + W_+-),
      (3) directed V4 program (W_-+ + W_--), (4) 93.75% non-bridge screen complement.
    """
    valid = [r for r in species_summaries if r["species"] in matrix_cache]
    S = len(valid)
    if S < 3:
        return None
    names = [r["species"] for r in valid]
    clades = [r.get("clade", NAME2CLADE.get(r["species"], "unknown")) for r in valid]
    gc3 = np.array([float(r.get("gc3_test_mean", 0.5)) for r in valid], dtype=np.float64)
    tables = [int(r.get("transl_table", 1)) for r in valid]
    DS = _layer_avg_matrix(SUF16, 16)
    DO = _layer_avg_matrix(OUT16, 16)

    Ps = [matrix_cache[nm]["P1_test"] for nm in names]
    # Precompute orthogonal projections for each species
    P_stat, P_dir, P_nonbr = [], [], []
    for P in Ps:
        TP = P.T
        RP = P.T[np.ix_(MIRROR, MIRROR)]
        SP = P[np.ix_(MIRROR, MIRROR)]
        P_stat.append(0.5 * (P + TP))         # W_++ + W_+-
        P_dir.append(0.5 * (P - TP))          # W_-+ + W_--
        P_nonbr.append(P - DS @ P @ DO)       # 93.75% H3b-invisible complement

    D_gc3 = np.abs(gc3[:, None] - gc3[None, :])
    D_stat = np.zeros((S, S))
    D_dir = np.zeros((S, S))
    D_nonbr = np.zeros((S, S))
    D_dir_cos = np.zeros((S, S))
    D_nonbr_cos = np.zeros((S, S))
    P_dir_u = [M / max(float(np.linalg.norm(M, ord="fro")), TINY) for M in P_dir]
    P_nonbr_u = [M / max(float(np.linalg.norm(M, ord="fro")), TINY) for M in P_nonbr]
    for i in range(S):
        for j in range(i + 1, S):
            ds = float(np.linalg.norm(P_stat[i] - P_stat[j], ord="fro"))
            dd = float(np.linalg.norm(P_dir[i] - P_dir[j], ord="fro"))
            dn = float(np.linalg.norm(P_nonbr[i] - P_nonbr[j], ord="fro"))
            ddc = 1.0 - float(np.sum(P_dir_u[i] * P_dir_u[j]))
            dnc = 1.0 - float(np.sum(P_nonbr_u[i] * P_nonbr_u[j]))
            D_stat[i, j] = D_stat[j, i] = ds
            D_dir[i, j] = D_dir[j, i] = dd
            D_nonbr[i, j] = D_nonbr[j, i] = dn
            D_dir_cos[i, j] = D_dir_cos[j, i] = ddc
            D_nonbr_cos[i, j] = D_nonbr_cos[j, i] = dnc

    def _loo_1nn(D):
        hits = elig = 0
        clade_counts = {}
        for c in clades:
            clade_counts[c] = clade_counts.get(c, 0) + 1
        nn_list = []
        for i in range(S):
            row = D[i].copy()
            row[i] = np.inf
            j = int(np.argmin(row))
            nn_list.append(names[j])
            if clade_counts[clades[i]] >= 2:
                elig += 1
                if clades[j] == clades[i]:
                    hits += 1
        return {"hits": hits, "eligible": elig,
                "accuracy": (hits / elig) if elig > 0 else None,
                "nearest_neighbors": nn_list}

    # Exact cross-species held-out d2b_transfer (zero RNG!) for all same-table pairs
    d2b_mat = np.full((S, S), np.nan)
    for i in range(S):
        lp_src = matrix_cache[names[i]]["logP_train"]
        lq2b_src = matrix_cache[names[i]]["logQ2b_train"]
        diff_src = lp_src - lq2b_src
        for j in range(S):
            if tables[i] != tables[j]:
                continue
            d2b_mat[i, j] = float(np.sum(Ps[j] * diff_src))

    # Multi-table genetic code fingerprint summary (Theorem Q.2)
    fp_groups = {}
    for r in valid:
        q = (r.get("theorem_Q") or {}).get("AA", {})
        c_sq = round(float(q.get("commutator_R_fro", 0.0) ** 2), 6)
        hs_sq = round(float(q.get("superoperator_commutator_HS_sq", 0.0)), 6)
        key = f"F2={c_sq:.6f}|HS2={hs_sq:.6f}"
        fp_groups.setdefault(key, []).append(
            {"species": r["species"], "table": r["transl_table"], "clade": r.get("clade")})

    return {
        "n_species": S,
        "species": names,
        "clades": clades,
        "tables": tables,
        "gc3": gc3.tolist(),
        "loo_1nn_clade_purity": {
            "GC3_scalar": _loo_1nn(D_gc3),
            "Static_Form_Wpp_Wpm": _loo_1nn(D_stat),
            "Directed_Program_Wmp_Wmm": _loo_1nn(D_dir),
            "Directed_Program_Cosine": _loo_1nn(D_dir_cos),
            "NonBridge_93p75_Screen": _loo_1nn(D_nonbr),
            "NonBridge_Cosine": _loo_1nn(D_nonbr_cos),
        },
        "code_fingerprint_groups": fp_groups,
        "d2b_transfer_matrix": d2b_mat.tolist(),
        "dist_static_matrix": D_stat.tolist(),
        "dist_directed_matrix": D_dir.tolist(),
    }


def self_test_v30(ensG=None, aa1=None):
    ok = {}
    if aa1 is None:
        aa1 = AA_TABLE1
    # 1. Verify BIOSPHERE_100 catalog integrity (exact count = 100, 5 tables)
    assert len(BIOSPHERE_100) == 100, len(BIOSPHERE_100)
    tabs_present = sorted({v[1] for v in BIOSPHERE_100.values()})
    assert tabs_present == [1, 4, 6, 11, 12], tabs_present
    ok["biosphere_100_count"] = len(BIOSPHERE_100)
    ok["biosphere_100_tables"] = tabs_present

    # 2. [THM-AD] Exact forward/reverse 3-term screen KL chain rule + duality + P_ind reduction
    if ensG is not None:
        Bk = lag_counts(ensG, NC, 2)
        Ek = het_expected_counts(ensG, NC, 2)
        skl = screen_kl_chain_rule(Bk, Ek)
        assert skl["all_nonnegative"]
        assert skl["fwd_residual"] < 1e-14, skl["fwd_residual"]
        assert skl["rev_residual"] < 1e-14, skl["rev_residual"]
        # Antiparallel duplex duality under R
        Bka = Bk.T[np.ix_(MIRROR, MIRROR)]
        Eka = Ek.T[np.ix_(MIRROR, MIRROR)]
        skla = screen_kl_chain_rule(Bka, Eka)
        e_ad_dual = max(
            abs(skl["K_bridge"] - skla["K_bridge"]),
            abs(skl["K_lead_given_SO"] - skla["K_wob_given_SO"]),
            abs(skl["K_wob_given_SO"] - skla["K_lead_given_SO"]),
            abs(skl["K_wob_given_XO"] - skla["K_lead_given_SY"]),
            abs(skl["K_lead_given_SY"] - skla["K_wob_given_XO"]),
        )
        assert e_ad_dual < 1e-13, e_ad_dual
        # Reduction to Theorem T when Q = P_ind = p_X x p_Y
        Pind = (Bk.sum(1)[:, None] * Bk.sum(0)[None, :]) / Bk.sum()
        skl_ind = screen_kl_chain_rule(Bk, Pind)
        spT = two_sided_screen_split(Bk, aa1)
        e_ad_red = max(
            abs(skl_ind["K_bridge"] - spT["I_bridge"]),
            abs(skl_ind["K_lead_given_SO"] - spT["D_5prime_lead"]),
            abs(skl_ind["K_wob_given_SO"] - spT["D_3prime_wobble"]),
            abs(skl_ind["K_lead_given_SY"] - spT["Delta_x1"]),
            abs(skl_ind["K_wob_given_XO"] - (spT["D_3prime_wobble"] + spT["D_1_3"])),
        )
        assert e_ad_red < 1e-13, e_ad_red
        ad_rows = thm_AD_isochore_screen_audit(ensG, [1, 5, 30], aa1)
        assert len(ad_rows) == 3
        ok["thmAD_kl_chain_max_err"] = max(skl["fwd_residual"], skl["rev_residual"])
        ok["thmAD_duplex_dual_max_err"] = e_ad_dual
        ok["thmAD_Pind_reduction_max_err"] = e_ad_red

        # 3. [THM-AB] Deterministic ZNF C2H2 periodic repeat detector & stratification test
        # Build 4 synthetic ZNF genes with 28-aa tandem C2H2 arrays + 12 background genes
        znf_aa = "CPYCGKRFTRSDNLQRHTRTHTGEKPYK"  # exact 28-aa C2H2 zinc-finger repeat!
        aa2cod = {AA_STRING[a]: int(np.flatnonzero(aa1 == a)[0]) for a in range(NAA)}
        znf_unit = np.array([aa2cod[ch] for ch in znf_aa], dtype=np.int64)
        znf_gene = np.tile(znf_unit, 8)  # 224 codons = 8 tandem 28-aa fingers
        bg_genes = [ensG[0][i * 200:(i + 1) * 200] for i in range(12)]
        ab = thm_AB_periodicity_znf_audit(bg_genes + [znf_gene] * 4, aa1, kmax=60)
        assert ab["n_znf_c2h2_genes"] == 4, ab["n_znf_c2h2_genes"]
        assert ab["stratified"] is not None
        assert ab["stratified"]["znf_only_mp"]["dip14_log10"] is not None
        ok["thmAB_znf_detected"] = ab["n_znf_c2h2_genes"]
        ok["thmAB_verdict_fixture"] = ab["verdict"]
        ok["thmAB_znf_dft_peak_period"] = ab["stratified"]["znf_only_mp"]["dft_peak_period"]

    ok["v30_all_passed"] = True
    return ok


def self_test(verbose=False):
    ok = {}
    assert len(CODONS) == 64 and np.all(MIRROR[MIRROR] == np.arange(64))
    assert np.all(MIRROR != np.arange(64)) and len(set(MIRROR_PAIR.tolist())) == 32
    assert AA_STRING[aa_array(2)[CODON_IDX["AGA"]]] == "*"
    assert AA_STRING[aa_array(12)[CODON_IDX["CTG"]]] == "S"
    assert np.array_equal(aa_array(11), aa_array(1))
    aa1 = aa_array(1)
    fib1 = fibers_of(aa1)
    layD, FD, szD, d8 = duplex_layers(aa1)
    assert int(d8.sum()) == 16 and int((szD == 1).sum()) == 48 and FD == 56
    ok["genetic_code_mirror_duplex"] = True

    cpg = cpg_quartet_partition(aa1)
    assert cpg["n_cpg_doublets"] == 4 and cpg["n_leu_doublets"] == 4
    assert cpg["n_mixed_doublets"] == 0
    assert cpg["cpg_quartet_rho_invariant"] and cpg["leu_quartet_rho_invariant"]
    assert {CODONS[i] for i in np.flatnonzero(cpg["cpg_mask"])} == {"TCG", "TCT", "CCG", "CCT", "CGA", "AGA", "CGG", "AGG"}
    ok["thmJ_cpg_leu_quartets"] = True

    dfF = duplex_df_table(aa1)
    assert (dfF["df_AA"], dfF["df_d2b"], dfF["df_anti"], dfF["df_dup"]) == (400, 860, 2205, 504)
    assert dfF["df_sum_check"] == 3969
    dfS = duplex_df_table(aa1, support=(aa1 != STOP))
    assert (dfS["df_AA"], dfS["df_d2b"], dfS["df_anti"], dfS["df_dup"]) == (361, 779, 2040, 420)
    assert dfS["df_sum_check"] == 3600
    ok["thmG2_df_sense"] = {k: dfS[k] for k in ("df_AA", "df_d2b", "df_anti", "df_dup")}

    pP = {}
    for tid in sorted(SUPPORTED_TABLES):
        rP = theorem_P_transversality(aa_array(tid))
        pP[str(tid)] = rP["n_joint_AA_mirror_classes"]
        if tid == 5:
            assert rP["n_joint_AA_mirror_classes"] == 63, rP
            assert rP["n_self_complementary_aa_mirror_pairs"] == 1, rP
            assert set(rP["self_complementary_pairs"][0][:2]) == {"TCT", "AGA"}, rP
        else:
            assert rP["join_is_identity"], (tid, rP)
            assert rP["relic_diagonal_total"] == 0, (tid, rP)
    ok["thmP_transversality_by_table"] = pP

    qQ = klein_commutator_audit(aa1)
    assert qQ["AA"]["commutator_R_fro"] > 1e-6
    assert qQ["dup"]["commutes_with_V4"] and qQ["mirror"]["commutes_with_V4"]
    for nm in ("AA", "dup", "mirror"):
        assert qQ[nm]["commutator_identity_abs_err"] < 1e-10, qQ[nm]
        assert qQ[nm]["commutator_T_fro"] < 1e-12, qQ[nm]
        assert qQ[nm]["quadrant_pythagoras_max_err"] < 1e-12, qQ[nm]
        C_ = qQ[nm]["n_layers"]
        S2_ = qQ[nm]["S2_trace_DRDR"]
        assert abs(qQ[nm]["superoperator_commutator_HS_sq"]
                   - (2 * C_ * C_ - 2 * S2_ * S2_)) < 1e-9
    assert abs(qQ["AA"]["S2_trace_DRDR"] - 569.0 / 72.0) < 1e-12
    assert abs(qQ["AA"]["commutator_R_fro"] ** 2 - 943.0 / 36.0) < 1e-11
    assert abs(qQ["AA"]["superoperator_commutator_HS_sq"] - 1962383.0 / 2592.0) < 1e-8
    assert abs(qQ["AA"]["S1_linear_weight_sum"] - 509.0 / 72.0) < 1e-12
    assert qQ["duplex_rho_well_defined"]
    assert qQ["dims_image_sum"] == 56 * 56
    assert qQ["dims_kernel_sum"] == 64 * 64 - 56 * 56
    assert qQ["dims_duplex_image"] == {"W_++": 812, "W_+-": 784, "W_-+": 784, "W_--": 756}
    assert qQ["dims_D8_kernel"] == {"W_++": 244, "W_+-": 240, "W_-+": 240, "W_--": 236}
    ok["thmQ_commutators"] = {nm: qQ[nm]["commutator_R_fro"] for nm in ("AA", "dup", "mirror")}
    ok["thmQ2_AA_invariants"] = {
        "C": qQ["AA"]["n_layers"], "S2": qQ["AA"]["S2_trace_DRDR"],
        "S1": qQ["AA"]["S1_linear_weight_sum"],
        "comm_F_sq": qQ["AA"]["commutator_R_fro"] ** 2,
        "superop_HS_sq": qQ["AA"]["superoperator_commutator_HS_sq"]}
    ok["thmQ_dims"] = {"image": qQ["dims_duplex_image"], "kernel": qQ["dims_D8_kernel"]}

    dims, tr = klein_dims()
    assert tr == {"id": 4096, "T": 64, "R": 64, "S": 0}, tr
    assert (dims["W_++"], dims["W_+-"], dims["W_-+"], dims["W_--"]) == (1056, 1024, 1024, 992)
    assert dims["sum"] == 4096
    ok["thmM_klein_traces"] = tr
    ok["thmM_klein_dims"] = dims
    ix = np.arange(64)
    raw = (1 + ((ix[:, None] * 17 + ix[None, :] * 31 + (ix[:, None] ^ ix[None, :])) % 97)).astype(float)
    P = raw / raw.sum(axis=1, keepdims=True)
    logP = np.log(P)
    Bfx = bigram_counts([np.array([(7 * i + 3 * i * i + 11) % 64 for i in range(400)],
                                  dtype=np.int64)], 64)
    qM = theorem_M_klein_quadrants(Bfx)
    assert qM["pythagorean_residual"] < 1e-15
    assert qM["time_flux_split_residual"] < 1e-15
    assert qM["mirror_flux_split_residual"] < 1e-15
    assert qM["reconstruction_max_abs_err"] < 1e-15
    ok["thmM_pythagorean_residual"] = qM["pythagorean_residual"]

    # ---- Theorem A exhaustive permutations ----
    n = 23
    s = np.array([(7 * i + 3 * i * i + 11) % 64 for i in range(n)], dtype=np.int64)
    max_err, ncase, ebl = 0.0, 0, 0.0
    for L in (2, 3, 4, 5, 7):
        m = n // L
        if m < 2 or m > 6:
            continue
        ana = analytic_h4([s], logP, L)["ll_shuf_expected"] * (n - 1)
        tot, cnt = 0.0, 0
        for perm in itertools.permutations(range(m)):
            out = s[:m * L].reshape(m, L)[list(perm)].reshape(-1)
            if n > m * L:
                out = np.concatenate([out, s[m * L:]])
            tot += float(np.sum(logP[out[:-1], out[1:]]))
            cnt += 1
        max_err = max(max_err, abs(ana - tot / cnt))
        ncase += 1
        ends, starts = s[L - 1:m * L:L], s[0:m * L:L]
        Mr = logP[np.ix_(ends, starts)]
        vm = (float(Mr.sum()) - float(np.trace(Mr))) / m
        ce = np.bincount(ends, minlength=64).astype(float)
        cs = np.bincount(starts, minlength=64).astype(float)
        vb = (float(ce @ logP @ cs) - float(np.sum(logP[ends, starts]))) / m
        ebl = max(ebl, abs(vm - vb))
    assert ncase >= 3 and max_err < 1e-9 and ebl < 1e-12
    ok["theorem_A_exhaustive_permutations_max_abs_err"] = max_err

    def brute_j(nn, LL):
        mm = nn // LL
        if mm < 2:
            return 0
        idx = [k * LL - 1 for k in range(1, mm)]
        if nn - mm * LL > 0:
            idx.append(mm * LL - 1)
        return len(idx)

    assert max(abs(junction_count(nn, LL) - brute_j(nn, LL))
               for nn in range(2, 200) for LL in range(2, 40)) == 0
    pf = psi_constant([s], logP)
    assert pf["identity_abs_err"] < 1e-12 and pf["edge_abs_err"] < 1e-12
    Bf = bigram_counts([s], 64)
    Ij, Jrev, J = jeffreys_of_counts(Bf)
    pf_emp = psi_constant([s], np.log(np.maximum(Bf / Bf.sum(), TINY)))
    assert abs(pf_emp["Psi_LR"] - J) < 1e-10
    f3 = np.array([1.0, 1.0, -2.0])
    Pc = (1.0 + 0.1 * np.outer(f3, f3)) / 9.0
    Ic, _, Jc = jeffreys_of_counts(Pc * 1e6)
    assert Jc < 2 * Ic
    ok["thmA4_J_ge_2I_is_FALSE"] = {"J": Jc, "two_I": 2 * Ic}

    cdf = np.cumsum(P, axis=1); cdf[:, -1] = 1.0
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    u = (np.arange(SELFTEST_LONG_N) * phi) % 1.0
    sL = np.empty(SELFTEST_LONG_N, dtype=np.int64); sL[0] = 0
    for t in range(1, SELFTEST_LONG_N):
        sL[t] = int(np.searchsorted(cdf[sL[t - 1]], u[t], side="right"))
    Bl = bigram_counts([sL], 64)
    lpL = np.log((Bl + 0.5) / (Bl.sum(axis=1, keepdims=True) + 0.5 * 64))
    psiL = psi_constant([sL], lpL)["Psi_LR"]
    worst = 0.0
    for LL in (2, 3, 5, 10, 20, 29, 30, 50, 100):
        r = analytic_h4([sL], lpL, LL)
        worst = max(worst, abs(r["delta_H4_analytic"] / (r["gamma_L"] * psiL) - 1.0))
    assert worst < 0.02
    ok["theorem_A2_Rjunc_max_dev_from_1"] = worst

    ens = [sL[:400], sL[400:1100], s, sL[1100:1300]]
    psiE = psi_constant(ens, logP)
    eE = eW = 0.0
    for LL in (2, 3, 5, 10, 30):
        rec = analytic_h4(ens, logP, LL, collect_terms=True)
        tE = theorem_E(ens, logP, LL, psiE, rec)
        if tE is None:
            continue
        eE = max(eE, tE["identity_abs_err"]); eW = max(eW, tE["within_split_abs_err"])
    assert eE < 1e-12 and eW < 1e-12
    ok["thmE_identity_max_abs_err"] = eE

    q4 = 4
    s4 = np.array([0, 0, 1, 2, 2], dtype=np.int64)
    eh = 0.0
    for k in (1, 2, 3):
        B = np.zeros((q4, q4))
        for perm in itertools.permutations(range(s4.size)):
            p = s4[list(perm)]
            B += np.bincount(p[:-k] * q4 + p[k:], minlength=q4 * q4).reshape(q4, q4)
        B /= math.factorial(s4.size)
        eh = max(eh, float(np.abs(B - het_expected_counts([s4], q4, k)).max()))
    assert eh < 1e-12
    ok["theorem_C2_het_null_vs_all_permutations_max_err"] = eh

    ensG = [sL[:3000], sL[3000:7000]]
    eLsym = 0.0
    for k in (1, 2, 5, 17, 40):
        Bk = lag_counts(ensG, 64, k)
        Ek = het_expected_counts(ensG, 64, k)
        tl = theorem_L_diagnostics(Bk, Ek)
        eLsym = max(eLsym, tl["null_symmetry_max_abs_err"])
        assert tl["null_is_symmetric"] and tl["JS_ge_half_A2"]
        qm = theorem_M_klein_quadrants(Ek)
        assert qm["norm2_W_mp_sym_time_flux"] < 1e-24
        assert qm["norm2_W_mm_duplex_time_flux"] < 1e-24
    assert eLsym < 1e-9
    ok["thmL_null_symmetry_max_abs_err"] = eLsym
    Bq = lag_counts(ensG, 64, 2)
    Eq = het_expected_counts(ensG, 64, 2)
    kr = klein_quadrants_of_matrix(Bq / Bq.sum() - Eq / Eq.sum())
    assert kr["pythagorean_residual"] < 1e-15
    ok["fix30_signed_klein_ok"] = True

    Pm = np.array([[0.5, 0.3, 0.2], [0.1, 0.6, 0.3], [0.4, 0.2, 0.4]])
    pim = stationary(Pm)
    J3 = np.einsum("x,xy,yz->xyz", pim, Pm, Pm); J3 /= J3.sum()
    I2 = 0.0
    for y in range(3):
        py = J3[:, y, :].sum()
        for x in range(3):
            for z in range(3):
                p = J3[x, y, z]
                if p > 0:
                    I2 += p * math.log(p / (py * (J3[x, y, :].sum() / py) * (J3[:, y, z].sum() / py)))
    assert abs(I2) < 1e-14
    perm = np.array([2, 0, 3, 1]); T = np.eye(4)[perm]
    # [FIX-41] T @ diag(d) @ T.T is DIAGONAL -> row-normalise gives I.
    # Reversible mixing chain: P0 = (1-eps)*1*pi^T + eps*I, spec {1,eps,eps,eps}.
    pi0 = np.array([0.4, 0.3, 0.2, 0.1]); eps0 = 0.4
    P0 = (1.0 - eps0) * np.ones((4, 1)) * pi0[None, :] + eps0 * np.eye(4)
    Pperm = T @ P0 @ T.T
    Pperm = Pperm / Pperm.sum(axis=1, keepdims=True)
    prp, specP = markov_lag_prediction(Pperm, 30)
    assert np.max(np.abs(Pperm - np.eye(4))) > 0.1, 'degenerate fixture'
    assert specP['is_reversible'], specP['reversibility_defect_L1']
    assert abs(specP['lambda2_modulus'] - eps0) < 1e-9, specP['lambda2_modulus']
    assert prp[29]["I_reversible_spectral_sum"] < 1e-6
    Jsym = np.array([[9., 2., 1., .5], [2., 7., 1.5, 1.], [1., 1.5, 8., 2.], [.5, 1., 2., 6.]])
    Prev = Jsym / Jsym.sum(axis=1, keepdims=True)
    predR, specR = markov_lag_prediction(Prev, 40)
    assert specR["is_reversible"] and specR["kappa_min"] > 1 - 1e-8
    assert abs(predR[19]["chi2_ratio"] - 1) < 1e-4
    predN, specN = markov_lag_prediction(P, 12)
    assert not specN["is_reversible"] and specN["kappa_min"] > 1 - 1e-7
    assert max(r["biorthogonal_rel_err"] for r in predN) < 1e-6
    ok["thmF2_kappa_min_nonreversible"] = specN["kappa_min"]

    eG_tri = eG_dup = eG_cpg = eG_mir = 0.0
    for k in (1, 2, 5, 17):
        g = lag_duplex_decomposition(ensG, 64, k, aa1, d8, df=dfF, cpg=cpg)
        eG_tri = max(eG_tri, abs(g["tripartite_residual"]))
        eG_dup = max(eG_dup, abs(g["duplex_residual"]))
        eG_cpg = max(eG_cpg, g["S_dup_split_abs_err"])
        eG_mir = max(eG_mir, abs(g["mirror_capture_residual"]), abs(g["mirror_cond_residual"]))
        assert g["D8_bound_ok"] and g["cpg_bound_ok"] and g["leu_bound_ok"]
        assert g["df_anti_eff"] <= g["df_anti_full"]
    assert eG_tri < 1e-13 and eG_dup < 1e-13 and eG_cpg < 1e-13 and eG_mir < 1e-13
    ok["thmG_tri_max_err"] = eG_tri
    ok["thmP_mirror_capture_max_err"] = eG_mir

    lay0 = np.zeros(q4, dtype=np.int64)
    eH = 0.0
    for k in (1, 2, 3):
        E0, B0, _ = layer_tower_expected([s4], q4, k, lay0, 1)
        eH = max(eH, float(np.abs(E0 - het_expected_counts([s4], q4, k)).max()))
    assert eH < 1e-12
    ok["thmH_lambda0_equals_C2_max_abs_err"] = eH
    lay4 = np.array([0, 0, 1, 2])
    s6 = np.array([0, 1, 2, 0, 3], dtype=np.int64)
    groups = [np.flatnonzero(lay4[s6] == a) for a in range(3)]
    eHb = 0.0
    for k in (1, 2):
        Bacc = np.zeros((q4, q4)); cnt = 0
        for perms in itertools.product(*[list(itertools.permutations(g)) for g in groups]):
            out = s6.copy()
            for g, pm in zip(groups, perms):
                out[g] = s6[list(pm)]
            Bacc += np.bincount(out[:-k] * q4 + out[k:], minlength=q4 * q4).reshape(q4, q4)
            cnt += 1
        Bacc /= cnt
        Eh, _, _ = layer_tower_expected([s6], q4, k, lay4, 3)
        eHb = max(eHb, float(np.abs(Bacc - Eh).max()))
    assert eHb < 1e-12
    ok["thmH_layer_null_vs_all_permutations_max_abs_err"] = eHb

    e31 = 0.0
    for lay, F in ((aa1, NAA), (layD, FD), (MIRROR32, 32)):
        Ef, _, _ = layer_tower_expected(ensG[:2], 64, 2, lay, F)
        Er, _ = _layer_tower_reference(ensG[:2], 64, 2, lay, F)
        e31 = max(e31, float(np.abs(Ef - Er).max()))
    assert e31 < 1e-12, e31
    ok["fix31_fast_vs_reference_max_abs_err"] = e31

    proj_errs = {}
    for name, lay, F in (("aa", aa1, NAA), ("dup", layD, FD), ("mirror", MIRROR32, 32)):
        Ek, Bk, meta = layer_tower_expected(ensG, 64, 3, lay, F)
        ep, et = tower_projection_check(Ek, Bk, lay, F)
        rel = ep / max(1.0, float(Bk.max()))
        assert rel < 1e-12, f"projection rel err {rel} for {name}"
        assert meta["no_clipping_required"], f"negative E for {name}"
        h3 = lemma_H3_check(Ek, Bk)
        assert h3["absolute_continuity_holds"]
        proj_errs[name] = {"abs": ep, "rel": rel, "total": et}
        ok[f"lemH3_{name}"] = h3
        kw = kl_four_way(Bk, Ek, aa1)
        assert abs(kw["residual"]) < 1e-12 and kw["all_nonnegative"]
        kn = kl_symmetric_duplex_tower(Bk, Ek, aa1)
        assert abs(kn["residual"]) < 1e-14 and kn["all_nonnegative"]
        if name in ("aa", "dup"):
            assert abs(kw["D_AA"]) < 1e-12
            assert abs(kn["D_AA_sym"]) < 1e-12
            assert kn["pi_AA_max_abs_diff"] < 1e-14
        if name == "dup":
            assert abs(kn["D_anti_sym"]) < 1e-12
            assert kn["pi_dup_max_abs_diff"] < 1e-14
            assert abs(kn["KL_total"] - kn["D_D8_sym"]) < 1e-12
        loc = theorem_K_localisation(Ek, Bk, lay)
        assert loc["KL_localised_on_D8"]
        ok["thmN_dup_annihilation"] = {
            "D_AA_sym": kn["D_AA_sym"], "D_anti_sym": kn["D_anti_sym"],
            "D_D8_sym": kn["D_D8_sym"], "KL_total": kn["KL_total"]}
        ok["thmK_dup_localisation"] = loc
    ok["fix25_projection_errors"] = proj_errs

    Bk = lag_counts(ensG, 64, 2)
    Pind = (Bk.sum(1)[:, None] * Bk.sum(0)[None, :]) / Bk.sum()
    kwI = kl_four_way(Bk, Pind, aa1)
    gI = lag_duplex_decomposition(ensG, 64, 2, aa1, d8)
    errs = [abs(kwI["D_AA"] - gI["I_AA"]), abs(kwI["D_d2b"] - gI["d2b"]),
            abs(kwI["D_anti"] - gI["S_anti"]), abs(kwI["D_dup"] - gI["S_dup"])]
    assert max(errs) < 1e-12
    ok["thmI_independence_reduces_to_MI_max_err"] = max(errs)

    # ================= [FIX-34] de Bruijn THREE-period null =================
    cd_null, dbmeta = debruijn_phase_null(n_periods=3, order=6)
    assert dbmeta["period_len"] == 4096
    assert dbmeta["period_len_mod3"] == 1
    assert dbmeta["n_codons"] == 4098
    oO = theorem_O_frame_phase_audit([cd_null], 1, aa1, return_counts=True)
    assert oO is not None
    assert oO["frame_pair_counts_equal"], oO["frame_pair_counts"]
    assert oO["frame_pair_counts"][0] == 4096, oO["frame_pair_counts"]
    for Bfr in oO["_B_frames"]:
        assert np.all(Bfr == 1.0), "each codon bigram must occur EXACTLY once"
    assert oO["Omega_phase_JS_nats"] < 1e-14, oO["Omega_phase_JS_nats"]
    assert oO["Omega_marg_JS_nats"] < 1e-14
    assert oO["Omega_cond_JS_nats"] < 1e-14
    assert oO["chain_rule_abs_err"] < 1e-13
    assert oO["df_phase_full"] == oO["df_phase_marg"] + oO["df_phase_cond"]
    ok["fix34_debruijn_meta"] = dbmeta
    ok["fix34_debruijn_all_ones_verified"] = True
    ok["thmO_debruijn_null_Omega_phase"] = oO["Omega_phase_JS_nats"]
    cd_bad, _ = debruijn_phase_null(n_periods=2, order=6)
    oBad = theorem_O_frame_phase_audit([cd_bad], 1, aa1)
    assert oBad is not None and oBad["Omega_phase_JS_nats"] > 1e-4
    ok["fix34_two_period_fixture_is_NOT_null"] = oBad["Omega_phase_JS_nats"]
    oReal = theorem_O_frame_phase_audit(ensG, 1, aa1)
    assert oReal["Omega_phase_JS_nats"] > 1e-6
    assert oReal["chain_rule_abs_err"] < 1e-12
    assert oReal["frame_pair_counts_equal"]
    ok["thmO_real_Omega_phase"] = oReal["Omega_phase_JS_nats"]
    ok["thmO_real_Omega_cond"] = oReal["Omega_cond_JS_nats"]
    ok["thmO_real_bias_cond"] = oReal["bias_cond"]

    # ============ [FIX-38] cached frame cohorts == per-call frames ============
    fc = build_frame_cohorts(ensG, aa1)
    e38 = 0.0
    for k in (1, 5, 17):
        o1 = theorem_O_from_cohorts(fc, k, aa1)
        o2 = theorem_O_frame_phase_audit(ensG, k, aa1)
        assert o1["frame_pair_counts_equal"]
        e38 = max(e38, abs(o1["Omega_phase_JS_nats"] - o2["Omega_phase_JS_nats"]),
                  abs(o1["Omega_cond_JS_nats"] - o2["Omega_cond_JS_nats"]))
    assert e38 < 1e-15, e38
    ml = theorem_O_multi_lag(ensG, [1, 5, 17], aa1)
    assert [o["lag"] for o in ml] == [1, 5, 17]
    ok["fix38_cached_frames_match_max_err"] = e38

    # ============ [FIX-37] EXACT nt projection vs naive unpacking ============
    Ca_n, Cp_n = _nt_order2_counts(ensG)
    Ca_p, Cp_p = _nt_order2_counts_projected(ensG)
    e37 = max(float(np.abs(Ca_n - Ca_p).max()), float(np.abs(Cp_n - Cp_p).max()))
    assert e37 == 0.0, e37
    assert abs(float(Cp_p[2].sum()) - sum(s.size for s in ensG)) < 1e-9
    assert abs(float(Cp_p[0].sum()) - sum(s.size - 1 for s in ensG)) < 1e-9
    assert abs(float(Cp_p[1].sum()) - sum(s.size - 1 for s in ensG)) < 1e-9
    assert abs(float(Ca_p.sum()) - sum(3 * s.size - 2 for s in ensG)) < 1e-9
    ok["fix37_nt_projection_exact_max_err"] = e37
    ok["fix37_phase_totals"] = {"phase2": float(Cp_p[2].sum()),
                               "phase0": float(Cp_p[0].sum()),
                               "phase1": float(Cp_p[1].sum())}

    # ============ [FIX-35 / THM-S] semigroup + rank-16 ============
    LAGS = [1, 2, 3, 5, 8]
    kerR, metaR = build_h3b_lag_kernels(ensG, LAGS, alpha=0.1)
    Pa = metaR["P_all"]
    Q1_ref = (Pa[C_N2[:, None], C_N3[:, None], C_N1[None, :]]
              * Pa[C_N3[:, None], C_N1[None, :], C_N2[None, :]]
              * Pa[C_N1[None, :], C_N2[None, :], C_N3[None, :]])
    e_r1 = float(np.abs(kerR[1]["Q"] - Q1_ref).max())
    assert e_r1 < 1e-15, e_r1
    Q2bf = np.zeros((NC, NC))
    for x in range(NC):
        uu, vv = int(C_N2[x]), int(C_N3[x])
        for b4 in range(4):
            p4 = Pa[uu, vv, b4]
            for b5 in range(4):
                p5 = p4 * Pa[vv, b4, b5]
                for b6 in range(4):
                    p6 = p5 * Pa[b4, b5, b6]
                    for y1 in range(4):
                        p7 = p6 * Pa[b5, b6, y1]
                        for y2 in range(4):
                            p8 = p7 * Pa[b6, y1, y2]
                            base = 16 * y1 + 4 * y2
                            for y3 in range(4):
                                Q2bf[x, base + y3] += p8 * Pa[y1, y2, y3]
    e_r2 = float(np.abs(kerR[2]["Q"] - Q2bf).max())
    assert e_r2 < 1e-14, e_r2
    sg_max = sg_max_p = 0.0
    for k in LAGS:
        kk = kerR[k]
        assert kk["row_sum_err"] < 1e-13
        assert kk["row_sum_err_phase"] < 1e-13
        assert kk["bridge_row_sum_err"] < 1e-13
        assert kk["bridge_row_sum_err_phase"] < 1e-13
        assert kk["semigroup_closure_err"] < 1e-13, (k, kk["semigroup_closure_err"])
        assert kk["semigroup_closure_err_phase"] < 1e-13
        assert kk["suffix_row_max_dev"] < 1e-15
        assert kk["rank_Q"] <= 16 and kk["rank_Q_phase"] <= 16
        sg_max = max(sg_max, kk["semigroup_closure_err"])
        sg_max_p = max(sg_max_p, kk["semigroup_closure_err_phase"])
    assert metaR["n_zero_eigenvalues_Q1"] >= 48
    assert metaR["rank_Q1"] <= 16
    ok["thmR_k1_factorisation_abs_err"] = e_r1
    ok["thmR_k2_bruteforce_abs_err"] = e_r2
    ok["fix35_semigroup_max_err"] = sg_max
    ok["fix35_semigroup_phase_max_err"] = sg_max_p
    ok["thmS_rank_Q1"] = metaR["rank_Q1"]
    ok["thmS_n_zero_eigenvalues_Q1"] = metaR["n_zero_eigenvalues_Q1"]
    ok["thmS_ranks_by_lag"] = {str(k): kerR[k]["rank_Q"] for k in LAGS}

    # ================= [THM-T] four-term screen + df + duplex dual =================
    eT4 = eT_dual = 0.0
    for k in (1, 2, 5, 17):
        Bkk = lag_counts(ensG, 64, k)
        f = two_sided_screen_split(Bkk, aa1)
        eT4 = max(eT4, abs(f["four_term_residual"]))
        assert f["all_nonnegative"]
        assert abs(f["Delta_x1"] - (f["D_5prime_lead"] + f["D_1_3"])) < 1e-13
        assert abs(f["I_suffix16"] - (f["I_bridge"] + f["D_3prime_wobble"])) < 1e-13
        assert f["df_full"]["sum"] == 3969 == 63 ** 2
        assert f["df_eff"]["sum"] <= 3969
        sp2 = suffix16_split(Bkk, aa1)
        assert abs(sp2["delta_x1"] - f["Delta_x1"]) < 1e-13
        assert abs(sp2["I_suffix16"] - f["I_suffix16"]) < 1e-13
        d = screen_antiparallel_dual(Bkk, aa1)
        eT_dual = max(eT_dual, d["bridge_abs_err"], d["d13_abs_err"],
                      d["lead_vs_wobble_abs_err"], d["wobble_vs_lead_abs_err"],
                      d["I_CC_abs_err"])
    assert eT4 < 1e-13 and eT_dual < 1e-13, (eT4, eT_dual)
    ok["thmT_four_term_max_abs_err"] = eT4
    ok["thmT_duplex_dual_max_abs_err"] = eT_dual
    ok["thmT_H3b_visible_df_fraction"] = 225.0 / 3969.0

    # H3b model joint: ALL THREE non-bridge layers must vanish identically
    p_left = np.bincount(ensG[0], minlength=NC).astype(float)
    p_left /= p_left.sum()
    for k in (1, 3):
        Bmod = p_left[:, None] * kerR[k]["Q"] * 1e6
        fm = two_sided_screen_split(Bmod, aa1)
        assert abs(fm["D_5prime_lead"]) < 1e-12, fm["D_5prime_lead"]
        assert abs(fm["D_3prime_wobble"]) < 1e-12, fm["D_3prime_wobble"]
        assert abs(fm["D_1_3"]) < 1e-12, fm["D_1_3"]
        Bmodp = p_left[:, None] * kerR[k]["Q_phase"] * 1e6
        fp = two_sided_screen_split(Bmodp, aa1)
        assert abs(fp["D_5prime_lead"]) < 1e-12
        assert abs(fp["D_3prime_wobble"]) < 1e-12
        assert abs(fp["D_1_3"]) < 1e-12
    ok["thmT_H3b_three_layers_vanish"] = True
    ok["thmT_model_dims"] = {
        "H3b_homogeneous": metaR["model_free_params_homogeneous"],
        "H3b_phase": metaR["model_free_params_phase"],
        "unlinked_screen": metaR["unlinked_screen_free_params"],
        "full_kernel": metaR["full_kernel_free_params"]}
    assert metaR["model_free_params_homogeneous"] == 48
    assert metaR["model_free_params_phase"] == 144
    assert metaR["unlinked_screen_free_params"] == 288
    assert 720 + 720 + 2304 + 240 + 48 == 4032 == 64 * 63

    # ============ [FIX-40] EXACT SENSE df OF THE SCREEN (all tables) ============
    dff = screen_df_bounds(None)
    assert (dff["bridge"], dff["lead"], dff["wobble"], dff["d13"]) == (225, 720, 720, 2304)
    assert dff["sum"] == 3969 and dff["identity_ok"]
    assert dff["delta_x1"] == 3024 and dff["suffix16"] == 945
    dfs1 = screen_df_bounds(aa1 != STOP)
    assert (dfs1["bridge"], dfs1["lead"], dfs1["wobble"], dfs1["d13"]) == (225, 675, 675, 2025)
    assert dfs1["sum"] == 3600 == 60 ** 2 and dfs1["identity_ok"]
    assert dfs1["sum_left_df"] == 45 and dfs1["sum_right_df"] == 45
    assert dfs1["delta_x1"] == 2700 and dfs1["suffix16"] == 900
    # sense kernel parameter accounting: lead+wobble+d13 + A(240) + R(45) = 61*60
    assert dfs1["kernel_total_df"] == 61 * 60 and dfs1["kernel_identity_ok"]
    assert dfs1["kernel_A_df"] == 240 and dfs1["kernel_R_df"] == 45
    df_by_table = {}
    for tid in sorted(SUPPORTED_TABLES):
        aat = aa_array(tid)
        sm = (aat != STOP)
        qs = int(sm.sum())
        ds = screen_df_bounds(sm)
        dd = duplex_df_table(aat, support=sm)
        assert ds["identity_ok"] and ds["sum"] == (qs - 1) ** 2, (tid, ds)
        assert ds["sum"] == dd["df_sum_check"], (tid, ds["sum"], dd["df_sum_check"])
        assert ds["kernel_total_df"] == qs * (qs - 1) and ds["kernel_identity_ok"]
        bp = box_partition(aat)
        assert bp["rho_consistent"], tid
        assert bp["n_unsplit_prefix"] + bp["n_split_prefix"] == 16
        assert bp["n_unsplit_suffix"] + bp["n_split_suffix"] == 16
        df_by_table[str(tid)] = {"q_sense": qs, "screen": ds,
                                 "duplex_sum": dd["df_sum_check"],
                                 "n_unsplit_prefix": bp["n_unsplit_prefix"]}
    assert df_by_table["2"]["screen"]["sum"] == 3481 == 59 ** 2
    assert df_by_table["4"]["screen"]["sum"] == 3721 == 61 ** 2
    assert df_by_table["6"]["screen"]["sum"] == 3844 == 62 ** 2
    ok["fix40_screen_df_by_table"] = df_by_table
    spS = suffix16_split(lag_counts(ensG, 64, 1), aa1)
    assert spS["df_delta_x1_full_torus"] == 3024
    assert spS["df_delta_x1_sense"] == 2700
    assert spS["df_suffix16_sense"] == 900
    ok["thmS_df_full_torus"] = 3024
    ok["thmS_df_sense"] = 2700
    ok["thmS_df_sense_partition"] = {k: dfs1[k] for k in ("bridge", "lead", "wobble", "d13", "sum")}

    # ================= [THM-U] 8+8 box partition and componentwise duality =======
    bp1 = box_partition(aa1)
    assert bp1["n_unsplit_prefix"] == 8 and bp1["n_unsplit_suffix"] == 8
    assert set(bp1["B4_boxes"]) == {"CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"}
    assert set(bp1["S4_boxes"]) == {"AG", "AC", "GA", "GG", "GT", "GC", "CG", "CC"}
    assert bp1["n_codons_B4"] == 32 and bp1["n_codons_Bsplit"] == 32
    eU = eUd = 0.0
    for k in (1, 2, 5, 17):
        Bkk = lag_counts(ensG, 64, k)
        f = two_sided_screen_split(Bkk, aa1)
        eU = max(eU, abs(f["box_split_lead_residual"]), abs(f["box_split_wobble_residual"]))
        d = screen_antiparallel_dual(Bkk, aa1)
        assert d["df_bridge_abs_err"] == 0 and d["df_d13_abs_err"] == 0
        assert d["df_lead_vs_wobble_abs_err"] == 0 and d["df_wobble_vs_lead_abs_err"] == 0
        eUd = max(eUd, d["U_leadS4_vs_wobB4_abs_err"], d["U_leadSsplit_vs_wobBsplit_abs_err"],
                  d["U_wobB4_vs_leadS4_abs_err"], d["U_wobBsplit_vs_leadSsplit_abs_err"])
    assert eU < 1e-13 and eUd < 1e-12, (eU, eUd)
    ok["thmU_box_split_max_abs_err"] = eU
    ok["thmU_box_duplex_dual_max_abs_err"] = eUd
    ok["thmU_B4_boxes"] = bp1["B4_boxes"]
    ok["thmU_S4_boxes"] = bp1["S4_boxes"]

    # ===== [FIX-39 / THM-T.5] EXACT held-out identities =====
    profT, _ = h3b_lag_profile(ensG[:1], ensG[1:], LAGS, alpha=0.1, aa=aa1)
    assert [r["lag"] for r in profT] == LAGS
    e_ho = e_t4 = eS = e_t5 = e_t5p = e_t5d = e_hop = 0.0
    for r in profT:
        assert r["D_gen_train_test"] >= -1e-12
        assert r["D_suffix_test"] >= -1e-12
        assert r["Delta_x1_test"] >= -1e-12
        assert r["rank_Q"] <= 16
        e_ho = max(e_ho, r["thmT_heldout_identity_abs_err"])
        e_t4 = max(e_t4, r["thmT_suffix_split_abs_err"])
        eS = max(eS, r["thmS_identity_abs_err"])
        e_t5 = max(e_t5, r["thmT5_five_term_abs_err"])
        e_t5p = max(e_t5p, r["thmT5_phase_five_term_abs_err"])
        e_t5d = max(e_t5d, r["thmT5_phase_diff_abs_err"])
        e_hop = max(e_hop, r["thmT_heldout_identity_phase_abs_err"])
    assert e_ho < 1e-13, e_ho
    assert e_t4 < 1e-13, e_t4
    assert eS < 1e-13, eS
    assert e_t5 < 1e-13 and e_t5p < 1e-13 and e_t5d < 1e-13 and e_hop < 1e-13
    ok["fix39_heldout_identity_max_abs_err"] = e_ho
    ok["thmT_suffix_split_max_abs_err"] = e_t4
    ok["thmS_insample_identity_max_abs_err"] = eS
    ok["thmT5_five_term_max_abs_err"] = e_t5
    ok["thmT5_phase_five_term_max_abs_err"] = e_t5p
    ok["thmT5_phase_diff_max_abs_err"] = e_t5d
    ok["thmT5_heldout_phase_identity_max_abs_err"] = e_hop
    ok["thmT_delta_x1_test_by_lag"] = {str(r["lag"]): r["Delta_x1_test"] for r in profT}
    ok["thmT_D_gen_by_lag"] = {str(r["lag"]): r["D_gen_train_test"] for r in profT}

    # [THM-T.5.2] R_wobble_phase == P_train(Y3|Y1,Y2): on train==test the phase
    # wobble KL vanishes up to Dirichlet smoothing, while the homogeneous one does not.
    profS, metaS = h3b_lag_profile(ensG, ensG, [1, 3], alpha=1e-9, aa=aa1)
    Rw2_chk = metaS["R_wobble_phase"]
    marg_chk = np.zeros(NC)
    for s_ in ensG:
        marg_chk += np.bincount(s_, minlength=NC)
    M3 = marg_chk.reshape(4, 4, 4)
    P3 = (M3 + 1e-9) / (M3.sum(axis=2, keepdims=True) + 4e-9)
    e_rw = max(float(np.abs(Rw2_chk[o] - P3[o // 4, o % 4, :]).max()) for o in range(16))
    assert e_rw < 1e-15, e_rw
    for r in profS:
        assert r["D_wobble_kl_phase_aligned"] < 1e-9, r["D_wobble_kl_phase_aligned"]
        assert abs(r["D_wobble_phase_5prime_edge"]
                   - (r["D_wobble_kl_phase"] - r["D_wobble_kl_phase_aligned"])) < 1e-15
        assert r["D_wobble_kl"] >= r["D_wobble_kl_phase_aligned"] - 1e-12
    ok["thmT5_R_wobble_phase_is_train_conditional_max_err"] = e_rw
    ok["thmT5_D_wobble_phase_on_train_eq_test"] = {
        str(r["lag"]): r["D_wobble_kl_phase"] for r in profS}
    ok["thmT5_D_wobble_homog_on_train_eq_test"] = {
        str(r["lag"]): r["D_wobble_kl"] for r in profS}

    # ---- Theorem D vs exhaustive enumeration ----
    raw4 = (1 + ((np.arange(q4)[:, None] * 7 + np.arange(q4)[None, :] * 13) % 11)).astype(float)
    lP4 = np.log(raw4 / raw4.sum(axis=1, keepdims=True))
    tot, cnt = 0.0, 0
    for perms in itertools.product(*[list(itertools.permutations(g)) for g in groups]):
        out = s6.copy()
        for g, pm in zip(groups, perms):
            out[g] = s6[list(pm)]
        tot += float(np.sum(lP4[out[:-1], out[1:]])); cnt += 1
    a6 = layer_shuffle_analytic_gene(s6, q4, lay4, 3, lP4, np.diag(lP4))
    assert abs(tot / cnt - a6) < 1e-12
    i0, i1, i2 = CODON_IDX["ATG"], CODON_IDX["GCT"], CODON_IDX["TTA"]
    s7 = np.array([i0, MIRROR[i0], i1, MIRROR[i1], i2], dtype=np.int64)
    gr7 = [np.flatnonzero(MIRROR32[s7] == a) for a in np.unique(MIRROR32[s7])]
    tot, cnt = 0.0, 0
    for perms in itertools.product(*[list(itertools.permutations(g)) for g in gr7]):
        out = s7.copy()
        for g, pm in zip(gr7, perms):
            out[g] = s7[list(pm)]
        tot += float(np.sum(logP[out[:-1], out[1:]])); cnt += 1
    a7 = layer_shuffle_analytic_gene(s7, NC, MIRROR32, 32, logP, np.diag(logP))
    assert abs(tot / cnt - a7) < 1e-12
    Btrc = ((ix[:, None] * 13 + ix[None, :] * 7) % 41).astype(np.float64)
    Qs = build_h2b_sampler_from_counts(Btrc, fib1, alpha=0.1)
    Qt, Rt = h5_pair_tables(Qs, logP, fib1)
    s5 = np.array([CODON_IDX["ATG"], CODON_IDX["TGG"], CODON_IDX["TTT"], CODON_IDX["TAT"]],
                  dtype=np.int64)
    tot = 0.0
    for traj in itertools.product(*[fib1[int(aa1[c])] for c in s5]):
        if traj[0] != s5[0]:
            continue
        pr = 1.0
        for t in range(1, len(traj)):
            pr *= Qs[traj[t - 1], int(aa1[s5[t]]), traj[t]]
        tot += pr * sum(logP[traj[t - 1], traj[t]] for t in range(1, len(traj)))
    a5, _ = h5_analytic_gene(s5, aa1, fib1, Qt, Rt)
    assert abs(tot - a5) < 1e-12
    okgene = np.array([fib1[int(a)][0] for a in aa1[sL[:5000]]], dtype=np.int64)
    _, devL = h5_analytic_gene(okgene, aa1, fib1, Qt, Rt)
    assert devL < MASS_TOL
    ok["thmD_H5_mass_dev_5000_steps"] = devL

    with tempfile.TemporaryDirectory() as td:
        fp = Path(td) / "synthetic_cds_from_genomic.fna.gz"
        with gzip.open(fp, "wt", encoding="utf-8") as fh:
            for gi in range(30):
                pfx = "ZNF" if gi >= 15 else "AAA"
                fh.write(f">seq{gi} [gene={pfx}{gi:03d}] [transl_table=1]\n")
                fh.write("ATG" + "GCT" * 20 + "TAA\n")
        r29, t29, s29 = parse_cds(fp, max_genes=10, trim=True, default_table=1, min_codons=8, gene_order="hash")
        assert any(g.startswith("ZNF") for g, _ in r29), "hash cap must not amputate Z-prefixed genes"
        assert len(r29) == 10, s29
        assert isinstance(s29["n_atg_start_kept"], int), s29
        assert s29["n_atg_start_kept"] == 10 and s29["n_stop_end_kept"] == 10, s29
        assert s29["n_atg_start_raw"] == 30 and s29["n_stop_end_raw"] == 30, s29
        ok["fix29_counters_after_slice"] = {
            "kept": s29["kept"], "atg_kept": s29["n_atg_start_kept"],
            "stop_kept": s29["n_stop_end_kept"], "atg_raw": s29["n_atg_start_raw"]}

    Blap = bigram_counts([sL[:2000]], 64)
    ref = np.log((Blap + 0.1) / (Blap.sum(axis=1, keepdims=True) + 64 * 0.1))
    assert np.max(np.abs(laplace_logp(Blap, 0.1) - ref)) < 1e-15
    sa = [f"SAOUHSC_{i:05d}" for i in range(1, 301)]
    assert len({family_root(g) for g in sa}) == 300
    assert len({family_root(g) for g in ["TP53-RA", "TP53-RB", "TP53"]}) == 1
    mask, info = assign_split(sa, 80)
    assert info["mode"] == "family_root" and min(info["n_train"], info["n_test"]) > 0
    sm = np.array([(2 * i) % 4 for i in range(120)], dtype=np.int64)
    mm = build_models([sm], 4, 1, kappa=8.0, min_ctx_count=1)
    lp2 = full_first_order_logp(mm, 4)
    assert lp2.shape == (4, 4) and np.allclose(lp2[1], lp2[3]) and not np.allclose(lp2[0], lp2[1])
    for bad in (None, float("nan"), float("inf"), "x", {}):
        assert fnum(bad) == "n/a"
    ok["all_kernel_and_split_guards"] = True
    ok.update(self_test_v22(ensG, aa1))
    ok.update(self_test_v23(ensG, aa1))
    ok.update(self_test_v24(ensG, aa1))
    ok.update(self_test_v30(ensG, aa1))

    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


# =====================================================================
# 18. Species driver
# =====================================================================
def analyze_species(name, path, args):
    t0 = time.time()
    print(f"\n[{name}] {path.name}")
    macc = re.search(r"(GCF_\d+)", path.name)
    default_tb = ACC2TABLE.get(macc.group(1), NAME2TABLE.get(name, 1)) if macc else NAME2TABLE.get(name, 1)
    recs, table, stats = parse_cds(path, max_genes=args.max_genes, trim=not args.no_trim,
                                   default_table=default_tb, min_codons=args.min_codons,
                                   gene_order=getattr(args, "gene_order", "hash"))
    if len(recs) < args.min_genes:
        print(f"  [skip] too few valid CDS ({len(recs)})")
        return None
    aa = aa_array(table)
    fibers = fibers_of(aa)
    layD, FD, szD, d8 = duplex_layers(aa)
    cpg = cpg_quartet_partition(aa)
    sup = np.zeros(NC, dtype=bool)
    for _, cd in recs:
        sup[cd] = True
    dfF = duplex_df_table(aa)
    dfS = duplex_df_table(aa, support=sup)
    dfS_theory = duplex_df_table(aa, support=(aa != STOP))
    mask, sinfo = assign_split([g for g, _ in recs], args.train_pct)
    train_cd = [cd for (_, cd), m in zip(recs, mask) if m]
    test_cd = [cd for (_, cd), m in zip(recs, mask) if not m]
    if len(train_cd) < 20 or len(test_cd) < 20:
        print(f"  [skip] degenerate split ({sinfo})")
        return None
    gc3_test = np.array([float(np.mean(IS_GC3[cd])) for cd in test_cd])
    print(f"  table={table} train={len(train_cd)} test={len(test_cd)} "
          f"codons={stats['n_codons_kept']} median_gene={stats['median_gene_codons']}")
    print(f"  [split] mode={sinfo['mode']} roots={sinfo['n_family_roots']} | "
          f"ATG_kept={stats['n_atg_start_kept']} stop_kept={stats['n_stop_end_kept']} "
          f"(raw {stats['n_atg_start_raw']}/{stats['n_stop_end_raw']}) | "
          f"rej nonACGT={stats['rej_nonacgt']} frame={stats['rej_frame']} "
          f"intStop={stats['rej_internal_stop']}")
    bp = box_partition(aa)
    dfs = screen_df_bounds(aa != STOP)
    dfs_emp = screen_df_bounds(sup)
    dff = screen_df_bounds(None)
    print(f"  [FIX-40] screen df full={dff['sum']} sense={dfs['sum']} "
          f"(bridge {dfs['bridge']}, lead {dfs['lead']}, "
          f"wobble {dfs['wobble']}, d13 {dfs['d13']}) | G.2 sense={dfS['df_sum_check']}")
    clade_nm = NAME2CLADE.get(name, ACC2CLADE.get(macc.group(1), "unknown") if macc else "unknown")
    res = {"species": name, "clade": clade_nm, "file": path.name, "transl_table": int(table),
           "parse_stats": stats, "split_info": sinfo,
           "gc3_test_mean": float(gc3_test.mean()),
           "df_full": dfF, "df_sense_support": dfS, "df_sense_theory": dfS_theory,
           "screen_df_full_torus": dff, "screen_df_sense": dfs,
           "screen_df_sense_empirical": dfs_emp,
           "df_support_convention_ok": bool(dfs["sum"] == dfS_theory["df_sum_check"]
                                            and dfs_emp["sum"] == dfS["df_sum_check"]),
           "box_partition": {k: bp[k] for k in ("B4_boxes", "Bsplit_boxes",
                                                "S4_boxes", "Ssplit_boxes",
                                                "n_unsplit_prefix", "n_unsplit_suffix",
                                                "n_codons_B4", "n_codons_Bsplit",
                                                "rho_consistent")},
           "thmJ_partition": {k: v for k, v in cpg.items()
                              if k not in ("cpg_mask", "leu_mask")}}
    logP64 = build_kernel(train_cd, 64, args.kernel_mode, args.kappa, args.laplace_alpha)
    row_ent = float(np.mean(-np.sum(np.exp(logP64) * logP64, axis=1)))
    if row_ent > math.log(64) - 1e-6:
        raise RuntimeError("degenerate uniform first-order kernel")
    res["kernel_diagnostics"] = {"mode": args.kernel_mode,
                                 "mean_row_entropy_nats": row_ent}
    lengths = sorted({int(x) for x in args.h4_lengths.split(",") if x.strip() and int(x) >= 2})
    psi = psi_constant(test_cd, logP64)
    pv = psi["Psi_LR"] if psi else None
    h4 = []
    for L in lengths:
        r = analytic_h4(test_cd, logP64, L, collect_terms=True)
        if r is None:
            continue
        r["L_times_delta"] = L * r["delta_H4_analytic"]
        r["L_times_gamma"] = L * r["gamma_L"]
        okp = _finite(pv) and pv != 0
        r["ratio_to_Psi_naive"] = (r["L_times_delta"] / pv) if okp else None
        r["R_junc"] = (r["delta_H4_analytic"] / (r["gamma_L"] * pv)) \
            if (okp and r["gamma_L"] > 0) else None
        tE = theorem_E(test_cd, logP64, L, psi, r)
        if tE:
            r["theorem_E"] = tE
        c = cohort_decomposition(test_cd, logP64, L)
        if c:
            r["cohort"] = c
        if args.mc_check > 0 and L in (3, 10, 30):
            mc = mc_h4_check(test_cd, logP64, L, args.mc_check, args.seed + L)
            if mc:
                r["mc_abs_err"] = abs(mc["delta_H4_montecarlo"] - r["delta_H4_analytic"])
        for kk in ("_Dsum", "_Nsum", "_NcompSum", "_NisoSum"):
            r.pop(kk, None)
        h4.append(r)
    rj = [r["R_junc"] for r in h4 if _finite(r["R_junc"])]
    rn = [r["ratio_to_Psi_naive"] for r in h4 if _finite(r["ratio_to_Psi_naive"])]
    if not _finite(pv) or pv <= 0:
        verdict = "PSI_NONPOSITIVE_KERNEL_MISMATCH"
    elif not rj:
        verdict = "DEGENERATE_NO_RATIOS"
    elif (max(rj) - min(rj)) < args.dilution_tol:
        verdict = "PURE_JUNCTION_DILUTION"
    else:
        verdict = "POSITIONAL_STRUCTURE_DETECTED"
    res["theorem_A"] = {"Psi": psi, "by_length": h4,
                        "R_junc_min": float(min(rj)) if rj else None,
                        "R_junc_max": float(max(rj)) if rj else None,
                        "R_junc_spread": float(max(rj) - min(rj)) if rj else None,
                        "naive_spread": float(max(rn) - min(rn)) if rn else None,
                        "verdict": verdict}
    kmax = int(args.max_lag)
    prof = lagged_mi_profile(test_cd, 64, kmax)
    pred, spec = markov_lag_prediction(np.exp(logP64), kmax)
    for o, p in zip(prof, pred):
        o.update({k: v for k, v in p.items() if k != "lag"})
        o["excess_over_markov"] = (o["I_MM"] - p["I_markov"]) if _finite(o.get("I_MM")) else None
    res["theorem_C"] = {"spectral": spec, "lag_profile": prof}
    print(f"  [F.2] |l2|={fnum(spec['lambda2_modulus'],4)} xi={fnum(spec['xi_codons'],2)} "
          f"kappa_mean={fnum(spec.get('kappa_mean'),4)}")
    if args.het_windows:
        wl = [int(x) for x in args.het_window_lags.split(",") if x.strip()]
        wr = {}
        for W in [int(x) for x in args.het_windows.split(",") if x.strip()]:
            wc = window_cache(test_cd, 64, W, args.het_max_windows)
            if wc is None:
                continue
            rows = []
            for k in wl:
                if k >= W:
                    continue
                o, e = window_pairs_cached(wc, 64, k)
                Io, *_ = _mi_from_counts(o, mm=True)
                Ie, *_ = _mi_from_counts(e, mm=False)
                rows.append({"lag": k, "N": int(o.sum()), "I_obs_MM": Io, "I_het_window": Ie,
                             "het_fraction": (Ie / Io) if Io > 0 else None,
                             "KL_obs_vs_het_window": kl_of_joint(o, e),
                             "theorem_L": theorem_L_diagnostics(o, e, f"het_W{W}")})
            wr[str(W)] = {"meta": {k: v for k, v in wc.items()
                                   if k not in ("S_full", "C_full", "tails", "C_tail", "L_tail")},
                          "rows": rows}
        res["theorem_C3_windowed"] = wr
    if not args.no_duplex:
        gl = [int(x) for x in args.duplex_lags.split(",") if x.strip()]
        grows = [g for g in (lag_duplex_decomposition(test_cd, 64, k, aa, d8, df=dfS, cpg=cpg)
                             for k in gl) if g is not None]
        res["theorem_G"] = {"lags": grows}
        for g in grows:
            if g["lag"] in (1, 5, 10, 30):
                print(f"  [G/J/P] k={g['lag']:2d} raw_anti%={fnum(100*(g['antisense_fraction_raw'] or 0),1)}"
                      f" mirror%={fnum(100*(g.get('mirror_capture_fraction') or 0),1)}"
                      f" enrich_raw={fnum(g.get('enrich_per_df_raw'),3)}"
                      f" CpG%={fnum(100*(g.get('cpg_share_of_S_dup') or 0),1)}")
    if not args.no_tower:
        tT = time.time()
        parts = {"het": (np.zeros(NC, dtype=np.int64), 1),
                 "aa": (aa, NAA), "dup": (layD, FD), "mirror": (MIRROR32, 32)}
        want = [x.strip() for x in args.tower_partitions.split(",") if x.strip() in parts]
        tl = [int(x) for x in args.tower_lags.split(",") if x.strip()]
        tower = {}
        for pname in want:
            lay, F = parts[pname]
            rows = []
            for k in tl:
                Ek, Bk, meta = layer_tower_expected(test_cd, 64, k, lay, F,
                                                    max_genes=args.tower_max_genes)
                if Bk.sum() <= 0:
                    continue
                ep, et = tower_projection_check(Ek, Bk, lay, F)
                kw = kl_four_way(Bk, Ek, aa)
                kn = kl_symmetric_duplex_tower(Bk, Ek, aa)
                if kw is None or kn is None:
                    continue
                Nb, Ne = float(Bk.sum()), float(Ek.sum())
                row = {"lag": k, "N": int(Nb), "n_genes": meta["n_genes"],
                       "projection_block_max_err": ep,
                       "projection_rel_err": ep / max(1.0, float(Bk.max())),
                       "no_clipping_required": meta.get("no_clipping_required"),
                       "lemma_H3": lemma_H3_check(Ek, Bk),
                       "klein_null": klein_quadrants_of_matrix(Ek / max(Ne, TINY)),
                       "klein_residual": klein_quadrants_of_matrix(
                           Bk / max(Nb, TINY) - Ek / max(Ne, TINY)),
                       "one_sided": kw, "symmetric": kn}
                if pname == "dup":
                    row["theorem_K"] = theorem_K_localisation(Ek, Bk, lay)
                rows.append(row)
            tower[pname] = rows
            if pname in ("aa", "dup"):
                for r in rows:
                    if r["lag"] in (1, 5, 10, 30):
                        sy = r["symmetric"]
                        print(f"  [N/{pname}] k={r['lag']:2d} D_AA_sym={sy['D_AA_sym']:.1e} "
                              f"D_anti_sym={sy['D_anti_sym']:.1e} "
                              f"D_D8_sym={fnum(sy['D_D8_sym'],6,True)} "
                              f"proj_rel={r['projection_rel_err']:.1e}")
        res["theorem_tower"] = {"partitions": tower, "seconds": time.time() - tT}
    if not args.no_phase:
        pl = [int(x) for x in args.phase_lags.split(",") if x.strip()]
        prow = theorem_O_multi_lag(test_cd, pl, aa)      # [FIX-38] one build
        res["theorem_O"] = {"lags": prow}
        for o in prow:
            if o["lag"] in (1, 10, 30):
                fr = o["frames"]
                print(f"  [O] k={o['lag']:2d} Omega_phase={fnum(o['Omega_phase_JS_nats'],6)} "
                      f"Omega_cond={fnum(o['Omega_cond_JS_nats'],6)} "
                      f"bias_cond={o['bias_cond']:.2e} "
                      f"pairs_equal={o['frame_pair_counts_equal']} "
                      f"S_syn(0/1/2)={fnum(fr['frame_plus0']['S_syn'],4)}/"
                      f"{fnum(fr['frame_plus1']['S_syn'],4)}/"
                      f"{fnum(fr['frame_plus2']['S_syn'],4)}")
    # ---------- [THM-R/S/T/T.5/U] order-2 nt propagator, screens, exact splits ----------
    if not args.no_h3b:
        tR = time.time()
        kl = [int(x) for x in args.h3b_lags.split(",") if x.strip()]
        rowsR, metaR = h3b_lag_profile(train_cd, test_cd, kl,
                                       alpha=args.laplace_alpha, aa=aa)
        res["theorem_R"] = {"lags": rowsR,
                            "rank_Q1": metaR["rank_Q1"],
                            "n_zero_eigenvalues_Q1": metaR["n_zero_eigenvalues_Q1"],
                            "model_dims": {
                                "H3b_homogeneous": metaR["model_free_params_homogeneous"],
                                "H3b_phase": metaR["model_free_params_phase"],
                                "unlinked_screen": metaR["unlinked_screen_free_params"],
                                "full_kernel": metaR["full_kernel_free_params"]},
                            "seconds": time.time() - tR}
        for r in rowsR:
            print(f"  [R/S] k={r['lag']:2d} d3b={fnum(r['d3b'],6,True)} "
                  f"d3b_ph={fnum(r['d3b_phase'],6,True)} "
                  f"dx1={fnum(r['Delta_x1_test'],6,True)} "
                  f"Dgen={fnum(r['D_gen_train_test'],6,True)} "
                  f"rank={r['rank_Q']:2d} ho_id={r['thmT_heldout_identity_abs_err']:.1e}")
            print(f"  [T]   k={r['lag']:2d} I_br={fnum(r['I_bridge'],6,True)} "
                  f"lead={fnum(r['D_5prime_lead'],6,True)} "
                  f"wob={fnum(r['D_3prime_wobble'],6,True)} "
                  f"13={fnum(r['D_1_3'],6,True)} "
                  f"T4={r['thmT_four_term_abs_err']:.1e} "
                  f"Tkl={r['thmT_suffix_split_abs_err']:.1e}")
            print(f"  [T.5] k={r['lag']:2d} 5term={r['thmT5_five_term_abs_err']:.1e} "
                  f"ph_diff={r['thmT5_phase_diff_abs_err']:.1e} "
                  f"Dbr={fnum(r['D_bridge_kl'],6,True)}/{fnum(r['D_bridge_kl_phase'],6,True)} "
                  f"Dwob={fnum(r['D_wobble_kl'],6,True)}/{fnum(r['D_wobble_kl_phase'],6,True)}")
            print(f"  [U]   k={r['lag']:2d} wobB4={fnum(r['D_3prime_wobble_B4'],6,True)} "
                  f"wobSplit={fnum(r['D_3prime_wobble_Bsplit'],6,True)} "
                  f"leadS4={fnum(r['D_5prime_lead_S4'],6,True)} "
                  f"leadSplit={fnum(r['D_5prime_lead_Ssplit'],6,True)} "
                  f"df_eff={r['df_eff']['sum']}/{dfs['sum']}")
    res["theorem_X"] = thm_X_screen_v4()
    res["theorem_AA_left"] = {"sense": left_screen_structure(aa, True),
                              "full": left_screen_structure(aa, False)}
    _ls = res["theorem_AA_left"]["sense"]
    print(f"  [AA] left |S,A|={_ls['n_SA']} comp={_ls['n_components']} "
          f"cycles={_ls['cycle_rank']} forest={_ls['is_forest']} "
          f"tower={_ls['tower']}")
    res["duplex_multi"] = duplex_multi_partition(aa)
    Bt1 = lag_counts(test_cd, NC, 1)
    res["theorem_Y"] = {"sense": thm_Y_diamond(Bt1, aa, True),
                        "full": thm_Y_diamond(Bt1, aa, False)}
    res["theorem_Z"] = {"sense": thm_Z_tower(Bt1, aa, True),
                        "full": thm_Z_tower(Bt1, aa, False)}
    ys, zs = res["theorem_Y"]["sense"], res["theorem_Z"]["sense"]
    print(f"  [Y] O={ys['n_O']} A={ys['n_A']} OA={ys['n_OA']} comp={ys['n_components']} "
          f"forest={ys['is_forest']} balance={ys['balance_abs_err']:.1e}")
    print(f"  [Z] tower={ys['n_O']}->{zs['tower'][1]}->{zs['tower'][2]}->{zs['tower'][3]} "
          f"Sdup={fnum(zs['S_dup'],6,True)} chain={zs['chain_abs_err']:.1e} "
          f"diamond={zs['diamond_abs_err']:.1e}")
    print(f"  [X] comm_DS_R^2={fnum(res['theorem_X']['comm_DS_R_sq'],4)} "
          f"C2 dims={res['theorem_X']['C2_quadrant_dims']}")
    res["theorem_P"] = theorem_P_transversality(aa)
    res["theorem_Q"] = klein_commutator_audit(aa)
    print(f"  [P] (alpha,M) classes={res['theorem_P']['n_joint_AA_mirror_classes']}/64 "
          f"relic_diag={res['theorem_P']['relic_diagonal_total']}")
    print(f"  [Q] ||[Pi_AA,R]||_F={fnum(res['theorem_Q']['AA']['commutator_R_fro'],4)} "
          f"HS_superop={fnum(res['theorem_Q']['AA']['superoperator_commutator_HS'],4)} "
          f"dup={res['theorem_Q']['dup']['commutator_R_fro']:.1e}")
    if not args.no_h567:
        Qs = build_h2b_sampler(train_cd, aa, fibers, alpha=args.h5_alpha)
        res["theorem_D"] = h567_audit(test_cd, logP64, aa, fibers, Qs,
                                      h5_max_genes=args.h5_max_genes)
        for kk, vv in res["theorem_D"].items():
            print(f"  [D] {kk:18s} n={vv['n_genes_common']:5d} "
                  f"delta={fnum(vv['delta_per_transition_gene_mean'],6,True)} "
                  f"t={fnum(vv['t_stat'],2,True)}")
    res["theorem_B"] = {}
    for alb in [a.strip() for a in args.alphabets.split(",") if a.strip()]:
        if alb not in ALPHABETS:
            continue
        q, _, dm = ALPHABETS[alb]
        mo = min(int(args.max_order), dm) if args.max_order > 0 else dm
        tr = [encode(c, alb, aa) for c in train_cd]
        te = [encode(c, alb, aa) for c in test_cd]
        lad = memory_ladder(tr, te, q, mo, args.kappa, args.min_ctx_count,
                            args.max_ctx, gene_gc3=gc3_test)
        if lad is None:
            continue
        lad["alphabet"] = alb
        res["theorem_B"][alb] = lad
        print(f"  [B] {alb:9s} p*={lad['p_star_argmax']:2d} "
              f"range={lad['max_significant_range']:2d} "
              f"gain={lad['memory_beyond_first_order_nats']:+.5f} "
              f"cov={lad['coverage_at_p_star']:.3f}")
    # ---------- [THM-AD] Isochore-free screen decomposition against E_het^(k) ----------
    if not args.no_h3b:
        kl_ad = [int(x) for x in args.h3b_lags.split(",") if x.strip()]
        ad_rows = thm_AD_isochore_screen_audit(test_cd, kl_ad, aa)
        res["theorem_AD"] = {"lags": ad_rows}
        for r_ad in ad_rows:
            if r_ad["lag"] in (1, 30):
                sk = r_ad["screen_kl"]
                print(f"  [AD]  k={r_ad['lag']:2d} wob_MM={fnum(r_ad['D_3prime_wobble_MM'],5)} "
                      f"wob_het={fnum(r_ad['D_3prime_wobble_het'],5)} "
                      f"wob_intra={fnum(r_ad['D_3prime_wobble_intra_MM'],5,True)} | "
                      f"dx1_intra={fnum(r_ad['Delta_x1_intra_MM'],5,True)} "
                      f"K_nonbr={fnum(sk['K_nonbridge'],5)} K_l|SY={fnum(sk['K_lead_given_SY'],5)}")
    # ---------- [THM-AB] Periodicity & C2H2 Zinc-Finger stratification ----------
    if not args.no_znf:
        ab_res = thm_AB_periodicity_znf_audit(test_cd, aa, lag_profile=prof, kmax=kmax)
        res["theorem_AB"] = ab_res
        print(f"  [AB]  ZNF_C2H2={ab_res['n_znf_c2h2_genes']}/{ab_res['n_genes_total']} "
              f"dip14(all)={fnum(ab_res['all_mp']['dip14_log10'],3,True)} "
              f"peakP={fnum(ab_res['all_mp']['dft_peak_period'],1)} "
              f"SNR14={fnum(ab_res['all_mp']['dft_snr_at_14'],2)} -> {ab_res['verdict']}")
    # ---------- Compact matrices for [THM-AC] cross-species phylogeny (kept in memory only) ----------
    B1_tr = lag_counts(train_cd, NC, 1)
    P1_te = Bt1 / max(float(Bt1.sum()), TINY)
    # Build exact H2b conditional log-kernel for zero-RNG cross-species d2b_transfer
    bg_a = B1_tr + args.laplace_alpha
    P_aa_tr = np.zeros((NAA, NAA), dtype=np.float64)
    for a_ in range(NAA):
        fa_ = fibers[a_]
        if fa_.size:
            for b_ in range(NAA):
                fb_ = fibers[b_]
                if fb_.size:
                    P_aa_tr[a_, b_] = bg_a[np.ix_(fa_, fb_)].sum()
    P_aa_tr /= np.maximum(P_aa_tr.sum(axis=1, keepdims=True), TINY)
    Q2b_tr = np.zeros((NC, NC), dtype=np.float64)
    for b_ in range(NAA):
        fb_ = fibers[b_]
        if fb_.size:
            sub = bg_a[:, fb_]
            sub = sub / np.maximum(sub.sum(axis=1, keepdims=True), TINY)
            Q2b_tr[:, fb_] = P_aa_tr[aa, b_][:, None] * sub
    res["_matrix_cache"] = {
        "P1_test": P1_te,
        "logP_train": logP64,
        "logQ2b_train": np.log(np.maximum(Q2b_tr, TINY)),
    }
    res["seconds_total"] = time.time() - t0
    A = res["theorem_A"]
    print(f"  [A] R_junc=[{fnum(A['R_junc_min'])},{fnum(A['R_junc_max'])}] "
          f"spread={fnum(A['R_junc_spread'])} -> {A['verdict']}")
    return res


# =====================================================================
# 19. Reporting
# =====================================================================
def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    return o


def atomic_write(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, p)


def write_md(payload, path):
    R = payload["results"]
    L = [f"# Codon Memory Spectrum v{VERSION} — real NCBI RefSeq audit", "",
         "ZERO RNG in every proof path. Theorems A–AD.", "",
         f"kernel `{payload['config']['kernel_mode']}`; taxa **{len(R)}**.", "",
         "## 1. Theorem M — Klein V_4 quadrants (1056/1024/1024/992)", "",
         "| Species | k | ‖W_++‖² | ‖W_+-‖² static | ‖W_-+‖² | ‖W_--‖² | dyn share of mirror |",
         "|---|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        d = {x["lag"]: x for x in r["theorem_C"]["lag_profile"]}
        for k in (1, 10, 30, 60):
            x = d.get(k)
            if not x or not x.get("theorem_M"):
                continue
            m = x["theorem_M"]
            L.append(f"| {r['species']} | {k} | {fnum(m['norm2_W_pp'],8)} | "
                     f"{fnum(m['norm2_W_pm_static_strand'],8)} | "
                     f"{fnum(m['norm2_W_mp_sym_time_flux'],8)} | "
                     f"{fnum(m['norm2_W_mm_duplex_time_flux'],8)} | "
                     f"{fnum(m['dynamical_share_of_mirror_flux_F2'],4)} |")
    L += ["", "## 1b. Theorem Q / Q.2 — both commutator invariants", "",
          "| Species | C | S₂ | S₁ | ‖[D,R]‖_F² | ‖[Π,R]‖_HS² | dup comm | mirror comm |",
          "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        q = r.get("theorem_Q") or {}
        if not q:
            continue
        a = q["AA"]
        L.append(f"| {r['species']} | {a['n_layers']} | {fnum(a['S2_trace_DRDR'],6)} | "
                 f"{fnum(a['S1_linear_weight_sum'],6)} | "
                 f"{fnum(a['commutator_R_fro']**2,6)} | "
                 f"{fnum(a['superoperator_commutator_HS_sq'],6)} | "
                 f"{q['dup']['commutator_R_fro']:.1e} | {q['mirror']['commutator_R_fro']:.1e} |")
    L += ["", "Table 1 AA: C=21, S₂=569/72, S₁=509/72 ⟹ 943/36 ≈ 26.1944 and",
          "2C²−2S₂² = 1962383/2592 ≈ 757.0922.  All computed, never hardcoded.", "",
          "## 1c. Theorem P — (α,M) transversality and mirror capture", "",
          "| Species | classes/64 | relic diag | k | S_syn | S_mirror | capture % | residual |",
          "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        p = r.get("theorem_P") or {}
        for g in (r.get("theorem_G") or {}).get("lags", []):
            if g["lag"] not in (1, 5, 10, 30):
                continue
            L.append(f"| {r['species']} | {p.get('n_joint_AA_mirror_classes','n/a')} | "
                     f"{p.get('relic_diagonal_total','n/a')} | {g['lag']} | "
                     f"{fnum(g['S_syn'],6,True)} | {fnum(g.get('S_mirror_orbit'),6,True)} | "
                     f"{fnum(100*(g.get('mirror_capture_fraction') or 0),2)} | "
                     f"{g.get('mirror_capture_residual',0):.1e} |")
    L += ["", "## 2. Theorem N — symmetric duplex tower", ""]
    for r in R:
        tw = (r.get("theorem_tower") or {}).get("partitions") or {}
        if not tw:
            continue
        L += [f"#### {r['species']}", "",
              "| Λ | k | KL | D_AA_sym | D_anti_sym | D_D8_sym | proj rel | H.3 min E |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
        for pname, rows in tw.items():
            for x in rows:
                sy = x["symmetric"]
                h3 = x.get("lemma_H3") or {}
                L.append(f"| `{pname}` | {x['lag']} | {fnum(sy['KL_total'],6,True)} | "
                         f"{sy['D_AA_sym']:.1e} | {sy['D_anti_sym']:.1e} | "
                         f"{fnum(sy['D_D8_sym'],6,True)} | "
                         f"{x['projection_rel_err']:.1e} | "
                         f"{h3.get('min_E_on_support',0):.2e} |")
        L.append("")
    L += ["## 3. Theorem T — two-sided 16×16 screen: four-term MI split", "",
          "I(C_t;C_{t+k}) = I_bridge + D_5'lead + D_3'wobble + D_1↔3.",
          "Full-torus df: **225 + 720 + 720 + 2304 = 3969 = 63²**.",
          "Sense df [FIX-40]: **225 + 675 + 675 + 2025 = 3600 = 60²** (Table 1),",
          "the SAME 3600 as G.2's 361+779+2040+420. H3b sees ONLY the 225-df bridge",
          "(5.669% full / 6.25% sense); all three other layers vanish under it.", "",
          "| Species | k | I_CC | I_bridge | D_5'lead | D_3'wobble | D_1↔3 | 4-term err | df_eff sum | df_sense |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            L.append(f"| {r['species']} | {x['lag']} | {fnum(x['I_CC_test'],6,True)} | "
                     f"{fnum(x['I_bridge'],6,True)} | {fnum(x['D_5prime_lead'],6,True)} | "
                     f"{fnum(x['D_3prime_wobble'],6,True)} | {fnum(x['D_1_3'],6,True)} | "
                     f"{x['thmT_four_term_abs_err']:.1e} | {x['df_eff']['sum']} | "
                     f"{(x.get('df_sense') or {}).get('sum','n/a')} |")
    L += ["", "### 3b. Theorem T.3 / U — duplex duality (lead ↔ wobble under ρ)", "",
          "| Species | k | lead(sense) | wobble(anti) | swap err | bridge err | D_1↔3 err | df lead↔wob err |",
          "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            d = x.get("antiparallel_dual") or {}
            if not d:
                continue
            L.append(f"| {r['species']} | {x['lag']} | "
                     f"{fnum(d['sense']['D_5prime_lead'],6,True)} | "
                     f"{fnum(d['antisense']['D_3prime_wobble'],6,True)} | "
                     f"{d['lead_vs_wobble_abs_err']:.1e} | "
                     f"{d['bridge_abs_err']:.1e} | {d['d13_abs_err']:.1e} | "
                     f"{d.get('df_lead_vs_wobble_abs_err','n/a')} |")
    L += ["", "Leading-base (amino-acid-core) memory on the sense strand IS wobble",
          "(synonymous) memory on the antisense strand — one object rotated by ρ.", "",
          "## 3c. Theorem R / FIX-39 — EXACT held-out three-term identity", "",
          "**d3b = Δ_x₁(test) + D_suffix(test) − D_gen**, all three terms ≥ 0.",
          "D_gen is the generalisation gap (estimation error + train/test shift).",
          "It explains negative held-out d3b at large lag even when Δ_x₁ > 0.", "",
          "| Species | k | d3b | d3b_phase | Δ_x₁ | D_suffix | D_gen | identity err | rank | semigroup err |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            L.append(f"| {r['species']} | {x['lag']} | {fnum(x['d3b'],6,True)} | "
                     f"{fnum(x['d3b_phase'],6,True)} | {fnum(x['Delta_x1_test'],6,True)} | "
                     f"{fnum(x['D_suffix_test'],6,True)} | {fnum(x['D_gen_train_test'],6,True)} | "
                     f"{x['thmT_heldout_identity_abs_err']:.1e} | {x['rank_Q']} | "
                     f"{x['semigroup_closure_err_k']:.1e} |")
    L += ["", "### 3d. Theorem T.4 — three-term split of D_suffix", "",
          "D_suffix = D_3'wobble + D_bridge(KL) + D_wobble(KL).",
          "Full-torus kernel accounting: 720+720+2304+240+48 = **4032 = 64·63**.",
          "Sense accounting [FIX-40]: 675+675+2025+240+45 = **3660 = 61·60**.",
          "**288 = 240+48 is the UNLINKED screen family**; H3b itself is 48-dim",
          "(homogeneous) / 144-dim (phase), since A_k and R share the same P_nt.", "",
          "| Species | k | D_suffix | D_3'wobble | D_bridge | D_wobble | split err |",
          "|---|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            L.append(f"| {r['species']} | {x['lag']} | {fnum(x['D_suffix_test'],6,True)} | "
                     f"{fnum(x['D_3prime_wobble'],6,True)} | {fnum(x['D_bridge_kl'],6,True)} | "
                     f"{fnum(x['D_wobble_kl'],6,True)} | "
                     f"{x['thmT_suffix_split_abs_err']:.1e} |")
    L += ["", "### 3e. Theorem T.5 — EXACT five-term held-out decomposition", "",
          "**d3b = D_5'lead + D_3'wobble + D_1↔3 + D_bridge(KL) + D_wobble(KL) − D_gen.**",
          "The three MI layers and D_gen are model-free, hence",
          "**d3b − d3b_phase = (D_bridge − D_bridge_ph) + (D_wobble − D_wobble_ph)**.",
          "By FIX-37, R_wobble_phase = P_train(Y₃|Y₁,Y₂), so D_wobble_ph measures only",
          "the train→test wobble-marginal shift (zero on train = test).", "",
          "| Species | k | 3 MI layers | D_bridge | D_bridge_ph | D_wobble | D_wobble_ph | D_gen | 5-term err | phase-diff err |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            L.append(f"| {r['species']} | {x['lag']} | {fnum(x.get('I_nonbridge_test'),6,True)} | "
                     f"{fnum(x['D_bridge_kl'],6,True)} | {fnum(x.get('D_bridge_kl_phase'),6,True)} | "
                     f"{fnum(x['D_wobble_kl'],6,True)} | {fnum(x.get('D_wobble_kl_phase'),6,True)} | "
                     f"{fnum(x['D_gen_train_test'],6,True)} | "
                     f"{x.get('thmT5_five_term_abs_err',0):.1e} | "
                     f"{x.get('thmT5_phase_diff_abs_err',0):.1e} |")
    L += ["", "### 3f. Theorem U — 8 unsplit 4-boxes vs 8 split boxes (duplex duality)", "",
          "B4 = {CT,GT,TC,CC,AC,GC,CG,GG} (sense wobble 100% synonymous);",
          "S4 = ρ̃(B4) = {AG,AC,GA,GG,GT,GC,CG,CC} (antisense wobble 100% synonymous).",
          "(D_5'lead^S4)^anti = D_3'wobble^B4 and (D_5'lead^Ssplit)^anti = D_3'wobble^Bsplit.", "",
          "| Species | k | D_wob(B4) | D_wob(split) | D_lead(S4) | D_lead(Ssplit) | box err | anti-box err |",
          "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            d = x.get("antiparallel_dual") or {}
            eU = max(d.get("U_leadS4_vs_wobB4_abs_err", 0.0),
                     d.get("U_leadSsplit_vs_wobBsplit_abs_err", 0.0),
                     d.get("U_wobB4_vs_leadS4_abs_err", 0.0),
                     d.get("U_wobBsplit_vs_leadSsplit_abs_err", 0.0)) if d else None
            L.append(f"| {r['species']} | {x['lag']} | {fnum(x.get('D_3prime_wobble_B4'),6,True)} | "
                     f"{fnum(x.get('D_3prime_wobble_Bsplit'),6,True)} | "
                     f"{fnum(x.get('D_5prime_lead_S4'),6,True)} | "
                     f"{fnum(x.get('D_5prime_lead_Ssplit'),6,True)} | "
                     f"{max(x.get('thmU_box_split_lead_abs_err',0), x.get('thmU_box_split_wobble_abs_err',0)):.1e} | "
                     f"{('n/a' if eU is None else format(eU, '.1e'))} |")
    L += ["", "### 3g. [FIX-40] Screen degrees of freedom: full torus vs sense", "",
          "| Species | table | q_sense | bridge | lead | wobble | d13 | screen theory | G.2 theory | screen emp | G.2 emp | B4 boxes |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in R:
        ds = r.get("screen_df_sense") or {}
        bpx = r.get("box_partition") or {}
        dse = r.get("screen_df_sense_empirical") or {}
        dst = r.get("df_sense_theory") or {}
        L.append(f"| {r['species']} | {r['transl_table']} | {ds.get('q_codons','n/a')} | "
                 f"{ds.get('bridge','n/a')} | {ds.get('lead','n/a')} | {ds.get('wobble','n/a')} | "
                 f"{ds.get('d13','n/a')} | {ds.get('sum','n/a')} | "
                 f"{dst.get('df_sum_check','n/a')} | {dse.get('sum','n/a')} | "
                 f"{r['df_sense_support']['df_sum_check']} | "
                 f"{','.join(bpx.get('B4_boxes',[]))} |")
    L += ["", "## 3h. Theorems Y / Z — exact bridge between G.2 and T (sense convention)", "",
    "Y: bipartite O–A graph is a FOREST; |O|+|A| = |O∨A|+|O∧A| (modular diamond).",
    "Balance: D_3'wobble + D_1↔3 − D_split-AA = S_syn − S_syn^hex-box.",
    "Z: tower σ(O) ⊂ σ(O,A) ⊂ σ(O,A,A^ρ) ⊂ σ(Y); S_dup = S_dup^lead + S_dup^wob.",
    "Corrected diamond: S_syn^hex-box + I(X;A^ρ|A,O) = S_anti + S_dup^lead.", "",
    "| Species | k | O | A | O∨A | comp | forest | Y balance | Z tower | Z chain | Z diamond |",
    "|---|---:|---:|---:|---:|---:|---|---:|---|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            yv, zv = x.get("thmY"), x.get("thmZ")
            if not yv or not zv:
                continue
            L.append(f"| {r['species']} | {x['lag']} | {yv['n_O']} | {yv['n_A']} | "
                     f"{yv['n_OA']} | {yv['n_components']} | `{yv['is_forest']}` | "
                     f"{yv['balance_abs_err']:.1e} | {'→'.join(map(str, zv['tower']))} | "
                     f"{zv['chain_abs_err']:.1e} | {zv['diamond_abs_err']:.1e} |")
    L += ["", "### 3j. Theorem AA — left tower, d2b = screen + lead-doublets, 3x3 self-dual matrix", "",
    "Left sigma(S,A) is NOT a forest: 4 components / 25 cycles (sense), 3 / 26 (full).",
    "d2b = I(S;A'|A) + I(X;A'|A,S); sense df 703 + 76 = 779, full 780 + 80 = 860.",
    "M_ij = I(l_i ; r_j | L_<i, R_<j); M(B_anti) = M(B_sense)^T exactly.", "",
    "| Species | k | d2b | d2b_screen | d2b_lead | AA err | M22 | M22 share of D13 | dual err |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            av, dv = x.get("thmAA"), x.get("thmAA_d2b")
            if not av or not dv:
                continue
            _m = av["M"]; _d13 = sum(_m[i][j] for i in (1, 2) for j in (1, 2))
            _sh = (_m[1][1] / _d13) if _d13 > 0 else float("nan")
            L.append(f"| {r['species']} | {x['lag']} | {fnum(dv['d2b'],6,True)} | "
                     f"{fnum(dv['d2b_screen'],6,True)} | {fnum(dv['d2b_lead'],6,True)} | "
                     f"{max(av['T_split_abs_err'], dv['residual']):.1e} | "
                     f"{fnum(_m[1][1],6,True)} | {fnum(100*_sh,2)}% | "
                     f"{av['dual_transpose_abs_err']:.1e} |")
    L += ["", "### 3i. Theorem W — stop-leakage / AA-split / pure synonymous", "",
    "| Species | k | D_stop | D_AA_split | D_syn | W err | D_stop(phase) | kernel shift ≥0 |",
    "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_R") or {}).get("lags", []):
            wv, wp = x.get("thmW"), x.get("thmW_phase")
            if not wv:
                continue
            L.append(f"| {r['species']} | {x['lag']} | {fnum(wv['D_stop_leakage'],6,True)} | "
                     f"{fnum(wv['D_AA_split'],6,True)} | {fnum(wv['D_syn'],6,True)} | "
                     f"{wv['residual']:.1e} | "
                     f"{fnum((wp or {}).get('D_stop_leakage'),6,True)} | "
                     f"{fnum(x.get('D_wobble_kernel_shift_nonneg'),6,True)} |")
    L += ["", "## 4. Theorem O — reading-frame phase locking (frames aligned)", "",
          "| Species | k | Ω_phase | Ω_marg | Ω_cond | bias_cond | pairs equal | S_syn(+0/+1/+2) |",
          "|---|---:|---:|---:|---:|---:|---|---|"]
    for r in R:
        for o in (r.get("theorem_O") or {}).get("lags", []):
            fr = o["frames"]
            L.append(f"| {r['species']} | {o['lag']} | {fnum(o['Omega_phase_JS_nats'],6)} | "
                     f"{fnum(o['Omega_marg_JS_nats'],6)} | {fnum(o['Omega_cond_JS_nats'],6)} | "
                     f"{o['bias_cond']:.2e} | `{o['frame_pair_counts_equal']}` | "
                     f"{fnum(fr['frame_plus0']['S_syn'],4)}/"
                     f"{fnum(fr['frame_plus1']['S_syn'],4)}/"
                     f"{fnum(fr['frame_plus2']['S_syn'],4)} |")
    L += ["", "χ² df is a BIAS SCALE only (frames are dependent overlapping readings).",
          "Null fixture: B(4,6) needs THREE periods since 4⁶ ≡ 1 (mod 3).", "",
          "## 5. Theorem G.2 / J — duplex split, gated enrichment, CpG quartets", "",
          "| Species | k | raw anti % | enrich_raw | enrich_eff | gate | CpG % of S_dup | R_cap |",
          "|---|---:|---:|---:|---:|---|---:|---:|"]
    for r in R:
        for g in (r.get("theorem_G") or {}).get("lags", []):
            L.append(f"| {r['species']} | {g['lag']} | "
                     f"{fnum(100*(g['antisense_fraction_raw'] or 0),1)} | "
                     f"{fnum(g.get('enrich_per_df_raw'),3)} | "
                     f"{fnum(g.get('enrich_per_df_eff'),3)} | "
                     f"`{g.get('enrich_gate_ok_eff')}` | "
                     f"{fnum(100*(g.get('cpg_share_of_S_dup') or 0),1)} | "
                     f"{fnum(g.get('R_cap'),3)} |")
    L += ["", "## 6. Theorem D — exact H5/H6/H7 on the COMMON cohort", "",
          "| Species | control | n common | Δ per transition | 95% CI | t |",
          "|---|---|---:|---:|---|---:|"]
    for r in R:
        for kk, v in (r.get("theorem_D") or {}).items():
            L.append(f"| {r['species']} | `{kk}` | {v['n_genes_common']} | "
                     f"{fnum(v['delta_per_transition_gene_mean'],6,True)} | "
                     f"[{fnum(v['ci95_low'],6,True)}, {fnum(v['ci95_high'],6,True)}] | "
                     f"{fnum(v['t_stat'],2,True)} |")
    L += ["", "## 7. Parse / split audit", "",
          "| Species | records | kept | ATG raw | ATG kept | stop raw | stop kept | mode | roots |",
          "|---|---:|---:|---:|---:|---:|---:|---|---:|"]
    for r in R:
        s, si = r["parse_stats"], r["split_info"]
        L.append(f"| {r['species']} | {s['n_records']} | {s['kept']} | "
                 f"{s['n_atg_start_raw']} | {s['n_atg_start_kept']} | "
                 f"{s['n_stop_end_raw']} | {s['n_stop_end_kept']} | "
                 f"`{si['mode']}` | {si['n_family_roots']} |")
    L += ["", "## 8. Theorem AD — Isochore-free screen decomposition against E_het^(k)", "",
          "Separates between-gene GC3 heterogeneity (isochores) from true within-gene memory across the 16×16 screen.",
          "Exact KL chain rule: **D_KL(P‖P_het) = K_bridge + K_wob|SO + K_lead|SY** (all ≥ 0).", "",
          "| Species | k | wob_MM | wob_het | wob_intra | dx1_MM | dx1_het | dx1_intra | K_nonbridge | K_nonbr_MM | K_lead|SY | K_l|SY_MM |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        for x in (r.get("theorem_AD") or {}).get("lags", []):
            if x["lag"] not in (1, 10, 30):
                continue
            sk = x.get("screen_kl") or {}
            L.append(f"| {r['species']} | {x['lag']} | {fnum(x['D_3prime_wobble_MM'],6,True)} | "
                     f"{fnum(x['D_3prime_wobble_het'],6,True)} | {fnum(x['D_3prime_wobble_intra_MM'],6,True)} | "
                     f"{fnum(x['Delta_x1_MM'],6,True)} | {fnum(x['Delta_x1_het'],6,True)} | "
                     f"{fnum(x['Delta_x1_intra_MM'],6,True)} | {fnum(sk.get('K_nonbridge'),6,True)} | "
                     f"{fnum(sk.get('K_nonbridge_MM'),6,True)} | {fnum(sk.get('K_lead_given_SY'),6,True)} | "
                     f"{fnum(sk.get('K_lead_given_SY_MM'),6,True)} |")
    L += ["", "## 9. Theorem AB — Periodicity & C2H2 Zinc-Finger stratification (Policy 10.1)", "",
          "| Species | Clade | ZNF genes | dip14 (all) | dip14 (non-ZNF) | dip14 (ZNF-only) | DFT peak P | SNR@14 | Verdict |",
          "|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for r in R:
        ab = r.get("theorem_AB") or {}
        if not ab:
            continue
        st = ab.get("stratified") or {}
        d_non = (st.get("non_znf_mp") or {}).get("dip14_log10")
        d_znf = (st.get("znf_only_mp") or {}).get("dip14_log10")
        L.append(f"| {r['species']} | {r.get('clade','n/a')} | {ab['n_znf_c2h2_genes']}/{ab['n_genes_total']} | "
                 f"{fnum(ab['all_mp'].get('dip14_log10'),3,True)} | {fnum(d_non,3,True)} | {fnum(d_znf,3,True)} | "
                 f"{fnum(ab['all_mp'].get('dft_peak_period'),1)} | {fnum(ab['all_mp'].get('dft_snr_at_14'),2)} | "
                 f"`{ab['verdict']}` |")
    ac = payload.get("theorem_AC") or {}
    if ac:
        L += ["", "## 10. Theorem AC — Zero-RNG Biosphere Phylogeny & Multi-Code Fingerprints", "",
              "| Channel | LOO 1-NN Hits | Eligible | Clade Purity % |",
              "|---|---:|---:|---:|"]
        for ch, v in (ac.get("loo_1nn_clade_purity") or {}).items():
            acc_pct = (100.0 * v["accuracy"]) if v.get("accuracy") is not None else None
            L.append(f"| `{ch}` | {v['hits']} | {v['eligible']} | {fnum(acc_pct,2)}% |")
        L += ["", "**Genetic-code commutator fingerprint clusters (zero species labels):**", ""]
        for fp_k, members in (ac.get("code_fingerprint_groups") or {}).items():
            tabs_set = sorted({m["table"] for m in members})
            sample_sp = ", ".join(m["species"] for m in members[:6]) + ("..." if len(members) > 6 else "")
            L.append(f"- `{fp_k}`: **{len(members)} taxa** (tables {tabs_set}): {sample_sp}")
    L.append("")
    atomic_write(path, "\n".join(L) + "\n")


# =====================================================================
# 19b. Unified Publication-Grade Plotter (merged replot_memory_v24.py
#      + 100-Taxon Clade-Median/IQR Mode + New Discovery Figs 5, 6, 7)
# =====================================================================
CLADE_COLORS = {
    "mammal": "#1f77b4", "bird": "#ff7f0e", "reptile": "#8c564b",
    "amphibian": "#e377c2", "fish": "#17becf", "insect": "#2ca02c",
    "nematode": "#bcbd22", "invert": "#7f7f7f", "plant": "#1b9e77",
    "fungus": "#9467bd", "protist": "#d62728", "bacterium": "#393b79",
    "archaeon": "#e6ab02", "unknown": "#555555",
}


def build_species_palette(n):
    """Combines tab20 + tab20b + tab20c (60 distinct colours); turbo beyond 60."""
    c1 = plt.cm.tab20(np.linspace(0, 1, 20))
    c2 = plt.cm.tab20b(np.linspace(0, 1, 20))
    c3 = plt.cm.tab20c(np.linspace(0, 1, 20))
    joined = np.vstack([c1, c2, c3])
    if n <= len(joined):
        return [joined[i] for i in range(n)]
    return [plt.cm.turbo(i / max(1, n - 1)) for i in range(n)]


def _add_bottom_legend(fig, labels, colors_map, ncol=9, fontsize=7.0):
    handles = [
        Line2D([0], [0], color=colors_map[lb], lw=1.6, marker="o", markersize=3.2, label=lb)
        for lb in labels
    ]
    n_rows = math.ceil(len(labels) / max(1, ncol))
    bottom_pad = min(0.42, 0.07 + 0.025 * n_rows)
    fig.subplots_adjust(top=0.87, bottom=bottom_pad, left=0.055, right=0.985, wspace=0.22)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=ncol, fontsize=fontsize, frameon=True, framealpha=0.95,
               edgecolor="#cccccc", handlelength=1.4, columnspacing=0.9, labelspacing=0.3)


def _plot_clade_bands(ax, clade_series_dict, ls="-", marker="o", alpha_band=0.16, lw=1.8):
    """Plots clade medians with 25-75% IQR bands (clean for 50-100+ taxa)."""
    for cl, curves in sorted(clade_series_dict.items()):
        if not curves:
            continue
        common_k = sorted(set.intersection(*(set(c.keys()) for c in curves)))
        if not common_k:
            continue
        M = np.array([[c[k] for k in common_k] for c in curves], dtype=np.float64)
        med = np.nanmedian(M, axis=0)
        q25 = np.nanpercentile(M, 25, axis=0)
        q75 = np.nanpercentile(M, 75, axis=0)
        col = CLADE_COLORS.get(cl, "#444444")
        ax.plot(common_k, med, color=col, ls=ls, marker=marker if len(common_k) <= 15 else None,
                ms=3.2, lw=lw, label=f"{cl} (n={len(curves)})")
        if len(curves) >= 3:
            ax.fill_between(common_k, q25, q75, color=col, alpha=alpha_band, lw=0)


def make_figures(payload, out_dir, dpi=300, fmt="png", plot_mode="auto"):
    if not HAVE_MPL or not payload.get("results"):
        return []
    fd = Path(out_dir) / "figures"
    fd.mkdir(parents=True, exist_ok=True)
    res = payload["results"]
    saved = []
    n_sp = len(res)
    species_list = [r["species"] for r in res]
    clade_of = {r["species"]: r.get("clade", NAME2CLADE.get(r["species"], "unknown")) for r in res}
    pal = build_species_palette(n_sp)
    sp_colors = {sp: pal[i] for i, sp in enumerate(species_list)}
    modes = (["lines", "clade"] if (plot_mode == "both" or (plot_mode == "auto" and n_sp > 20))
             else ([plot_mode] if plot_mode in ("lines", "clade") else ["lines"]))

    def _sv(fig, stem):
        p_out = fd / f"{stem}.{fmt}"
        fig.savefig(p_out, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        saved.append(str(p_out))
        print(f"  [FIG] {p_out}")

    def _finish_leg(fig, use_clade):
        if use_clade:
            cl_present = sorted({clade_of[sp] for sp in species_list})
            cl_map = {c: CLADE_COLORS.get(c, "#444444") for c in cl_present}
            _add_bottom_legend(fig, cl_present, cl_map, ncol=min(7, max(1, len(cl_present))), fontsize=8.5)
        else:
            ncol = 9 if n_sp > 30 else (5 if n_sp > 15 else max(1, min(5, n_sp)))
            fs = 6.2 if n_sp > 60 else (7.3 if n_sp > 25 else 8.2)
            _add_bottom_legend(fig, species_list, sp_colors, ncol=ncol, fontsize=fs)

    # ---- FIGS 1-4 (per mode: lines and/or clade median+IQR) ----
    for md in modes:
        use_cl = (md == "clade")
        sfx = "_clades" if use_cl else ""

        # Fig 1: Theorem M — Klein V_4 quadrants
        fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
        fig.suptitle(f"Fig 1. Thm M — static form vs directed program ({n_sp} taxa{', clade medians+IQR' if use_cl else ''})",
                     fontsize=13, fontweight="bold", y=0.96)
        panels1 = [
            ("norm2_W_pm_static_strand", "‖W_+-‖² static strand (het null CAN reach)", False),
            ("norm2_W_mp_sym_time_flux", "‖W_-+‖² directed (het null = 0)", True),
            ("norm2_W_mm_duplex_time_flux", "‖W_--‖² directed duplex (het null = 0)", True),
        ]
        for ax, (key, title, v30) in zip(axes, panels1):
            if use_cl:
                by_cl = {}
                for r in res:
                    pr = r.get("theorem_C", {}).get("lag_profile", [])
                    d = {row["lag"]: max(row["theorem_M"].get(key, 0.0), 1e-12) for row in pr if row.get("theorem_M")}
                    if d:
                        by_cl.setdefault(clade_of[r["species"]], []).append(d)
                _plot_clade_bands(ax, by_cl)
            else:
                for r in res:
                    pr = r.get("theorem_C", {}).get("lag_profile", [])
                    ks = [row["lag"] for row in pr if row.get("theorem_M")]
                    vs = [max(row["theorem_M"].get(key, 0.0), 1e-12) for row in pr if row.get("theorem_M")]
                    if ks:
                        ax.plot(ks, vs, color=sp_colors[r["species"]], lw=0.95, alpha=0.82)
            ax.set_yscale("log"); ax.set_title(title, fontsize=10.5); ax.set_xlabel("lag k (codons)")
            ax.grid(True, which="both", ls=":", alpha=0.35)
            if v30:
                ax.axvline(30, color="crimson", ls=":", lw=1.2)
        _finish_leg(fig, use_cl)
        _sv(fig, f"fig1_klein_quadrants{sfx}")

        # Fig 2: Theorem T — two-sided 16x16 screen layers
        fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
        fig.suptitle(f"Fig 2. Thm T — two-sided 16×16 screen (225 of 3600 sense df visible to H3b; {n_sp} taxa)",
                     fontsize=13, fontweight="bold", y=0.96)
        panels2 = [
            ("I_bridge", "I_bridge (225 df) — ALL H3b can see"),
            ("D_5prime_lead", "D_5'lead (675 sense df) — H3b ≡ 0"),
            ("D_3prime_wobble", "D_3'wobble (675 sense df) — H3b ≡ 0"),
        ]
        for ax, (key, title) in zip(axes, panels2):
            if use_cl:
                by_cl = {}
                for r in res:
                    rl = (r.get("theorem_R") or {}).get("lags", [])
                    d = {x["lag"]: x[key] for x in rl if _finite(x.get(key))}
                    if d:
                        by_cl.setdefault(clade_of[r["species"]], []).append(d)
                _plot_clade_bands(ax, by_cl)
            else:
                for r in res:
                    rl = (r.get("theorem_R") or {}).get("lags", [])
                    if rl:
                        ax.plot([x["lag"] for x in rl], [x.get(key, np.nan) for x in rl],
                                color=sp_colors[r["species"]], marker="o", ms=2.8, lw=1.0, alpha=0.82)
            ax.axhline(0.0, color="k", ls="--", lw=0.9); ax.axvline(30, color="crimson", ls=":", lw=1.2)
            ax.set_title(title, fontsize=10.5); ax.set_xlabel("lag k (codons)"); ax.grid(True, ls=":", alpha=0.35)
        _finish_leg(fig, use_cl)
        _sv(fig, f"fig2_thmT_screen_layers{sfx}")

        # Fig 3: FIX-39 / T.5 — Exact held-out identity
        fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
        fig.suptitle(f"Fig 3. FIX-39 / T.5 — d3b = Δ_x₁ + D_suffix − D_gen (EXACT, held-out; {n_sp} taxa)",
                     fontsize=13, fontweight="bold", y=0.96)
        for ax, key, title in ((axes[0], "d3b", "Held-out d3b(k)"),
                               (axes[2], "D_gen_train_test", "D_gen — generalisation gap (subtracted)")):
            if use_cl:
                by_cl = {}
                for r in res:
                    rl = (r.get("theorem_R") or {}).get("lags", [])
                    d = {x["lag"]: x[key] for x in rl if _finite(x.get(key))}
                    if d:
                        by_cl.setdefault(clade_of[r["species"]], []).append(d)
                _plot_clade_bands(ax, by_cl)
            else:
                for r in res:
                    rl = (r.get("theorem_R") or {}).get("lags", [])
                    if rl:
                        ax.plot([x["lag"] for x in rl], [x.get(key, np.nan) for x in rl],
                                color=sp_colors[r["species"]], marker="o", ms=2.8, lw=1.0, alpha=0.82)
            ax.axhline(0.0, color="k", ls="--", lw=0.9); ax.axvline(30, color="crimson", ls=":", lw=1.2)
            ax.set_title(title, fontsize=10.5); ax.set_xlabel("lag k (codons)"); ax.grid(True, ls=":", alpha=0.35)
        # Middle panel: Delta_x1 (solid) vs D_suffix (dashed)
        ax1 = axes[1]
        if use_cl:
            b1, b2 = {}, {}
            for r in res:
                rl = (r.get("theorem_R") or {}).get("lags", [])
                d1 = {x["lag"]: x["Delta_x1_test"] for x in rl if _finite(x.get("Delta_x1_test"))}
                d2 = {x["lag"]: x["D_suffix_test"] for x in rl if _finite(x.get("D_suffix_test"))}
                if d1:
                    b1.setdefault(clade_of[r["species"]], []).append(d1)
                if d2:
                    b2.setdefault(clade_of[r["species"]], []).append(d2)
            _plot_clade_bands(ax1, b1, ls="-", alpha_band=0.14)
            _plot_clade_bands(ax1, b2, ls="--", marker="s", alpha_band=0.08)
        else:
            for r in res:
                rl = (r.get("theorem_R") or {}).get("lags", [])
                if rl:
                    ks = [x["lag"] for x in rl]
                    ax1.plot(ks, [x.get("Delta_x1_test", np.nan) for x in rl], color=sp_colors[r["species"]], ls="-", marker="o", ms=2.6, lw=1.0, alpha=0.85)
                    ax1.plot(ks, [x.get("D_suffix_test", np.nan) for x in rl], color=sp_colors[r["species"]], ls="--", marker="s", ms=2.3, lw=0.9, alpha=0.72)
        ax1.axhline(0.0, color="k", ls="--", lw=0.9); ax1.axvline(30, color="crimson", ls=":", lw=1.2)
        ax1.set_title("Δ_x₁ (solid) and D_suffix (dashed)", fontsize=10.5); ax1.set_xlabel("lag k (codons)"); ax1.grid(True, ls=":", alpha=0.35)
        ax1.legend(handles=[Line2D([0], [0], color="#333", lw=1.5, ls="-", marker="o", ms=3.5, label="Δ_x₁ (solid)"),
                            Line2D([0], [0], color="#333", lw=1.5, ls="--", marker="s", ms=3.5, label="D_suffix (dashed)")],
                   loc="upper right", fontsize=8.5, framealpha=0.9)
        _finish_leg(fig, use_cl)
        _sv(fig, f"fig3_heldout_identity{sfx}")

        # Fig 4: Theorems N / O / U
        fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
        fig.suptitle(f"Fig 4. Thm N / O / U ({n_sp} taxa)", fontsize=13, fontweight="bold", y=0.96)
        if use_cl:
            n_aa, n_dp, o_cd, u_b4, u_sp = {}, {}, {}, {}, {}
            for r in res:
                cl = clade_of[r["species"]]
                pts = (r.get("theorem_tower") or {}).get("partitions") or {}
                if pts.get("aa"):
                    n_aa.setdefault(cl, []).append({x["lag"]: x["symmetric"]["D_D8_sym"] for x in pts["aa"] if "symmetric" in x})
                if pts.get("dup"):
                    n_dp.setdefault(cl, []).append({x["lag"]: x["symmetric"]["D_D8_sym"] for x in pts["dup"] if "symmetric" in x})
                ol = (r.get("theorem_O") or {}).get("lags", [])
                if ol:
                    o_cd.setdefault(cl, []).append({o["lag"]: o["Omega_cond_JS_nats"] for o in ol})
                rl = (r.get("theorem_R") or {}).get("lags", [])
                if rl:
                    u_b4.setdefault(cl, []).append({x["lag"]: x.get("D_3prime_wobble_B4", np.nan) for x in rl})
                    u_sp.setdefault(cl, []).append({x["lag"]: x.get("D_3prime_wobble_Bsplit", np.nan) for x in rl})
            _plot_clade_bands(axes[0], n_aa, ls="-"); _plot_clade_bands(axes[0], n_dp, ls="--", marker="s", alpha_band=0.08)
            _plot_clade_bands(axes[1], o_cd, ls="-")
            _plot_clade_bands(axes[2], u_b4, ls="-"); _plot_clade_bands(axes[2], u_sp, ls="--", marker="s", alpha_band=0.08)
        else:
            for r in res:
                col = sp_colors[r["species"]]
                pts = (r.get("theorem_tower") or {}).get("partitions") or {}
                if pts.get("aa"):
                    axes[0].plot([x["lag"] for x in pts["aa"]], [x["symmetric"]["D_D8_sym"] for x in pts["aa"]], color=col, ls="-", marker="o", ms=2.6, lw=1.0, alpha=0.82)
                if pts.get("dup"):
                    axes[0].plot([x["lag"] for x in pts["dup"]], [x["symmetric"]["D_D8_sym"] for x in pts["dup"]], color=col, ls="--", marker="s", ms=2.3, lw=0.9, alpha=0.72)
                ol = (r.get("theorem_O") or {}).get("lags", [])
                if ol:
                    axes[1].plot([o["lag"] for o in ol], [o["Omega_cond_JS_nats"] for o in ol], color=col, ls="-", marker="o", ms=2.6, lw=1.0, alpha=0.82)
                    axes[1].plot([o["lag"] for o in ol], [o["bias_cond"] for o in ol], color=col, ls=":", lw=0.8, alpha=0.6)
                rl = (r.get("theorem_R") or {}).get("lags", [])
                if rl:
                    axes[2].plot([x["lag"] for x in rl], [x.get("D_3prime_wobble_B4", np.nan) for x in rl], color=col, ls="-", marker="o", ms=2.6, lw=1.0, alpha=0.82)
                    axes[2].plot([x["lag"] for x in rl], [x.get("D_3prime_wobble_Bsplit", np.nan) for x in rl], color=col, ls="--", marker="s", ms=2.3, lw=0.9, alpha=0.72)
        for ax, ttl in ((axes[0], "Thm N: pure 8-doublet residual"),
                        (axes[1], "Thm O: Ω_cond vs bias scale"),
                        (axes[2], "Thm U: wobble memory in B4 (solid) vs split (dashed)")):
            ax.axhline(0.0, color="k", ls="--", lw=0.9); ax.axvline(30, color="crimson", ls=":", lw=1.2)
            ax.set_title(ttl, fontsize=10.5); ax.set_xlabel("lag k (codons)"); ax.grid(True, ls=":", alpha=0.35)
        axes[0].legend(handles=[Line2D([0], [0], color="#333", lw=1.5, ls="-", label="Λ = aa"),
                                Line2D([0], [0], color="#333", lw=1.5, ls="--", label="Λ = dup")], loc="upper right", fontsize=8.5)
        axes[2].legend(handles=[Line2D([0], [0], color="#333", lw=1.5, ls="-", label="B4 unsplit"),
                                Line2D([0], [0], color="#333", lw=1.5, ls="--", label="Split boxes")], loc="upper right", fontsize=8.5)
        _finish_leg(fig, use_cl)
        _sv(fig, f"fig4_tower_phase_boxes{sfx}")

    # ---- FIG 5 (NEW): Theorem AD — Isochore-Free Screen Decomposition ----
    has_ad = any((r.get("theorem_AD") or {}).get("lags") for r in res)
    if has_ad:
        use_cl = (n_sp > 20)
        fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
        fig.suptitle(f"Fig 5. Thm AD — Isochore-free within-gene screen memory vs E_het^(k) ({n_sp} taxa)",
                     fontsize=13, fontweight="bold", y=0.96)
        has_knb_mm = any(
            "K_nonbridge_MM" in (x.get("screen_kl") or {})
            for r in res for x in (r.get("theorem_AD") or {}).get("lags", [])
        )
        knb_key = "K_nonbridge_MM" if has_knb_mm else "K_nonbridge"
        knb_ttl = ("Bias-corrected non-bridge KL (K_nonbr − df/(2N))" if has_knb_mm
                   else "Exact isochore-free non-bridge KL (D_KL − K_bridge ≥ 0)")
        specs5 = [
            ("D_3prime_wobble_intra_MM", "D_3'wobble intra-gene excess (MM − E_het)"),
            ("Delta_x1_intra_MM", "Δ_x₁ intra-gene excess (MM − E_het)"),
            (knb_key, knb_ttl),
        ]
        for ax, (key, title) in zip(axes, specs5):
            if use_cl:
                by_cl = {}
                for r in res:
                    al = (r.get("theorem_AD") or {}).get("lags", [])
                    d = {}
                    for x in al:
                        v = x["screen_kl"].get(key) if key.startswith("K_nonbridge") else x.get(key)
                        if _finite(v):
                            d[x["lag"]] = v
                    if d:
                        by_cl.setdefault(clade_of[r["species"]], []).append(d)
                _plot_clade_bands(ax, by_cl)
            else:
                for r in res:
                    al = (r.get("theorem_AD") or {}).get("lags", [])
                    if al:
                        ks = [x["lag"] for x in al]
                        vs = [x["screen_kl"].get(key, np.nan) if key.startswith("K_nonbridge") else x.get(key, np.nan) for x in al]
                        ax.plot(ks, vs, color=sp_colors[r["species"]], marker="o", ms=2.8, lw=1.1, alpha=0.85)
            ax.axhline(0.0, color="k", ls="--", lw=0.9); ax.axvline(30, color="crimson", ls=":", lw=1.2)
            ax.set_title(title, fontsize=10.5); ax.set_xlabel("lag k (codons)"); ax.grid(True, ls=":", alpha=0.35)
        _finish_leg(fig, use_cl)
        _sv(fig, "fig5_thmAD_isochore_free_screen")

    # ---- FIG 6 (NEW): Theorem AB — ZNF C2H2 Stratification & Delta_k=14 Dips ----
    ab_taxa = [r for r in res if (r.get("theorem_AB") or {}).get("all_mp", {}).get("dip14_log10") is not None]
    if ab_taxa:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.4))
        fig.suptitle(f"Fig 6. Thm AB — C2H2 Zinc-Finger Stratification & Δk=14 Periodicity Audit ({len(ab_taxa)} taxa)",
                     fontsize=13, fontweight="bold", y=0.98)
        # Panel A: dip14_log10 vs ZNF gene fraction across biosphere
        for r in ab_taxa:
            ab = r["theorem_AB"]
            cl = clade_of[r["species"]]
            axes[0].scatter(100.0 * ab["znf_gene_fraction"], ab["all_mp"]["dip14_log10"],
                            color=CLADE_COLORS.get(cl, "#444"), s=42, alpha=0.85, edgecolors="k", lw=0.4)
            if (ab["all_mp"]["dip14_log10"] < -0.08 or ab["znf_gene_fraction"] > 0.028
                    or r["species"] in ("Human_T2T", "M_musculus")):
                axes[0].annotate(r["species"], (100.0 * ab["znf_gene_fraction"], ab["all_mp"]["dip14_log10"]),
                                 fontsize=7, xytext=(4, 2), textcoords="offset points")
        axes[0].axhline(0.0, color="k", ls="--", lw=0.9)
        axes[0].set_xlabel("C2H2 ZNF genes (% of test cohort)"); axes[0].set_ylabel("Δk=14 dip depth (log₁₀)")
        axes[0].set_title("Biosphere Δk=14 dip vs C2H2 ZNF abundance"); axes[0].grid(True, ls=":", alpha=0.35)
        # Panel B: Strongest dip14 taxon with ZNF stratification: All vs Non-ZNF vs ZNF-only
        strat_taxa = [r for r in ab_taxa if (r["theorem_AB"].get("stratified") or {}).get("mp_all_series")]
        if strat_taxa:
            best = min(strat_taxa, key=lambda r: r["theorem_AB"]["all_mp"]["dip14_log10"])
            st = best["theorem_AB"]["stratified"]
            ks = np.arange(1, len(st["mp_all_series"]) + 1)
            axes[1].plot(ks, st["mp_all_series"], color="#1f77b4", lw=1.4, label="All genes")
            axes[1].plot(ks, st["mp_non_znf_series"], color="#2ca02c", lw=1.4, ls="--", label="Non-ZNF genes")
            axes[1].plot(ks, st["mp_znf_only_series"], color="#d62728", lw=1.2, alpha=0.8, label="ZNF-only genes")
            for h in (14, 28, 42, 56):
                if h <= ks[-1]:
                    axes[1].axvline(h, color="grey", ls=":", lw=0.9)
            axes[1].set_yscale("log"); axes[1].set_xlabel("lag k (codons)"); axes[1].set_ylabel("‖W_-+‖²")
            axes[1].set_title(f"Stratification showcase: {best['species']}"); axes[1].legend(fontsize=8); axes[1].grid(True, ls=":", alpha=0.35)
        else:
            axes[1].text(0.5, 0.5, "No taxon with >=3 ZNF genes", ha="center", va="center")
        # Panel C: DFT SNR at period 14 across taxa
        snrs = [(r["species"], r["theorem_AB"]["all_mp"].get("dft_snr_at_14") or 0.0, clade_of[r["species"]]) for r in ab_taxa]
        snrs.sort(key=lambda x: -x[1])
        top = snrs[:min(15, len(snrs))]
        axes[2].barh([t[0] for t in top][::-1], [t[1] for t in top][::-1],
                     color=[CLADE_COLORS.get(t[2], "#444") for t in top][::-1])
        axes[2].set_xlabel("DFT power / median (period 14)"); axes[2].set_title("Top Δk=14 spectral peaks")
        axes[2].grid(True, axis="x", ls=":", alpha=0.35)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        _sv(fig, "fig6_thmAB_znf_periodicity")

    # ---- FIG 7 (NEW): Theorem AC — Directed vs Static Phylogeny & Code Fingerprints ----
    ac = payload.get("theorem_AC")
    if ac and ac.get("loo_1nn_clade_purity"):
        fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
        fig.suptitle(f"Fig 7. Thm AC — Zero-RNG Phylogeny Channels & Genetic-Code Fingerprints ({ac['n_species']} taxa)",
                     fontsize=13, fontweight="bold", y=0.98)
        pur = ac["loo_1nn_clade_purity"]
        ch_meta = {
            "GC3_scalar": ("GC3\nscalar", "#7f7f7f"),
            "Static_Form_Wpp_Wpm": ("Static Form\n(W_++ + W_+-)", "#1f77b4"),
            "Directed_Program_Wmp_Wmm": ("Directed\n(Frob.)", "#d62728"),
            "Directed_Program_Cosine": ("Directed\n(Cosine)", "#ff7f0e"),
            "NonBridge_93p75_Screen": ("Non-Bridge\n(Frob.)", "#2ca02c"),
            "NonBridge_Cosine": ("Non-Bridge\n(Cosine)", "#17becf"),
        }
        ch_names = list(pur.keys())
        accs = [100.0 * (pur[c]["accuracy"] or 0.0) for c in ch_names]
        cols = [ch_meta.get(c, (c, "#444444"))[1] for c in ch_names]
        lbls = [ch_meta.get(c, (c, "#444444"))[0] for c in ch_names]
        axes[0].bar(range(len(ch_names)), accs, color=cols)
        axes[0].set_xticks(range(len(ch_names)))
        axes[0].set_xticklabels(lbls, fontsize=8 if len(ch_names) > 4 else 9)
        axes[0].set_ylabel("LOO 1-NN Clade Purity (%)"); axes[0].set_ylim(0, 108)
        for i, v in enumerate(accs):
            axes[0].text(i, v + 2.0, f"{v:.1f}%", ha="center", fontsize=9, fontweight="bold")
        axes[0].set_title("Clade classification without sequence alignment"); axes[0].grid(True, axis="y", ls=":", alpha=0.35)
        # Panel B: Genetic-code commutator fingerprint scatter (Tables 1/11, 4, 6, 12)
        for r in res:
            q = (r.get("theorem_Q") or {}).get("AA", {})
            if not q:
                continue
            f2 = q["commutator_R_fro"] ** 2
            hs2 = q["superoperator_commutator_HS_sq"]
            tb = r.get("transl_table", 1)
            col = {1: "#1f77b4", 11: "#17becf", 4: "#ff7f0e", 6: "#d62728", 12: "#9467bd"}.get(tb, "#333")
            axes[1].scatter(f2, hs2, color=col, s=65, alpha=0.85, edgecolors="k", lw=0.5)
        # Annotate unique algebraic commutator clusters (merges identical Tables 1 & 11)
        by_coord = {}
        for r in res:
            q = (r.get("theorem_Q") or {}).get("AA", {})
            if q:
                f2_val = q["commutator_R_fro"] ** 2
                hs2_val = q["superoperator_commutator_HS_sq"]
                ckey = (round(f2_val, 4), round(hs2_val, 4))
                by_coord.setdefault(ckey, []).append((int(r.get("transl_table", 1)), f2_val, hs2_val))
        offsets = {
            (1, 11): (8, 5),
            (1,): (8, 5),
            (4,): (-82, -14),
            (6,): (8, 5),
            (12,): (-82, -14),
        }
        for ckey, items in sorted(by_coord.items()):
            tbs = tuple(sorted({it[0] for it in items}))
            mx = float(np.mean([it[1] for it in items]))
            my = float(np.mean([it[2] for it in items]))
            if len(tbs) == 1:
                lbl = f"Table {tbs[0]} (n={len(items)})"
            else:
                lbl = "Tables " + " & ".join(map(str, tbs)) + f" (n={len(items)})"
            dx, dy = offsets.get(tbs, (8, 5))
            axes[1].annotate(lbl, (mx, my), fontsize=8.5, fontweight="bold",
                             xytext=(dx, dy), textcoords="offset points",
                             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75))
        axes[1].margins(x=0.14, y=0.12)
        axes[1].set_xlabel("‖[D_AA, R]‖_F²"); axes[1].set_ylabel("‖[Π_AA, R]‖_HS²")
        axes[1].set_title("Blind genetic-code commutator fingerprints (Thm Q.2)"); axes[1].grid(True, ls=":", alpha=0.35)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        _sv(fig, "fig7_thmAC_phylogeny_codes")

    return saved


# =====================================================================
# 20. CLI
# =====================================================================
def main(argv=None):
    p = argparse.ArgumentParser(description=f"Codon Memory Spectrum v{VERSION} (Unified Zero-RNG Engine & Plotter)")
    p.add_argument("--cds-dir", default="cds_cache")
    p.add_argument("--out-dir", default="results_memory_v30")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--replot-only", type=str, default="", help="Regenerate Figs 1-7 from an existing JSON (v2.4 or v3.0) and exit")
    p.add_argument("--dpi", type=int, default=300, help="Figure resolution (default: 300)")
    p.add_argument("--plot-format", choices=["png", "pdf", "svg", "jpg"], default="png")
    p.add_argument("--plot-mode", choices=["auto", "lines", "clade", "both"], default="auto")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--auto-download", action="store_true")
    p.add_argument("--biosphere-100", action="store_true", help="Run full curated 100-species cellular panel (Tables 1,4,6,11,12)")
    p.add_argument("--resume", action="store_true", help="Resume from existing partial JSON checkpoint in --out-dir")
    p.add_argument("--panel", default=",".join(list(PANEL)[:9]))
    p.add_argument("--all-cached", action="store_true")
    p.add_argument("--limit-species", type=int, default=0)
    p.add_argument("--max-genes", type=int, default=6000, help="0 = all genes")
    p.add_argument("--gene-order", choices=["hash", "alpha"], default="hash",
                   help="Deterministic ordering before --max-genes cap (hash avoids ZNF truncation)")
    p.add_argument("--min-genes", type=int, default=40)
    p.add_argument("--min-codons", type=int, default=8)
    p.add_argument("--train-pct", type=int, default=80)
    p.add_argument("--no-trim", action="store_true")
    p.add_argument("--kernel-mode", choices=["backoff", "laplace"], default="backoff")
    p.add_argument("--laplace-alpha", type=float, default=0.1)
    p.add_argument("--alphabets", default="codon64,aa20,mirror32,wobble4,gc3bin2")
    p.add_argument("--max-order", type=int, default=0)
    p.add_argument("--kappa", type=float, default=8.0)
    p.add_argument("--min-ctx-count", type=int, default=3)
    p.add_argument("--max-ctx", type=int, default=200000)
    p.add_argument("--max-lag", type=int, default=60)
    p.add_argument("--h4-lengths", default="2,3,5,10,20,29,30,50,100,200")
    p.add_argument("--dilution-tol", type=float, default=0.15)
    p.add_argument("--het-windows", default="10,30,100")
    p.add_argument("--het-window-lags", default="1,2,5,10,20,30,60")
    p.add_argument("--het-max-windows", type=int, default=150000)
    p.add_argument("--duplex-lags", default="1,2,3,5,8,10,15,20,25,30,40,60")
    p.add_argument("--no-duplex", action="store_true")
    p.add_argument("--no-tower", action="store_true")
    p.add_argument("--tower-partitions", default="het,aa,dup,mirror")
    p.add_argument("--tower-lags", default="1,2,3,5,8,10,15,20,25,30,40,60")
    p.add_argument("--tower-max-genes", type=int, default=0)
    p.add_argument("--no-phase", action="store_true")
    p.add_argument("--phase-lags", default="1,2,3,5,10,20,30")
    p.add_argument("--no-h3b", action="store_true")
    p.add_argument("--h3b-lags", default="1,2,3,5,10,20,30")
    p.add_argument("--no-h567", action="store_true")
    p.add_argument("--no-znf", action="store_true", help="Skip Theorem AB ZNF/periodicity stratification")
    p.add_argument("--h5-alpha", type=float, default=0.1)
    p.add_argument("--h5-max-genes", type=int, default=0)
    p.add_argument("--mc-check", type=int, default=0)
    p.add_argument("--seed", type=int, default=20261114)
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args(argv)

    # ---- Mode 1: Standalone zero-RNG self-test ----
    if args.self_test:
        print("=" * 100)
        print(f"SELF-TEST v{VERSION} — Theorems A..AA + AB (ZNF/Periodicity) + AC (Phylogeny) + AD (Isochore-Free Screen)")
        print("=" * 100)
        self_test(verbose=True)
        print("ALL SELF-TESTS PASSED.")
        return 0

    # ---- Mode 2: Standalone publication re-plotter (replaces replot_memory_v24.py) ----
    if args.replot_only:
        jp = Path(args.replot_only)
        if not jp.exists():
            print(f"[ERROR] JSON not found: {jp}")
            return 1
        payload = json.loads(jp.read_text(encoding="utf-8"))
        out_p = Path(args.out_dir).resolve() if args.out_dir != "results_memory_v30" else jp.parent
        figs = make_figures(payload, out_p, dpi=args.dpi, fmt=args.plot_format, plot_mode=args.plot_mode)
        print(f"[REPLOT DONE] {len(figs)} figures saved in {out_p / 'figures'}")
        return 0

    st = self_test(verbose=False)
    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    cds_path = Path(args.cds_dir)
    found = discover_cache(cds_path)
    name2file = {ACC2NAME.get(a, a): pp for a, pp in found.items()}

    print("=" * 100)
    print(f"CODON MEMORY SPECTRUM v{VERSION} | {len(found)} cached assemblies | "
          f"catalog={len(BIOSPHERE_100)} taxa | kernel={args.kernel_mode} | max_genes={args.max_genes or 'ALL'}")
    print("=" * 100)

    if args.biosphere_100:
        targets = []
        for nm, (acc, tbl, cl, org) in BIOSPHERE_100.items():
            if acc and acc in found:
                targets.append((nm, found[acc]))
            elif nm in name2file:
                targets.append((nm, name2file[nm]))
            elif args.auto_download:
                pp = auto_download(acc, cds_path) if acc else None
                if pp is None:
                    pp = auto_download_by_organism(nm, org, tbl, cl, cds_path)
                if pp is not None:
                    targets.append((nm, pp))
            else:
                print(f"[skip] {nm:24s} ({org}, Table {tbl}): not in cache (use --auto-download)")
    elif args.all_cached:
        targets = [(ACC2NAME.get(a, a), pp) for a, pp in found.items()]
    else:
        targets = []
        for nm in [x.strip() for x in args.panel.split(",") if x.strip()]:
            if nm not in BIOSPHERE_100:
                continue
            acc, tbl, cl, org = BIOSPHERE_100[nm]
            if acc and acc in found:
                targets.append((nm, found[acc]))
            elif nm in name2file:
                targets.append((nm, name2file[nm]))
            elif args.auto_download:
                pp = auto_download(acc, cds_path) if acc else None
                if pp is None:
                    pp = auto_download_by_organism(nm, org, tbl, cl, cds_path)
                if pp is not None:
                    targets.append((nm, pp))
            else:
                print(f"[skip] {nm}: {acc or org} not in cache")

    if args.limit_species > 0:
        targets = targets[:args.limit_species]
    if args.dry_run:
        for nm, pp in targets:
            print(f"  [target] {nm:24s} (table={NAME2TABLE.get(nm,1):2d}, clade={NAME2CLADE.get(nm,'?'):10s}) {pp.name}")
        return 0
    if not targets:
        print("[ERROR] no target assemblies.")
        return 1

    # Optional --resume support (with persistent matrix_cache_v30.npz for Thm AC)
    results, failures = [], []
    done_names = set()
    matrix_cache = {}
    partial_path = out / "memory_spectrum_v30.partial.json"
    mcache_path = out / "matrix_cache_v30.npz"
    if args.resume and partial_path.exists():
        try:
            prev = json.loads(partial_path.read_text(encoding="utf-8"))
            results = prev.get("results", [])
            done_names = {r["species"] for r in results}
            if mcache_path.exists():
                npz = np.load(mcache_path)
                for nm_ in list(done_names):
                    if f"{nm_}__P1" in npz and f"{nm_}__LP" in npz and f"{nm_}__LQ" in npz:
                        matrix_cache[nm_] = {
                            "P1_test": np.array(npz[f"{nm_}__P1"]),
                            "logP_train": np.array(npz[f"{nm_}__LP"]),
                            "logQ2b_train": np.array(npz[f"{nm_}__LQ"]),
                        }
                npz.close()
            print(f"[RESUME] Loaded {len(results)} species ({len(matrix_cache)} matrices) from {partial_path.name}")
        except Exception as e:
            print(f"[RESUME] Could not load partial checkpoint: {e}")
    t_all = time.time()
    for i, (nm, pp) in enumerate(targets, 1):
        if nm in done_names:
            print(f"--- [{i}/{len(targets)}] {nm} (resumed, skipping) ---")
            continue
        print(f"\n--- [{i}/{len(targets)}] ---")
        try:
            r = analyze_species(nm, pp, args)
            if r:
                mc = r.pop("_matrix_cache", None)
                if mc is not None:
                    matrix_cache[nm] = mc
                    npz_dict = {}
                    for k_nm, k_mc in matrix_cache.items():
                        npz_dict[f"{k_nm}__P1"] = k_mc["P1_test"]
                        npz_dict[f"{k_nm}__LP"] = k_mc["logP_train"]
                        npz_dict[f"{k_nm}__LQ"] = k_mc["logQ2b_train"]
                    tmp_npz = out / "matrix_cache_v30.tmp.npz"
                    np.savez_compressed(tmp_npz, **npz_dict)
                    os.replace(tmp_npz, mcache_path)
                results.append(r)
                atomic_write(partial_path,
                             json.dumps(_clean({"version": VERSION,
                                                "n_species": len(results),
                                                "results": results}), indent=2))
        except Exception as e:
            import traceback
            print(f"[ERROR] {nm}: {e}")
            traceback.print_exc()
            failures.append({"species": nm, "error": repr(e)})

    # ---------- [THM-AC] Biosphere-wide phylogeny & multi-code matrix ----------
    thm_ac = thm_AC_biosphere_phylogeny(results, matrix_cache)
    if thm_ac:
        print("\n" + "=" * 100)
        print("[THM-AC] ZERO-RNG BIOSPHERE LOO 1-NN CLADE PURITY:")
        for ch, v in thm_ac["loo_1nn_clade_purity"].items():
            print(f"  {ch:28s}: {v['hits']:3d}/{v['eligible']:3d} ({fnum(100*(v['accuracy'] or 0),2)}%)")
        print(f"  Genetic-code fingerprint clusters: {len(thm_ac['code_fingerprint_groups'])}")
        print("=" * 100)

    payload = {"version": VERSION, "self_test": st, "config": vars(args),
               "n_species": len(results), "n_failed": len(failures),
               "failures": failures, "runtime_seconds": time.time() - t_all,
               "theorem_AC": thm_ac, "results": results}
    atomic_write(out / "memory_spectrum_v30.json", json.dumps(_clean(payload), indent=2))
    write_md(payload, out / "memory_spectrum_v30.md")
    figs = [] if args.no_plots else make_figures(payload, out, dpi=args.dpi, fmt=args.plot_format, plot_mode=args.plot_mode)
    print("\n" + "=" * 100)
    print(f"DONE {len(results)} OK, {len(failures)} failed, {(time.time()-t_all)/60:.1f} min")
    print(f"[JSON] {out/'memory_spectrum_v30.json'}")
    print(f"[MD]   {out/'memory_spectrum_v30.md'}")
    for f in figs:
        print(f"[FIG]  {f}")
    print("=" * 100)
    return 0


if __name__ == "__main__":
    sys.exit(main())
