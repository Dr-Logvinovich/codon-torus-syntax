# From Form to Program & The Double-Stranded Codon Torus ($\mathbb{T}_4^3 = (\mathbb{Z}/4\mathbb{Z})^3$): Conditional Codon Syntax, Markovian Lumpability & Synthetic mRNA Forensics

[![Codon Torus & Syntax CI](https://github.com/Dr-Logvinovich/codon-torus-syntax/actions/workflows/ci.yml/badge.svg)](https://github.com/Dr-Logvinovich/codon-torus-syntax/actions/workflows/ci.yml)
[![ORCID: V. Logvinovich](https://img.shields.io/badge/ORCID-0009--0007--2040--6955-a6ce39.svg)](https://orcid.org/0009-0007-2040-6955)
[![ORCID: J.-C. Perez](https://img.shields.io/badge/ORCID-0000--0001--6446--2042-a6ce39.svg)](https://orcid.org/0000-0001-6446-2042)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B_%7C_NumPy_%7C_SciPy-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Official Open-Science Verification Repository for the Genomic Information Theory & Double-Stranded Codon Torus Research Programme.**
>
> * **Authors:** **Victor Logvinovich** (Master of Physical and Mathematical Sciences; Researcher in Pedagogical Sciences (Mathematics), Kobrin, Belarus; [ORCID: 0009-0007-2040-6955](https://orcid.org/0009-0007-2040-6955)) and **Jean-Claude Perez** (PhD in Mathematics and Computer Science, Retired IBM AI European Research Centre, Montpellier; Luc Montagnier Foundation; [ORCID: 0000-0001-6446-2042](https://orcid.org/0000-0001-6446-2042)).

---

## Repository Contents (Core 4 Artifacts)

| File | Type | Description |
| :--- | :---: | :--- |
| **[`code_logos_perez_v14_8.tex.pdf`](code_logos_perez_v14_8.tex.pdf)** | **Preprint I (34 pp.)** | *"From Form to Program: Medium-Range Conditional Codon Syntax Inside Real Genes — Strict Null-Model Evidence Across 54 Cellular RefSeq Taxa, a Codon-Mirror-Symmetry Control, Two-Scale Phylogenetic Transfer, Intron-Exon Compartment Confinement, Synthetic Genome Forensics, and the Biotech 'Max-CAI' Codon-Optimization Crash"* (v14.8). |
| **[`super_codon_proof_v14_8.py`](super_codon_proof_v14_8.py)** | **Verification Engine I** | Autonomous self-testing pipeline (v14.6–v14.8) auditing 54 cellular RefSeq taxa + SARS-CoV-2 across a 10-level null model hierarchy (`H0`–`H7`), 14-gene human RefSeqGene intron-exon compartment confinement (`TP53`, `KRAS`, `BRCA1`, etc.), 13 frontier genomes (Giant Viruses, 2.8 km deep boreholes, $122^\circ\text{C}$ archaea, Venter's `JCVI-syn1.0` & `JCVI-syn3.0`), viral "Genetic Trojans" (`HPV-16 E7`, `Adenovirus E3 RID-α`), organelle OS rewriting, and the `Max-CAI` codon-optimization crash test. |
| **[`codon_torus_logos_perez_v3.tex.pdf`](codon_torus_logos_perez_v3.tex.pdf)** | **Preprint II (33 pp.)** | *"The Double-Stranded Codon Torus: Duplex Mirror Fibration of DNA Coding Space, Markovian Lumpability, and Information-Theoretic Consequences"* (v3). |
| **[`codon_torus_real_audit_v3_bnt162b2.py`](codon_torus_real_audit_v3_bnt162b2.py)** | **Verification Engine II** | Simulation-free, zero-RNG analytic audit suite (v3.2 / v3.5) verifying exact chain-rule decompositions on $\mathbb{T}_4^3 = (\mathbb{Z}/4\mathbb{Z})^3$, the two-strand duplex fibration $\alpha_{\text{dup}}: 64 \to 56$, Kemeny–Snell Markov lumpability ($G^2$), antisense program dominance ($74.1\%\text{–}88.5\%$), 6-frame duplex ORF forensics on **Pfizer BNT162b2** (3,825 nt) and **Moderna mRNA-1273** (3,828 nt / `OR134578.1`), and the constructive **Corollary 4.4 Duplex-Safe Viterbi Compiler**. |

---

## Executive Summary & Key Mathematical / Biological Discoveries

### Part I: From Static Compositional "Form" to Dynamic Sequential "Program" (`v14.5–v14.8`)

Traditional numerical genomics characterizes zero-order marginal codon counts, GC3 isochores, and Chargaff/Perez 32-mirror-pair strand parity (**static form**). We operationalize the distinction between static form and dynamic sequential **program** via held-out first-order Markovian conditional transitions $P(c_{i+1} \mid c_i)$ evaluated on 20% held-out genes against a **10-level null model hierarchy**:
* **`H0`–`H3b` (Analytical Compositional & Trigram Nulls):** Marginal codon usage (`H0`), boundary nucleotide $P(c_{i+1}\mid n_3(c_i))$ (`H1`), amino-acid bigram syntax + marginal synonymous usage (`H2`), amino-acid syntax + empirical codon-pair usage (`H2b`), folded dinucleotides (`H3`), and folded order-2 nucleotide trigrams (`H3b`).
* **`H4` (Medium-Range Block Shuffle):** Intra-gene block permutations at $L \in \{3, 10, 30, 100\}$ codons, isolating the ribosomal exit-tunnel scale $L = 30$ codons ($29$ codon transitions = 8th Lucas number; $89$ nucleotide steps = 11th Fibonacci number).
* **`H5`–`H7` (Surrogate, Synonymous & Chargaff Mirror-Pair Controls):** Amino-acid-preserving `H2b` surrogates (`H5`), per-gene synonymous multiset shuffles (`H6`), and exact 32-pair reverse-complement mirror shuffles (`H7`).

**Core Empirical Findings (`code_logos_perez_v14_8.tex.pdf`):**
1. **Universal Biosphere Invariant ($54$ Cellular RefSeq Taxa, $0.1416 \le \text{GC3} \le 0.9177$):** All 52 cellular assemblies meeting the within-species sample-size guard ($\ge 20$ units and $\ge 50,000$ transitions in both train and test splits) pass all 7 strict bootstrap criteria (`d2b`, `d3`, `d3b`, `ΔH4(30)`, `Δ5`, `Δ6`, `Δ7` $> 0$ at 95% CI), from *Dictyostelium discoideum* ($\text{GC3}=0.1416$) to *Thermus thermophilus* ($\text{GC3}=0.9177$).
2. **Two-Scale Phylogenetic Transfer & Euarchontoglires Directional Asymmetry:** Coarse conditional syntax ($d_{2b}^{\text{transfer}} > 0.02$) transfers symmetrically across $132/132$ directed mammalian pairs (12 species, 6 orders) and $20/20$ four-class vertebrate pairs, whereas fine `H5` codon-choice conservation exhibits a reproducible **Primate $\to$ Rodent directional asymmetry** at near-zero GC3 distance ($\Delta\text{GC3} = 0.000\text{–}0.008$).
3. **14/14 Intron–Exon Compartment Confinement & Oncogenic Hotspots (`v14.6`):** Across 14 human RefSeqGene loci (`TP53`, `TP63`, `TP73`, `KRAS`, `NRAS`, `HRAS`, `BRCA1`, `PTEN`, `RB1`, `GAPDH`, `ACTB`, `HBB`, `MYC`, `APOE`), $d_{3b}$ is strictly positive in concatenated coding exons (mean $+0.1721\text{ nats/codon}$) and strictly negative in both excised introns (mean $-0.3097$, $p = 2^{-14} = 6.1 \times 10^{-5}$) and unspliced pre-mRNA across all 3 reading frames, while 32-pair mirror asymmetry ($\text{Mirror}_{L_1}$) jumps from $0.1615$ (introns) to $0.4176$ (exons). In `TP53`, the mutable DNA-binding core (exons 5–8, codons 125–300, housing 86% of cancer mutations) carries $2.34\times$ higher $d_{2b}$ ($+0.0989$ vs. $+0.0423$) than flanking exons.
4. **Frontier Biosphere & Synthetic Genome Forensics on Venter's `JCVI-syn1.0` / `JCVI-syn3.0` (`v14.7`):** All 4 giant DNA viruses (*Mimivirus*, *Megavirus*, *Pandoravirus*, 30,000-yr permafrost *Pithovirus*), all 4 deep-subsurface/abyssal extremophiles (*Desulforudis audaxviator* at 2.8 km depth, *Methanopyrus kandleri* at $122^\circ\text{C}$), and the Asgard archaeon *Prometheoarchaeum syntrophicum* MK-D1 pass `STRICT 7/7`. Scoring Venter's synthetic cells under wild-type *Mycoplasma mycoides* (Table 4) confirms $0/412$ essential genes disrupted in `JCVI-syn3.0`, while blindly pinpointing the foreign engineered vector genes in `JCVI-syn1.0` (*E. coli* `lacZ` $d_{3b} = -0.4564$, yeast `HIS3` $-0.4286$, `bla` $-0.2353$) and ASCII watermark cassettes ($d_{3b} \in [-0.4825, -0.1964]$).
5. **Viral "Genetic Trojans", Organelle Rewriting & Biotech `Max-CAI` Crash (`v14.8`):** Standard industrial `Max-CAI` codon optimization crashes $d_{2b}$ from $+0.0427$ to $-0.0589\text{ nats/codon}$ across $14/14$ human genes ($p = 6.1 \times 10^{-5}$). Blind viral scoring on the human transition matrix isolates host-mimicking effectors: **HPV-16 `E7` oncoprotein** is the sole gene in HPV-16 with positive human synonymous syntax ($d_{2b} = +0.0327$), and **Adenovirus `E3 RID-α`** scores $1.94\times$ the native human mean ($d_{2b} = +0.0958$). In *Arabidopsis* chloroplasts, $65/67$ plastid genes ($97.0\%$, $p = 3.9 \times 10^{-17}$) have been rewritten from cyanobacterial syntax to host nuclear syntax.

---

### Part II: The Double-Stranded Codon Torus $\mathbb{T}_4^3 = (\mathbb{Z}/4\mathbb{Z})^3$ & Synthetic mRNA Duplex Forensics (`v3.2–v3.5`)

In `codon_torus_logos_perez_v3.tex.pdf` and `codon_torus_real_audit_v3_bnt162b2.py`, we construct an exact, simulation-free mathematical theory on the double-stranded codon state space $\mathcal{D} = \{(c, \rho(c)) : c \in \mathbb{T}_4^3\}$:

1. **Fixed-Point-Free Isometric Involution & Orthogonal Parity Split (Theorem 2.5 & Corollary 2.6):** Watson–Crick reverse complementation $\rho(c) = -Rc + 3\mathbf{1} \pmod 4$ is a fixed-point-free affine isometry in both Hamming ($d_H$) and Lee ($d_{\text{Lee}}$) metrics, partitioning $\mathbb{T}_4^3$ into **32 duplex mirror orbits** and decomposing $\mathbb{R}^{64} = V_+ \oplus V_-$ ($\dim V_+ = \dim V_- = 32$) into mirror-even form ($V_+$) and mirror-odd program ($V_-$).
2. **The Two-Strand Duplex Fibration $\alpha_{\text{dup}}: 64 \to 56$ & the 8-Doublet Theorem (Theorem 4.3):** While the single-strand genetic code $\alpha: \mathbb{T}_4^3 \to \mathcal{A}$ ($64 \to 21$) leaves a **43-dimensional** synonymous simplex bundle ($41$ on sense codons), the two-strand fibration $\alpha_{\text{dup}}(c) = (\alpha(c), \alpha(\rho(c)))$ maps the sense wobble base onto the antisense leading base, partitioning $\mathbb{T}_4^3$ into **48 singletons** and **8 binary doublets** (`{TAA,TAG}`, `{CAA,CAG}`, `{TCG,TCT}`, `{CCG,CCT}` and their $\rho$-duals `{TTA,CTA}`, `{TTG,CTG}`, `{CGA,AGA}`, `{CGG,AGG}`), reducing the duplex synonymous dimension from $43$ to $64 - 56 = \mathbf{8}$.
3. **Exact Chain-Rule Decompositions & Antisense Channel Dominance (Theorems 5.3, 5.4, 5.5):**
   $$I(C_{i-1}; C_i) = I(A_{i-1}; A_i) + \underbrace{I(C_{i-1}; A_i \mid A_{i-1})}_{d_{2b}} + \underbrace{I(A_{i-1}; C_i \mid A_i) + I(C_{i-1}; C_i \mid A_{i-1}, A_i)}_{S_{\text{syn}} = I(C_{i-1}; C_i \mid A_i)}$$
   Across 9 RefSeq genomes ($9,834,024$ held-out transitions), the two-strand split $S_{\text{syn}} = S_{\text{anti}} + S_{\text{dup}}$ proves that **$74.1\%\text{–}88.5\%$** ($82.5\%$ in `Human_T2T`) of the sense-strand synonymous program is carried directly by the antisense amino-acid channel $A_i^\rho = \alpha(\rho(C_i))$.
4. **Kemeny–Snell Markov Lumpability (Theorem 5.6):** The amino-acid projection $A_i = \alpha(C_i)$ is a lumped Markov chain if and only if $d_{2b} = 0$ (rejected with $\text{df} = 779$, $p < 10^{-300}$ in all 9 species).
5. **Six-Frame Duplex Forensics on Pfizer `BNT162b2` & Moderna `mRNA-1273` vs. Corollary 4.4 Duplex-Safe Viterbi Compiler (Section 12.1):**
   * By **Corollary 4.4**, the sense preimages of antisense frame `-0` stop codons are $\rho^{-1}(\text{Stop}) = \{\text{TTA}, \text{CTA}, \text{TCA}\}$, and the preimage of antisense `ATG` is $\rho^{-1}(\{\text{ATG}\}) = \{\text{CAT}\}$.
   * **Wild-type SARS-CoV-2 Spike (`NC_045512.2`):** Contains $63$ $\{\text{TTA}, \text{CTA}, \text{TCA}\}$ codons ($28/9/26$), placing $63$ stop barriers in antisense frame `-0` (max unbroken antisense span: $106$ codons; max `-0` `ATG->Stop` ORF: $79$ codons).
   * **Pfizer/BioNTech `BNT162b2` (3,825 nt):** GC3 optimization ($\text{GC3}=0.8105$) deletes all $63$ $\{\text{TTA}, \text{CTA}, \text{TCA}\}$ codons ($0/0/0$), wiping out all stop barriers across all $1,275$ triplets of antisense frame `-0`, while retaining $2$ sense `CAT` codons—opening a continuous **1,101-codon (3,303-nt) `ATG->End` antisense ORF** (`175..1275`) and incurring a negative human syntax score ($d_{2b} = -0.0549$).
   * **Moderna `mRNA-1273` (3,828 nt & `OR134578.1`):** Replaces all `CAT` with `CAC` ($\text{CAT}=0$, yielding $0$ `ATG->End` ORFs in antisense frame `-0`), but extreme wobble-T depletion ($\text{GC3}=0.9552$) collapses shifted sense frame `+2` stops from $145$ (wild-type) to $2$, opening an **820-codon (2,460-nt) unbroken span in sense frame `+2`** and crashing $d_{2b}$ to $-0.0616$.
   * **Constructive Solution (`Human_DuplexSafe_Viterbi_Spike_2P`):** Our deterministic zero-RNG Viterbi compiler enforces $\text{CAT}=0$ and periodic $\{\text{TTA}, \text{CTA}, \text{TCA}\}$ barriers while maximizing first-order human Markov transitions, achieving strongly positive human syntax ($d_{2b} = +0.2142$, $S_{\text{syn}} = +0.1557$, $d_{3b} = +0.3304$), **0 antisense `ATG->End` ORFs in frame `-0`**, $30$ antisense `-0` stop barriers (max span $174$ codons), and $39$ sense `+2` stop barriers (max span $144$ codons).

---

## Quick-Start & Complete Reproduction Commands

```bash
# 1. Install scientific dependencies:
python3 -m pip install numpy scipy matplotlib

# 2. Run built-in mathematical & genetic-code self-tests:
python3 super_codon_proof_v14_8.py --self-test
python3 -c "import codon_torus_real_audit_v3_bnt162b2 as ct; print(ct.self_test_torus_theorems())"

# 3. Run BNT162b2 & Moderna mRNA-1273 Duplex Torus & Viterbi Compiler Audit:
python3 codon_torus_real_audit_v3_bnt162b2.py --cds-dir cds_cache --out-dir results_torus_v3 --bnt162b2-only

# 4. Run 9-Genome RefSeq Analytic Torus Audit + 17q22 Hourglass + BNT162b2 Audit:
python3 codon_torus_real_audit_v3_bnt162b2.py --cds-dir cds_cache --out-dir results_torus_v3 --bnt162b2-audit

# 5. Run 14-Gene Human RefSeqGene Intron-Exon & TP53/RAS Oncogenic Hotspot Audit (v14.6):
python3 super_codon_proof_v14_8.py --intron-exon-audit --cds-dir cds_cache --out-dir results_intron_exon_v14_6

# 6. Run 13-Genome Frontier Biosphere & Venter JCVI-syn1.0/3.0 Forensics Audit (v14.7):
python3 super_codon_proof_v14_8.py --frontier-audit --cds-dir cds_cache --out-dir results_frontier_v14_7

# 7. Run Viral "Genetic Trojan" Scanner, Organelle Rewriting & Biotech Max-CAI Crash Test (v14.8):
python3 super_codon_proof_v14_8.py --chimera-audit --cds-dir cds_cache --out-dir results_chimera_v14_8
```

---

## Citation (BibTeX)

```bibtex
@misc{LogvinovichPerez2026FormToProgram,
  author       = {Logvinovich, Victor and Perez, Jean-Claude},
  title        = {{From Form to Program: Medium-Range Conditional Codon Syntax Inside Real Genes}},
  year         = {2026},
  month        = oct,
  version      = {v14.8},
  note         = {GitHub Repository: [https://github.com/Dr-Logvinovich/codon-torus-syntax](https://github.com/Dr-Logvinovich/codon-torus-syntax), ORCID: 0009-0007-2040-6955, 0000-0001-6446-2042}
}

@misc{LogvinovichPerez2026CodonTorus,
  author       = {Logvinovich, Victor and Perez, Jean-Claude},
  title        = {{The Double-Stranded Codon Torus: Duplex Mirror Fibration of DNA Coding Space, Markovian Lumpability, and Information-Theoretic Consequences}},
  year         = {2026},
  month        = oct,
  version      = {v3.5},
  note         = {GitHub Repository: [https://github.com/Dr-Logvinovich/codon-torus-syntax](https://github.com/Dr-Logvinovich/codon-torus-syntax), ORCID: 0009-0007-2040-6955, 0000-0001-6446-2042}
}
```
