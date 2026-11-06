# From Form to Program & The Double-Stranded Codon Torus ($\mathbb{T}_4^3 = (\mathbb{Z}/4\mathbb{Z})^3$): Conditional Codon Syntax, Markovian Lumpability & Synthetic mRNA Forensics

[![Codon Torus & Syntax CI](https://github.com/Dr-Logvinovich/codon-torus-syntax/actions/workflows/ci.yml/badge.svg)](https://github.com/Dr-Logvinovich/codon-torus-syntax/actions/workflows/ci.yml)
[![ORCID: V. Logvinovich](https://img.shields.io/badge/ORCID-0009--0007--2040--6955-a6ce39.svg)](https://orcid.org/0009-0007-2040-6955)
[![ORCID: J.-C. Perez](https://img.shields.io/badge/ORCID-0000--0001--6446--2042-a6ce39.svg)](https://orcid.org/0000-0001-6446-2042)
[![Offline Cache Included](https://img.shields.io/badge/NCBI_Cache_Included-568_MB_(103_Files)-008080.svg)](cds_cache)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B_%7C_NumPy_%7C_SciPy-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Official Open-Science Verification Repository for the Genomic Information Theory & Double-Stranded Codon Torus Research Programme (Complete with 568 MB Pre-Cached NCBI RefSeq/GenBank Archives & Pre-Computed Audit Reports).**
>
> * **Authors:** **Victor Logvinovich** (Master of Physical and Mathematical Sciences; Researcher in Pedagogical Sciences (Mathematics), Kobrin, Belarus; [ORCID: 0009-0007-2040-6955](https://orcid.org/0009-0007-2040-6955)) and **Jean-Claude Perez** (PhD in Mathematics and Computer Science, Retired IBM AI European Research Centre, Montpellier; Luc Montagnier Foundation; [ORCID: 0000-0001-6446-2042](https://orcid.org/0000-0001-6446-2042)).

---

## Complete Repository Architecture (100% Self-Contained & Offline-Ready)

Unlike standard bioinformatics repositories that require hours of external FTP downloads before a single test can run, this repository includes the **entire 568 MB local NCBI cache (`cds_cache/`, 103 verified genomic files)** and **all pre-computed verification outputs (`results_verified/`)**:

| Path / Artifact | Category | Contents & Verification Scope |
| :--- | :---: | :--- |
| **[`code_logos_perez_v14_8.tex.pdf`](code_logos_perez_v14_8.tex.pdf)** | **Preprint I (34 pp.)** | *"From Form to Program: Medium-Range Conditional Codon Syntax Inside Real Genes — Strict Null-Model Evidence Across 54 Cellular RefSeq Taxa, a Codon-Mirror-Symmetry Control, Two-Scale Phylogenetic Transfer, Intron-Exon Compartment Confinement, Synthetic Genome Forensics, and the Biotech 'Max-CAI' Codon-Optimization Crash"* (`v14.8`). |
| **[`codon_torus_logos_perez_v3.tex.pdf`](codon_torus_logos_perez_v3.tex.pdf)** | **Preprint II (33 pp.)** | *"The Double-Stranded Codon Torus: Duplex Mirror Fibration of DNA Coding Space, Markovian Lumpability, and Information-Theoretic Consequences"* (`v3.5`). |
| **[`super_codon_proof_v14_8.py`](super_codon_proof_v14_8.py)** | **Verification Engine I** | Autonomous self-testing pipeline (`v14.6–v14.8`): 10-level null model hierarchy (`H0`–`H7`) across 54 cellular RefSeq taxa + SARS-CoV-2, 14-gene human RefSeqGene intron-exon compartment confinement (`TP53`, `KRAS`, `BRCA1`, etc.), 13 frontier genomes (Giant Viruses, 2.8 km boreholes, $122^\circ\text{C}$ archaea, Venter's `JCVI-syn1.0` & `JCVI-syn3.0`), viral "Genetic Trojans" (`HPV-16 E7`, `Adenovirus E3 RID-α`), organelle OS rewriting, and the `Max-CAI` crash test. |
| **[`codon_torus_real_audit_v3_bnt162b2.py`](codon_torus_real_audit_v3_bnt162b2.py)** | **Verification Engine II** | Simulation-free, **zero-RNG** analytic audit suite (`v3.2 / v3.5`): verifies exact chain-rule decompositions on $\mathbb{T}_4^3 = (\mathbb{Z}/4\mathbb{Z})^3$ to $< 10^{-15}\text{ nats/codon}$, the two-strand duplex fibration $\alpha_{\text{dup}}: 64 \to 56$, Kemeny–Snell Markov lumpability ($G^2$), antisense program dominance ($74.1\%\text{–}88.5\%$), 6-frame duplex ORF forensics on **Pfizer BNT162b2** (3,825 nt) and **Moderna mRNA-1273** (3,828 nt / `OR134578.1`), and the constructive **Corollary 4.4 Duplex-Safe Viterbi Compiler**. |
| **[`run_all_audits.sh`](run_all_audits.sh)** | **1-Click Master Runner** | Executes all built-in mathematical unit tests, runs both verification engines across all cached datasets, and generates all JSON/Markdown reports and publication figures in one command. |
| **[`cds_cache/`](cds_cache)** | **Local NCBI Cache (568 MB, 103 Files)** | **55** official RefSeq `*_cds_from_genomic.fna.gz` assemblies, **16** human RefSeqGene GenBank loci (`NG_*.gb` + `JCVI-syn1.0/3.0`), **13** Frontier Biosphere CDS archives, **19** viral/chimera/vaccine FASTA records (`NC_045512.2`, `OR134577.1`, `OR134578.1`, `mRNA-1273`), and the 100-kb human chromosome 17 window `NC_000017.11:55154472-55254472`. |
| **[`results_verified/`](results_verified)** | **Verified Outputs & Figures** | Complete JSON, Markdown, and log outputs for all 5 audit stages (`01_codon_torus_and_bnt162b2_v3_5` through `05_biosphere_55taxa_smoke`) plus all 8 publication figures (`PNG` & `PDF`). |

---

## Methodological & Mathematical Protocol (Read First)

1. **Strict Separation of Static "Form" vs. Dynamic Sequential "Program":**
   Zero-order nucleotide frequencies, GC3 isochores, marginal codon usage (`CAI`), and single-stranded Chargaff/Perez 32-mirror-pair balances describe the **static compositional form** of a sequence (the 32-dimensional mirror-even subspace $V_+ \subset \mathbb{R}^{64}$). This programme measures the **dynamic sequential program**—the first-order and medium-range ($L = 30$ codons) conditional transition structure $P(C_i \mid C_{i-1})$ living on the mirror-odd subspace $V_-$ and the synonymous fiber bundle after all compositional, amino-acid, and strand-parity marginals are absorbed.
2. **Zero Monte Carlo Surrogates in the Torus Proofs (`codon_torus_real_audit_v3_bnt162b2.py`):**
   Every theorem in Preprint II (`codon_torus_logos_perez_v3.tex.pdf`) is proved deductively in exact algebra and information theory and audited analytically via exact plug-in conditional mutual information, held-out Dirichlet-smoothed log-likelihoods, and asymptotic/exact contingency-table likelihood-ratio ($G^2$) statistics with closed-form degrees of freedom ($\text{df} = 2709, 2460, 1849, 1681, 860, 779, 400, 361$). Zero random number generators are used (`zero_random_generators_used = True`).
3. **Every Integer Is Derived (Anti-Numerology Policy):**
   Every integer appearing in the torus framework ($2, 3, 4, 6, 7, 8, 16, 19, 20, 21, 24, 32, 40, 41, 42, 43, 48, 54, 56, 60, 61, 62, 63, 64, 361, 400, 779, 860, 1681, 1849, 2460, 2709$) is derived explicitly from a group order, free orbit count, simplex dimension $\sum_a (k_a - 1)$, or contingency-table degree-of-freedom formula.

---

## Part I: The Double-Stranded Codon Torus $\mathbb{T}_4^3 = (\mathbb{Z}/4\mathbb{Z})^3$ & Exact Decompositions (`v3.2–v3.5`)

### 1. Duplex Geometry, Mirror Involution $\rho$, and the Two-Strand Fibration $\alpha_{\text{dup}}: 64 \to 56$
* **Discrete Torus & Fixed-Point-Free Involution (Theorem 2.5 & Corollary 2.6):** Under $e: \{\text{A,C,G,T}\} \to \mathbb{Z}/4\mathbb{Z}$, the 64 codons form the discrete torus $\mathbb{T}_4^3 = (\mathbb{Z}/4\mathbb{Z})^3$. Watson–Crick reverse complementation $\rho(c) = -Rc + 3\mathbf{1} \pmod 4$ is a fixed-point-free affine isometry in both Hamming ($d_H$) and Lee ($d_{\text{Lee}}$) metrics, partitioning $\mathbb{T}_4^3$ into **32 mirror orbits** and splitting $\mathbb{R}^{64} = V_+ \oplus V_-$ ($\dim V_+ = \dim V_- = 32$).
* **The 8-Doublet Theorem (Theorem 4.3):** While the single-strand genetic code $\alpha: \mathbb{T}_4^3 \to \mathcal{A}$ ($64 \to 21$) leaves a **43-dimensional** synonymous simplex bundle ($41$ on sense codons), evaluating the code on both strands via $\alpha_{\text{dup}}(c) = (\alpha(c), \alpha(\rho(c)))$ maps the sense wobble base onto the antisense leading base, partitioning $\mathbb{T}_4^3$ into **48 singletons** and **8 binary doublets** (`{TAA,TAG}`, `{CAA,CAG}`, `{TCG,TCT}`, `{CCG,CCT}` and their $\rho$-duals `{TTA,CTA}`, `{TTG,CTG}`, `{CGA,AGA}`, `{CGG,AGG}`) across **56 joint channels**, collapsing the duplex synonymous dimension from $43$ to $64 - 56 = \mathbf{8}$.

### 2. Exact Tripartite, Quadripartite, and Two-Strand Decompositions Across 9 RefSeq Genomes (9,834,024 Held-Out Transitions)

For any joint distribution of consecutive codons $(C_{i-1}, C_i)$ with $A_i = \alpha(C_i)$ and $A_i^\rho = \alpha(\rho(C_i))$, total bigram mutual information decomposes with **zero stationarity assumptions** (Theorems 5.3, 5.4, 5.5):

$$I(C_{i-1}; C_i) = \underbrace{I(A_{i-1}; A_i)}_{\text{Protein dipeptide syntax}} + \underbrace{I(C_{i-1}; A_i \mid A_{i-1})}_{d_{2b}\text{ (Codon}\to\text{AA context)}} + \underbrace{I(A_{i-1}; C_i \mid A_i) + I(C_{i-1}; C_i \mid A_{i-1}, A_i)}_{S_{\text{syn}} = I(C_{i-1}; C_i \mid A_i) = S_{\text{anti}} + S_{\text{dup}}}$$

![Fig 1. Exact Tripartite Decomposition](results_verified/01_codon_torus_and_bnt162b2_v3_5/figures/fig1_tripartite.png)
![Fig 2. Exact Quadripartite Decomposition](results_verified/01_codon_torus_and_bnt162b2_v3_5/figures/fig2_quadripartite.png)

| Species (Held-Out Test Split) | Table | Test Genes | Test Trans | $I(C;C)$ | $I(A;A)$ | $d_{2b}$ | $S_{\text{syn}}$ | $S_{\text{syn}}^{(2)}$ | $S_{\text{anti}}$ (Antisense %) | Chain Residual | $\text{df}(G^2 / G_{(2)}^2 / \text{lump})$ | Gene $t$-stat (Pos %) |
| :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: | ---: |
| **`Human_T2T`** (`GCF_009914755.1`) | 1 | 3,999 | 2,342,893 | `+0.0842` | `+0.0099` | `+0.0483` | `+0.0260` | `+0.0222` | `+0.0215` (**82.5%**) | `4.2e-17` | `2460 / 1681 / 779` | `+52.58` (83.0%) |
| **`M_musculus`** (`GCF_000001635.27`) | 1 | 4,509 | 2,411,478 | `+0.0770` | `+0.0091` | `+0.0543` | `+0.0136` | `+0.0104` | `+0.0110` (**81.1%**) | `2.2e-16` | `2460 / 1681 / 779` | `+49.46` (84.5%) |
| **`D_melanogaster`** (`GCF_000001215.4`) | 1 | 2,698 | 1,493,418 | `+0.0632` | `+0.0125` | `+0.0343` | `+0.0164` | `+0.0103` | `+0.0134` (**81.8%**) | `1.9e-16` | `2460 / 1681 / 779` | `+54.63` (87.7%) |
| **`A_thaliana`** (`GCF_000001735.4`) | 1 | 5,426 | 2,186,312 | `+0.0458` | `+0.0082` | `+0.0253` | `+0.0123` | `+0.0083` | `+0.0099` (**80.4%**) | `2.8e-17` | `2460 / 1681 / 779` | `+64.57` (85.1%) |
| **`S_cerevisiae`** (`GCF_000146045.2`) | 1 | 1,171 | 563,029 | `+0.0254` | `+0.0065` | `+0.0111` | `+0.0078` | `+0.0053` | `+0.0062` (**79.0%**) | `3.5e-17` | `2460 / 1681 / 779` | `+17.43` (73.6%) |
| **`E_coli`** (`GCF_000005845.2`) | 11 | 875 | 272,548 | `+0.0943` | `+0.0103` | `+0.0475` | `+0.0366` | `+0.0204` | `+0.0304` (**83.2%**) | `2.8e-17` | `2460 / 1680 / 779` | `+33.85` (90.1%) |
| **`B_subtilis`** (`GCF_000009045.1`) | 11 | 897 | 251,034 | `+0.0811` | `+0.0102` | `+0.0439` | `+0.0269` | `+0.0179` | `+0.0200` (**74.1%**) | `1.9e-16` | `2460 / 1681 / 779` | `+27.98` (86.7%) |
| **`M_jannaschii`** (`GCF_000007745.1`) | 11 | 439 | 147,216 | `+0.1152` | `+0.0101` | `+0.0660` | `+0.0392` | `+0.0255` | `+0.0338` (**86.4%**) | `3.1e-16` | `2456 / 1676 / 779` | `+30.51` (91.8%) |
| **`H_salinarum`** (`GCF_000015585.1`) | 11 | 480 | 166,096 | `+0.0721` | `+0.0085` | `+0.0387` | `+0.0249` | `+0.0125` | `+0.0221` (**88.5%**) | `6.9e-17` | `2455 / 1671 / 779` | `+17.09` (85.4%) |

![Fig 3. Deterministic Max-CAI Collapse & KL Cross-Scoring Crash](results_verified/01_codon_torus_and_bnt162b2_v3_5/figures/fig3_maxcai_collapse.png)
![Fig 4. Non-Equilibrium Dynamics and Mirror Symmetry Breaking](results_verified/01_codon_torus_and_bnt162b2_v3_5/figures/fig4_dynamics.png)
![Fig 5. Cluster-Robust Gene-Level Evidence on Held-Out Test Genes](results_verified/01_codon_torus_and_bnt162b2_v3_5/figures/fig5_per_gene.png)
![Fig 6. Empirical Hourglass Profile on 17q22](results_verified/01_codon_torus_and_bnt162b2_v3_5/figures/fig6_hourglass_17q22.png)

---

## Part II: Six-Frame Duplex Forensics on Pfizer `BNT162b2` & Moderna `mRNA-1273` vs. Corollary 4.4 Duplex-Safe Viterbi Compiler (`v3.5`)

By **Corollary 4.4**, under $\rho(c) = -Rc + 3\mathbf{1} \pmod 4$, the sense preimages of the termination codons $\text{Stop} = \{\text{TAA, TAG, TGA}\}$ in antiparallel antisense frame `-0` are exclusively the three wobble-A codons $\rho^{-1}(\text{Stop}) = F_{\text{Leu,Stop}} \cup F_{\text{Ser,Stop}} = \{\text{TTA, CTA, TCA}\}$, while the sense preimage of the canonical antisense start codon `ATG` is the wobble-T Histidine singleton $\rho^{-1}(\{\text{ATG}\}) = F_{\text{His,Met}} = \{\text{CAT}\}$.

![Fig 8. BNT162b2 & Moderna mRNA-1273 Duplex Torus Forensics vs. Duplex-Safe Viterbi Compiler](results_verified/01_codon_torus_and_bnt162b2_v3_5/figures/fig8_bnt162b2_ambigram_duplex.png)

| Spike-2P Construct (`SHA-256` Verified) | Codons (Len nt) | GC3 | CpG / 100nt | $H_{\text{cod}}$ (bits) | Cross $d_{2b}$ (`Human_T2T`) | Cross $S_{\text{syn}}$ | $\sigma_{\text{dup}}$ (nats) | Sense `{TTA,CTA,TCA}` (`-0` Stops) | Sense `CAT` (`-0 ATG`) | Sense Stops (`+0/+1/+2`) | Antisense Stops (`-0/-1/-2`) | Max Unbroken Antisense Span | Frame `-0` `ATG->End` ORF |
| :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1. `WT_SARS2_Spike_NC045512`** (`sha256=7d5621cd3b3e`) | 1,274 (3,822) | `0.2665` | `0.76` | `5.5173` | **`+0.0222`** | `+0.0619` | `17.8925` | `28/9/26` (**63**) | `13` | `1 / 55 / 145` | `63 / 64 / 134` | `106c` (`-1`) | `79c` (`650..728`) |
| **2. `BNT162b2_Full_Spike_2P`** (Pfizer 3,825 nt, `sha256=4a4a21a52bba`) | 1,275 (3,825) | `0.8105` | `5.52` | `4.9006` | **`-0.0549`** | `+0.0300` | `21.2378` | `0/0/0` (**0**) | **`2`** | `2 / 59 / 20` | **`0`** `/ 43 / 72` | **`1,275c` (`-0`)** | **`1,101c` (`175..1275`)** |
| **3. `mRNA1273_Full_Spike_2P`** (Moderna 3,828 nt, `sha256=e895ebf67081`) | 1,276 (3,828) | `0.9552` | `7.21` | `4.4700` | **`-0.0616`** | `+0.1154` | `25.0785` | `1/1/2` (**4**) | **`0`** | `3 / 76 /` **`2`** | `4 / 36 / 54` | `458c` (`-0`) | **`0c` (`0..0`)** |
| **4. `mRNA1273_NCBI_OR134578`** (NCBI Vector 3,819 nt, `sha256=b249652999cb`) | 1,273 (3,819) | `0.9425` | `7.18` | `4.5336` | **`-0.0635`** | `+0.1059` | `24.3469` | `1/1/2` (**4**) | **`0`** | `2 / 76 /` **`5`** | `4 / 36 / 54` | `458c` (`-0`) | **`0c` (`0..0`)** |
| **5. `Human_MaxCAI_Spike_2P`** (Theorem 6.1 Deterministic $\phi$) | 1,274 (3,822) | `0.9670` | `6.60` | `4.1578` | **`-0.0647`** | `+0.1213` | `26.4177` | `0/0/0` (**0**) | `0` | `1 / 79 /` **`0`** | **`0`** `/ 32 / 61` | `1,274c` (`-0`) | `0c` (`0..0`) |
| **6. `Human_MaxGC3_Spike_2P`** (Deterministic $\text{GC3}=1.0000$) | 1,274 (3,822) | `1.0000` | `6.60` | `4.1578` | **`-0.0666`** | `+0.1339` | `26.4180` | `0/0/0` (**0**) | `0` | `1 / 79 /` **`0`** | **`0`** `/ 32 / 61` | `1,274c` (`-0`) | `0c` (`0..0`) |
| **7. `Human_Markov_H2b_Spike_2P`** (Unconstrained Viterbi, Prop. 6.3) | 1,274 (3,822) | `0.6588` | `3.25` | `5.6133` | **`+0.2194`** | `+0.1553` | `19.7486` | `22/0/5` (**27**) | `9` | `1 / 41 / 40` | `27 / 54 / 86` | `197c` (`-0`) | `103c` (`620..722`) |
| **8. `Human_DuplexSafe_Viterbi_Spike_2P`** (**Corollary 4.4 Compiler**) | 1,274 (3,822) | `0.6769` | `3.51` | `5.5910` | **`+0.2142`** | **`+0.1557`** | **`19.7983`** | `23/1/6` (**30**) | **`0`** | `1 / 41 /` **`39`** | **`30`** `/ 47 / 87` | **`174c` (`-0`)** | **`0c` (`0..0`)** |

---

## Part III: Multi-Scale Biosphere, Intron–Exon Splicing, Synthetic Cells & Viral Trojans (`v14.5–v14.8`)

1. **54 Cellular RefSeq Taxa Across the Biosphere ($0.1416 \le \text{GC3} \le 0.9177$):** All 52 guard-passing cellular assemblies pass all 7 strict bootstrap criteria (`d2b`, `d3`, `d3b`, `ΔH4(30)`, `Δ5`, `Δ6`, `Δ7` $> 0$ at 95% CI), and all 55 cached assemblies are verified in [`results_verified/05_biosphere_55taxa_smoke/report_smoke.md`](results_verified/05_biosphere_55taxa_smoke/report_smoke.md).
2. **14/14 Human RefSeqGene Intron–Exon Confinement & `TP53` / `RAS` Oncogenic Hotspots (`v14.6`):** Full report in [`results_verified/02_intron_exon_p53_ras_v14_6/p53_intron_exon_audit.md`](results_verified/02_intron_exon_p53_ras_v14_6/p53_intron_exon_audit.md). Across all 14 human genes, concatenated coding exons have $d_{3b} > 0$ (mean $+0.1721\text{ nats/codon}$, $0.00$ internal stops/100 codons) while excised introns ($-0.3097$, contrast $+0.48181\text{ nats/codon}$) and whole unspliced pre-mRNAs across all 3 reading frames (contrast $+0.41616\text{ nats/codon}$) are strictly negative ($14/14$, $p = 2^{-14} = 6.1 \times 10^{-5}$).
3. **13-Genome Frontier Biosphere & Venter's `JCVI-syn1.0` / `JCVI-syn3.0` Forensics (`v14.7`):** Full report in [`results_verified/03_frontier_and_jcvi_syn_v14_7/frontier_audit.md`](results_verified/03_frontier_and_jcvi_syn_v14_7/frontier_audit.md). All 4 giant DNA viruses (*Mimivirus*, *Megavirus*, *Pandoravirus*, *Pithovirus*), all 4 deep-biosphere extremophiles (*Desulforudis audaxviator* at 2.8 km depth, *Methanopyrus kandleri* at $122^\circ\text{C}$, *Aquifex aeolicus*, *Pyrococcus furiosus*), and the Asgard archaeon *Prometheoarchaeum syntrophicum* MK-D1 pass `STRICT_7/7_PASS`. Blind cross-scoring under wild-type *Mycoplasma mycoides* isolates the foreign engineered vector genes (`MMSYN1_0915 lacZ` $d_{3b} = -0.4564$, `MMSYN1_0918 HIS3` $-0.4286$, `MMSYN1_0914 bla` $-0.2353$) in `JCVI-syn1.0` and confirms $0/412$ disrupted essential genes in `JCVI-syn3.0` (while all 10 intergenic watermark cassettes $\ge 600\text{ nt}$ collapse to $d_{3b} \in [-0.4559, -0.1964]$).
4. **Viral "Genetic Trojans", Organelle OS Rewriting & Biotech `Max-CAI` Crash Test (`v14.8`):** Full report in [`results_verified/04_viral_trojan_and_maxcai_crash_v14_8/chimera_hgt_trojan_audit.md`](results_verified/04_viral_trojan_and_maxcai_crash_v14_8/chimera_hgt_trojan_audit.md). Standard industrial `Max-CAI` optimization crashes $d_{2b}$ from $+0.0427$ to $-0.0589\text{ nats/codon}$ and drops 64-codon Shannon entropy from $5.322$ to $4.091\text{ bits}$ across $14/14$ human genes, whereas first-order `Markov_H2b` compilation restores $d_{3b} = +0.1129\text{ nats/codon}$. Blind viral scoring on `Human_T2T` isolates **HPV-16 `E7`** ($d_{2b} = +0.0327$, sole positive gene in HPV-16) and **Adenovirus `E3 RID-α`** ($d_{2b} = +0.0958$).

---

## One-Click Reproduction (Zero Downloads Required — Uses Pre-Populated `cds_cache/`)

Clone the repository and run the master verification script `run_all_audits.sh`:

```bash
git clone [https://github.com/Dr-Logvinovich/codon-torus-syntax.git](https://github.com/Dr-Logvinovich/codon-torus-syntax.git)
cd codon-torus-syntax
bash run_all_audits.sh
```

Or run individual modules directly:

```bash
# 1. Built-in unit tests (8 NCBI genetic code tables, 32 mirror orbits, 56 duplex channels, chain-rule < 1e-15):
python3 super_codon_proof_v14_8.py --self-test
python3 -c "import codon_torus_real_audit_v3_bnt162b2 as ct; print(ct.self_test_torus_theorems())"

# 2. Double-Stranded Codon Torus Audit (9 RefSeq genomes + 17q22 locus + BNT162b2 & mRNA-1273 forensics):
python3 codon_torus_real_audit_v3_bnt162b2.py --cds-dir cds_cache --out-dir results_verified/01_codon_torus_and_bnt162b2_v3_5 --bnt162b2-audit

# 3. 14-Gene Human RefSeqGene Intron-Exon, Splicing & TP53/RAS Oncogenic Hotspot Audit (v14.6):
python3 super_codon_proof_v14_8.py --intron-exon-audit --intron-exon-genes TP53,TP63,TP73,BRCA1,PTEN,RB1,GAPDH,ACTB,HBB,KRAS,NRAS,HRAS,MYC,APOE --cds-dir cds_cache --out-dir results_verified/02_intron_exon_p53_ras_v14_6

# 4. 13-Genome Frontier Biosphere & Venter JCVI-syn1.0/3.0 Forensics Audit (v14.7):
python3 super_codon_proof_v14_8.py --frontier-audit --cds-dir cds_cache --out-dir results_verified/03_frontier_and_jcvi_syn_v14_7

# 5. Viral "Genetic Trojan" Scanner, Organelle OS Rewriting & Biotech Max-CAI Crash Test (v14.8):
python3 super_codon_proof_v14_8.py --chimera-audit --cds-dir cds_cache --out-dir results_verified/04_viral_trojan_and_maxcai_crash_v14_8
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
  url          = {[https://github.com/Dr-Logvinovich/codon-torus-syntax](https://github.com/Dr-Logvinovich/codon-torus-syntax)},
  note         = {GitHub Repository: [https://github.com/Dr-Logvinovich/codon-torus-syntax](https://github.com/Dr-Logvinovich/codon-torus-syntax), ORCID: 0009-0007-2040-6955, 0000-0001-6446-2042}
}

@misc{LogvinovichPerez2026CodonTorus,
  author       = {Logvinovich, Victor and Perez, Jean-Claude},
  title        = {{The Double-Stranded Codon Torus: Duplex Mirror Fibration of DNA Coding Space, Markovian Lumpability, and Information-Theoretic Consequences}},
  year         = {2026},
  month        = oct,
  version      = {v3.5},
  url          = {[https://github.com/Dr-Logvinovich/codon-torus-syntax](https://github.com/Dr-Logvinovich/codon-torus-syntax)},
  note         = {GitHub Repository: [https://github.com/Dr-Logvinovich/codon-torus-syntax](https://github.com/Dr-Logvinovich/codon-torus-syntax), ORCID: 0009-0007-2040-6955, 0000-0001-6446-2042}
}
```

## License

This project is licensed under the **MIT License** (code & scripts) and **CC BY 4.0** (manuscripts & figures) — see the [`LICENSE`](LICENSE) file for details.
