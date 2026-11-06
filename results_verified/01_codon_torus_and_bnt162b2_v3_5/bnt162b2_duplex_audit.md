# BNT162b2 & Moderna mRNA-1273 Duplex Torus, Ribosomal Pausing & Ambigrammatic ORF Forensics (v3.5)

## 1. Cross-Scoring on the `Human_T2T` Reference Model and Two-Strand Torus Decomposition

| Construct | Codons | GC3 | CpG/100nt | H_codon (bits) | Cross d2b | Cross d3b | Cross S_syn | Internal S_syn | sigma_dup | Slippery (UUUU) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **WT_SARS2_Spike_NC045512** | 1274 | 0.2665 | 0.76 | 5.5173 | +0.0222 | +0.1029 | +0.0619 | +0.5134 | 17.8925 | 67 (39) |
| **BNT162b2_Full_Spike_2P** | 1275 | 0.8105 | 5.52 | 4.9006 | -0.0549 | +0.2015 | +0.0300 | +0.2250 | 21.2378 | 27 (1) |
| **mRNA1273_Full_Spike_2P** | 1276 | 0.9552 | 7.21 | 4.4700 | -0.0616 | +0.3175 | +0.1154 | +0.1263 | 25.0785 | 34 (0) |
| **mRNA1273_NCBI_OR134578** | 1273 | 0.9425 | 7.18 | 4.5336 | -0.0635 | +0.3064 | +0.1059 | +0.1468 | 24.3469 | 34 (0) |
| **Human_MaxCAI_Spike_2P** | 1274 | 0.9670 | 6.60 | 4.1578 | -0.0647 | +0.3055 | +0.1213 | +0.0000 | 26.4177 | 37 (0) |
| **Human_MaxGC3_Spike_2P** | 1274 | 1.0000 | 6.60 | 4.1578 | -0.0666 | +0.3248 | +0.1339 | +0.0000 | 26.4180 | 40 (0) |
| **Human_Markov_H2b_Spike_2P** | 1274 | 0.6588 | 3.25 | 5.6133 | +0.2194 | +0.3250 | +0.1553 | +0.6342 | 19.7486 | 50 (8) |
| **Human_DuplexSafe_Viterbi_Spike_2P** | 1274 | 0.6769 | 3.51 | 5.5910 | +0.2142 | +0.3304 | +0.1557 | +0.6190 | 19.7983 | 49 (8) |

## 2. Corollary 4.4 Ambigrammatic Overlapping ORF Forensics (`{TTA, CTA, TCA}` vs. `CAT` Antisense `ATG`)

| Construct | Sense `{TTA, CTA, TCA}` | Sense `CAT` (`-0 ATG`) | Sense Stops (`+0/+1/+2`) | Antisense Stops (`-0/-1/-2`) | Max Unbroken Antisense Span | Frame `-0` `ATG->End` ORF |
|---|---:|---:|---|---|---:|---:|
| **WT_SARS2_Spike_NC045512** | `28/9/26 (63)` | `13` | `1/55/145` | `63/64/134` | `106 cod (-1)` | `79 cod (650..728)` |
| **BNT162b2_Full_Spike_2P** | `0/0/0 (0)` | `2` | `2/59/20` | `0/43/72` | `1275 cod (-0)` | `1101 cod (175..1275)` |
| **mRNA1273_Full_Spike_2P** | `1/1/2 (4)` | `0` | `3/76/2` | `4/36/54` | `458 cod (-0)` | `0 cod (0..0)` |
| **mRNA1273_NCBI_OR134578** | `1/1/2 (4)` | `0` | `2/76/5` | `4/36/54` | `458 cod (-0)` | `0 cod (0..0)` |
| **Human_MaxCAI_Spike_2P** | `0/0/0 (0)` | `0` | `1/79/0` | `0/32/61` | `1274 cod (-0)` | `0 cod (0..0)` |
| **Human_MaxGC3_Spike_2P** | `0/0/0 (0)` | `0` | `1/79/0` | `0/32/61` | `1274 cod (-0)` | `0 cod (0..0)` |
| **Human_Markov_H2b_Spike_2P** | `22/0/5 (27)` | `9` | `1/41/40` | `27/54/86` | `197 cod (-0)` | `103 cod (620..722)` |
| **Human_DuplexSafe_Viterbi_Spike_2P** | `23/1/6 (30)` | `0` | `1/41/39` | `30/47/87` | `174 cod (-0)` | `0 cod (0..0)` |
