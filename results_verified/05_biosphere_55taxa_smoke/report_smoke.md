# Super Codon Conditional Syntax Report — Stage: `smoke` (v14.6)

## Within-Species Conditional Syntax Summary

| Species | Clade | Table | Train/Test | GC3 | d0 | d2b (AA-pair) | d3b (nt-tri) | dH4(L=3) | dH5 (Surr) | dH6 (Syn) | dH7 (Mirror) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Drosophila | insect | 1 | 321/79 | 0.634 | +0.05514 | +0.03261 | +0.16153 | +0.04083 | +0.05996 | +0.09456 | +0.07388 |
| Mouse | mammal | 1 | 320/79 | 0.564 | +0.05934 | +0.04721 | +0.08742 | +0.04824 | +0.08296 | +0.12653 | +0.10967 |
| Arabidopsis | plant | 1 | 320/79 | 0.412 | +0.03051 | +0.02278 | +0.08816 | +0.02608 | +0.03617 | +0.07046 | +0.05375 |
| Zebrafish | vertebrate | 1 | 321/79 | 0.536 | +0.05400 | +0.04427 | +0.09214 | +0.04813 | +0.08028 | +0.12828 | +0.10871 |
| Canis_lupus_familiaris | mammal | 1 | 320/80 | 0.670 | +0.08060 | +0.04128 | +0.14429 | +0.05304 | +0.11004 | +0.12767 | +0.11115 |
| Trypanosoma_brucei | protist | 1 | 384/16 | 0.506 | +0.03972 | +0.01437 | +0.07983 | +0.02606 | +0.02785 | +0.05483 | +0.04528 |
| Chlamydomonas_reinhardtii | plant-like protist | 1 | 320/80 | 0.804 | +0.10569 | +0.04335 | +0.19510 | +0.07420 | +0.06582 | +0.14018 | +0.11320 |
| Leishmania_major | protist | 1 | 397/3 | 0.743 | +0.05114 | +0.03897 | +0.14960 | +0.03386 | +0.03039 | +0.09136 | +0.06722 |
| Plasmodium_falciparum | protist | 1 | 363/37 | 0.164 | +0.03257 | +0.00948 | +0.20663 | +0.02911 | +0.02909 | +0.03586 | +0.04341 |
| Populus_trichocarpa | plant | 1 | 320/80 | 0.386 | +0.03117 | +0.03035 | +0.08347 | +0.03006 | +0.05307 | +0.07750 | +0.06462 |
| Schizosaccharomyces_pombe | fungus | 1 | 320/80 | 0.320 | +0.01236 | +0.00867 | +0.12929 | +0.01454 | +0.01973 | +0.03104 | +0.02396 |
| C.elegans | nematode | 1 | 317/81 | 0.396 | +0.04014 | +0.02553 | +0.07489 | +0.03518 | +0.03533 | +0.08836 | +0.07210 |
| Sus_scrofa | mammal | 1 | 320/80 | 0.659 | +0.07255 | +0.04277 | +0.13057 | +0.04558 | +0.09270 | +0.12199 | +0.10189 |
| Oryctolagus_cuniculus | mammal | 1 | 320/79 | 0.609 | +0.06174 | +0.03454 | +0.10855 | +0.04464 | +0.05546 | +0.11215 | +0.09836 |
| Xenopus_tropicalis | vertebrate | 1 | 320/80 | 0.479 | +0.04165 | +0.03601 | +0.07532 | +0.03956 | +0.07207 | +0.10162 | +0.08513 |
| Glycine_max | plant | 1 | 320/80 | 0.437 | +0.03548 | +0.02634 | +0.09027 | +0.03487 | +0.03653 | +0.08583 | +0.07347 |
| Caenorhabditis_briggsae | nematode | 1 | 316/79 | 0.428 | +0.03677 | +0.02563 | +0.09158 | +0.03527 | +0.03931 | +0.08325 | +0.07269 |
| Dictyostelium_discoideum | protist | 1 | 320/80 | 0.144 | +0.08452 | +0.02283 | +0.28707 | +0.05580 | +0.03745 | +0.08543 | +0.09393 |
| Zea_mays | plant | 1 | 320/80 | 0.605 | +0.06908 | +0.01407 | +0.14166 | +0.04094 | +0.00684 | +0.09504 | +0.06526 |
| Anopheles_gambiae | insect | 1 | 320/80 | 0.701 | +0.07340 | +0.03854 | +0.17857 | +0.05563 | +0.05802 | +0.12949 | +0.09789 |
| E.coli | bacterium | 1 | 312/78 | 0.546 | +0.04876 | +0.03397 | +0.16832 | +0.05447 | +0.03153 | +0.14349 | +0.09777 |
| Toxoplasma_gondii | protist | 1 | 393/7 | 0.655 | +0.04771 | +0.01242 | +0.08240 | +0.03638 | +0.06585 | +0.08179 | +0.08231 |
| Pseudomonas_aeruginosa | bacterium | 1 | 320/80 | 0.880 | +0.08700 | +0.06376 | +0.33443 | +0.07567 | +0.14424 | +0.20582 | +0.13113 |
| Streptococcus_pyogenes | bacterium | 1 | 320/79 | 0.321 | +0.02908 | +0.02564 | +0.14698 | +0.03630 | +0.02129 | +0.08646 | +0.06466 |
| Lactococcus_lactis | bacterium | 1 | 320/80 | 0.249 | +0.03073 | +0.02419 | +0.19318 | +0.03903 | +0.04239 | +0.09718 | +0.06599 |
| Salmonella_enterica | bacterium | 1 | 320/80 | 0.600 | +0.06192 | +0.04596 | +0.18413 | +0.06056 | +0.06881 | +0.15021 | +0.09788 |
| M.jannaschii | archaeon | 1 | 320/80 | 0.416 | +0.07938 | +0.05697 | +0.19876 | +0.07147 | +0.09850 | +0.19139 | +0.13865 |
| Thermus_thermophilus | bacterium | 1 | 320/80 | 0.914 | +0.17378 | +0.11127 | +0.35430 | +0.15583 | +0.26718 | +0.38585 | +0.29196 |
| Helicobacter_pylori | bacterium | 1 | 320/80 | 0.427 | +0.16078 | +0.08505 | +0.25605 | +0.14340 | +0.19354 | +0.37893 | +0.28496 |
| B.subtilis | bacterium | 1 | 321/79 | 0.448 | +0.05177 | +0.03549 | +0.13342 | +0.05270 | +0.05334 | +0.13078 | +0.10244 |
| Clostridioides_difficile | bacterium | 1 | 320/80 | 0.168 | +0.00887 | +0.01324 | +0.19548 | +0.02384 | +0.01531 | +0.05468 | +0.03934 |
| Synechocystis_sp_PCC6803 | cyanobacterium | 1 | 320/80 | 0.483 | +0.08310 | +0.05629 | +0.18800 | +0.07907 | +0.10532 | +0.20677 | +0.16971 |
| Staphylococcus_aureus | bacterium | 1 | 393/7 | 0.198 | +0.00555 | +0.01385 | +0.24141 | +0.01831 | +0.06810 | +0.03682 | +0.01925 |
| H.salinarum | archaeon | 1 | 320/80 | 0.877 | +0.03932 | +0.03168 | +0.28991 | +0.04058 | +0.08120 | +0.10428 | +0.07823 |
| Bacteroides_fragilis | bacterium | 1 | 320/80 | 0.462 | +0.04566 | +0.03378 | +0.14722 | +0.04368 | +0.05307 | +0.11164 | +0.08753 |
| Yeast | fungus | 1 | 320/79 | 0.374 | +0.00953 | +0.00748 | +0.10695 | +0.01367 | +0.00220 | +0.02678 | +0.02437 |
| Aspergillus_nidulans | fungus | 1 | 320/80 | 0.591 | +0.01207 | +0.00958 | +0.08605 | +0.01556 | +0.01582 | +0.03611 | +0.03064 |
| Gorilla_gorilla | mammal | 1 | 320/80 | 0.604 | +0.07357 | +0.04836 | +0.11574 | +0.05009 | +0.09571 | +0.12871 | +0.11512 |
| Felis_catus | mammal | 1 | 320/80 | 0.638 | +0.06347 | +0.03793 | +0.10167 | +0.04203 | +0.07449 | +0.11129 | +0.09485 |
| Neurospora_crassa | fungus | 1 | 320/80 | 0.651 | +0.03181 | +0.01548 | +0.11567 | +0.02575 | +0.02211 | +0.06178 | +0.05293 |
| Candida_albicans | fungus | 1 | 320/80 | 0.266 | +0.03979 | +0.02565 | +0.16052 | +0.03413 | +0.05235 | +0.08000 | +0.06459 |
| Solanum_lycopersicum | plant | 1 | 320/80 | 0.370 | +0.02541 | +0.02148 | +0.08655 | +0.02610 | +0.03178 | +0.06219 | +0.05408 |
| Mycobacterium_tuberculosis | bacterium | 1 | 320/80 | 0.800 | +0.04995 | +0.04308 | +0.17526 | +0.04669 | +0.08856 | +0.11771 | +0.09457 |
| Oryza_sativa | plant | 1 | 320/80 | 0.624 | +0.10566 | +0.01636 | +0.18358 | +0.05313 | +0.04800 | +0.13392 | +0.09031 |
| Bos_taurus | mammal | 1 | 320/80 | 0.592 | +0.06191 | +0.03480 | +0.11683 | +0.04251 | +0.04970 | +0.11087 | +0.09758 |
| Equus_caballus | mammal | 1 | 320/80 | 0.526 | +0.06611 | +0.03734 | +0.08805 | +0.04430 | +0.04381 | +0.11427 | +0.09731 |
| Chimp | mammal | 1 | 320/80 | 0.606 | +0.07321 | +0.05042 | +0.11082 | +0.05232 | +0.09301 | +0.13356 | +0.11728 |
| Apis_mellifera | insect | 1 | 320/80 | 0.368 | +0.05749 | +0.01242 | +0.13506 | +0.02650 | -0.00721 | +0.06124 | +0.05120 |
| Macaca_mulatta | mammal | 1 | 320/80 | 0.604 | +0.07394 | +0.04761 | +0.10762 | +0.05279 | +0.08892 | +0.13570 | +0.11778 |
| SARS-CoV-2 | virus | 1 | 9/2 | 0.285 | -0.22112 | -0.01654 | -0.14077 | +0.03319 | -0.24423 | +0.11229 | +0.07422 |
| Human_T2T | mammal | 1 | 316/78 | 0.602 | +0.06690 | +0.04077 | +0.11148 | +0.04990 | +0.07728 | +0.12807 | +0.11009 |
| Bombyx_mori | insect | 1 | 320/80 | 0.541 | +0.04343 | +0.00947 | +0.11748 | +0.02026 | +0.01978 | +0.04728 | +0.03205 |
| Rattus_norvegicus | mammal | 1 | 320/80 | 0.586 | +0.05880 | +0.04886 | +0.09374 | +0.04847 | +0.08435 | +0.12548 | +0.10633 |
| Chicken | vertebrate | 1 | 320/80 | 0.574 | +0.07612 | +0.04057 | +0.11069 | +0.05019 | +0.07347 | +0.13282 | +0.11633 |
| Takifugu_rubripes | vertebrate | 1 | 320/80 | 0.638 | +0.05514 | +0.03944 | +0.12397 | +0.04656 | +0.06980 | +0.11815 | +0.10287 |

## Cross-Species Transfer Summary

| Pair (Train->Test) | Tier | d0 | d2b | d3b | dH5_transfer | Gap (d2b - dH5_tr) |
|---|---|---:|---:|---:|---:|---:|
| `Human_T2T:Chimp` | strict_clean | +0.07671 | +0.05208 | +0.11598 | +0.11137 | -0.05929 |
| `Chimp:Human_T2T` | strict_clean | +0.06579 | +0.03957 | +0.10834 | +0.10711 | -0.06755 |
| `Human_T2T:Mouse` | strict_clean | +0.05501 | +0.04725 | +0.08364 | +0.09801 | -0.05076 |
| `Mouse:Human_T2T` | strict_clean | +0.06179 | +0.03863 | +0.10818 | +0.11797 | -0.07935 |
| `Chimp:Mouse` | strict_clean | +0.04704 | +0.04471 | +0.07622 | +0.09923 | -0.05452 |
| `Mouse:Chimp` | strict_clean | +0.07236 | +0.04977 | +0.11271 | +0.13257 | -0.08280 |
| `Zebrafish:Human_T2T` | strict_clean | +0.02531 | +0.02398 | +0.07639 | +0.15350 | -0.12952 |
| `Human_T2T:Zebrafish` | strict_clean | +0.02991 | +0.03431 | +0.07208 | +0.11750 | -0.08319 |
| `Chicken:Human_T2T` | strict_clean | +0.05768 | +0.03725 | +0.10196 | +0.10950 | -0.07224 |
| `Human_T2T:Chicken` | strict_clean | +0.06069 | +0.03449 | +0.09617 | +0.09417 | -0.05968 |
| `Zebrafish:Chicken` | strict_clean | +0.03157 | +0.02294 | +0.07183 | +0.13220 | -0.10925 |
| `Chicken:Zebrafish` | strict_clean | +0.03363 | +0.03591 | +0.07649 | +0.11538 | -0.07947 |
| `Yeast:Human_T2T` | low_signal_d2b_le_0p02 | -0.00854 | +0.00194 | -0.00233 | +0.26920 | -0.26726 |
| `Human_T2T:Yeast` | low_signal_d2b_le_0p02 | -0.02878 | -0.03886 | +0.04654 | +0.20041 | -0.23927 |
| `C.elegans:Human_T2T` | low_signal_d2b_le_0p02 | -0.03489 | -0.02334 | -0.00499 | +0.25046 | -0.27380 |
| `Human_T2T:C.elegans` | low_signal_d2b_le_0p02 | -0.03389 | -0.03784 | -0.01496 | +0.22570 | -0.26354 |
| `Drosophila:Human_T2T` | low_signal_d2b_le_0p02 | -0.02625 | -0.00914 | +0.02803 | +0.18573 | -0.19488 |
| `Human_T2T:Drosophila` | low_signal_d2b_le_0p02 | -0.04368 | -0.02165 | +0.07892 | +0.13825 | -0.15990 |
| `Arabidopsis:Human_T2T` | low_signal_d2b_le_0p02 | -0.03838 | -0.00853 | -0.00363 | +0.25441 | -0.26294 |
| `Human_T2T:Arabidopsis` | low_signal_d2b_le_0p02 | -0.02206 | -0.02138 | +0.01740 | +0.19754 | -0.21892 |
| `E.coli:Human_T2T` | low_signal_d2b_le_0p02 | -0.14767 | -0.05440 | -0.12036 | +0.27060 | -0.32501 |
| `Human_T2T:E.coli` | low_signal_d2b_le_0p02 | -0.12009 | -0.06730 | +0.01725 | +0.37315 | -0.44046 |
| `B.subtilis:Human_T2T` | low_signal_d2b_le_0p02 | -0.09611 | -0.03764 | -0.05944 | +0.27908 | -0.31672 |
| `Human_T2T:B.subtilis` | low_signal_d2b_le_0p02 | -0.06607 | -0.04085 | +0.01449 | +0.24649 | -0.28734 |
| `M.jannaschii:Human_T2T` | low_signal_d2b_le_0p02 | -0.25825 | -0.12880 | -0.22586 | +0.35838 | -0.48718 |
| `Human_T2T:M.jannaschii` | low_signal_d2b_le_0p02 | -0.13208 | -0.08605 | -0.00884 | +0.60243 | -0.68849 |
| `H.salinarum:Human_T2T` | low_signal_d2b_le_0p02 | -0.10570 | -0.05583 | -0.13126 | +0.04988 | -0.10571 |
| `Human_T2T:H.salinarum` | low_signal_d2b_le_0p02 | -0.05815 | -0.09988 | +0.19558 | +0.82515 | -0.92503 |
| `B.subtilis:E.coli` | low_signal_d2b_le_0p02 | -0.05967 | -0.02755 | +0.04967 | +0.30312 | -0.33067 |
| `E.coli:B.subtilis` | low_signal_d2b_le_0p02 | -0.07052 | -0.02535 | -0.01233 | +0.14123 | -0.16658 |
| `M.jannaschii:E.coli` | low_signal_d2b_le_0p02 | -0.04212 | -0.02400 | +0.07918 | +0.16826 | -0.19227 |
| `E.coli:M.jannaschii` | low_signal_d2b_le_0p02 | -0.00483 | +0.00417 | +0.10710 | +0.28619 | -0.28202 |
| `H.salinarum:M.jannaschii` | low_signal_d2b_le_0p02 | -0.07274 | -0.06355 | +0.00551 | +0.63703 | -0.70058 |
| `M.jannaschii:H.salinarum` | low_signal_d2b_le_0p02 | -0.15161 | -0.01764 | -0.09308 | +1.08739 | -1.10503 |
| `Drosophila:Yeast` | low_signal_d2b_le_0p02 | -0.01232 | -0.02607 | +0.02022 | +0.22887 | -0.25494 |
| `Yeast:Drosophila` | low_signal_d2b_le_0p02 | +0.00617 | -0.00050 | +0.04095 | +0.27398 | -0.27448 |
| `Arabidopsis:Drosophila` | low_signal_d2b_le_0p02 | -0.02697 | -0.00799 | +0.02094 | +0.26217 | -0.27016 |
| `Drosophila:Arabidopsis` | low_signal_d2b_le_0p02 | -0.01101 | -0.01135 | -0.01983 | +0.23103 | -0.24237 |
| `C.elegans:Drosophila` | low_signal_d2b_le_0p02 | -0.00787 | -0.00946 | +0.02641 | +0.23047 | -0.23993 |
| `Drosophila:C.elegans` | low_signal_d2b_le_0p02 | +0.00563 | -0.00370 | -0.03115 | +0.22922 | -0.23292 |
| `SARS-CoV-2:Human_T2T` | small_sample | -0.37941 | -0.10682 | -0.42047 | +0.25244 | -0.35925 |
| `Human_T2T:SARS-CoV-2` | small_sample | +0.05185 | +0.01200 | +0.06884 | +0.61658 | -0.60458 |
| `SARS-CoV-2:Chimp` | small_sample | -0.36823 | -0.09647 | -0.41456 | +0.24203 | -0.33850 |
| `Chimp:SARS-CoV-2` | small_sample | +0.04620 | +0.01457 | +0.07505 | +0.59079 | -0.57622 |
| `SARS-CoV-2:Mouse` | small_sample | -0.36292 | -0.08089 | -0.40809 | +0.23830 | -0.31919 |
| `Mouse:SARS-CoV-2` | small_sample | +0.02449 | +0.01428 | +0.04827 | +0.59480 | -0.58052 |
| `SARS-CoV-2:Yeast` | small_sample | -0.35013 | -0.11424 | -0.25858 | +0.00685 | -0.12110 |
| `Yeast:SARS-CoV-2` | small_sample | -0.00733 | +0.00980 | +0.06847 | +0.27883 | -0.26902 |
| `Human_T2T:Rattus_norvegicus` | strict_clean | +0.05394 | +0.04700 | +0.08971 | +0.09973 | -0.05272 |
| `Rattus_norvegicus:Human_T2T` | strict_clean | +0.05819 | +0.03714 | +0.10433 | +0.11625 | -0.07911 |
| `Human_T2T:Macaca_mulatta` | strict_clean | +0.07456 | +0.04729 | +0.10818 | +0.10842 | -0.06113 |
| `Macaca_mulatta:Human_T2T` | strict_clean | +0.06549 | +0.03957 | +0.10974 | +0.09722 | -0.05765 |
| `Human_T2T:Gorilla_gorilla` | strict_clean | +0.07254 | +0.04846 | +0.11426 | +0.11169 | -0.06323 |
| `Gorilla_gorilla:Human_T2T` | strict_clean | +0.06792 | +0.04064 | +0.11243 | +0.09969 | -0.05905 |
| `Macaca_mulatta:Rattus_norvegicus` | strict_clean | +0.05416 | +0.04697 | +0.09107 | +0.09915 | -0.05217 |
| `Rattus_norvegicus:Macaca_mulatta` | strict_clean | +0.06781 | +0.04595 | +0.10231 | +0.13099 | -0.08503 |
| `Chimp:Rattus_norvegicus` | strict_clean | +0.04922 | +0.04518 | +0.08534 | +0.10536 | -0.06018 |
| `Rattus_norvegicus:Chimp` | strict_clean | +0.07164 | +0.05076 | +0.11122 | +0.13263 | -0.08187 |
| `Mouse:Rattus_norvegicus` | strict_clean | +0.05958 | +0.04793 | +0.09513 | +0.10694 | -0.05901 |
| `Rattus_norvegicus:Mouse` | strict_clean | +0.05897 | +0.04774 | +0.08671 | +0.10219 | -0.05445 |
| `Human_T2T:Bos_taurus` | strict_clean | +0.06018 | +0.03409 | +0.11465 | +0.07603 | -0.04193 |
| `Bos_taurus:Human_T2T` | strict_clean | +0.06618 | +0.04035 | +0.11032 | +0.09271 | -0.05236 |
| `Human_T2T:Sus_scrofa` | strict_clean | +0.07256 | +0.04243 | +0.12936 | +0.12807 | -0.08564 |
| `Sus_scrofa:Human_T2T` | strict_clean | +0.06938 | +0.04066 | +0.11241 | +0.08496 | -0.04430 |
| `Bos_taurus:Sus_scrofa` | strict_clean | +0.07530 | +0.04267 | +0.13260 | +0.12264 | -0.07998 |
| `Sus_scrofa:Bos_taurus` | strict_clean | +0.06470 | +0.03514 | +0.11883 | +0.06072 | -0.02558 |
| `Human_T2T:Canis_lupus_familiaris` | strict_clean | +0.07695 | +0.04158 | +0.14200 | +0.12786 | -0.08628 |
| `Canis_lupus_familiaris:Human_T2T` | strict_clean | +0.06828 | +0.04038 | +0.11224 | +0.10225 | -0.06187 |
| `Human_T2T:Felis_catus` | strict_clean | +0.06391 | +0.03697 | +0.10310 | +0.09895 | -0.06198 |
| `Felis_catus:Human_T2T` | strict_clean | +0.06697 | +0.04030 | +0.11056 | +0.09449 | -0.05418 |
| `Canis_lupus_familiaris:Felis_catus` | strict_clean | +0.06637 | +0.03626 | +0.10480 | +0.09634 | -0.06008 |
| `Felis_catus:Canis_lupus_familiaris` | strict_clean | +0.07744 | +0.04197 | +0.14420 | +0.12452 | -0.08255 |
| `Human_T2T:Equus_caballus` | strict_clean | +0.06193 | +0.03594 | +0.08313 | +0.05460 | -0.01866 |
| `Equus_caballus:Human_T2T` | strict_clean | +0.06939 | +0.04096 | +0.11334 | +0.09852 | -0.05756 |
| `Human_T2T:Oryctolagus_cuniculus` | strict_clean | +0.05737 | +0.03208 | +0.10603 | +0.08870 | -0.05663 |
| `Oryctolagus_cuniculus:Human_T2T` | strict_clean | +0.06603 | +0.04023 | +0.10760 | +0.09127 | -0.05103 |
| `Zebrafish:Xenopus_tropicalis` | strict_clean | +0.01564 | +0.02112 | +0.05103 | +0.10747 | -0.08634 |
| `Xenopus_tropicalis:Zebrafish` | strict_clean | +0.04153 | +0.03484 | +0.08078 | +0.15168 | -0.11684 |
| `Human_T2T:Xenopus_tropicalis` | strict_clean | +0.03097 | +0.02740 | +0.05789 | +0.10477 | -0.07737 |
| `Xenopus_tropicalis:Human_T2T` | strict_clean | +0.05207 | +0.03199 | +0.08746 | +0.17306 | -0.14107 |
| `Chicken:Xenopus_tropicalis` | strict_clean | +0.03280 | +0.02778 | +0.06228 | +0.10151 | -0.07373 |
| `Xenopus_tropicalis:Chicken` | strict_clean | +0.05639 | +0.03126 | +0.08759 | +0.15672 | -0.12546 |
| `Zebrafish:Takifugu_rubripes` | strict_clean | +0.04114 | +0.03143 | +0.10460 | +0.14742 | -0.11599 |
| `Takifugu_rubripes:Zebrafish` | strict_clean | +0.05328 | +0.03871 | +0.08447 | +0.09836 | -0.05965 |
| `Human_T2T:Takifugu_rubripes` | strict_clean | +0.02715 | +0.02566 | +0.09924 | +0.12931 | -0.10365 |
| `Takifugu_rubripes:Human_T2T` | strict_clean | +0.04329 | +0.02769 | +0.09154 | +0.11293 | -0.08524 |
| `Drosophila:Anopheles_gambiae` | low_signal_d2b_le_0p02 | +0.00274 | -0.00136 | +0.10790 | +0.23394 | -0.23530 |
| `Anopheles_gambiae:Drosophila` | low_signal_d2b_le_0p02 | -0.00967 | -0.00830 | +0.09908 | +0.06682 | -0.07511 |
| `Drosophila:Apis_mellifera` | low_signal_d2b_le_0p02 | +0.00481 | -0.02882 | +0.01498 | +0.25868 | -0.28751 |
| `Apis_mellifera:Drosophila` | low_signal_d2b_le_0p02 | +0.04820 | -0.00487 | +0.03757 | +0.26554 | -0.27041 |
| `Drosophila:Bombyx_mori` | low_signal_d2b_le_0p02 | -0.02874 | -0.02686 | +0.03414 | +0.06696 | -0.09382 |
| `Bombyx_mori:Drosophila` | low_signal_d2b_le_0p02 | +0.00392 | -0.00108 | +0.09170 | +0.20857 | -0.20965 |
| `Anopheles_gambiae:Apis_mellifera` | low_signal_d2b_le_0p02 | +0.03572 | -0.04653 | +0.01032 | +0.25817 | -0.30470 |
| `Apis_mellifera:Anopheles_gambiae` | low_signal_d2b_le_0p02 | +0.10374 | +0.00093 | +0.06271 | +0.47070 | -0.46978 |
| `C.elegans:Caenorhabditis_briggsae` | strict_clean | +0.03092 | +0.02318 | +0.08459 | +0.06592 | -0.04274 |
| `Caenorhabditis_briggsae:C.elegans` | strict_clean | +0.03295 | +0.02616 | +0.06124 | +0.07702 | -0.05087 |
| `Caenorhabditis_briggsae:Human_T2T` | low_signal_d2b_le_0p02 | -0.03994 | -0.02085 | +0.00221 | +0.25613 | -0.27699 |
| `Human_T2T:Caenorhabditis_briggsae` | low_signal_d2b_le_0p02 | -0.04196 | -0.03225 | +0.01716 | +0.21220 | -0.24444 |
| `Arabidopsis:Oryza_sativa` | low_signal_d2b_le_0p02 | +0.00341 | -0.01159 | +0.06948 | +0.21667 | -0.22827 |
| `Oryza_sativa:Arabidopsis` | low_signal_d2b_le_0p02 | +0.00663 | +0.00708 | +0.04558 | +0.16696 | -0.15988 |
| `Arabidopsis:Zea_mays` | low_signal_d2b_le_0p02 | -0.01909 | -0.01251 | +0.04658 | +0.17889 | -0.19140 |
| `Zea_mays:Arabidopsis` | low_signal_d2b_le_0p02 | +0.02047 | +0.00770 | +0.05371 | +0.19047 | -0.18277 |
| `Arabidopsis:Solanum_lycopersicum` | low_signal_d2b_le_0p02 | +0.01117 | +0.01342 | +0.07646 | +0.10135 | -0.08793 |
| `Solanum_lycopersicum:Arabidopsis` | low_signal_d2b_le_0p02 | +0.01588 | +0.01268 | +0.07858 | +0.05179 | -0.03911 |
| `Arabidopsis:Glycine_max` | low_signal_d2b_le_0p02 | +0.01562 | +0.01556 | +0.07412 | +0.09012 | -0.07456 |
| `Glycine_max:Arabidopsis` | low_signal_d2b_le_0p02 | +0.00980 | +0.01086 | +0.07378 | +0.07023 | -0.05937 |
| `Arabidopsis:Populus_trichocarpa` | low_signal_d2b_le_0p02 | +0.01376 | +0.01829 | +0.07155 | +0.10606 | -0.08777 |
| `Populus_trichocarpa:Arabidopsis` | low_signal_d2b_le_0p02 | +0.01646 | +0.01073 | +0.08106 | +0.06226 | -0.05153 |
| `Oryza_sativa:Zea_mays` | low_signal_d2b_le_0p02 | +0.06440 | +0.01368 | +0.13633 | +0.02622 | -0.01254 |
| `Zea_mays:Oryza_sativa` | low_signal_d2b_le_0p02 | +0.09651 | +0.01534 | +0.17557 | +0.08012 | -0.06478 |
| `Solanum_lycopersicum:Glycine_max` | strict_clean | +0.02820 | +0.02253 | +0.07953 | +0.05174 | -0.02921 |
| `Glycine_max:Solanum_lycopersicum` | low_signal_d2b_le_0p02 | +0.02179 | +0.01885 | +0.08198 | +0.08518 | -0.06633 |
| `Populus_trichocarpa:Solanum_lycopersicum` | low_signal_d2b_le_0p02 | +0.02115 | +0.01903 | +0.08319 | +0.07393 | -0.05490 |
| `Solanum_lycopersicum:Populus_trichocarpa` | strict_clean | +0.02883 | +0.02880 | +0.08222 | +0.06831 | -0.03952 |
| `Yeast:Schizosaccharomyces_pombe` | low_signal_d2b_le_0p02 | +0.00169 | +0.00325 | +0.08696 | +0.11078 | -0.10752 |
| `Schizosaccharomyces_pombe:Yeast` | low_signal_d2b_le_0p02 | -0.00035 | +0.00346 | +0.07657 | +0.02496 | -0.02150 |
| `Yeast:Candida_albicans` | low_signal_d2b_le_0p02 | +0.01666 | +0.00930 | +0.12448 | +0.27362 | -0.26431 |
| `Candida_albicans:Yeast` | low_signal_d2b_le_0p02 | -0.03195 | -0.01352 | +0.05350 | -0.05953 | +0.04602 |
| `Yeast:Aspergillus_nidulans` | low_signal_d2b_le_0p02 | -0.01745 | -0.00670 | +0.01899 | +0.14644 | -0.15313 |
| `Aspergillus_nidulans:Yeast` | low_signal_d2b_le_0p02 | -0.01415 | -0.01074 | +0.05992 | +0.22452 | -0.23525 |
| `Yeast:Neurospora_crassa` | low_signal_d2b_le_0p02 | -0.00469 | -0.00516 | +0.03705 | +0.25923 | -0.26439 |
| `Neurospora_crassa:Yeast` | low_signal_d2b_le_0p02 | -0.02848 | -0.02946 | +0.04757 | +0.23730 | -0.26676 |
| `Candida_albicans:Aspergillus_nidulans` | low_signal_d2b_le_0p02 | -0.09912 | -0.05481 | -0.08494 | +0.17327 | -0.22808 |
| `Aspergillus_nidulans:Candida_albicans` | low_signal_d2b_le_0p02 | -0.00813 | -0.01517 | +0.06488 | +0.62661 | -0.64178 |
| `Neurospora_crassa:Aspergillus_nidulans` | low_signal_d2b_le_0p02 | -0.00523 | +0.00212 | +0.06562 | +0.00164 | +0.00048 |
| `Aspergillus_nidulans:Neurospora_crassa` | low_signal_d2b_le_0p02 | +0.01966 | +0.00975 | +0.09842 | +0.11329 | -0.10354 |
| `Schizosaccharomyces_pombe:Candida_albicans` | low_signal_d2b_le_0p02 | +0.00373 | +0.00469 | +0.12348 | +0.24922 | -0.24453 |
| `Candida_albicans:Schizosaccharomyces_pombe` | low_signal_d2b_le_0p02 | -0.04493 | -0.01528 | +0.07226 | -0.00184 | -0.01343 |
| `Plasmodium_falciparum:Human_T2T` | low_signal_d2b_le_0p02 | +0.03473 | -0.01733 | -0.03397 | +0.24430 | -0.26163 |
| `Human_T2T:Plasmodium_falciparum` | small_sample | +0.10192 | -0.04097 | +0.23063 | +0.79922 | -0.84019 |
| `Trypanosoma_brucei:Human_T2T` | low_signal_d2b_le_0p02 | -0.01436 | +0.00270 | +0.03391 | +0.26562 | -0.26292 |
| `Human_T2T:Trypanosoma_brucei` | small_sample | -0.04149 | -0.02429 | +0.00163 | +0.07572 | -0.10001 |
| `Leishmania_major:Human_T2T` | low_signal_d2b_le_0p02 | -0.02281 | -0.01680 | +0.04682 | +0.14855 | -0.16535 |
| `Human_T2T:Leishmania_major` | small_sample | -0.06073 | -0.06252 | +0.08041 | +0.32136 | -0.38387 |
| `Toxoplasma_gondii:Human_T2T` | low_signal_d2b_le_0p02 | -0.01328 | -0.01137 | +0.06951 | +0.21788 | -0.22925 |
| `Human_T2T:Toxoplasma_gondii` | small_sample | -0.05924 | -0.05953 | +0.00971 | +0.13262 | -0.19215 |
| `Dictyostelium_discoideum:Human_T2T` | low_signal_d2b_le_0p02 | -0.08158 | -0.06796 | -0.29015 | +0.23277 | -0.30073 |
| `Human_T2T:Dictyostelium_discoideum` | low_signal_d2b_le_0p02 | +0.09733 | -0.05055 | +0.19791 | +1.18631 | -1.23686 |
| `Plasmodium_falciparum:Toxoplasma_gondii` | small_sample | +0.05574 | -0.04537 | +0.04972 | +0.37961 | -0.42498 |
| `Toxoplasma_gondii:Plasmodium_falciparum` | small_sample | +0.00386 | -0.09456 | +0.27429 | +0.90414 | -0.99870 |
| `Trypanosoma_brucei:Leishmania_major` | small_sample | +0.02963 | +0.00775 | +0.06859 | +0.38576 | -0.37801 |
| `Leishmania_major:Trypanosoma_brucei` | small_sample | +0.00941 | -0.02117 | +0.02516 | +0.03469 | -0.05586 |
| `Dictyostelium_discoideum:Trypanosoma_brucei` | small_sample | -0.10984 | -0.06347 | -0.22788 | +0.09837 | -0.16184 |
| `Trypanosoma_brucei:Dictyostelium_discoideum` | low_signal_d2b_le_0p02 | +0.01169 | -0.01503 | +0.15390 | +1.08507 | -1.10010 |
| `Chlamydomonas_reinhardtii:Arabidopsis` | low_signal_d2b_le_0p02 | -0.01036 | -0.07060 | -0.00022 | +0.33587 | -0.40647 |
| `Arabidopsis:Chlamydomonas_reinhardtii` | low_signal_d2b_le_0p02 | -0.04627 | -0.05473 | +0.04327 | +0.83030 | -0.88502 |
| `Synechocystis_sp_PCC6803:Arabidopsis` | low_signal_d2b_le_0p02 | -0.13429 | -0.07347 | -0.12307 | +0.18936 | -0.26282 |
| `Arabidopsis:Synechocystis_sp_PCC6803` | low_signal_d2b_le_0p02 | -0.05083 | -0.02494 | +0.01915 | +0.39751 | -0.42244 |
| `Synechocystis_sp_PCC6803:Chlamydomonas_reinhardtii` | low_signal_d2b_le_0p02 | -0.21515 | -0.10677 | -0.15529 | +0.70019 | -0.80697 |
| `Chlamydomonas_reinhardtii:Synechocystis_sp_PCC6803` | low_signal_d2b_le_0p02 | -0.05985 | -0.08012 | +0.02388 | +0.46944 | -0.54956 |
| `E.coli:Pseudomonas_aeruginosa` | low_signal_d2b_le_0p02 | -0.02750 | +0.01612 | +0.13589 | +1.01925 | -1.00313 |
| `Pseudomonas_aeruginosa:E.coli` | low_signal_d2b_le_0p02 | -0.14461 | -0.11118 | -0.10178 | +0.09886 | -0.21005 |
| `E.coli:Salmonella_enterica` | strict_clean | +0.04045 | +0.03408 | +0.15997 | +0.11940 | -0.08532 |
| `Salmonella_enterica:E.coli` | strict_clean | +0.03675 | +0.02660 | +0.15278 | +0.08306 | -0.05646 |
| `Pseudomonas_aeruginosa:Salmonella_enterica` | low_signal_d2b_le_0p02 | -0.15793 | -0.11462 | -0.06276 | +0.06538 | -0.18000 |
| `Salmonella_enterica:Pseudomonas_aeruginosa` | low_signal_d2b_le_0p02 | +0.00399 | +0.01211 | +0.18609 | +0.94331 | -0.93120 |
| `B.subtilis:Staphylococcus_aureus` | small_sample | -0.05713 | -0.03922 | +0.10500 | +0.52081 | -0.56003 |
| `Staphylococcus_aureus:B.subtilis` | low_signal_d2b_le_0p02 | -0.07151 | -0.01730 | -0.05111 | +0.11914 | -0.13643 |
| `B.subtilis:Streptococcus_pyogenes` | low_signal_d2b_le_0p02 | -0.05234 | -0.03119 | +0.06053 | +0.25191 | -0.28310 |
| `Streptococcus_pyogenes:B.subtilis` | low_signal_d2b_le_0p02 | -0.06857 | -0.02865 | +0.00367 | +0.15041 | -0.17906 |
| `Staphylococcus_aureus:Streptococcus_pyogenes` | low_signal_d2b_le_0p02 | -0.03692 | -0.01171 | +0.05279 | +0.01851 | -0.03022 |
| `Streptococcus_pyogenes:Staphylococcus_aureus` | small_sample | -0.04032 | -0.02702 | +0.17580 | +0.29477 | -0.32179 |
| `Mycobacterium_tuberculosis:Human_T2T` | low_signal_d2b_le_0p02 | -0.15741 | -0.07360 | -0.06978 | +0.14853 | -0.22213 |
| `Human_T2T:Mycobacterium_tuberculosis` | low_signal_d2b_le_0p02 | -0.12049 | -0.11935 | +0.09672 | +0.54239 | -0.66174 |
| `Mycobacterium_tuberculosis:B.subtilis` | low_signal_d2b_le_0p02 | -0.17239 | -0.07705 | -0.09078 | +0.27713 | -0.35418 |
| `B.subtilis:Mycobacterium_tuberculosis` | low_signal_d2b_le_0p02 | -0.08888 | -0.02378 | -0.01119 | +0.73824 | -0.76202 |
| `Mycobacterium_tuberculosis:Streptococcus_pyogenes` | low_signal_d2b_le_0p02 | -0.24836 | -0.10574 | -0.18001 | +0.66856 | -0.77430 |
| `Streptococcus_pyogenes:Mycobacterium_tuberculosis` | low_signal_d2b_le_0p02 | -0.26128 | -0.10740 | -0.22404 | +0.97718 | -1.08459 |
| `Helicobacter_pylori:Bacteroides_fragilis` | low_signal_d2b_le_0p02 | -0.33664 | -0.16938 | -0.21235 | +0.18601 | -0.35539 |
| `Bacteroides_fragilis:Helicobacter_pylori` | low_signal_d2b_le_0p02 | -0.07751 | -0.05211 | +0.03061 | +0.72001 | -0.77211 |
| `Clostridioides_difficile:Bacteroides_fragilis` | low_signal_d2b_le_0p02 | -0.15463 | -0.05933 | -0.16260 | +0.14345 | -0.20278 |
| `Bacteroides_fragilis:Clostridioides_difficile` | low_signal_d2b_le_0p02 | -0.06547 | -0.03705 | +0.03201 | +0.79849 | -0.83554 |
| `Lactococcus_lactis:Streptococcus_pyogenes` | low_signal_d2b_le_0p02 | -0.01701 | +0.00408 | +0.09649 | +0.01539 | -0.01131 |
| `Streptococcus_pyogenes:Lactococcus_lactis` | low_signal_d2b_le_0p02 | -0.00325 | +0.00680 | +0.15536 | +0.19156 | -0.18476 |
| `Thermus_thermophilus:B.subtilis` | low_signal_d2b_le_0p02 | -0.35522 | -0.31245 | -0.34228 | +0.38531 | -0.69775 |
| `B.subtilis:Thermus_thermophilus` | low_signal_d2b_le_0p02 | -0.15405 | -0.07034 | -0.05524 | +1.85482 | -1.92516 |
| `Thermus_thermophilus:M.jannaschii` | low_signal_d2b_le_0p02 | -0.31504 | -0.28390 | -0.26803 | +0.80125 | -1.08516 |
| `M.jannaschii:Thermus_thermophilus` | low_signal_d2b_le_0p02 | -0.32475 | -0.09543 | -0.32155 | +1.94730 | -2.04273 |
| `Thermus_thermophilus:H.salinarum` | low_signal_d2b_le_0p02 | -0.32743 | -0.22023 | -0.02892 | +0.15696 | -0.37719 |
| `H.salinarum:Thermus_thermophilus` | low_signal_d2b_le_0p02 | -0.04343 | -0.03909 | +0.20175 | +0.73360 | -0.77268 |
