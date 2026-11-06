#!/usr/bin/env bash
# ==============================================================================
# MASTER AUTONOMOUS VERIFICATION & NCBI CACHE DOWNLOAD SUITE (1-CLICK RUNNER)
# Repository: https://github.com/Dr-Logvinovich/codon-torus-syntax
# Authors:    Victor Logvinovich & Jean-Claude Perez (2026)
#
# Usage:
#   bash run_all_audits.sh               # Download all caches + run all audits (default)
#   bash run_all_audits.sh --cache-only  # Download all 55 RefSeq + viral/locus caches & smoke test
#   bash run_all_audits.sh --full        # Also run full 120-rep bootstrap across multi-species panels
# ==============================================================================

set -euo pipefail

# Переходим в папку, где лежит сам скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MODE="${1:-standard}"
CDS_DIR="$SCRIPT_DIR/cds_cache"
OUT_ROOT="$SCRIPT_DIR/results_verified"
TIMEOUT=600

echo "================================================================================"
echo " CODON TORUS (v3.5) & CONDITIONAL CODON SYNTAX (v14.8) — MASTER PIPELINE"
echo " Working directory : $SCRIPT_DIR"
echo " Local NCBI cache  : $CDS_DIR"
echo " Output directory  : $OUT_ROOT"
echo " Execution mode    : $MODE"
echo "================================================================================"

mkdir -p "$CDS_DIR" "$OUT_ROOT"

# ------------------------------------------------------------------------------
# ШАГ 0. Проверка Python 3 и научных библиотек (NumPy, SciPy, Matplotlib)
# ------------------------------------------------------------------------------
echo ""
echo "[STEP 0/6] Checking Python 3 environment and scientific packages..."
python3 -c "import numpy, scipy, matplotlib" 2>/dev/null || {
    echo "  -> Installing numpy, scipy, matplotlib..."
    python3 -m pip install --upgrade pip
    python3 -m pip install numpy scipy matplotlib
}

echo "  -> Running built-in mathematical & genetic-code self-tests..."
python3 super_codon_proof_v14_8.py --self-test
python3 -c "import codon_torus_real_audit_v3_bnt162b2 as ct; res = ct.self_test_torus_theorems(); print('  [OK] Torus self-test passed:', res)"

# ------------------------------------------------------------------------------
# ШАГ 1. Аналитический аудит Кодонного Тора T_4^3 (9 геномов RefSeq + локус 17q22)
#        + 6-рамочная криминалистика Pfizer BNT162b2 и Moderna mRNA-1273 (v3.5)
# ------------------------------------------------------------------------------
echo ""
echo "[STEP 1/6] Running Double-Stranded Codon Torus Audit (9 RefSeq genomes + chr17 17q22 + BNT162b2 & mRNA-1273)..."
python3 codon_torus_real_audit_v3_bnt162b2.py \
    --cds-dir "$CDS_DIR" \
    --out-dir "$OUT_ROOT/01_codon_torus_and_bnt162b2_v3_5" \
    --bnt162b2-audit \
    --timeout "$TIMEOUT"

# ------------------------------------------------------------------------------
# ШАГ 2. Автономный аудит интронов/экзонов, фаз сплайсинга и онко-хотспотов
#        по 14 генам человека RefSeqGene (семейства TP53, RAS, BRCA1 и др. — v14.6)
# ------------------------------------------------------------------------------
echo ""
echo "[STEP 2/6] Running 14-Gene Human RefSeqGene Intron-Exon, Splicing & TP53/RAS Hotspot Audit (v14.6)..."
python3 super_codon_proof_v14_8.py \
    --intron-exon-audit \
    --intron-exon-genes TP53,TP63,TP73,BRCA1,PTEN,RB1,GAPDH,ACTB,HBB,KRAS,NRAS,HRAS,MYC,APOE \
    --cds-dir "$CDS_DIR" \
    --out-dir "$OUT_ROOT/02_intron_exon_p53_ras_v14_6" \
    --download-timeout "$TIMEOUT"

# ------------------------------------------------------------------------------
# ШАГ 3. Аудит экстремальной биосферы (Гигантские вирусы, скважины 2.8 км, 122C археи)
#        и синтетических клеток Крейга Вентера JCVI-syn1.0 / JCVI-syn3.0 (v14.7)
# ------------------------------------------------------------------------------
echo ""
echo "[STEP 3/6] Running 13-Genome Frontier Biosphere & Venter Synthetic Cell Forensics Audit (v14.7)..."
python3 super_codon_proof_v14_8.py \
    --frontier-audit \
    --cds-dir "$CDS_DIR" \
    --out-dir "$OUT_ROOT/03_frontier_and_jcvi_syn_v14_7" \
    --download-timeout "$TIMEOUT"

# ------------------------------------------------------------------------------
# ШАГ 4. Вирусные «Генетические трояны» (HPV-16 E7, Adenovirus E3 RID-a, EBV),
#        перезапись органелл, HGT-сканер (1500 генов) и краш-тест Max-CAI (v14.8)
# ------------------------------------------------------------------------------
echo ""
echo "[STEP 4/6] Running Viral Trojan, Organelle OS Rewriting, HGT Scan & Biotech Max-CAI Crash Audit (v14.8)..."
python3 super_codon_proof_v14_8.py \
    --chimera-audit \
    --cds-dir "$CDS_DIR" \
    --out-dir "$OUT_ROOT/04_viral_trojan_and_maxcai_crash_v14_8" \
    --download-timeout "$TIMEOUT"

# ------------------------------------------------------------------------------
# ШАГ 5. Полная выкачка всех 55 архивов RefSeq CDS (13 базовых + 41 расширенный
#        таксон биосферы + SARS-CoV-2 + митохондрия) и быстрый Smoke-аудит кэша
# ------------------------------------------------------------------------------
echo ""
echo "[STEP 5/6] Downloading & verifying all 55 RefSeq CDS archives across the biosphere (Smoke Stage)..."
python3 super_codon_proof_v14_8.py \
    --cds-dir "$CDS_DIR" \
    --out-dir "$OUT_ROOT/05_biosphere_55taxa_smoke" \
    --include-extended \
    --include-virus \
    --include-mito \
    --smoke-only \
    --download-timeout "$TIMEOUT"

# Опциональный полный бутстрап-прогон (--full)
if [[ "$MODE" == "--full" ]]; then
    echo ""
    echo "[OPTIONAL FULL MODE] Running full 13-taxon + SARS-CoV-2 bootstrap audit..."
    python3 super_codon_proof_v14_8.py \
        --cds-dir "$CDS_DIR" \
        --out-dir "$OUT_ROOT/06_primary_13taxa_full_bootstrap" \
        --include-virus \
        --full-only \
        --max-records 8000 \
        --min-trans 100 \
        --bootstrap-reps 120 \
        --h4-reps 8 \
        --h4-bootstrap-reps 40 \
        --h5-reps 3 \
        --h5-bootstrap-reps 20 \
        --h5-bootstrap-inner-reps 2 \
        --h6-reps 4 \
        --h6-bootstrap-reps 20 \
        --h7-reps 4 \
        --h7-bootstrap-reps 20 \
        --h5-transfer-reps 12 \
        --min-within-units 20 \
        --min-within-trans 50000 \
        --min-transfer-units 20 \
        --min-transfer-trans 50000 \
        --download-timeout "$TIMEOUT"
fi

# ------------------------------------------------------------------------------
# ШАГ 6. Сводка по скачанному кэшу cds_cache (проверка размеров файлов для GitHub)
# ------------------------------------------------------------------------------
echo ""
echo "================================================================================"
echo " [STEP 6/6] ALL AUDITS & DOWNLOADS COMPLETED SUCCESSFULLY!"
echo "================================================================================"
echo "Total size of local cache ($CDS_DIR):"
du -sh "$CDS_DIR" || true
echo ""
echo "File counts by category in cds_cache:"
echo "  - RefSeq CDS (*.fna.gz)          : $(find "$CDS_DIR" -maxdepth 1 -name "*_cds_from_genomic.fna.gz" | wc -l | tr -d ' ')"
echo "  - RefSeqGene Intron/Exon (*.gb)  : $(find "$CDS_DIR/refseqgene_gb" -name "*.gb" 2>/dev/null | wc -l | tr -d ' ')"
echo "  - Frontier Genomes (*.fna.gz)    : $(find "$CDS_DIR/frontier_cds" -name "*.fna.gz" 2>/dev/null | wc -l | tr -d ' ')"
echo "  - Chimera & Viral FASTA          : $(find "$CDS_DIR/chimera_entrez_cache" "$CDS_DIR/ncbi_live_viral" -type f 2>/dev/null | wc -l | tr -d ' ')"
echo ""
echo "Checking for files larger than 25 MB (GitHub web browser drag-and-drop limit):"
LARGE_FILES=$(find "$CDS_DIR" -type f -size +25M)
if [[ -z "$LARGE_FILES" ]]; then
    echo "  [GREAT] Every single file in cds_cache is under 25 MB! You can drag-and-drop all of them in the browser."
else
    echo "  [NOTE] The following files exceed 25 MB (upload via git CLI or compress before web upload):"
    find "$CDS_DIR" -type f -size +25M -exec ls -lh {} +
fi
echo "================================================================================"
