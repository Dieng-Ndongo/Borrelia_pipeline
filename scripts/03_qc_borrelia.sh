#!/bin/bash
# =============================================================================
# 03_qc_borrelia.sh
# CONTRÔLE QUALITÉ DES READS BORRELIA EXTRAITS PAR KRAKEN2
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}[INFO]${NC}    $1"
}

log_success() {
    echo -e "${GREEN}[OK]${NC}      $1"
}

log_error() {
    echo -e "${RED}[ERREUR]${NC}  $1"
    exit 1
}

log_step() {
    echo -e "\n${CYAN}============================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
}

# =============================================================================
# CONFIGURATION
# =============================================================================

ENV_NAME="borrelia_pipeline"

KRAKEN_DIR="kraken"
QC_DIR="qc/borrelia"

THREADS=8

# =============================================================================
# VÉRIFICATION ENVIRONNEMENT
# =============================================================================

log_step "VÉRIFICATION DE L'ENVIRONNEMENT"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
    log_error "Activez d'abord l'environnement :
    conda activate $ENV_NAME"
fi

log_success "Environnement : $ENV_NAME"

# =============================================================================
# VÉRIFICATION DES OUTILS
# =============================================================================

for tool in fastqc multiqc; do
    command -v "$tool" >/dev/null 2>&1 || \
        log_error "Outil '$tool' introuvable."

    log_success "$tool : $(command -v "$tool")"
done

# =============================================================================
# CRÉATION DES RÉPERTOIRES
# =============================================================================

mkdir -p "$QC_DIR"

# =============================================================================
# DÉTECTION DES READS BORRELIA
# =============================================================================

R1_FILES=()

while IFS= read -r file; do
    R1_FILES+=("$file")
done < <(
    find "$KRAKEN_DIR" \
        -type f \
        -name "*_borrelia_R1.fastq" \
        | sort
)

[[ ${#R1_FILES[@]} -gt 0 ]] || \
    log_error "Aucun read Borrelia R1 trouvé dans $KRAKEN_DIR"

log_info "Échantillons détectés : ${#R1_FILES[@]}"

# =============================================================================
# FASTQC
# =============================================================================

log_step "FASTQC — READS BORRELIA"

FASTQC_DIR="${QC_DIR}/fastqc"

mkdir -p "$FASTQC_DIR"

for R1 in "${R1_FILES[@]}"; do

    SAMPLE=$(basename "$R1" _borrelia_R1.fastq)

    R2=$(dirname "$R1")/${SAMPLE}_borrelia_R2.fastq

    [[ -f "$R2" ]] || \
        log_error "R2 introuvable pour $SAMPLE : $R2"

    log_info "Analyse : $SAMPLE"

    fastqc \
        "$R1" \
        "$R2" \
        --outdir "$FASTQC_DIR" \
        --threads "$THREADS"

    log_success "FastQC terminé : $SAMPLE"

done

# =============================================================================
# MULTIQC
# =============================================================================

log_step "MULTIQC — RAPPORT GLOBAL"

MULTIQC_DIR="${QC_DIR}/multiqc"

mkdir -p "$MULTIQC_DIR"

multiqc \
    "$FASTQC_DIR" \
    --outdir "$MULTIQC_DIR" \
    --filename "multiqc_borrelia.html" \
    --force

log_success "Rapport MultiQC généré"

# =============================================================================
# COMPTAGE DES READS
# =============================================================================

log_step "RÉSUMÉ DES READS BORRELIA"

SUMMARY="${QC_DIR}/borrelia_reads_summary.tsv"

echo -e "Sample\tR1_reads\tR2_reads\tStatus" > "$SUMMARY"

for R1 in "${R1_FILES[@]}"; do

    SAMPLE=$(basename "$R1" _borrelia_R1.fastq)

    R2=$(dirname "$R1")/${SAMPLE}_borrelia_R2.fastq

    R1_READS=$(awk 'END {print NR/4}' "$R1")
    R2_READS=$(awk 'END {print NR/4}' "$R2")

    if [[ "$R1_READS" -eq "$R2_READS" ]]; then
        STATUS="OK"
    else
        STATUS="ERROR"
    fi

    printf "%s\t%s\t%s\t%s\n" \
        "$SAMPLE" \
        "$R1_READS" \
        "$R2_READS" \
        "$STATUS" >> "$SUMMARY"

done

column -t -s $'\t' "$SUMMARY"

# =============================================================================
# RÉSUMÉ FINAL
# =============================================================================

log_step "QC BORRELIA TERMINÉ"

echo
echo "Résultats FastQC :"
echo "  └── $FASTQC_DIR/"

echo
echo "Rapport MultiQC :"
echo "  └── $MULTIQC_DIR/multiqc_borrelia.html"

echo
echo "Résumé des reads :"
echo "  └── $SUMMARY"

echo

if grep -q $'\tERROR$' "$SUMMARY"; then
    log_error "Une ou plusieurs paires R1/R2 présentent un nombre différent de reads."
else
    log_success "Toutes les paires R1/R2 sont correctement appariées."
fi
