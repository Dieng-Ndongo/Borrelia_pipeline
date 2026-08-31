#!/bin/bash
# =============================================================================
# 01_qc.sh
# CONTRÔLE QUALITÉ — Tous les échantillons
#
# Entrée :
#   data/raw/*_L001_R1_001.fastq.gz
#   data/raw/*_L001_R2_001.fastq.gz
#   data/raw/*_L002_R1_001.fastq.gz
#   data/raw/*_L002_R2_001.fastq.gz
#
# Sorties :
#   data/raw/*_merged_R1.fastq.gz
#   data/raw/*_merged_R2.fastq.gz
#   data/clean/*_R1_clean.fastq.gz
#   data/clean/*_R2_clean.fastq.gz
#   qc/*
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log_info()    { echo -e "${BLUE}[INFO]${NC}    $1"; }
log_success() { echo -e "${GREEN}[OK]${NC}      $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error()   { echo -e "${RED}[ERREUR]${NC}  $1"; exit 1; }

log_step() {
    echo -e "\n${CYAN}============================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
}

# -----------------------------------------------------------------------------
# Vérifications
# -----------------------------------------------------------------------------

ENV_NAME="borrelia_pipeline"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
    log_error "Activez l'environnement : conda activate $ENV_NAME"
fi

mkdir -p data/raw data/clean qc results

if [[ $(find data/raw -maxdepth 1 -name "*.fastq.gz" | wc -l) -eq 0 ]]; then
    log_error "Aucun fichier .fastq.gz trouvé dans data/raw/"
fi

for tool in fastqc fastp multiqc; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        log_error "Outil '$tool' introuvable."
    fi
done

# -----------------------------------------------------------------------------
# Détection des échantillons
# -----------------------------------------------------------------------------

log_step "Détection des échantillons"

SAMPLES=()

while IFS= read -r file; do
    SAMPLE=$(basename "$file" "_L001_R1_001.fastq.gz")
    SAMPLES+=("$SAMPLE")
done < <(find data/raw -maxdepth 1 -type f \
    -name "*_L001_R1_001.fastq.gz" | sort)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    log_error "Aucun échantillon détecté dans data/raw/"
fi

log_info "${#SAMPLES[@]} échantillon(s) détecté(s) :"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "    - $SAMPLE"
done

# =============================================================================
# TRAITEMENT DE TOUS LES ÉCHANTILLONS
# =============================================================================

for SAMPLE in "${SAMPLES[@]}"; do

    log_step "QC — $SAMPLE"

    L001_R1="data/raw/${SAMPLE}_L001_R1_001.fastq.gz"
    L001_R2="data/raw/${SAMPLE}_L001_R2_001.fastq.gz"
    L002_R1="data/raw/${SAMPLE}_L002_R1_001.fastq.gz"
    L002_R2="data/raw/${SAMPLE}_L002_R2_001.fastq.gz"

    for FILE in "$L001_R1" "$L001_R2" "$L002_R1" "$L002_R2"; do
        [[ -f "$FILE" ]] || log_error "Fichier introuvable : $FILE"
    done

    R1_RAW="data/raw/${SAMPLE}_merged_R1.fastq.gz"
    R2_RAW="data/raw/${SAMPLE}_merged_R2.fastq.gz"

    # -------------------------------------------------------------------------
    # Fusion des lanes
    # -------------------------------------------------------------------------

    if [[ ! -f "$R1_RAW" ]]; then
        log_info "Fusion des lanes R1..."
        cat "$L001_R1" "$L002_R1" > "$R1_RAW"
    else
        log_info "R1 fusionné déjà présent."
    fi

    if [[ ! -f "$R2_RAW" ]]; then
        log_info "Fusion des lanes R2..."
        cat "$L001_R2" "$L002_R2" > "$R2_RAW"
    else
        log_info "R2 fusionné déjà présent."
    fi

    # -------------------------------------------------------------------------
    # Vérification du nombre de reads
    # -------------------------------------------------------------------------

    N_R1=$(zcat "$R1_RAW" | wc -l)
    N_R2=$(zcat "$R2_RAW" | wc -l)

    N_R1=$((N_R1 / 4))
    N_R2=$((N_R2 / 4))

    log_info "Reads R1 : $N_R1"
    log_info "Reads R2 : $N_R2"

    if [[ "$N_R1" -ne "$N_R2" ]]; then
        log_error "Déséquilibre R1/R2 pour $SAMPLE"
    fi

    # -------------------------------------------------------------------------
    # FastQC brut
    # -------------------------------------------------------------------------

    log_info "FastQC des reads bruts..."

    fastqc \
        "$R1_RAW" \
        "$R2_RAW" \
        -o qc \
        -t 8 \
        --quiet

    # -------------------------------------------------------------------------
    # Fastp
    # -------------------------------------------------------------------------

    R1_CLEAN="data/clean/${SAMPLE}_R1_clean.fastq.gz"
    R2_CLEAN="data/clean/${SAMPLE}_R2_clean.fastq.gz"

    log_info "Nettoyage avec Fastp..."

    fastp \
        -i "$R1_RAW" \
        -I "$R2_RAW" \
        -o "$R1_CLEAN" \
        -O "$R2_CLEAN" \
        --html "qc/${SAMPLE}_fastp_report.html" \
        --json "qc/${SAMPLE}_fastp_report.json" \
        --thread 8 \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --detect_adapter_for_pe \
        2> "results/${SAMPLE}_fastp.log"

    # -------------------------------------------------------------------------
    # FastQC nettoyé
    # -------------------------------------------------------------------------

    log_info "FastQC des reads nettoyés..."

    fastqc \
        "$R1_CLEAN" \
        "$R2_CLEAN" \
        -o qc \
        -t 8 \
        --quiet

    log_success "$SAMPLE : QC terminé"

done

# =============================================================================
# MultiQC GLOBAL
# =============================================================================

log_step "MultiQC — Tous les échantillons"

multiqc qc/ \
    -o qc/ \
    --quiet \
    --filename "multiqc_all_samples"

log_success "QC terminé pour tous les échantillons"
log_success "Rapport : qc/multiqc_all_samples.html"
