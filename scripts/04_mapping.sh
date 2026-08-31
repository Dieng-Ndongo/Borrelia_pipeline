#!/bin/bash
# =============================================================================
# 03_mapping.sh
# ALIGNEMENT BWA + BAM + STATISTIQUES
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
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

# Référence principale
REF_FASTA="data/reference/Borrelia_crocidurae_Achema_CP003426.fasta"

# Nombre de threads
THREADS=8

# Répertoires
KRAKEN_DIR="kraken"
MAPPING_DIR="mapping"

# =============================================================================
# VÉRIFICATION ENVIRONNEMENT
# =============================================================================

log_step "VÉRIFICATION DE L'ENVIRONNEMENT"

[[ "${CONDA_DEFAULT_ENV:-}" == "$ENV_NAME" ]] || \
    log_error "Activez : conda activate $ENV_NAME"

log_success "Environnement : $ENV_NAME"

# =============================================================================
# VÉRIFICATION DE LA RÉFÉRENCE
# =============================================================================

[[ -f "$REF_FASTA" ]] || \
    log_error "Référence absente : $REF_FASTA"

log_success "Référence : $REF_FASTA"

# =============================================================================
# VÉRIFICATION DES OUTILS
# =============================================================================

for tool in bwa samtools; do
    command -v "$tool" >/dev/null 2>&1 || \
        log_error "Outil '$tool' introuvable"

    log_success "$tool : $(command -v "$tool")"
done

# =============================================================================
# CRÉATION DU DOSSIER DE SORTIE
# =============================================================================

mkdir -p "$MAPPING_DIR"

# =============================================================================
# INDEXATION DE LA RÉFÉRENCE
# =============================================================================

log_step "PRÉPARATION DE LA RÉFÉRENCE"

if [[ ! -f "${REF_FASTA}.bwt" ]]; then
    log_info "Indexation BWA..."
    bwa index "$REF_FASTA"
    log_success "Index BWA créé"
else
    log_info "Index BWA déjà présent"
fi

if [[ ! -f "${REF_FASTA}.fai" ]]; then
    log_info "Création de l'index FASTA..."
    samtools faidx "$REF_FASTA"
    log_success "Index FASTA créé"
else
    log_info "Index FASTA déjà présent"
fi

# =============================================================================
# DÉTECTION DES ÉCHANTILLONS
# =============================================================================

log_step "DÉTECTION DES READS BORRELIA"

SAMPLES=()

while IFS= read -r file; do
    SAMPLE=$(basename "$file" "_borrelia_R1.fastq")
    SAMPLES+=("$SAMPLE")
done < <(
    find "$KRAKEN_DIR" \
        -type f \
        -name "*_borrelia_R1.fastq" \
        | sort
)

[[ ${#SAMPLES[@]} -gt 0 ]] || \
    log_error "Aucun read Borrelia trouvé dans $KRAKEN_DIR"

log_info "Échantillons détectés : ${#SAMPLES[@]}"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "  - $SAMPLE"
done

# =============================================================================
# MAPPING
# =============================================================================

for SAMPLE in "${SAMPLES[@]}"; do

    log_step "MAPPING — $SAMPLE"

    # -------------------------------------------------------------------------
    # FICHIERS D'ENTRÉE
    # -------------------------------------------------------------------------

    R1="${KRAKEN_DIR}/${SAMPLE}/${SAMPLE}_borrelia_R1.fastq"
    R2="${KRAKEN_DIR}/${SAMPLE}/${SAMPLE}_borrelia_R2.fastq"

    [[ -f "$R1" ]] || \
        log_error "R1 introuvable : $R1"

    [[ -f "$R2" ]] || \
        log_error "R2 introuvable : $R2"

    # -------------------------------------------------------------------------
    # VÉRIFICATION DU NOMBRE DE READS
    # -------------------------------------------------------------------------

    R1_READS=$(awk 'END {print NR/4}' "$R1")
    R2_READS=$(awk 'END {print NR/4}' "$R2")

    log_info "Reads R1 : $R1_READS"
    log_info "Reads R2 : $R2_READS"

    if [[ "$R1_READS" -ne "$R2_READS" ]]; then
        log_error "$SAMPLE : R1 ($R1_READS) ≠ R2 ($R2_READS)"
    fi

    log_success "R1/R2 correctement appariés"

    # -------------------------------------------------------------------------
    # FICHIERS DE SORTIE
    # -------------------------------------------------------------------------

    BAM="${MAPPING_DIR}/${SAMPLE}_aligned.bam"
    FLAGSTAT="${MAPPING_DIR}/${SAMPLE}_flagstat.txt"
    COVERAGE="${MAPPING_DIR}/${SAMPLE}_coverage.txt"
    DEPTH="${MAPPING_DIR}/${SAMPLE}_depth.txt"
    LOG="${MAPPING_DIR}/${SAMPLE}_mapping.log"

    # -------------------------------------------------------------------------
    # BWA MEM + SAMTOOLS SORT
    # -------------------------------------------------------------------------

    log_info "Mapping avec BWA-MEM..."

    bwa mem \
        -t "$THREADS" \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        "$REF_FASTA" \
        "$R1" \
        "$R2" \
        2> "$LOG" \
        | samtools sort \
            -@ "$THREADS" \
            -o "$BAM" -

    log_success "Mapping BWA terminé"

    # -------------------------------------------------------------------------
    # INDEXATION BAM
    # -------------------------------------------------------------------------

    log_info "Indexation du BAM..."

    samtools index "$BAM"

    log_success "BAM indexé"

    # -------------------------------------------------------------------------
    # STATISTIQUES FLAGSTAT
    # -------------------------------------------------------------------------

    log_info "Calcul des statistiques flagstat..."

    samtools flagstat "$BAM" \
        > "$FLAGSTAT"

    # -------------------------------------------------------------------------
    # COUVERTURE
    # -------------------------------------------------------------------------

    log_info "Calcul de la couverture..."

    samtools coverage "$BAM" \
        > "$COVERAGE"

    # -------------------------------------------------------------------------
    # PROFONDEUR
    # -------------------------------------------------------------------------

    log_info "Calcul de la profondeur..."

    samtools depth -a "$BAM" \
        > "$DEPTH"

    # -------------------------------------------------------------------------
    # RÉSUMÉ
    # -------------------------------------------------------------------------

    log_success "$SAMPLE : mapping terminé"

    echo
    echo "  BAM       : $BAM"
    echo "  Flagstat  : $FLAGSTAT"
    echo "  Coverage  : $COVERAGE"
    echo "  Depth     : $DEPTH"

done

# =============================================================================
# RÉSUMÉ FINAL
# =============================================================================

log_step "MAPPING TERMINÉ POUR TOUS LES ÉCHANTILLONS"

echo
echo "Référence utilisée :"
echo "  $REF_FASTA"

echo
echo "Échantillons traités : ${#SAMPLES[@]}"

echo
echo "Résultats :"
echo "  └── $MAPPING_DIR/"
