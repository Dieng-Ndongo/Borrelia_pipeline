#!/bin/bash
# =============================================================================
# 05_variants.sh
# APPEL DE VARIANTS + FILTRAGE + CONSENSUS
#
# Entrées :
#   mapping/*_aligned.bam
#   data/reference/Borrelia_crocidurae_Achema_CP003426.fasta  ← même nom que 04_mapping.sh
#
# Sorties :
#   variants/SAMPLE_variants.vcf
#   variants/SAMPLE_variants_filtered.vcf.gz
#   variants/SAMPLE_consensus.fasta
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log_info()    { echo -e "${BLUE}[INFO]${NC}    $1"; }
log_success() { echo -e "${GREEN}[OK]${NC}      $1"; }
log_error()   { echo -e "${RED}[ERREUR]${NC}  $1"; exit 1; }

log_step() {
    echo -e "\n${CYAN}============================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
}

# =============================================================================
# CONFIGURATION
# =============================================================================

ENV_NAME="borrelia_pipeline"

# ⚠️  Même chemin que dans 04_mapping.sh
REF_FASTA="data/reference/Borrelia_crocidurae_Achema_CP003426.fasta"

MAPPING_DIR="mapping"
VARIANTS_DIR="variants"

MIN_QUAL=20
MIN_DEPTH=10

KRAKEN_DIR="kraken"

B_CROCIDURAE_TAXIDS=("29520" "1155096" "1293575")

# =============================================================================
# VÉRIFICATIONS
# =============================================================================

[[ "${CONDA_DEFAULT_ENV:-}" == "$ENV_NAME" ]] || \
    log_error "Activez : conda activate $ENV_NAME"

[[ -f "$REF_FASTA" ]] || \
    log_error "Référence absente : $REF_FASTA
    Ce chemin doit être identique à celui utilisé dans 04_mapping.sh"

[[ -d "$MAPPING_DIR" ]] || \
    log_error "Dossier mapping/ introuvable. Lancez d'abord 04_mapping.sh"

for tool in bcftools bgzip tabix; do
    command -v "$tool" >/dev/null 2>&1 || \
        log_error "Outil '$tool' introuvable"
done

mkdir -p "$VARIANTS_DIR" results

# =============================================================================
# DÉTECTION DES BAM
# =============================================================================

log_step "DÉTECTION DES ÉCHANTILLONS"

SAMPLES=()

while IFS= read -r bam; do
    SAMPLE=$(basename "$bam" "_aligned.bam")
    SAMPLES+=("$SAMPLE")
done < <(find "$MAPPING_DIR" -maxdepth 1 -type f -name "*_aligned.bam" | sort)

[[ ${#SAMPLES[@]} -gt 0 ]] || \
    log_error "Aucun BAM trouvé dans $MAPPING_DIR/"

log_info "${#SAMPLES[@]} échantillon(s) détecté(s) :"
for SAMPLE in "${SAMPLES[@]}"; do
    echo "    - $SAMPLE"
done

# =============================================================================
# TRAITEMENT
# =============================================================================

TOTAL=${#SAMPLES[@]}
SUCCESS=0
FAILED=0

for SAMPLE in "${SAMPLES[@]}"; do

    log_step "Variants — $SAMPLE"

    # -------------------------------------------------------------------------
    # FILTRE TAXONOMIQUE KRAKEN
    # -------------------------------------------------------------------------

    KRAKEN_REPORT="${KRAKEN_DIR}/${SAMPLE}/${SAMPLE}_report.txt"

    if [[ ! -r "${KRAKEN_REPORT}" ]]; then
        log_info "${SAMPLE} : rapport Kraken absent — IGNORÉ."
        FAILED=$((FAILED + 1))
        continue
    fi

    BC_READS=0
    for TAXID in "${B_CROCIDURAE_TAXIDS[@]}"; do
        N="$(awk -v t="${TAXID}" '$5 == t {print $2; exit}' "${KRAKEN_REPORT}")"
        BC_READS=$(( BC_READS + ${N:-0} ))
    done

    if (( BC_READS == 0 )); then
        log_info "${SAMPLE} : aucun read B. crocidurae (Kraken) — IGNORÉ."
        continue
    fi

    log_info "${SAMPLE} : ${BC_READS} reads B. crocidurae détectés."

    BAM="${MAPPING_DIR}/${SAMPLE}_aligned.bam"
    VCF="${VARIANTS_DIR}/${SAMPLE}_variants.vcf"
    FILTERED="${VARIANTS_DIR}/${SAMPLE}_variants_filtered.vcf"
    CONSENSUS="${VARIANTS_DIR}/${SAMPLE}_consensus.fasta"
    LOG="results/${SAMPLE}_variants.log"

    if [[ ! -f "$BAM" ]]; then
        log_info "${SAMPLE} : BAM introuvable — IGNORÉ."
        FAILED=$((FAILED + 1))
        continue
    fi

    # -------------------------------------------------------------------------
    # APPEL DE VARIANTS
    # -------------------------------------------------------------------------

    log_info "Appel de variants avec bcftools..."

    bcftools mpileup \
        -f "$REF_FASTA" \
        --min-MQ 20 \
        --min-BQ 20 \
        "$BAM" \
        2> "$LOG" \
        | bcftools call \
            -mv \
            -Ov \
            -o "$VCF"

    log_success "Variants appelés"

    # -------------------------------------------------------------------------
    # FILTRAGE
    # -------------------------------------------------------------------------

    log_info "Filtrage (QUAL>=$MIN_QUAL, DP>=$MIN_DEPTH)..."

    bcftools filter \
        -e "QUAL<${MIN_QUAL} || DP<${MIN_DEPTH}" \
        "$VCF" \
        > "$FILTERED"

    bgzip -f "$FILTERED"
    tabix -f "${FILTERED}.gz"

    N_SNPS=$(bcftools stats "${FILTERED}.gz" \
        | awk '/^SN.*number of SNPs/{print $NF}')

    log_info "SNPs filtrés : ${N_SNPS:-0}"

    # -------------------------------------------------------------------------
    # CONSENSUS
    # -------------------------------------------------------------------------

    log_info "Génération du consensus FASTA..."

    bcftools consensus \
        -f "$REF_FASTA" \
        "${FILTERED}.gz" \
        > "$CONSENSUS"

    # Renommer l'en-tête avec le nom de l'échantillon
    sed -i "s/>.*/>Borrelia_${SAMPLE}/" "$CONSENSUS"

    log_success "$SAMPLE — consensus : $CONSENSUS"

    echo
    echo "  VCF brut    : $VCF"
    echo "  VCF filtré  : ${FILTERED}.gz"
    echo "  Consensus   : $CONSENSUS"
    echo "  SNPs        : ${N_SNPS:-0}"

    SUCCESS=$((SUCCESS + 1))

done

# =============================================================================
# RÉSUMÉ
# =============================================================================

log_step "VARIANTS TERMINÉS"

echo
echo -e "  Réussis : ${GREEN}$SUCCESS / $TOTAL${NC}"
echo -e "  Échoués : ${RED}$FAILED / $TOTAL${NC}"
echo
echo "  Résultats : $VARIANTS_DIR/"
echo
