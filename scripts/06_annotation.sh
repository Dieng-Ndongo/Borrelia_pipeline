#!/bin/bash
# =============================================================================
# 06_annotation.sh
# ANNOTATION FONCTIONNELLE DES VARIANTS AVEC SNPEFF
#
# Pipeline Borrelia
#
# Entrées :
#   data/reference/Borrelia_crocidurae_Achema_CP003426.fasta
#   data/reference/annotation/CP003426.1.gb
#   variants/*_variants_filtered.vcf.gz
#
# Sorties :
#   data/reference/annotation/snpeff/
#       ├── snpEff.config
#       └── data/
#           └── Borrelia_crocidurae_CP003426/
#               ├── genes.gbk
#               ├── sequence.bin
#               └── snpEffectPredictor.bin
#
#   variants/annotated/
#       ├── SAMPLE_annotated.vcf.gz
#       └── SAMPLE_annotated.tsv
#
#   results/annotation/
#       ├── snpeff_build.log
#       └── SAMPLE_snpeff.log
#
# IMPORTANT :
#   - Aucun chemin absolu n'est enregistré dans le projet.
#   - Le chemin absolu utilisé par SnpEff est généré dynamiquement
#     avec $(pwd) au moment de l'exécution.
#   - SnpEff utilise -nodownload pour éviter toute tentative de
#     téléchargement d'une base personnalisée.
# =============================================================================

set -euo pipefail

# =============================================================================
# COULEURS
# =============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# =============================================================================
# FONCTIONS
# =============================================================================

log_info() {
    echo -e "${BLUE}[INFO]${NC}    $1"
}

log_success() {
    echo -e "${GREEN}[OK]${NC}      $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERREUR]${NC}  $1"
    exit 1
}

log_step() {
    echo
    echo -e "${CYAN}============================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
}

# =============================================================================
# 0 — RACINE DU PROJET
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

log_info "Racine du projet : $PROJECT_ROOT"

# =============================================================================
# CONFIGURATION
# =============================================================================

ENV_NAME="borrelia_pipeline"

# Référence
REF_FASTA="data/reference/Borrelia_crocidurae_Achema_CP003426.fasta"

# GenBank de référence
GENBANK="data/reference/annotation/CP003426.1.gb"

# Répertoire SnpEff
SNPEFF_DIR="data/reference/annotation/snpeff"

# Répertoire de données SnpEff
SNPEFF_DATA_DIR="${SNPEFF_DIR}/data"

# Génome SnpEff
GENOME_ID="Borrelia_crocidurae_CP003426"

# Nom du génome dans la configuration
GENOME_NAME="Borrelia_crocidurae_Achema"

# Dossier contenant les VCF
VARIANTS_DIR="variants"

# Dossier contenant les VCF annotés
ANNOTATED_DIR="${VARIANTS_DIR}/annotated"

# Logs
RESULTS_DIR="results/annotation"

BUILD_LOG="${RESULTS_DIR}/snpeff_build.log"

# Seuils / options
SNPEFF_FORMAT="-genbank"

# =============================================================================
# 1 — VÉRIFICATION DE L'ENVIRONNEMENT
# =============================================================================

log_step "1 — VÉRIFICATION DE L'ENVIRONNEMENT"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
    log_error "Environnement incorrect.

Environnement actuel :
    ${CONDA_DEFAULT_ENV:-aucun}

Environnement attendu :
    $ENV_NAME

Activez-le avec :
    conda activate $ENV_NAME"
fi

log_success "Environnement actif : $ENV_NAME"

# =============================================================================
# 2 — VÉRIFICATION DES OUTILS
# =============================================================================

log_step "2 — VÉRIFICATION DES OUTILS"

TOOLS=(
    snpEff
    bcftools
    bgzip
    tabix
    curl
)

for TOOL in "${TOOLS[@]}"; do

    if ! command -v "$TOOL" >/dev/null 2>&1; then
        log_error "Outil '$TOOL' introuvable dans l'environnement."
    fi

    log_success "$TOOL : $(command -v "$TOOL")"

done

echo
log_info "Version SnpEff :"
snpEff -version

# =============================================================================
# 3 — VÉRIFICATION DE LA RÉFÉRENCE
# =============================================================================

log_step "3 — VÉRIFICATION DE LA RÉFÉRENCE"

if [[ ! -f "$REF_FASTA" ]]; then
    log_error "Référence FASTA absente :

    $REF_FASTA

Cette référence doit être présente avant l'annotation."
fi

log_success "Référence FASTA trouvée :"
echo "    $REF_FASTA"

REF_HEADER=$(head -n 1 "$REF_FASTA")

log_info "En-tête FASTA : $REF_HEADER"

# =============================================================================
# 4 — PRÉPARATION DES DOSSIERS
# =============================================================================

log_step "4 — PRÉPARATION DES DOSSIERS"

mkdir -p \
    "$(dirname "$GENBANK")" \
    "$SNPEFF_DIR" \
    "$SNPEFF_DATA_DIR" \
    "$ANNOTATED_DIR" \
    "$RESULTS_DIR"

log_success "Dossiers préparés"

# =============================================================================
# 5 — TÉLÉCHARGEMENT AUTOMATIQUE DU GENBANK
# =============================================================================

log_step "5 — VÉRIFICATION / TÉLÉCHARGEMENT DU GENBANK"

if [[ -f "$GENBANK" ]]; then

    log_info "Annotation GenBank déjà présente."
    log_info "Fichier : $GENBANK"

else

    log_warning "Fichier GenBank absent."
    log_info "Téléchargement automatique depuis NCBI..."

    # -------------------------------------------------------------------------
    # Téléchargement du GenBank correspondant à CP003426.1
    # -------------------------------------------------------------------------

    NCBI_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP003426.1&rettype=gbwithparts&retmode=text"

    TEMP_GENBANK="${GENBANK}.tmp"

    rm -f "$TEMP_GENBANK"

    curl \
        -L \
        --fail \
        --retry 3 \
        --connect-timeout 30 \
        --max-time 300 \
        "$NCBI_URL" \
        -o "$TEMP_GENBANK"

    # Vérification du téléchargement
    if [[ ! -s "$TEMP_GENBANK" ]]; then
        rm -f "$TEMP_GENBANK"

        log_error "Échec du téléchargement du GenBank CP003426.1."
    fi

    # Vérification minimale du format GenBank
    if ! grep -q "^LOCUS" "$TEMP_GENBANK"; then

        rm -f "$TEMP_GENBANK"

        log_error "Le fichier téléchargé ne semble pas être un fichier GenBank valide."
    fi

    mv "$TEMP_GENBANK" "$GENBANK"

    log_success "GenBank téléchargé automatiquement :"
    echo "    $GENBANK"

fi

# =============================================================================
# 6 — CONTRÔLE DU GENBANK
# =============================================================================

log_step "6 — CONTRÔLE DE L'ANNOTATION GENBANK"

if [[ ! -s "$GENBANK" ]]; then
    log_error "Fichier GenBank vide ou illisible :
    $GENBANK"
fi

if ! grep -q "^LOCUS" "$GENBANK"; then
    log_error "Le fichier GenBank ne contient pas de ligne LOCUS."
fi

log_success "Fichier GenBank valide et lisible."

echo
log_info "Annotations ribosomales / gènes recherchés :"

grep -n -E "flaB|16S|23S|5S|ribosomal RNA" "$GENBANK" \
    | head -20 \
    || log_warning "Aucune annotation recherchée détectée."

# =============================================================================
# 7 — PRÉPARATION DE LA BASE SNPEFF
# =============================================================================

log_step "7 — PRÉPARATION DE LA BASE SNPEFF"

GENOME_DIR="${SNPEFF_DATA_DIR}/${GENOME_ID}"

mkdir -p "$GENOME_DIR"

GENES_GBK="${GENOME_DIR}/genes.gbk"

# -------------------------------------------------------------------------
# IMPORTANT :
#
# SnpEff attend précisément :
#
# data/
# └── GenomeID/
#       └── genes.gbk
#
# Nous recréons donc systématiquement genes.gbk à partir du GenBank source.
# -------------------------------------------------------------------------

log_info "Création / mise à jour de genes.gbk..."

cp "$GENBANK" "$GENES_GBK"

if [[ ! -s "$GENES_GBK" ]]; then
    log_error "Impossible de créer :
    $GENES_GBK"
fi

log_success "GenBank SnpEff créé :"
echo "    $GENES_GBK"

# Vérification
if ! grep -q "^LOCUS" "$GENES_GBK"; then
    log_error "genes.gbk n'est pas un GenBank valide."
fi

log_success "genes.gbk valide."

# =============================================================================
# 8 — CONFIGURATION SNPEFF
# =============================================================================

log_step "8 — CONFIGURATION SNPEFF"

SNPEFF_CONFIG="${SNPEFF_DIR}/snpEff.config"

cat > "$SNPEFF_CONFIG" << EOF
# =============================================================================
# Configuration SnpEff
# Générée automatiquement par 06_annotation.sh
#
# Projet : Borrelia pipeline
# Génome : ${GENOME_ID}
#
# IMPORTANT :
# Le chemin data.dir n'est volontairement PAS défini ici.
#
# Il est fourni dynamiquement à SnpEff avec :
#
#     -dataDir "\$(pwd)/data/reference/annotation/snpeff/data"
#
# Cela permet de conserver une configuration portable et reproductible.
# =============================================================================

genome : ${GENOME_ID}
${GENOME_ID}.genome : ${GENOME_NAME}
EOF

log_success "Configuration SnpEff créée :"
echo "    $SNPEFF_CONFIG"

echo
cat "$SNPEFF_CONFIG"

# =============================================================================
# 9 — CHEMIN ABSOLU DYNAMIQUE POUR SNPEFF
# =============================================================================

log_step "9 — CONFIGURATION DU DATA DIRECTORY"

# IMPORTANT :
# Ce chemin absolu n'est PAS écrit dans le projet.
# Il est uniquement construit pendant l'exécution.

SNPEFF_DATA_ABS="${PROJECT_ROOT}/${SNPEFF_DATA_DIR}"

log_info "Data directory utilisé pendant l'exécution :"
echo "    $SNPEFF_DATA_ABS"

# =============================================================================
# 10 — CONSTRUCTION / VÉRIFICATION DE LA BASE SNPEFF
# =============================================================================

log_step "10 — CONSTRUCTION DE LA BASE SNPEFF"

SNPEFF_BIN="${GENOME_DIR}/snpEffectPredictor.bin"
SNPEFF_SEQ="${GENOME_DIR}/sequence.bin"

# -------------------------------------------------------------------------
# Vérifier si la base est déjà complète
# -------------------------------------------------------------------------

if [[ -s "$SNPEFF_BIN" && -s "$SNPEFF_SEQ" ]]; then

    log_success "Base SnpEff déjà présente."

    echo
    echo "    Genome ID : $GENOME_ID"
    echo "    Base      : $SNPEFF_BIN"
    echo "    Sequence  : $SNPEFF_SEQ"

else

    log_info "Base SnpEff absente ou incomplète."
    log_info "Construction automatique avec GenBank..."
    log_info "Cette étape peut prendre quelques minutes."

    rm -f \
        "$SNPEFF_BIN" \
        "$SNPEFF_SEQ"

    if ! snpEff build \
        -c "$SNPEFF_CONFIG" \
        -dataDir "$SNPEFF_DATA_ABS" \
        -genbank \
        -noCheckCds \
        -noCheckProtein \
        -v \
        "$GENOME_ID" \
        > "$BUILD_LOG" \
        2>&1
    then

        echo
        log_error "Échec de construction de la base SnpEff.

Consultez :
    $BUILD_LOG"

    fi

    log_success "Commande de construction SnpEff terminée."

fi

# =============================================================================
# 11 — VÉRIFICATION DE LA BASE SNPEFF
# =============================================================================

log_step "11 — VÉRIFICATION DE LA BASE SNPEFF"

if [[ ! -s "$SNPEFF_BIN" ]]; then
    log_error "snpEffectPredictor.bin absent après construction :
    $SNPEFF_BIN"
fi

if [[ ! -s "$SNPEFF_SEQ" ]]; then
    log_error "sequence.bin absent après construction :
    $SNPEFF_SEQ"
fi

log_success "Base SnpEff opérationnelle."

echo
echo "    $SNPEFF_BIN"
echo "    $SNPEFF_SEQ"

# =============================================================================
# 12 — DÉTECTION DES VCF
# =============================================================================

log_step "12 — DÉTECTION DES VCF À ANNOTER"

VCF_FILES=()

while IFS= read -r VCF; do
    VCF_FILES+=("$VCF")
done < <(
    find "$VARIANTS_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*_variants_filtered.vcf.gz" \
        | sort
)

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    log_error "Aucun VCF filtré trouvé dans :
    $VARIANTS_DIR/

Attendu :
    *_variants_filtered.vcf.gz"
fi

log_success "${#VCF_FILES[@]} VCF trouvé(s)."

for VCF in "${VCF_FILES[@]}"; do
    echo "    - $VCF"
done

# =============================================================================
# 13 — ANNOTATION DES VARIANTS
# =============================================================================

log_step "13 — ANNOTATION DES VARIANTS"

TOTAL=${#VCF_FILES[@]}
SUCCESS=0
FAILED=0

for VCF in "${VCF_FILES[@]}"; do

    # -------------------------------------------------------------------------
    # Nom de l'échantillon
    # -------------------------------------------------------------------------

    BASENAME=$(basename "$VCF")

    SAMPLE="${BASENAME%_variants_filtered.vcf.gz}"

    OUTPUT_VCF="${ANNOTATED_DIR}/${SAMPLE}_annotated.vcf"

    OUTPUT_VCF_GZ="${OUTPUT_VCF}.gz"

    OUTPUT_TSV="${ANNOTATED_DIR}/${SAMPLE}_annotated.tsv"

    SAMPLE_LOG="${RESULTS_DIR}/${SAMPLE}_snpeff.log"

    # -------------------------------------------------------------------------
    # Affichage
    # -------------------------------------------------------------------------

    echo
    echo "------------------------------------------------------------"
    log_info "Échantillon : $SAMPLE"
    echo "------------------------------------------------------------"

    # -------------------------------------------------------------------------
    # Vérification du VCF
    # -------------------------------------------------------------------------

    if [[ ! -s "$VCF" ]]; then

        log_warning "VCF vide : $VCF"

        FAILED=$((FAILED + 1))

        continue

    fi

    # Vérifier que le VCF peut être lu
    if ! bcftools view -h "$VCF" >/dev/null 2>&1; then

        log_warning "VCF invalide ou illisible : $VCF"

        FAILED=$((FAILED + 1))

        continue

    fi

    N_INPUT=$(bcftools view -H "$VCF" | wc -l)

    log_info "Variants en entrée : $N_INPUT"

    if [[ "$N_INPUT" -eq 0 ]]; then

        log_warning "Le VCF ne contient aucun variant."

        FAILED=$((FAILED + 1))

        continue

    fi

    # -------------------------------------------------------------------------
    # Nettoyage d'anciennes sorties
    # -------------------------------------------------------------------------

    rm -f \
        "$OUTPUT_VCF" \
        "$OUTPUT_VCF_GZ" \
        "${OUTPUT_VCF_GZ}.tbi" \
        "$OUTPUT_TSV"

    # -------------------------------------------------------------------------
    # ANNOTATION SNPEFF
    # -------------------------------------------------------------------------

    log_info "Annotation avec SnpEff..."

    if ! snpEff ann \
        -c "$SNPEFF_CONFIG" \
        -dataDir "$SNPEFF_DATA_ABS" \
        -nodownload \
        -v \
        "$GENOME_ID" \
        "$VCF" \
        > "$OUTPUT_VCF" \
        2> "$SAMPLE_LOG"
    then

        log_warning "SnpEff a échoué pour $SAMPLE."

        echo
        echo "Consultez :"
        echo "    $SAMPLE_LOG"

        FAILED=$((FAILED + 1))

        rm -f "$OUTPUT_VCF"

        continue

    fi

    # -------------------------------------------------------------------------
    # VÉRIFICATION DU VCF ANNOTÉ
    # -------------------------------------------------------------------------

    if [[ ! -s "$OUTPUT_VCF" ]]; then

        log_warning "SnpEff a produit un fichier vide pour $SAMPLE."

        FAILED=$((FAILED + 1))

        rm -f "$OUTPUT_VCF"

        continue

    fi

    # Nombre de variants
    N_OUTPUT=$(grep -v '^#' "$OUTPUT_VCF" | wc -l)

    log_info "Variants en sortie : $N_OUTPUT"

    if [[ "$N_OUTPUT" -eq 0 ]]; then

        log_warning "Le VCF annoté ne contient aucun variant."

        FAILED=$((FAILED + 1))

        rm -f "$OUTPUT_VCF"

        continue

    fi

    # -------------------------------------------------------------------------
    # VÉRIFICATION DE L'ANNOTATION ANN
    # -------------------------------------------------------------------------

    N_ANN=$(grep -c "ANN=" "$OUTPUT_VCF" || true)

    if [[ "$N_ANN" -eq 0 ]]; then

        log_warning "Aucune annotation ANN détectée pour $SAMPLE."

        echo
        echo "Le VCF existe mais SnpEff n'a produit aucune annotation ANN."
        echo "Consultez :"
        echo "    $SAMPLE_LOG"

        FAILED=$((FAILED + 1))

        continue

    fi

    log_success "Annotation ANN détectée : $N_ANN variant(s)."

    # -------------------------------------------------------------------------
    # COMPRESSION
    # -------------------------------------------------------------------------

    log_info "Compression du VCF annoté..."

    bgzip -f "$OUTPUT_VCF"

    # -------------------------------------------------------------------------
    # INDEXATION
    # -------------------------------------------------------------------------

    log_info "Indexation avec tabix..."

    tabix -f -p vcf "$OUTPUT_VCF_GZ"

    # -------------------------------------------------------------------------
    # CRÉATION DU TSV
    #
    # On conserve :
    # CHROM
    # POS
    # REF
    # ALT
    # QUAL
    # FILTER
    # INFO
    #
    # L'information ANN reste disponible dans INFO.
    # -------------------------------------------------------------------------

    log_info "Création du TSV..."

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' \
        "$OUTPUT_VCF_GZ" \
        > "$OUTPUT_TSV"

    # Ajouter l'en-tête
    TMP_TSV="${OUTPUT_TSV}.tmp"

    {
        echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO"
        cat "$OUTPUT_TSV"
    } > "$TMP_TSV"

    mv "$TMP_TSV" "$OUTPUT_TSV"

    # -------------------------------------------------------------------------
    # VÉRIFICATION FINALE
    # -------------------------------------------------------------------------

    if [[ ! -s "$OUTPUT_VCF_GZ" ]]; then

        log_warning "VCF annoté compressé absent."

        FAILED=$((FAILED + 1))

        continue

    fi

    if [[ ! -s "${OUTPUT_VCF_GZ}.tbi" ]]; then

        log_warning "Index Tabix absent."

        FAILED=$((FAILED + 1))

        continue

    fi

    if [[ ! -s "$OUTPUT_TSV" ]]; then

        log_warning "TSV absent."

        FAILED=$((FAILED + 1))

        continue

    fi

    log_success "$SAMPLE — annotation terminée."

    echo
    echo "  VCF annoté :"
    echo "      $OUTPUT_VCF_GZ"

    echo
    echo "  Index :"
    echo "      ${OUTPUT_VCF_GZ}.tbi"

    echo
    echo "  TSV :"
    echo "      $OUTPUT_TSV"

    echo
    echo "  Variants : $N_OUTPUT"
    echo "  ANN       : $N_ANN"

    SUCCESS=$((SUCCESS + 1))

done

# =============================================================================
# 14 — RÉSUMÉ
# =============================================================================

log_step "14 — RÉSUMÉ"

echo
echo "VCF détectés          : $TOTAL"
echo -e "Annotations réussies  : ${GREEN}$SUCCESS${NC}"
echo -e "Annotations échouées  : ${RED}$FAILED${NC}"

echo
echo "Résultats :"
echo "    $ANNOTATED_DIR"

echo
echo "Logs :"
echo "    $RESULTS_DIR"

# =============================================================================
# 15 — VÉRIFICATION FINALE
# =============================================================================

log_step "15 — VÉRIFICATION FINALE"

N_ANNOTATED=$(find "$ANNOTATED_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*_annotated.vcf.gz" \
    | wc -l)

N_TSV=$(find "$ANNOTATED_DIR" \
    -maxdepth 1 \
    -type f \
    -name "*_annotated.tsv" \
    | wc -l)

echo
echo "  VCF annotés : $N_ANNOTATED"
echo "  TSV         : $N_TSV"

# -------------------------------------------------------------------------
# Vérification de cohérence
# -------------------------------------------------------------------------

if [[ "$SUCCESS" -eq "$TOTAL" ]]; then

    echo
    log_success "Toutes les annotations sont terminées avec succès."

elif [[ "$SUCCESS" -gt 0 ]]; then

    echo
    log_warning "Certaines annotations ont réussi et d'autres ont échoué."

else

    echo
    log_error "Aucune annotation n'a réussi.

Consultez les logs :
    $RESULTS_DIR"

fi

# =============================================================================
# FIN
# =============================================================================

log_step "06_annotation.sh terminé"

echo
