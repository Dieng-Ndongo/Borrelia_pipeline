#!/bin/bash
# =============================================================================
# install.sh
# INSTALLATION / MISE À JOUR — Pipeline Borrelia
#
# Utilisation recommandée :
#
#     cd ~/borrelia_project
#     source scripts/install.sh
#
# Le script :
#   1. Vérifie Conda
#   2. Détecte la racine du projet
#   3. Crée l'environnement s'il n'existe pas
#   4. Met à jour l'environnement depuis environment.yml s'il existe déjà
#   5. Installe KrakenTools
#   6. Vérifie les outils nécessaires
#   7. Active borrelia_pipeline dans le shell courant
#
# IMPORTANT :
# Pour que l'activation reste dans le terminal courant, utiliser :
#
#     source scripts/install.sh
#
# et non :
#
#     bash scripts/install.sh
#
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
    return 1
}

log_step() {
    echo
    echo -e "${CYAN}============================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
}

# =============================================================================
# CONFIGURATION
# =============================================================================

ENV_NAME="borrelia_pipeline"

# =============================================================================
# RACINE DU PROJET
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

ENV_FILE="${PROJECT_ROOT}/environment.yml"

log_step "DÉTECTION DU PROJET"

log_info "Racine du projet    : $PROJECT_ROOT"
log_info "Fichier environment : $ENV_FILE"

if [[ ! -f "$ENV_FILE" ]]; then
    log_error "environment.yml introuvable :
$ENV_FILE"
    exit 1
fi

log_success "environment.yml trouvé"

# =============================================================================
# VÉRIFICATION CONDA
# =============================================================================

log_step "VÉRIFICATION DE CONDA"

if ! command -v conda >/dev/null 2>&1; then

    log_error "Conda n'est pas installé ou n'est pas disponible dans le PATH."

    echo
    echo "Installez Miniconda ou Anaconda avant de continuer."
    echo

    return 1 2>/dev/null || exit 1

fi

CONDA_BASE="$(conda info --base)"

log_success "Conda : $(conda --version)"
log_success "Base  : $CONDA_BASE"

# =============================================================================
# INITIALISATION DE CONDA
# =============================================================================

log_step "INITIALISATION DE CONDA"

# Permet à conda activate de fonctionner dans ce script
eval "$(conda shell.bash hook)"

log_success "Conda initialisé"

# =============================================================================
# CONFIGURATION DE CONDA
# =============================================================================

log_step "CONFIGURATION DE CONDA"

conda config --set channel_priority flexible
conda config --set notify_outdated_conda false

log_success "channel_priority : flexible"
log_success "notify_outdated_conda : false"

# =============================================================================
# CRÉATION OU MISE À JOUR DE L'ENVIRONNEMENT
# =============================================================================

log_step "ENVIRONNEMENT : $ENV_NAME"

ENV_BIN="${CONDA_BASE}/envs/${ENV_NAME}/bin"

# Détecter si l'environnement existe réellement
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then

    # =========================================================================
    # ENVIRONNEMENT EXISTANT
    # =========================================================================

    log_info "L'environnement '$ENV_NAME' existe déjà."

    log_info "Mise à jour depuis :"
    echo "    $ENV_FILE"

    echo
    log_info "Conda va installer les nouveaux outils présents dans environment.yml."
    log_info "Les paquets déjà installés seront conservés."

    conda env update \
        -n "$ENV_NAME" \
        -f "$ENV_FILE" \
        --prune

    log_success "Environnement '$ENV_NAME' mis à jour."

else

    # =========================================================================
    # CRÉATION DE L'ENVIRONNEMENT
    # =========================================================================

    log_info "L'environnement '$ENV_NAME' n'existe pas."

    log_info "Création depuis environment.yml..."
    log_info "Cette opération peut prendre plusieurs minutes."

    conda env create \
        -n "$ENV_NAME" \
        -f "$ENV_FILE"

    log_success "Environnement '$ENV_NAME' créé."

fi

# =============================================================================
# VÉRIFICATION DU DOSSIER BIN
# =============================================================================

if [[ ! -d "$ENV_BIN" ]]; then

    log_error "Le dossier des exécutables de l'environnement est introuvable :

$ENV_BIN"

    return 1 2>/dev/null || exit 1

fi

log_success "Binaires de l'environnement :"
echo "    $ENV_BIN"

# =============================================================================
# KRAKENTOOLS
# =============================================================================

log_step "INSTALLATION DE KRAKENTOOLS"

KRAKENTOOLS_DIR="${PROJECT_ROOT}/scripts/krakentools"

if [[ -d "$KRAKENTOOLS_DIR" ]]; then

    log_info "KrakenTools est déjà présent."

else

    if ! command -v git >/dev/null 2>&1; then

        log_error "Git est nécessaire pour télécharger KrakenTools."

        return 1 2>/dev/null || exit 1

    fi

    log_info "Téléchargement de KrakenTools depuis GitHub..."

    git clone \
        https://github.com/jenniferlu717/KrakenTools.git \
        "$KRAKENTOOLS_DIR"

    log_success "KrakenTools téléchargé."

fi

# =============================================================================
# PERMISSIONS KRAKENTOOLS
# =============================================================================

chmod +x "${KRAKENTOOLS_DIR}"/*.py

log_success "Scripts KrakenTools rendus exécutables."

# =============================================================================
# INTÉGRATION KRAKENTOOLS
# =============================================================================

log_step "INTÉGRATION DE KRAKENTOOLS"

KRAKEN_TOOLS=(
    extract_kraken_reads.py
    kreport2mpa.py
    combine_kreports.py
)

for tool in "${KRAKEN_TOOLS[@]}"; do

    SOURCE_TOOL="${KRAKENTOOLS_DIR}/${tool}"
    TARGET_TOOL="${ENV_BIN}/${tool}"

    if [[ -f "$SOURCE_TOOL" ]]; then

        ln -sf "$SOURCE_TOOL" "$TARGET_TOOL"

        log_success "$tool → environnement"

    else

        log_warning "$tool introuvable dans KrakenTools"

    fi

done

# =============================================================================
# ACTIVATION DE L'ENVIRONNEMENT
# =============================================================================

log_step "ACTIVATION DE L'ENVIRONNEMENT"

conda activate "$ENV_NAME"

log_success "Environnement activé : $CONDA_DEFAULT_ENV"

# Vérification
if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then

    log_error "Impossible d'activer l'environnement '$ENV_NAME'."

    return 1 2>/dev/null || exit 1

fi

# =============================================================================
# VÉRIFICATION DES OUTILS
# =============================================================================

log_step "VÉRIFICATION DES OUTILS"

# IMPORTANT :
#
# Dans ton environnement actuel :
#
#     which iqtree
#
# donne :
#
#     /home/ndongo/anaconda3/envs/borrelia_pipeline/bin/iqtree
#
# Donc on vérifie "iqtree" et non "iqtree2".
#
# =============================================================================

TOOLS=(
    fastqc
    fastp
    multiqc
    kraken2
    bracken
    bwa
    samtools
    bcftools
    bgzip
    tabix
    mafft
    iqtree
    efetch
    snpEff
    extract_kraken_reads.py
)

ALL_OK=true

for tool in "${TOOLS[@]}"; do

    if command -v "$tool" >/dev/null 2>&1; then

        TOOL_PATH="$(command -v "$tool")"

        log_success "$tool ✓"
        echo "              $TOOL_PATH"

    else

        log_warning "$tool ✗ — introuvable"

        ALL_OK=false

    fi

done

# =============================================================================
# INFORMATIONS DE VERSION
# =============================================================================

log_step "VERSIONS DES OUTILS PRINCIPAUX"

echo

if command -v fastqc >/dev/null 2>&1; then
    fastqc --version 2>&1 | head -1 || true
fi

if command -v fastp >/dev/null 2>&1; then
    fastp --version 2>&1 | head -1 || true
fi

if command -v kraken2 >/dev/null 2>&1; then
    kraken2 --version 2>&1 | head -1 || true
fi

if command -v bwa >/dev/null 2>&1; then
    bwa 2>&1 | head -1 || true
fi

if command -v samtools >/dev/null 2>&1; then
    samtools --version 2>&1 | head -1 || true
fi

if command -v bcftools >/dev/null 2>&1; then
    bcftools --version 2>&1 | head -1 || true
fi

if command -v mafft >/dev/null 2>&1; then
    mafft --version 2>&1 | head -1 || true
fi

if command -v iqtree >/dev/null 2>&1; then
    iqtree -version 2>&1 | head -1 || true
fi

if command -v snpEff >/dev/null 2>&1; then
    snpEff -version 2>&1 | head -1 || true
fi

# =============================================================================
# RÉSUMÉ
# =============================================================================

log_step "RÉSUMÉ DE L'INSTALLATION"

echo
echo "Projet :"
echo "    $PROJECT_ROOT"
echo

echo "Environnement :"
echo "    $ENV_NAME"
echo

echo "Environnement actif :"
echo "    ${CONDA_DEFAULT_ENV:-aucun}"
echo

echo "Environment.yml :"
echo "    $ENV_FILE"
echo

echo "Binaires :"
echo "    $ENV_BIN"
echo

# =============================================================================
# SUCCÈS
# =============================================================================

if [[ "$ALL_OK" == true ]]; then

    echo -e "${GREEN}============================================================${NC}"
    echo -e "${GREEN}  INSTALLATION / MISE À JOUR RÉUSSIE${NC}"
    echo -e "${GREEN}============================================================${NC}"

    echo
    echo -e "${GREEN}Tous les outils nécessaires sont disponibles.${NC}"
    echo

    echo "Environnement actif :"
    echo -e "    ${GREEN}$CONDA_DEFAULT_ENV${NC}"
    echo

    echo "Projet prêt à être exécuté."
    echo

    echo "Exemples :"
    echo
    echo "    bash scripts/06_annotation.sh"
    echo "    bash scripts/07_genome_phylogeny.sh"
    echo

else

    echo -e "${YELLOW}============================================================${NC}"
    echo -e "${YELLOW}  ATTENTION : OUTILS MANQUANTS${NC}"
    echo -e "${YELLOW}============================================================${NC}"

    echo
    echo "Certains outils ne sont pas disponibles."
    echo
    echo "Vérifiez environment.yml puis relancez :"
    echo
    echo "    source scripts/install.sh"
    echo

fi

# =============================================================================
# FIN
# =============================================================================

echo
log_success "Installation terminée."
echo

# =============================================================================
# IMPORTANT
# =============================================================================
#
# Si le script est lancé avec :
#
#     source scripts/install.sh
#
# l'environnement reste activé dans le terminal courant.
#
# Si lancé avec :
#
#     bash scripts/install.sh
#
# l'environnement sera activé uniquement dans le sous-shell du script
# et ne restera PAS actif après sa fermeture.
#
# =============================================================================
