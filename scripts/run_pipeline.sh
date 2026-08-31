#!/bin/bash
# =============================================================================
# run_all_samples.sh
# PIPELINE COMPLET — Lance tous les scripts dans l'ordre
#
# ⚠️  Lancer depuis la racine du projet :
#     cd borrelia_project/
#     bash scripts/run_all_samples.sh
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

# =============================================================================
# VÉRIFICATIONS
# =============================================================================

ENV_NAME="borrelia_pipeline"

[[ "${CONDA_DEFAULT_ENV:-}" == "$ENV_NAME" ]] || \
    log_error "Activez l'environnement : conda activate $ENV_NAME"

# Vérifier qu'on est bien à la racine du projet
[[ -d "data" && -d "scripts" ]] || \
    log_error "Lancez depuis la racine du projet :
    cd borrelia_project/
    bash scripts/run_all_samples.sh"

# Vérifier que tous les scripts existent
SCRIPTS=(
    "scripts/01_qc.sh"
    "scripts/02_kraken.sh"
    "scripts/03_qc_borrelia.sh"
    "scripts/04_mapping.sh"
    "scripts/05_variants.sh"
    "scripts/06_annotation.sh"
    "scripts/07_genome_phylogeny.sh"
    "scripts/08_locus_detection.sh"
    "scripts/09_locus_phylogeny.sh"
)

for script in "${SCRIPTS[@]}"; do
    [[ -f "$script" ]] || \
        log_error "Script introuvable : $script"
done

log_success "Tous les scripts sont présents"

# =============================================================================
# LANCEMENT DU PIPELINE
# =============================================================================

START_GLOBAL=$(date +%s)

echo ""
echo -e "${CYAN}  Pipeline Borrelia crocidurae${NC}"
echo -e "${CYAN}  Démarré le : $(date)${NC}"
echo ""

run_step() {
    local NUM=$1
    local NOM=$2
    local SCRIPT=$3

    log_step "ÉTAPE $NUM — $NOM"
    log_info "Lancement de $SCRIPT..."

    bash "$SCRIPT"

    if [[ $? -eq 0 ]]; then
        log_success "Étape $NUM terminée : $NOM"
    else
        log_error "Échec à l'étape $NUM : $NOM
        Consultez les logs dans results/"
    fi
}

run_step 1 "Contrôle qualité + fusion des lanes"   "scripts/01_qc.sh"
run_step 2 "Classification Kraken2 + extraction"   "scripts/02_kraken.sh"
run_step 3 "QC des reads Borrelia"                 "scripts/03_qc_borrelia.sh"
run_step 4 "Alignement BWA"                        "scripts/04_mapping.sh"
run_step 5 "Variants + consensus"                  "scripts/05_variants.sh"
run_step 6 "Annotation des génomes"                  "scripts/06_annotation.sh"
run_step 7 "Phylogénie MAFFT + IQ-TREE2"           "scripts/07_genome_phylogeny.sh"
run_step 8 "Détection des locus"                   "scripts/08_locus_detection.sh"
run_step 9 "Phylogénie des locus"                  "scripts/09_locus_phylogeny.sh"

# =============================================================================
# RÉSUMÉ FINAL
# =============================================================================

END_GLOBAL=$(date +%s)
DURATION=$(( (END_GLOBAL - START_GLOBAL) / 60 ))

log_step "PIPELINE TERMINÉ"

echo ""
echo -e "  ${GREEN}✔ Pipeline Borrelia crocidurae terminé avec succès !${NC}"
echo -e "  Durée totale : ${DURATION} minutes"
echo ""
echo "  Résultats :"
echo "  ├── qc/                              → rapports qualité"
echo "  ├── kraken/                          → classification + reads Borrelia"
echo "  ├── mapping/                         → fichiers BAM alignés"
echo "  ├── variants/                        → VCF + consensus FASTA"
echo "  ├── annotation/                      → annotations Prokka"
echo "  ├── phylogeny/borrelia_tree.treefile → arbre phylogénétique génome"
echo "  └── phylogeny/locus_v2/             → arbres phylogénétiques par locus"
echo ""
echo "  Visualiser les arbres :"
echo "  → FigTree local : figtree phylogeny/borrelia_tree.treefile"
echo "  → iTOL en ligne : https://itol.embl.de"
echo ""
