#!/bin/bash
# =============================================================================
# 08_locus_detection.sh
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
log_step()    { echo; echo -e "${CYAN}============================================================${NC}"; echo -e "${CYAN}  $1${NC}"; echo -e "${CYAN}============================================================${NC}"; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

ENV_NAME="borrelia_pipeline"
GENBANK="data/reference/annotation/CP003426.1.gb"
OUTPUT_DIR="results/locus_detection"
GLOBAL_TSV="${OUTPUT_DIR}/loci_detected.tsv"
LOG_FILE="${OUTPUT_DIR}/locus_detection.log"

mkdir -p "$OUTPUT_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

log_info "Racine du projet : $PROJECT_ROOT"

log_step "1 — VÉRIFICATION DE L'ENVIRONNEMENT"
if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
    log_error "L'environnement '$ENV_NAME' n'est pas actif. Activez-le avec : conda activate $ENV_NAME"
fi
log_success "Environnement actif : $ENV_NAME"

log_step "2 — VÉRIFICATION DES OUTILS"
for tool in awk sed grep; do
    if command -v "$tool" >/dev/null 2>&1; then
        log_success "$tool : $(command -v "$tool")"
    else
        log_error "Outil '$tool' introuvable."
    fi
done

log_step "3 — VÉRIFICATION DU GENBANK"
[[ ! -f "$GENBANK" ]] && log_error "GenBank introuvable : $GENBANK"
[[ ! -r "$GENBANK" ]] && log_error "GenBank non lisible : $GENBANK"
log_success "GenBank trouvé : $GENBANK"

log_step "4 — IDENTIFICATION DE LA RÉFÉRENCE"
REFERENCE_ID=$(grep -m1 '^LOCUS' "$GENBANK" | awk '{print $2}')
[[ -z "$REFERENCE_ID" ]] && log_error "Impossible de déterminer l'identifiant LOCUS."
log_success "Identifiant GenBank : $REFERENCE_ID"

log_step "5 — PRÉPARATION DES FICHIERS"
rm -f "$GLOBAL_TSV" "$OUTPUT_DIR/16S.tsv" "$OUTPUT_DIR/23S.tsv" "$OUTPUT_DIR/5S.tsv" "$OUTPUT_DIR/flaB.tsv"
printf "reference\tlocus\tfeature\tstart\tend\tstrand\tgene\tproduct\tlabel\n" > "$GLOBAL_TSV"
for locus in 16S 23S 5S flaB; do
    printf "reference\tlocus\tfeature\tstart\tend\tstrand\tgene\tproduct\tlabel\n" > "$OUTPUT_DIR/${locus}.tsv"
done
log_success "Fichiers de sortie initialisés."

log_step "6 — ANALYSE DES FEATURES GENBANK"
log_info "Lecture de l'annotation GenBank..."

# Teste d'abord que la détection fonctionne avec l'awk corrigé
awk '
BEGIN { in_features=0; feature=""; location=""; gene=""; product="" }
/^FEATURES/ { in_features=1; next }
/^ORIGIN/   { in_features=0 }
!in_features { next }
/^[[:space:]]{5}[A-Za-z0-9_]+[[:space:]]+/ {
    if (feature != "") print feature "\t" location "\t" gene "\t" product
    feature=$1; location=""
    line=$0; sub(/^[[:space:]]{5}[A-Za-z0-9_]+[[:space:]]+/,"",line); location=line
    gene=""; product=""; next
}
/\/gene="/    { line=$0; sub(/^.*\/gene="/,"",line);    sub(/".*$/,"",line); gene=line;    next }
/\/product="/  { line=$0; sub(/^.*\/product="/,"",line); sub(/".*$/,"",line); product=line; next }
END { if (feature != "") print feature "\t" location "\t" gene "\t" product }
' data/reference/annotation/CP003426.1.gb > results/locus_detection/.features.tmp

wc -l results/locus_detection/.features.tmp
grep -i 'rRNA\|flaB' results/locus_detection/.features.tmp

log_step "7 — IDENTIFICATION DES LOCI"
TOTAL=0; COUNT_16S=0; COUNT_23S=0; COUNT_5S=0; COUNT_FLAB=0

while IFS=$'\t' read -r FEATURE LOCATION GENE PRODUCT; do
    [[ -z "$FEATURE" ]] && continue
    SEARCH_TEXT="${FEATURE} ${GENE} ${PRODUCT}"
    LOCUS=""
    echo "$SEARCH_TEXT" | grep -Eiq '(^|[^0-9])16S([^0-9]|$)|16S rRNA|16S ribosomal RNA|small subunit ribosomal RNA' && LOCUS="16S"
    [[ -z "$LOCUS" ]] && echo "$SEARCH_TEXT" | grep -Eiq '(^|[^0-9])23S([^0-9]|$)|23S rRNA|23S ribosomal RNA|large subunit ribosomal RNA' && LOCUS="23S"
    [[ -z "$LOCUS" ]] && echo "$SEARCH_TEXT" | grep -Eiq '(^|[^0-9])5S([^0-9]|$)|5S rRNA|5S ribosomal RNA' && LOCUS="5S"
    [[ -z "$LOCUS" ]] && echo "$SEARCH_TEXT" | grep -Eiq '(^|[^[:alnum:]])flaB([^[:alnum:]]|$)|flagellin' && LOCUS="flaB"
    [[ -z "$LOCUS" ]] && continue

    STRAND="+"
    echo "$LOCATION" | grep -q 'complement' && STRAND="-"

    START=$(echo "$LOCATION" | grep -oE '[0-9]+' | head -1 || true)
    END=$(echo   "$LOCATION" | grep -oE '[0-9]+' | tail -1 || true)
    [[ -z "$START" || -z "$END" ]] && { log_warning "Coordonnées impossibles : $LOCATION"; continue; }

    LABEL="${PRODUCT:-${GENE:-$LOCUS}}"
    LINE="${REFERENCE_ID}\t${LOCUS}\t${FEATURE}\t${START}\t${END}\t${STRAND}\t${GENE}\t${PRODUCT}\t${LABEL}"
    printf "%b\n" "$LINE" >> "$GLOBAL_TSV"
    printf "%b\n" "$LINE" >> "$OUTPUT_DIR/${LOCUS}.tsv"
    TOTAL=$((TOTAL+1))
    case "$LOCUS" in
        16S)  COUNT_16S=$((COUNT_16S+1))   ;;
        23S)  COUNT_23S=$((COUNT_23S+1))   ;;
        5S)   COUNT_5S=$((COUNT_5S+1))     ;;
        flaB) COUNT_FLAB=$((COUNT_FLAB+1)) ;;
    esac
done < "$OUTPUT_DIR/.features.tmp"
rm -f "$OUTPUT_DIR/.features.tmp"

log_step "8 — TRI DES COORDONNÉES"
{ head -n1 "$GLOBAL_TSV"; tail -n+2 "$GLOBAL_TSV" | sort -k2,2 -k4,4n; } > "${GLOBAL_TSV}.tmp" && mv "${GLOBAL_TSV}.tmp" "$GLOBAL_TSV"
for locus in 16S 23S 5S flaB; do
    FILE="$OUTPUT_DIR/${locus}.tsv"
    { head -n1 "$FILE"; tail -n+2 "$FILE" | sort -k4,4n; } > "${FILE}.tmp" && mv "${FILE}.tmp" "$FILE"
done
log_success "Coordonnées triées."

log_step "9 — RÉSUMÉ DE LA DÉTECTION"
echo
echo "Référence : $REFERENCE_ID"
echo
printf "  %-8s : %s annotation(s)\n" "16S"  "$COUNT_16S"
printf "  %-8s : %s annotation(s)\n" "23S"  "$COUNT_23S"
printf "  %-8s : %s annotation(s)\n" "5S"   "$COUNT_5S"
printf "  %-8s : %s annotation(s)\n" "flaB" "$COUNT_FLAB"
echo
echo "  Total   : $TOTAL annotation(s)"
echo

log_step "10 — COORDONNÉES DÉTECTÉES"
if [[ "$TOTAL" -gt 0 ]]; then
    column -t -s $'\t' "$GLOBAL_TSV"
else
    log_error "Aucun locus d'intérêt détecté dans le GenBank."
fi

log_step "11 — CONTRÔLE DES LOCI PRINCIPAUX"
[[ "$COUNT_16S"  -gt 0 ]] && log_success "16S détecté."  || log_warning "16S non détecté."
[[ "$COUNT_23S"  -gt 0 ]] && log_success "23S détecté."  || log_warning "23S non détecté."
[[ "$COUNT_5S"   -gt 0 ]] && log_success "5S détecté."   || log_warning "5S non détecté."
[[ "$COUNT_FLAB" -gt 0 ]] && log_success "flaB détecté." || log_warning "flaB non détecté dans l'annotation GenBank."

log_step "12 — FICHIERS PRODUITS"
echo
echo "  Fichier principal : $GLOBAL_TSV"
echo "  Fichiers par locus : $OUTPUT_DIR/{16S,23S,5S,flaB}.tsv"
echo "  Log : $LOG_FILE"
echo

log_step "08_locus_detection.sh terminé"
log_success "Détection des loci terminée avec succès."
echo
echo "Les coordonnées seront utilisées par :"
echo "    09_locus_phylogeny.sh"
echo
