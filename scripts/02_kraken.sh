#!/bin/bash
# =============================================================================
# 02_kraken.sh
#
# KRAKEN2 + BRACKEN
# BASE PERSONNALISÉE 100 % LOCALE
#
# OBJECTIFS
#   - aucun génome téléchargé
#   - aucune taxonomie NCBI téléchargée
#   - utilisation uniquement des FASTA locaux (souches_reference/)
#   - attribution explicite des TaxID avec kraken:taxid
#   - création d'une taxonomie minimale locale
#   - construction Kraken2
#   - préparation Bracken
#   - classification des reads
#   - extraction de Borrelia crocidurae
#
# =============================================================================

set -euo pipefail


# =============================================================================
# 0 — RACINE DU PROJET
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"


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

log_info()    { echo -e "${BLUE}[INFO]${NC}    $1"; }
log_success() { echo -e "${GREEN}[OK]${NC}      $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error()   { echo -e "${RED}[ERREUR]${NC}  $1"; exit 1; }

log_step() {
    echo ""
    echo -e "${CYAN}============================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}============================================================${NC}"
}


# =============================================================================
# 1 — CONFIGURATION
# =============================================================================

ENV_NAME="borrelia_pipeline"

FASTA_DIR="data/reference/souches_reference"
KRAKEN_DB="data/kraken_db/borrelia_custom"

CLEAN_DIR="data/clean"
KRAKEN_DIR="kraken"
RESULTS_DIR="results"

THREADS=8

KMER_LENGTH=35
MINIMIZER_LENGTH=31

# -------------------------------------------------------------------------
# TaxID
# -------------------------------------------------------------------------

ROOT_TAXID=1
BACTERIA_TAXID=2

CROCIDURAE_TAXID=29520
ACHEMA_TAXID=1155096
DOU_TAXID=1293575
DUTTONII_TAXID=29521
HERMSII_TAXID=29522
TURICATAE_TAXID=29523
PARKERI_TAXID=29524
RECURRENTIS_TAXID=29525
MIYAMOTOI_TAXID=47466

# -------------------------------------------------------------------------
# Association fichier FASTA → TaxID
# Les clés doivent correspondre EXACTEMENT aux noms de fichiers
# -------------------------------------------------------------------------

declare -A TAXIDS

TAXIDS["Borrelia_crocidurae_03-02_GCF_000825665.2.fasta"]="$CROCIDURAE_TAXID"
TAXIDS["Borrelia_crocidurae_Achema_CP003426.fasta"]="$ACHEMA_TAXID"
TAXIDS["Borrelia_crocidurae_DOU_CP004267.fasta"]="$DOU_TAXID"
TAXIDS["Borrelia duttonii CR2A Contig0001_Eest Africa.fasta"]="$DUTTONII_TAXID"
TAXIDS["Borrelia duttonii Ly_East Africa.fasta"]="$DUTTONII_TAXID"
TAXIDS["Borrelia hermsii DAH_North America.fasta"]="$HERMSII_TAXID"
TAXIDS["Borrelia miyamotoi FR64b_Japon Asia.fasta"]="$MIYAMOTOI_TAXID"
TAXIDS["Borrelia parkeri SLO_North America.fasta"]="$PARKERI_TAXID"
TAXIDS["Borrelia recurrentis A1_Africa.fasta"]="$RECURRENTIS_TAXID"
TAXIDS["Borrelia turicatae 91E135_North America.fasta"]="$TURICATAE_TAXID"

EXPECTED_FASTA=10

# -------------------------------------------------------------------------
# TaxID utilisés pour l'extraction de B. crocidurae
# -------------------------------------------------------------------------

CROCIDURAE_TAXIDS=(
    "$CROCIDURAE_TAXID"
    "$ACHEMA_TAXID"
    "$DOU_TAXID"
)

READ_LENGTH=""


# =============================================================================
# 2 — ENVIRONNEMENT
# =============================================================================

log_step "1 — VÉRIFICATION DE L'ENVIRONNEMENT"

[[ "${CONDA_DEFAULT_ENV:-}" == "$ENV_NAME" ]] || \
    log_error "Activez : conda activate $ENV_NAME"

log_success "Environnement actif : $ENV_NAME"


# =============================================================================
# 3 — OUTILS
# =============================================================================

log_step "2 — VÉRIFICATION DES OUTILS"

TOOLS=(kraken2 kraken2-build bracken bracken-build extract_kraken_reads.py)

for tool in "${TOOLS[@]}"; do
    command -v "$tool" >/dev/null 2>&1 || log_error "Outil introuvable : $tool"
    log_success "$tool : $(command -v "$tool")"
done


# =============================================================================
# 4 — VERSIONS
# =============================================================================

log_step "3 — VERSIONS"

echo "Kraken2 :"; kraken2 --version | head -1
echo ""; echo "Bracken :"; bracken -v 2>&1 | head -1 || true


# =============================================================================
# 5 — DOSSIERS
# =============================================================================

log_step "4 — VÉRIFICATION DES DOSSIERS"

[[ -d "$FASTA_DIR" ]] || log_error "Répertoire FASTA absent : $FASTA_DIR"

mkdir -p "$CLEAN_DIR" "$KRAKEN_DIR" "$RESULTS_DIR"

log_success "Répertoires vérifiés."


# =============================================================================
# 6 — DÉTECTION DES FASTA
# =============================================================================

log_step "5 — DÉTECTION DES GÉNOMES DE RÉFÉRENCE"

mapfile -t FASTA_FILES < <(
    find "$FASTA_DIR" \
        -maxdepth 1 \
        -type f \
        \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) \
        | sort
)

[[ ${#FASTA_FILES[@]} -eq $EXPECTED_FASTA ]] || \
    log_error "Le script attend $EXPECTED_FASTA FASTA. Trouvé : ${#FASTA_FILES[@]}"

log_info "${#FASTA_FILES[@]} FASTA détectés :"
for fasta in "${FASTA_FILES[@]}"; do
    echo "    - $(basename "$fasta")"
done


# =============================================================================
# 7 — VÉRIFICATION DES FASTA
# =============================================================================

log_step "6 — VÉRIFICATION DES HEADERS FASTA"

for fasta in "${FASTA_FILES[@]}"; do

    BASENAME=$(basename "$fasta")

    [[ -n "${TAXIDS[$BASENAME]+x}" ]] || \
        log_error "TaxID non défini pour : $BASENAME"

    HEADER=$(grep '^>' "$fasta" | head -1 || true)

    [[ -n "$HEADER" ]] || log_error "Aucun header FASTA dans : $fasta"

    echo "  $BASENAME"
    echo "  $HEADER"
    echo "  TaxID : ${TAXIDS[$BASENAME]}"
    echo ""

done

log_success "Les $EXPECTED_FASTA génomes sont valides."


# =============================================================================
# 8 — DÉTECTION DES ÉCHANTILLONS
# =============================================================================

log_step "7 — DÉTECTION DES ÉCHANTILLONS"

mapfile -t SAMPLES < <(
    find "$CLEAN_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*_R1_clean.fastq.gz" \
        -printf "%f\n" \
        | sed 's/_R1_clean.fastq.gz$//' \
        | sort
)

[[ ${#SAMPLES[@]} -gt 0 ]] || log_error "Aucun échantillon trouvé dans : $CLEAN_DIR"

log_info "${#SAMPLES[@]} échantillon(s) détecté(s)."
for sample in "${SAMPLES[@]}"; do echo "    - $sample"; done


# =============================================================================
# 9 — VÉRIFICATION R1/R2
# =============================================================================

log_step "8 — VÉRIFICATION DES PAIRES"

for sample in "${SAMPLES[@]}"; do
    [[ -f "$CLEAN_DIR/${sample}_R1_clean.fastq.gz" ]] || \
        log_error "R1 absent : $CLEAN_DIR/${sample}_R1_clean.fastq.gz"
    [[ -f "$CLEAN_DIR/${sample}_R2_clean.fastq.gz" ]] || \
        log_error "R2 absent : $CLEAN_DIR/${sample}_R2_clean.fastq.gz"
done

log_success "Toutes les paires R1/R2 sont présentes."


# =============================================================================
# 10 — LONGUEUR DES READS
# =============================================================================

log_step "9 — DÉTECTION AUTOMATIQUE DE LA LONGUEUR DES READS"

FIRST_SAMPLE="${SAMPLES[0]}"
FIRST_R1="$CLEAN_DIR/${FIRST_SAMPLE}_R1_clean.fastq.gz"
FIRST_R2="$CLEAN_DIR/${FIRST_SAMPLE}_R2_clean.fastq.gz"

detect_read_lengths() {
    local FASTQ="$1"
    python - "$FASTQ" <<'PY'
import sys, gzip
from collections import Counter

counter = Counter()
n = 0

with gzip.open(sys.argv[1], "rt") as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            counter[len(line.strip())] += 1
            n += 1
            if n >= 1000:
                break

for length, count in sorted(counter.items()):
    print(f"{length}\t{count}")
PY
}

log_info "Échantillon utilisé : $FIRST_SAMPLE"
log_info "Analyse de R1..."

R1_LENGTHS=$(detect_read_lengths "$FIRST_R1")
[[ -n "$R1_LENGTHS" ]] || log_error "Impossible de déterminer la longueur R1."

echo "Distribution R1 :"
echo "$R1_LENGTHS" | while IFS=$'\t' read -r length count; do
    echo "    ${length} bp : ${count} reads"
done

READ_LENGTH=$(echo "$R1_LENGTHS" | sort -k2,2nr | head -1 | cut -f1)
log_success "Longueur dominante R1 : ${READ_LENGTH} bp"

log_info "Analyse de R2..."
R2_LENGTHS=$(detect_read_lengths "$FIRST_R2")
R2_LENGTH=$(echo "$R2_LENGTHS" | sort -k2,2nr | head -1 | cut -f1)
log_info "Longueur dominante R2 : ${R2_LENGTH} bp"

[[ "$READ_LENGTH" == "$R2_LENGTH" ]] || \
    log_warning "Longueurs dominantes différentes. R1: ${READ_LENGTH} bp / R2: ${R2_LENGTH} bp"

(( READ_LENGTH >= KMER_LENGTH )) || \
    log_error "Reads trop courts : ${READ_LENGTH} bp < k-mer ${KMER_LENGTH}"

log_success "Longueur retenue pour Bracken : ${READ_LENGTH} bp"


# =============================================================================
# 11 — NETTOYAGE ANCIENNE BASE
# =============================================================================

log_step "10 — PRÉPARATION DE LA BASE KRAKEN2"

if [[ -d "$KRAKEN_DB" ]] && \
   [[ ! -f "$KRAKEN_DB/hash.k2d" || ! -f "$KRAKEN_DB/opts.k2d" || ! -f "$KRAKEN_DB/taxo.k2d" ]]; then
    log_warning "Ancienne base incomplète détectée. Suppression..."
    rm -rf "$KRAKEN_DB"
    log_success "Ancienne base supprimée."
fi

mkdir -p "$KRAKEN_DB"


# =============================================================================
# 12 — TAXONOMIE LOCALE
# =============================================================================

if [[ ! -f "$KRAKEN_DB/taxonomy/nodes.dmp" ]] || \
   [[ ! -f "$KRAKEN_DB/taxonomy/names.dmp" ]]; then

    log_step "10.1 — CRÉATION DE LA TAXONOMIE LOCALE"
    mkdir -p "$KRAKEN_DB/taxonomy"

    cat > "$KRAKEN_DB/taxonomy/nodes.dmp" <<EOF
1	|	1	|	no rank	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
2	|	1	|	superkingdom	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
29520	|	2	|	species	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
1155096	|	29520	|	strain	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
1293575	|	29520	|	strain	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
29521	|	2	|	species	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
29522	|	2	|	species	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
29523	|	2	|	species	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
29524	|	2	|	species	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
29525	|	2	|	species	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
47466	|	2	|	species	|	-	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|	0	|
EOF

    cat > "$KRAKEN_DB/taxonomy/names.dmp" <<EOF
1	|	root	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
29520	|	Borrelia crocidurae	|		|	scientific name	|
1155096	|	Borrelia crocidurae str. Achema	|		|	scientific name	|
1293575	|	Borrelia crocidurae DOU	|		|	scientific name	|
29521	|	Borrelia duttonii	|		|	scientific name	|
29522	|	Borrelia hermsii	|		|	scientific name	|
29523	|	Borrelia turicatae	|		|	scientific name	|
29524	|	Borrelia parkeri	|		|	scientific name	|
29525	|	Borrelia recurrentis	|		|	scientific name	|
47466	|	Borrelia miyamotoi	|		|	scientific name	|
EOF

    log_success "Taxonomie locale créée."

else
    log_success "Taxonomie locale déjà présente."
fi


# =============================================================================
# 13 — PRÉPARATION DES FASTA
# =============================================================================

log_step "10.2 — PRÉPARATION DES FASTA"

PREP_DIR="$KRAKEN_DB/prepared_fasta"
rm -rf "$PREP_DIR"
mkdir -p "$PREP_DIR"

for fasta in "${FASTA_FILES[@]}"; do

    BASENAME=$(basename "$fasta")
    TAXID="${TAXIDS[$BASENAME]}"
    OUTPUT="$PREP_DIR/$BASENAME"

    log_info "Préparation : $BASENAME (TaxID: $TAXID)"

    awk -v taxid="$TAXID" '
        /^>/ {
            header=$0
            sub(/^>/, "", header)
            split(header, fields, /[[:space:]]+/)
            accession=fields[1]
            rest=header
            sub(/^[^[:space:]]+[[:space:]]*/, "", rest)
            if (rest == header) { rest="" }
            if (rest != "") {
                print ">" accession "|kraken:taxid|" taxid " " rest
            } else {
                print ">" accession "|kraken:taxid|" taxid
            }
            next
        }
        { print }
    ' "$fasta" > "$OUTPUT"

    NEW_HEADER=$(grep '^>' "$OUTPUT" | head -1 || true)

    [[ "$NEW_HEADER" == *"kraken:taxid|${TAXID}"* ]] || \
        log_error "Annotation TaxID incorrecte : $NEW_HEADER"

    log_success "$BASENAME préparé."

done


# =============================================================================
# 14 — AJOUT DES FASTA À LA LIBRARY
# =============================================================================

log_step "10.3 — AJOUT DES FASTA À LA LIBRARY"

for fasta in "$PREP_DIR"/*; do
    [[ -f "$fasta" ]] || continue
    log_info "Ajout : $(basename "$fasta")"
    kraken2-build --add-to-library "$fasta" --db "$KRAKEN_DB"
    log_success "$(basename "$fasta") ajouté."
done


# =============================================================================
# 15 — VÉRIFICATION DU MAPPING
# =============================================================================

log_step "10.4 — VÉRIFICATION DU MAPPING TAXONOMIQUE"

for fasta in "$PREP_DIR"/*; do
    COUNT=$(grep -c '^>.*kraken:taxid|' "$fasta" || true)
    [[ "$COUNT" -gt 0 ]] || log_error "Aucun TaxID trouvé dans : $fasta"
    log_success "$(basename "$fasta") : $COUNT séquence(s) annotée(s)."
done


# =============================================================================
# 16 — CONSTRUCTION
# =============================================================================

log_step "10.5 — CONSTRUCTION DE KRAKEN2"

log_info "k-mer     : $KMER_LENGTH"
log_info "minimizer : $MINIMIZER_LENGTH"
log_info "threads   : $THREADS"

kraken2-build \
    --build \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --kmer-len "$KMER_LENGTH" \
    --minimizer-len "$MINIMIZER_LENGTH"

log_success "Base Kraken2 construite."


# =============================================================================
# 17 — VÉRIFICATION BASE
# =============================================================================

log_step "11 — VÉRIFICATION DE LA BASE KRAKEN2"

[[ -f "$KRAKEN_DB/hash.k2d" ]] || log_error "hash.k2d absent."
[[ -f "$KRAKEN_DB/opts.k2d" ]] || log_error "opts.k2d absent."
[[ -f "$KRAKEN_DB/taxo.k2d" ]] || log_error "taxo.k2d absent."

log_success "Base Kraken2 valide."


# =============================================================================
# 18 — PRÉPARATION BRACKEN
# =============================================================================

log_step "12 — PRÉPARATION DE BRACKEN"

BRACKEN_DISTRIBUTION="$KRAKEN_DB/database${READ_LENGTH}mers.kmer_distrib"

if [[ -f "$BRACKEN_DISTRIBUTION" ]]; then
    log_success "Bracken déjà préparé pour ${READ_LENGTH} bp."
else
    log_info "Préparation Bracken (read length: ${READ_LENGTH} bp, k-mer: ${KMER_LENGTH})..."
    bracken-build \
        -d "$KRAKEN_DB" \
        -t "$THREADS" \
        -k "$KMER_LENGTH" \
        -l "$READ_LENGTH"
    [[ -f "$BRACKEN_DISTRIBUTION" ]] || log_error "Bracken n'a pas créé : $BRACKEN_DISTRIBUTION"
    log_success "Bracken préparé."
fi


# =============================================================================
# 19 — CLASSIFICATION
# =============================================================================

log_step "13 — CLASSIFICATION KRAKEN2"

TOTAL=${#SAMPLES[@]}
SUCCESS=0
FAILED=0
START_GLOBAL=$(date +%s)

for SAMPLE in "${SAMPLES[@]}"; do

    log_step "ÉCHANTILLON — $SAMPLE"

    R1="$CLEAN_DIR/${SAMPLE}_R1_clean.fastq.gz"
    R2="$CLEAN_DIR/${SAMPLE}_R2_clean.fastq.gz"
    SAMPLE_DIR="$KRAKEN_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_DIR"

    KRAKEN_FILE="$SAMPLE_DIR/${SAMPLE}.kraken"
    KRAKEN_REPORT="$SAMPLE_DIR/${SAMPLE}_report.txt"
    BRACKEN_FILE="$SAMPLE_DIR/${SAMPLE}_bracken.txt"
    BORRELIA_R1="$SAMPLE_DIR/${SAMPLE}_borrelia_R1.fastq"
    BORRELIA_R2="$SAMPLE_DIR/${SAMPLE}_borrelia_R2.fastq"
    LOG_FILE="$RESULTS_DIR/${SAMPLE}_kraken.log"

    # Kraken2
    log_info "Classification Kraken2..."
    kraken2 \
        --db "$KRAKEN_DB" \
        --paired \
        --threads "$THREADS" \
        --output "$KRAKEN_FILE" \
        --report "$KRAKEN_REPORT" \
        "$R1" "$R2" \
        2> "$LOG_FILE"
    log_success "Kraken2 terminé."

    # Bracken
    log_info "Estimation d'abondance Bracken..."
    bracken \
        -d "$KRAKEN_DB" \
        -i "$KRAKEN_REPORT" \
        -o "$BRACKEN_FILE" \
        -r "$READ_LENGTH" \
        -l S \
        2>> "$LOG_FILE"
    log_success "Bracken terminé."

    # Extraction B. crocidurae
    log_info "Extraction de Borrelia crocidurae..."
    rm -f "$BORRELIA_R1" "$BORRELIA_R2"

    extract_kraken_reads.py \
        -k "$KRAKEN_FILE" \
        -r "$KRAKEN_REPORT" \
        -s "$R1" \
        -s2 "$R2" \
        -t "$CROCIDURAE_TAXID" \
        -t "$ACHEMA_TAXID" \
        -t "$DOU_TAXID" \
        --include-children \
        --fastq-output \
        -o "$BORRELIA_R1" \
        -o2 "$BORRELIA_R2" \
        2>> "$LOG_FILE"

    if [[ ! -f "$BORRELIA_R1" ]] || [[ ! -f "$BORRELIA_R2" ]]; then
        log_warning "Fichiers B. crocidurae non générés."
        FAILED=$((FAILED + 1))
        continue
    fi

    N_READS_R1=$(awk 'END {print int(NR/4)}' "$BORRELIA_R1")
    N_READS_R2=$(awk 'END {print int(NR/4)}' "$BORRELIA_R2")

    log_info "Reads B. crocidurae R1 : $N_READS_R1"
    log_info "Reads B. crocidurae R2 : $N_READS_R2"

    if [[ "$N_READS_R1" -ne "$N_READS_R2" ]]; then
        log_warning "Nombre de reads R1/R2 différent."
        FAILED=$((FAILED + 1))
        continue
    fi

    if [[ "$N_READS_R1" -eq 0 ]]; then
        log_warning "Aucun read B. crocidurae extrait."
        FAILED=$((FAILED + 1))
        continue
    fi

    log_success "$SAMPLE : $N_READS_R1 paires B. crocidurae extraites."
    SUCCESS=$((SUCCESS + 1))

done


# =============================================================================
# 20 — RÉSUMÉ
# =============================================================================

END_GLOBAL=$(date +%s)
DURATION=$((END_GLOBAL - START_GLOBAL))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

log_step "14 — RÉSUMÉ GLOBAL"

echo ""
echo "  Projet                 : $PROJECT_ROOT"
echo "  FASTA références      : $FASTA_DIR"
echo "  Base Kraken2          : $KRAKEN_DB"
echo "  Taxonomie             : locale"
echo "  Génomes téléchargés   : NON"
echo "  Taxonomie téléchargée : NON"
echo "  Longueur reads        : ${READ_LENGTH} bp"
echo "  Threads               : $THREADS"
echo ""
echo "  B. crocidurae extraite via TaxID :"
echo "      29520 (B. crocidurae)"
echo "      1155096 (str. Achema)"
echo "      1293575 (DOU)"
echo ""
echo "  Échantillons          : $TOTAL"
echo -e "  Réussis               : ${GREEN}$SUCCESS${NC}"
echo -e "  Échoués               : ${RED}$FAILED${NC}"
echo "  Durée                 : ${MINUTES} min ${SECONDS} sec"
echo ""
echo "Résultats Kraken2 : $KRAKEN_DIR/<échantillon>/"
echo "Reads B. crocidurae : *_borrelia_R1.fastq / *_borrelia_R2.fastq"
echo "Logs : $RESULTS_DIR/"
echo ""

if [[ "$FAILED" -gt 0 ]]; then
    log_warning "Certains échantillons ont échoué. Consulte : $RESULTS_DIR/"
else
    log_success "Tous les échantillons ont été traités avec succès."
fi
