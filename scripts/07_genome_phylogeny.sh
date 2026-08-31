#!/bin/bash
# =============================================================================
# 06_phylogeny.sh
#
# PHYLOGÉNIE GÉNOMIQUE GLOBALE DE BORRELIA
#
# Étapes :
#   1. Détection des consensus
#   2. Détection des références
#   3. Fusion des séquences
#   4. Contrôle qualité des FASTA
#   5. Contrôle des longueurs
#   6. Vérification des séquences identiques
#   7. Alignement multiple avec MAFFT
#   8. Contrôle de l'alignement
#   9. Construction de l'arbre avec IQ-TREE2
#  10. Enracinement avec B. miyamotoi FR64b
#  11. Génération des métadonnées
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
    exit 1
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

# -------------------------------------------------------------------------
# Entrées
# -------------------------------------------------------------------------

VARIANTS_DIR="variants"

REF_SOUCHES="data/reference/souches_reference"

# -------------------------------------------------------------------------
# Organisation des résultats
# -------------------------------------------------------------------------

PHYLOGENY_DIR="phylogeny"

GENOME_DIR="${PHYLOGENY_DIR}/genome"

INPUT_DIR="${GENOME_DIR}/input"

ALIGNMENT_DIR="${GENOME_DIR}/alignment"

QC_DIR="${GENOME_DIR}/qc"

METADATA_DIR="${GENOME_DIR}/metadata"

TREE_DIR="${GENOME_DIR}/tree"

# -------------------------------------------------------------------------
# Paramètres
# -------------------------------------------------------------------------

THREADS=8

BOOTSTRAP=1000

# -------------------------------------------------------------------------
# OUTGROUP
#
# Identifiant EXACT présent dans le fichier FASTA.
#
# >B_miyamotoi|Asie_Japon|FR64b_CP004217
#
# -------------------------------------------------------------------------

OUTGROUP_ID="B_miyamotoi|Asie_Japon|FR64b_CP004217"

# =============================================================================
# 1. VÉRIFICATION DE L'ENVIRONNEMENT
# =============================================================================

log_step "1 — VÉRIFICATION DE L'ENVIRONNEMENT"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then

    log_error "Environnement Conda incorrect.

Activez l'environnement avec :

    conda activate $ENV_NAME"

fi

log_success "Environnement actif : $ENV_NAME"

# =============================================================================
# 2. VÉRIFICATION DES OUTILS
# =============================================================================

log_step "2 — VÉRIFICATION DES OUTILS"

for tool in mafft iqtree awk grep sed sort uniq sha256sum column; do

    if command -v "$tool" >/dev/null 2>&1; then

        log_success "$tool : $(command -v "$tool")"

    else

        log_error "Outil '$tool' introuvable."

    fi

done

# =============================================================================
# 3. PRÉPARATION DES DOSSIERS
# =============================================================================

log_step "3 — PRÉPARATION DES DOSSIERS"

[[ -d "$VARIANTS_DIR" ]] || \
    log_error "Dossier introuvable : $VARIANTS_DIR"

[[ -d "$REF_SOUCHES" ]] || \
    log_error "Dossier introuvable : $REF_SOUCHES"

mkdir -p \
    "$INPUT_DIR" \
    "$ALIGNMENT_DIR" \
    "$QC_DIR" \
    "$METADATA_DIR" \
    "$TREE_DIR"

log_success "Structure phylogeny/genome/ créée"

# =============================================================================
# 4. DÉTECTION DES CONSENSUS
# =============================================================================

log_step "4 — DÉTECTION DES CONSENSUS"

CONSENSUS=()

while IFS= read -r file; do
    CONSENSUS+=("$file")
done < <(
    find "$VARIANTS_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*_consensus.fasta" \
        | sort
)

[[ ${#CONSENSUS[@]} -gt 0 ]] || \
    log_error "Aucun consensus trouvé dans $VARIANTS_DIR/"

log_success "${#CONSENSUS[@]} consensus détecté(s)"

for file in "${CONSENSUS[@]}"; do
    echo "    - $(basename "$file")"
done

# =============================================================================
# 5. DÉTECTION DES RÉFÉRENCES
# =============================================================================

log_step "5 — DÉTECTION DES SOUCHES DE RÉFÉRENCE"

REFERENCES=()

while IFS= read -r file; do
    REFERENCES+=("$file")
done < <(
    find "$REF_SOUCHES" \
        -maxdepth 1 \
        -type f \
        \( \
            -name "*.fasta" \
            -o -name "*.fa" \
            -o -name "*.fas" \
            -o -name "*.fna" \
        \) \
        | sort
)

[[ ${#REFERENCES[@]} -gt 0 ]] || \
    log_error "Aucune référence trouvée dans $REF_SOUCHES/"

log_success "${#REFERENCES[@]} référence(s) détectée(s)"

for file in "${REFERENCES[@]}"; do
    echo "    - $(basename "$file")"
done

# =============================================================================
# 6. FUSION DES SÉQUENCES
# =============================================================================

log_step "6 — CONSTRUCTION DU JEU DE DONNÉES"

ALL_SEQUENCES="${INPUT_DIR}/all_sequences.fasta"

rm -f "$ALL_SEQUENCES"

# Références d'abord
for file in "${REFERENCES[@]}"; do
    cat "$file" >> "$ALL_SEQUENCES"
done

# Consensus ensuite
for file in "${CONSENSUS[@]}"; do
    cat "$file" >> "$ALL_SEQUENCES"
done

[[ -s "$ALL_SEQUENCES" ]] || \
    log_error "Le fichier FASTA global est vide."

log_success "Jeu de données créé"

echo
echo "Fichier :"
echo "  $ALL_SEQUENCES"

# =============================================================================
# 7. COMPTAGE DES SÉQUENCES
# =============================================================================

log_step "7 — COMPTAGE DES SÉQUENCES"

N_SEQ=$(grep -c "^>" "$ALL_SEQUENCES" || true)

[[ "$N_SEQ" -ge 2 ]] || \
    log_error "Moins de 2 séquences détectées."

log_info "Nombre total de séquences : $N_SEQ"

# =============================================================================
# 8. AFFICHAGE DES SÉQUENCES
# =============================================================================

log_step "8 — SÉQUENCES INCLUSES"

grep "^>" "$ALL_SEQUENCES" |
    sed 's/^>/    /'

# =============================================================================
# 9. VÉRIFICATION DES IDENTIFIANTS
# =============================================================================

log_step "9 — VÉRIFICATION DES IDENTIFIANTS"

DUPLICATE_HEADERS=$(
    grep "^>" "$ALL_SEQUENCES" |
    sed 's/^>//' |
    sort |
    uniq -d
)

if [[ -n "$DUPLICATE_HEADERS" ]]; then

    log_error "Noms de séquences dupliqués :

$DUPLICATE_HEADERS"

fi

log_success "Aucun nom de séquence dupliqué"

# =============================================================================
# 10. CONTRÔLE DES SÉQUENCES
# =============================================================================

log_step "10 — CONTRÔLE QUALITÉ DES SÉQUENCES"

QC_FILE="${QC_DIR}/sequence_qc.tsv"

echo -e \
"Sequence\tLength\tN\tAmbiguous_IUPAC\tInvalid_characters\tStatus" \
> "$QC_FILE"


check_fasta_sequence() {

    local file="$1"

    local sequence_name
    local sequence
    local length
    local n_count
    local ambiguous_count
    local invalid_count
    local status

    sequence_name=$(
        grep "^>" "$file" |
        head -1 |
        sed 's/^>//'
    )

    sequence=$(
        grep -v "^>" "$file" |
        tr -d '\n\r'
    )

    length=${#sequence}

    n_count=$(
        printf '%s' "$sequence" |
        grep -o -i "N" |
        wc -l || true
    )

    ambiguous_count=$(
        printf '%s' "$sequence" |
        grep -o -i "[RYSWKMBDHV]" |
        wc -l || true
    )

    invalid_count=$(
        printf '%s' "$sequence" |
        grep -o -i "[^ACGTNRYWSKMBDHV]" |
        wc -l || true
    )

    status="OK"

    if [[ "$length" -eq 0 ]]; then
        status="ERROR"
    fi

    if [[ "$invalid_count" -gt 0 ]]; then
        status="ERROR"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sequence_name" \
        "$length" \
        "$n_count" \
        "$ambiguous_count" \
        "$invalid_count" \
        "$status" \
        >> "$QC_FILE"
}


# Contrôle références

for file in "${REFERENCES[@]}"; do

    log_info "Contrôle : $(basename "$file")"

    check_fasta_sequence "$file"

done


# Contrôle consensus

for file in "${CONSENSUS[@]}"; do

    log_info "Contrôle : $(basename "$file")"

    check_fasta_sequence "$file"

done


echo

column -t -s $'\t' "$QC_FILE"

# =============================================================================
# 11. VÉRIFICATION DES ERREURS
# =============================================================================

if awk -F'\t' \
    'NR > 1 && $6 == "ERROR" {found=1}
     END {exit !found}' \
    "$QC_FILE"; then

    log_error "Une ou plusieurs séquences contiennent des caractères invalides."

fi

log_success "Toutes les séquences FASTA sont valides"

# =============================================================================
# 12. VÉRIFICATION DES LONGUEURS
# =============================================================================

log_step "12 — VÉRIFICATION DES LONGUEURS"

REFERENCE_LENGTH=$(
    awk -F'\t' '
    NR == 2 {print $2}
    ' "$QC_FILE"
)

log_info "Longueur de référence attendue : $REFERENCE_LENGTH bp"

LENGTH_ERROR=0

while IFS=$'\t' read -r sequence length n ambiguous invalid status; do

    [[ "$sequence" == "Sequence" ]] && continue

    if [[ "$length" -ne "$REFERENCE_LENGTH" ]]; then

        log_warning "$sequence : $length bp"

        LENGTH_ERROR=1

    else

        log_success "$sequence : $length bp"

    fi

done < "$QC_FILE"


if [[ "$LENGTH_ERROR" -eq 1 ]]; then

    log_warning "Certaines séquences n'ont pas exactement la longueur de la référence."

    log_warning "Elles seront conservées : MAFFT peut gérer ces différences."

else

    log_success "Toutes les séquences ont la même longueur"

fi

# =============================================================================
# 13. RECHERCHE DES SÉQUENCES IDENTIQUES
# =============================================================================

log_step "13 — VÉRIFICATION DES SÉQUENCES DUPLIQUÉES"

HASH_FILE="${QC_DIR}/sequence_hashes.tsv"

rm -f "$HASH_FILE"

while IFS= read -r header; do

    sequence_name=$(echo "$header" | sed 's/^>//')

    sequence=$(
        awk -v name="$header" '
            $0 == name {found=1; next}
            /^>/ && found {exit}
            found {printf "%s", $0}
        ' "$ALL_SEQUENCES"
    )

    hash=$(printf '%s' "$sequence" |
        sha256sum |
        awk '{print $1}')

    echo -e "${hash}\t${sequence_name}" >> "$HASH_FILE"

done < <(grep "^>" "$ALL_SEQUENCES")


DUPLICATE_HASHES=$(
    cut -f1 "$HASH_FILE" |
    sort |
    uniq -d
)

if [[ -n "$DUPLICATE_HASHES" ]]; then

    log_warning "Des séquences identiques ont été détectées."

    while IFS= read -r hash; do

        echo
        grep "^${hash}" "$HASH_FILE"

    done <<< "$DUPLICATE_HASHES"

    log_warning "Les séquences identiques sont conservées."

else

    log_success "Aucune séquence identique détectée"

fi

# =============================================================================
# 14. ALIGNEMENT MULTIPLE — MAFFT
# =============================================================================

log_step "14 — ALIGNEMENT MULTIPLE AVEC MAFFT"

ALIGNED="${ALIGNMENT_DIR}/all_sequences_aligned.fasta"

rm -f "$ALIGNED"

mafft \
    --auto \
    --thread "$THREADS" \
    "$ALL_SEQUENCES" \
    > "$ALIGNED" \
    2> "${ALIGNMENT_DIR}/mafft.log"

[[ -s "$ALIGNED" ]] || \
    log_error "L'alignement MAFFT est vide."

log_success "Alignement MAFFT terminé"

echo
echo "Alignement :"
echo "  $ALIGNED"

echo
echo "Log MAFFT :"
echo "  ${ALIGNMENT_DIR}/mafft.log"

# =============================================================================
# 15. CONTRÔLE DE L'ALIGNEMENT
# =============================================================================

log_step "15 — CONTRÔLE DE L'ALIGNEMENT"

ALIGNMENT_QC="${QC_DIR}/alignment_qc.txt"

ALIGNMENT_SEQ=$(grep -c "^>" "$ALIGNED" || true)

ALIGNMENT_LENGTH=$(
    awk '
    /^>/ {
        if (seq != "") {
            print length(seq)
            exit
        }

        seq=""
        next
    }

    {
        seq=seq $0
    }

    END {
        if (seq != "")
            print length(seq)
    }
    ' "$ALIGNED"
)


TOTAL_GAPS=$(
    grep -v "^>" "$ALIGNED" |
    tr -d '\n\r' |
    grep -o "-" |
    wc -l || true
)


TOTAL_AMBIGUOUS=$(
    grep -v "^>" "$ALIGNED" |
    tr -d '\n\r' |
    grep -o -i "[NRYWSKMBDHV]" |
    wc -l || true
)


TOTAL_BASES=$(
    grep -v "^>" "$ALIGNED" |
    tr -d '\n\r' |
    wc -c
)


{
    echo "============================================================"
    echo "CONTRÔLE DE L'ALIGNEMENT MAFFT"
    echo "============================================================"
    echo
    echo "Nombre de séquences : $ALIGNMENT_SEQ"
    echo "Longueur alignement : $ALIGNMENT_LENGTH"
    echo "Gaps totaux         : $TOTAL_GAPS"
    echo "Bases ambiguës      : $TOTAL_AMBIGUOUS"
    echo "Caractères totaux   : $TOTAL_BASES"
    echo
} > "$ALIGNMENT_QC"


if [[ "$ALIGNMENT_SEQ" -ne "$N_SEQ" ]]; then

    log_error "Le nombre de séquences a changé après MAFFT."

fi

log_success "Nombre de séquences conservé : $ALIGNMENT_SEQ"

log_info "Longueur alignement : $ALIGNMENT_LENGTH bp"

log_info "Gaps totaux : $TOTAL_GAPS"

log_info "Bases ambiguës : $TOTAL_AMBIGUOUS"

# =============================================================================
# 16. VÉRIFICATION DE L'OUTGROUP
# =============================================================================

log_step "16 — VÉRIFICATION DE L'OUTGROUP"

log_info "Outgroup configuré :"
echo
echo "    $OUTGROUP_ID"
echo


# Vérification exacte de l'en-tête

OUTGROUP_COUNT=$(
    grep "^>" "$ALL_SEQUENCES" |
    sed 's/^>//' |
    awk -v id="$OUTGROUP_ID" '
        $0 == id {count++}
        END {print count+0}
    '
)


if [[ "$OUTGROUP_COUNT" -eq 0 ]]; then

    log_error "Outgroup introuvable dans :

$ALL_SEQUENCES

Identifiant recherché :

$OUTGROUP_ID

Vérifiez l'en-tête FASTA."

fi


if [[ "$OUTGROUP_COUNT" -gt 1 ]]; then

    log_error "L'outgroup apparaît plusieurs fois :

$OUTGROUP_ID"

fi


log_success "Outgroup trouvé exactement une fois"

# =============================================================================
# 17. CRÉATION DU FICHIER DE MÉTADONNÉES
# =============================================================================

log_step "17 — CRÉATION DES MÉTADONNÉES"

METADATA_FILE="${METADATA_DIR}/samples_metadata.tsv"

echo -e \
"ID\tSpecies\tLocation\tStrain_or_Sample\tType" \
> "$METADATA_FILE"


while IFS= read -r header; do

    ID="${header#>}"

    # ---------------------------------------------------------------------
    # Extraire les champs séparés par |
    #
    # Exemple :
    #
    # B_crocidurae|Senegal|MAD-G2-0155_S84
    #
    # ---------------------------------------------------------------------

    SPECIES=$(echo "$ID" | cut -d'|' -f1)

    LOCATION=$(echo "$ID" | cut -d'|' -f2)

    STRAIN=$(echo "$ID" | cut -d'|' -f3-)

    # Type
    if [[ "$ID" == *"MAD-G2"* ||
          "$ID" == *"DEG-G1"* ]]; then

        TYPE="Sample"

    elif [[ "$ID" == "$OUTGROUP_ID" ]]; then

        TYPE="Outgroup"

    else

        TYPE="Reference"

    fi

    echo -e \
    "${ID}\t${SPECIES}\t${LOCATION}\t${STRAIN}\t${TYPE}" \
    >> "$METADATA_FILE"

done < <(grep "^>" "$ALL_SEQUENCES")


log_success "Métadonnées créées"

echo
echo "Fichier :"
echo "  $METADATA_FILE"

# =============================================================================
# 18. CONSTRUCTION DE L'ARBRE — IQ-TREE
# =============================================================================

log_step "18 — CONSTRUCTION DE L'ARBRE AVEC IQ-TREE"

TREE_PREFIX="${TREE_DIR}/borrelia_tree"

log_info "Modèle : MFP"

log_info "Bootstrap ultrarapide : $BOOTSTRAP"

log_info "Threads : $THREADS"

log_info "Outgroup : $OUTGROUP_ID"

iqtree \
    -s "$ALIGNED" \
    -m MFP \
    -B "$BOOTSTRAP" \
    -T "$THREADS" \
    -o "$OUTGROUP_ID" \
    --prefix "$TREE_PREFIX" \
    --redo \
    > "${TREE_PREFIX}.log" \
    2>&1


# =============================================================================
# 19. VÉRIFICATION DE L'ARBRE
# =============================================================================

TREE_FILE="${TREE_PREFIX}.treefile"

[[ -s "$TREE_FILE" ]] || \
    log_error "L'arbre IQ-TREE n'a pas été généré."

log_success "Arbre IQ-TREE généré"

echo
echo "Arbre :"
echo "  $TREE_FILE"

# =============================================================================
# 20. EXTRACTION DU MODÈLE
# =============================================================================

BEST_MODEL=$(
    grep -i "Best-fit model" \
        "${TREE_PREFIX}.iqtree" |
    head -1 |
    sed 's/.*Best-fit model: *//' ||
    true
)


if [[ -z "$BEST_MODEL" ]]; then

    BEST_MODEL="Voir ${TREE_PREFIX}.iqtree"

fi

log_info "Modèle retenu : $BEST_MODEL"

# =============================================================================
# 21. RAPPORT FINAL
# =============================================================================

FINAL_REPORT="${TREE_DIR}/phylogeny_summary.txt"

{
    echo "============================================================"
    echo "PHYLOGÉNIE GÉNOMIQUE DE BORRELIA"
    echo "============================================================"
    echo
    echo "Nombre de consensus        : ${#CONSENSUS[@]}"
    echo "Nombre de références       : ${#REFERENCES[@]}"
    echo "Nombre total de séquences  : $N_SEQ"
    echo
    echo "Longueur référence         : $REFERENCE_LENGTH bp"
    echo "Longueur alignement        : $ALIGNMENT_LENGTH bp"
    echo "Gaps totaux                : $TOTAL_GAPS"
    echo "Bases ambiguës             : $TOTAL_AMBIGUOUS"
    echo
    echo "Modèle                     : $BEST_MODEL"
    echo "Bootstrap                  : $BOOTSTRAP"
    echo
    echo "Outgroup                   : $OUTGROUP_ID"
    echo
    echo "FASTA                      : $ALL_SEQUENCES"
    echo "Alignement                 : $ALIGNED"
    echo "Métadonnées                : $METADATA_FILE"
    echo "Arbre                      : $TREE_FILE"
    echo
} > "$FINAL_REPORT"


# =============================================================================
# 22. RÉSUMÉ FINAL
# =============================================================================

log_step "22 — PHYLOGÉNIE TERMINÉE"

echo
echo "============================================================"
echo "                 RÉSUMÉ FINAL"
echo "============================================================"
echo

echo "Séquences analysées :"
echo "  $N_SEQ"

echo

echo "FASTA global :"
echo "  $ALL_SEQUENCES"

echo

echo "Alignement MAFFT :"
echo "  $ALIGNED"

echo

echo "Contrôle qualité :"
echo "  $QC_DIR/"

echo

echo "Métadonnées :"
echo "  $METADATA_FILE"

echo

echo "Arbre phylogénétique :"
echo "  $TREE_FILE"

echo

echo "Rapport IQ-TREE :"
echo "  ${TREE_PREFIX}.iqtree"

echo

echo "Log IQ-TREE :"
echo "  ${TREE_PREFIX}.log"

echo

echo "Résumé :"
echo "  $FINAL_REPORT"

echo

echo "Modèle :"
echo "  $BEST_MODEL"

echo

echo "Bootstrap :"
echo "  $BOOTSTRAP"

echo

echo "Outgroup :"
echo "  $OUTGROUP_ID"

echo

echo "============================================================"

log_success "ARBRE GÉNOMIQUE ENRACINÉ AVEC SUCCÈS"

echo
echo "Pour afficher l'arbre :"
echo
echo "  cat $TREE_FILE"
echo
echo "Vous pouvez ensuite importer ce fichier .treefile dans iTOL."
echo
echo "============================================================"
