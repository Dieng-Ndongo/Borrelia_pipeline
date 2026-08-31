#!/usr/bin/env bash

# =============================================================================
# 10_filter_16S.sh
#
# OBJECTIF
# --------
# Filtrer l'alignement 16S puis reconstruire un arbre phylogénétique IQ-TREE.
#
# INPUT
# -----
#   phylogeny/locus/16S_alignment.fasta
#
# OUTPUT
# ------
#   phylogeny/locus/16S_alignment_filtered.fasta
#
#   phylogeny/locus/tree/borrelia_16S_filtered.treefile
#   phylogeny/locus/tree/borrelia_16S_filtered.contree
#   phylogeny/locus/tree/borrelia_16S_filtered.iqtree
#   phylogeny/locus/tree/borrelia_16S_filtered.log
#   etc.
#
# FILTRAGE
# --------
# Une colonne est supprimée si >= 50 % des séquences contiennent :
#   - N
#   - n
#   - -
#
# PHYLOGÉNIE
# ----------
# IQ-TREE :
#   - ModelFinder
#   - 1000 ultrafast bootstraps
#   - 1000 SH-aLRT
#   - outgroup = Borrelia miyamotoi FR64b
#
# REPRODUCTIBILITÉ
# ----------------
# L'alignement original n'est jamais modifié.
# L'arbre original n'est jamais écrasé.
# =============================================================================


set -Eeuo pipefail


# =============================================================================
# 0 — RACINE DU PROJET
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${PROJECT_ROOT}"

echo "[INFO]    Racine du projet : ${PROJECT_ROOT}"


# =============================================================================
# 1 — PARAMÈTRES
# =============================================================================

ENV_NAME="borrelia_pipeline"

ALIGNMENT="phylogeny/locus/16S_alignment.fasta"

OUTPUT_DIR="phylogeny/locus"

TREE_DIR="${OUTPUT_DIR}/tree"

LOG_DIR="results/locus_phylogeny"

FILTERED_ALIGNMENT="${OUTPUT_DIR}/16S_alignment_filtered.fasta"

TREE_PREFIX="${TREE_DIR}/borrelia_16S_filtered"

THREADS="${THREADS:-8}"

BOOTSTRAP="${BOOTSTRAP:-1000}"

SH_ALRT="${SH_ALRT:-1000}"

#
# Une colonne est supprimée si >= 50 % des séquences
# sont N ou gap.
#

MAX_MISSING="${MAX_MISSING:-0.50}"

#
# Outgroup utilisé dans l'arbre original
#

OUTGROUP_ID="Borrelia_miyamotoi_FR64b_Japon_Asia|16S"


mkdir -p \
    "${OUTPUT_DIR}" \
    "${TREE_DIR}" \
    "${LOG_DIR}"


# =============================================================================
# 2 — FONCTIONS
# =============================================================================

die() {
    echo "[ERREUR]  $*" >&2
    exit 1
}

info() {
    echo "[INFO]    $*"
}

ok() {
    echo "[OK]      $*"
}

warn() {
    echo "[WARNING] $*" >&2
}

section() {
    echo
    echo "============================================================"
    echo "  $*"
    echo "============================================================"
}


# =============================================================================
# 3 — ENVIRONNEMENT
# =============================================================================

section "1 — VÉRIFICATION DE L'ENVIRONNEMENT"

if [[ "${CONDA_DEFAULT_ENV:-}" == "${ENV_NAME}" ]]; then
    ok "Environnement actif : ${ENV_NAME}"
else
    die "L'environnement ${ENV_NAME} doit être actif."
fi


# =============================================================================
# 4 — OUTILS
# =============================================================================

section "2 — VÉRIFICATION DES OUTILS"

REQUIRED_TOOLS=(
    awk
    grep
    sort
    iqtree2
)

for tool in "${REQUIRED_TOOLS[@]}"; do

    if command -v "${tool}" >/dev/null 2>&1; then
        ok "${tool} : $(command -v "${tool}")"
    else
        die "Outil absent : ${tool}"
    fi

done


# =============================================================================
# 5 — ALIGNEMENT ORIGINAL
# =============================================================================

section "3 — VÉRIFICATION DE L'ALIGNEMENT"

[[ -r "${ALIGNMENT}" ]] \
    || die "Alignement absent : ${ALIGNMENT}"

ok "Alignement trouvé : ${ALIGNMENT}"

SEQUENCE_COUNT="$(
    grep -c '^>' "${ALIGNMENT}"
)"

[[ "${SEQUENCE_COUNT}" -ge 3 ]] \
    || die "Moins de 3 séquences dans l'alignement."

ORIGINAL_LENGTH="$(
    awk '
    /^>/ {
        if (seq != "")
            print length(seq)
        seq=""
        next
    }
    {
        seq=seq $0
    }
    END {
        if (seq != "")
            print length(seq)
    }' "${ALIGNMENT}" \
    | sort -nu
)"

[[ "$(echo "${ORIGINAL_LENGTH}" | wc -l)" -eq 1 ]] \
    || die "Les séquences n'ont pas toutes la même longueur."

ok "${SEQUENCE_COUNT} séquences."
ok "Longueur de l'alignement : ${ORIGINAL_LENGTH} bp"


# =============================================================================
# 6 — VÉRIFICATION DE L'OUTGROUP
# =============================================================================

section "4 — VÉRIFICATION DE L'OUTGROUP"

OUTGROUP_COUNT="$(
    grep -Fx ">${OUTGROUP_ID}" "${ALIGNMENT}" | wc -l
)"

if [[ "${OUTGROUP_COUNT}" -eq 0 ]]; then
    die "Outgroup introuvable :

${OUTGROUP_ID}"
fi

if [[ "${OUTGROUP_COUNT}" -gt 1 ]]; then
    die "L'outgroup apparaît plusieurs fois :
${OUTGROUP_ID}"
fi

ok "Outgroup trouvé exactement une fois."
echo "    ${OUTGROUP_ID}"


# =============================================================================
# 7 — FILTRAGE
# =============================================================================

section "5 — FILTRAGE DES COLONNES"

echo
echo "Seuil de données manquantes : ${MAX_MISSING}"
echo

awk -v max_missing="${MAX_MISSING}" '

/^>/ {

    if (name != "") {
        names[++n] = name
        seqs[n] = seq
    }

    name=$0
    sub(/^>/, "", name)

    seq=""

    next
}

{
    seq=seq $0
}

END {

    if (name != "") {
        names[++n] = name
        seqs[n] = seq
    }

    L=length(seqs[1])

    kept=0
    removed=0

    for (i=1; i<=L; i++) {

        missing=0

        for (j=1; j<=n; j++) {

            base=toupper(substr(seqs[j], i, 1))

            if (base == "N" || base == "-")
                missing++
        }

        fraction=missing/n

        if (fraction < max_missing) {

            keep[i]=1
            kept++

        } else {

            keep[i]=0
            removed++
        }
    }

    #
    # Écriture FASTA
    #

    for (j=1; j<=n; j++) {

        printf ">%s\n", names[j]

        line=""

        for (i=1; i<=L; i++) {

            if (keep[i] == 1)
                line=line substr(seqs[j], i, 1)

            if (length(line) >= 80) {
                printf "%s\n", line
                line=""
            }
        }

        if (length(line) > 0)
            printf "%s\n", line
    }

    #
    # Résumé envoyé sur stderr
    #

    printf "\nOriginal      : %d positions\n", L > "/dev/stderr"
    printf "Conservées    : %d positions\n", kept > "/dev/stderr"
    printf "Supprimées    : %d positions\n", removed > "/dev/stderr"
    printf "Séquences     : %d\n", n > "/dev/stderr"
    printf "Missing max   : %.2f%%\n", max_missing*100 > "/dev/stderr"
}
' "${ALIGNMENT}" > "${FILTERED_ALIGNMENT}"


[[ -s "${FILTERED_ALIGNMENT}" ]] \
    || die "Échec de création de l'alignement filtré."

ok "Alignement filtré créé : ${FILTERED_ALIGNMENT}"


# =============================================================================
# 8 — CONTRÔLE DU FASTA FILTRÉ
# =============================================================================

section "6 — CONTRÔLE DE L'ALIGNEMENT FILTRÉ"

FILTERED_COUNT="$(
    grep -c '^>' "${FILTERED_ALIGNMENT}"
)"

[[ "${FILTERED_COUNT}" -eq "${SEQUENCE_COUNT}" ]] \
    || die "Le nombre de séquences a changé après filtrage."

FILTERED_LENGTH="$(
    awk '
    /^>/ {
        if (seq != "")
            print length(seq)
        seq=""
        next
    }
    {
        seq=seq $0
    }
    END {
        if (seq != "")
            print length(seq)
    }' "${FILTERED_ALIGNMENT}" \
    | sort -nu
)"

[[ "$(echo "${FILTERED_LENGTH}" | wc -l)" -eq 1 ]] \
    || die "Les séquences filtrées n'ont pas toutes la même longueur."

ok "${FILTERED_COUNT} séquences conservées."
ok "Longueur filtrée : ${FILTERED_LENGTH} bp"


# =============================================================================
# 9 — SITES VARIABLES ET INFORMATIFS
# =============================================================================

section "7 — CALCUL DES SITES PHYLOGÉNÉTIQUES"

awk '

/^>/ {

    if (name != "")
        seqs[++n] = seq

    name=$0
    seq=""

    next
}

{
    seq=seq $0
}

END {

    if (name != "")
        seqs[++n] = seq

    L=length(seqs[1])

    variable=0
    informative=0

    for (i=1; i<=L; i++) {

        delete count

        for (j=1; j<=n; j++) {

            base=toupper(substr(seqs[j], i, 1))

            if (base ~ /^[ACGT]$/)
                count[base]++
        }

        distinct=0
        states2=0

        for (b in count) {

            if (count[b] > 0)
                distinct++

            if (count[b] >= 2)
                states2++
        }

        if (distinct >= 2)
            variable++

        if (states2 >= 2)
            informative++
    }

    printf "Sites variables     : %d\n", variable
    printf "Sites informatifs   : %d\n", informative
    printf "%% variables         : %.2f%%\n", 100*variable/L
    printf "%% informatifs       : %.2f%%\n", 100*informative/L
}' "${FILTERED_ALIGNMENT}"


# =============================================================================
# 10 — CONSTRUCTION DE L'ARBRE
# =============================================================================

section "8 — CONSTRUCTION DE L'ARBRE IQ-TREE"

info "Alignement : ${FILTERED_ALIGNMENT}"
info "Modèle     : MFP"
info "Bootstrap  : ${BOOTSTRAP}"
info "SH-aLRT    : ${SH_ALRT}"
info "Threads    : ${THREADS}"
info "Outgroup   : ${OUTGROUP_ID}"

iqtree2 \
    -s "${FILTERED_ALIGNMENT}" \
    -m MFP \
    -B "${BOOTSTRAP}" \
    --alrt "${SH_ALRT}" \
    -T "${THREADS}" \
    -o "${OUTGROUP_ID}" \
    --prefix "${TREE_PREFIX}" \
    --redo \
    > "${LOG_DIR}/iqtree_16S_filtered.log" \
    2>&1


# =============================================================================
# 11 — VÉRIFICATION DE L'ARBRE
# =============================================================================

section "9 — VÉRIFICATION DE L'ARBRE"

TREEFILE="${TREE_PREFIX}.treefile"

[[ -f "${TREEFILE}" ]] \
    || die "IQ-TREE n'a pas produit le treefile."

[[ -s "${TREEFILE}" ]] \
    || die "Le treefile est vide."

ok "Arbre construit : ${TREEFILE}"


# =============================================================================
# 12 — RÉSUMÉ
# =============================================================================

section "10 — RÉSUMÉ FINAL"

echo
echo "Alignement original :"
echo "  ${ALIGNMENT}"
echo "  ${SEQUENCE_COUNT} séquences"
echo "  ${ORIGINAL_LENGTH} bp"

echo
echo "Alignement filtré :"
echo "  ${FILTERED_ALIGNMENT}"
echo "  ${FILTERED_COUNT} séquences"
echo "  ${FILTERED_LENGTH} bp"

echo
echo "Arbre :"
echo "  ${TREEFILE}"

echo
echo "Outgroup :"
echo "  ${OUTGROUP_ID}"

echo
echo "Autres fichiers IQ-TREE :"
ls -1 "${TREE_PREFIX}".* 2>/dev/null || true

echo
echo "Logs :"
echo "  ${LOG_DIR}/iqtree_16S_filtered.log"

echo
echo "============================================================"
echo "  10_filter_16S.sh terminé"
echo "============================================================"
