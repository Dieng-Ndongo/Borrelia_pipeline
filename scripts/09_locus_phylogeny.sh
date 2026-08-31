#!/usr/bin/env bash

# =============================================================================
# 09_locus_phylogeny.sh
#
# OBJECTIF
# --------
# Construire une phylogénie ciblée de Borrelia crocidurae à partir :
#
#   1. des échantillons positifs B. crocidurae selon Kraken2
#   2. des génomes complets de référence B. crocidurae
#
# LOCI ÉVALUÉS
# ------------
#   - 16S
#   - 23S
#   - 5S
#   - flaB
#
# PIPELINE
# --------
# Kraken
#   ↓
# sélection B. crocidurae
#   ↓
# mapping BWA sur Achema
#   ↓
# évaluation des loci
#   ↓
# sélection du meilleur locus
#   ↓
# extraction références par BLAST
#   ↓
# extraction consensus échantillons
#   ↓
# MAFFT
#   ↓
# IQ-TREE
#
# IMPORTANT
# ---------
# Kraken est utilisé comme filtre taxonomique.
#
# Un échantillon identifié uniquement comme :
#   B. miyamotoi
#   B. hermsii
#   B. duttonii
#   B. parkeri
#   B. recurrentis
#   B. turicatae
#
# n'est PAS inclus dans la phylogénie B. crocidurae.
#
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

REFERENCE="data/reference/Borrelia_crocidurae_Achema_CP003426.fasta"

LOCUS_FILE="results/locus_detection/loci_detected.tsv"

RAW_DIR="data/raw"

KRAKEN_DIR="kraken"

SOUCHES_DIR="data/reference/souches_reference"

# Nouveau dossier pour ne pas écraser l'ancienne analyse
OUTPUT_DIR="phylogeny/locus"

WORK_DIR="${OUTPUT_DIR}/work"

TREE_DIR="${OUTPUT_DIR}/tree"

LOG_DIR="results/locus_phylogeny"

THREADS="${THREADS:-8}"

# TaxID B. crocidurae
BORRELIA_CROCIDURAE_TAXID="29520"

# Seuil de couverture du locus
MIN_COVERAGE="${MIN_COVERAGE:-40}"

# Fraction minimale d'échantillons exploitables
MIN_SAMPLE_FRACTION="${MIN_SAMPLE_FRACTION:-0.30}"

# Profondeur minimale
MIN_DEPTH="${MIN_DEPTH:-3}"

# Longueur minimale d'un locus
MIN_LENGTH="${MIN_LENGTH:-100}"


mkdir -p \
    "${OUTPUT_DIR}" \
    "${WORK_DIR}" \
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
    echo
}

progress() {
    echo "[PROGRESS] $*"
}

cleanup_tmp() {
    rm -f \
        "${WORK_DIR}"/*.tmp \
        "${WORK_DIR}"/*.sam \
        2>/dev/null || true
}

trap cleanup_tmp EXIT


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
    bwa
    samtools
    bcftools
    mafft
    iqtree2
    blastn
    makeblastdb
    seqkit
    awk
    grep
    sed
    sort
)

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "${tool}" >/dev/null 2>&1; then
        ok "${tool} : $(command -v "${tool}")"
    else
        die "Outil absent : ${tool}"
    fi
done


# =============================================================================
# 5 — VÉRIFICATION RÉFÉRENCE
# =============================================================================

section "3 — VÉRIFICATION DE LA RÉFÉRENCE"

[[ -r "${REFERENCE}" ]] \
    || die "Référence absente : ${REFERENCE}"

ok "Référence : ${REFERENCE}"

echo "[INFO]    $(head -n 1 "${REFERENCE}")"


# =============================================================================
# 6 — LOCI
# =============================================================================

section "4 — LECTURE DES LOCI DÉTECTÉS"

[[ -r "${LOCUS_FILE}" ]] \
    || die "Fichier absent : ${LOCUS_FILE}"

cat "${LOCUS_FILE}"


# =============================================================================
# 7 — EXTRACTION COORDONNÉES
# =============================================================================

section "5 — EXTRACTION DES COORDONNÉES DES LOCI"

declare -A LOCUS_START
declare -A LOCUS_END
declare -A LOCUS_STRAND

for locus in 16S 23S 5S flaB; do

    LINE="$(
        awk -F '\t' -v L="${locus}" '
            NR > 1 && $2 == L {
                print $4 "\t" $5 "\t" $6
                exit
            }
        ' "${LOCUS_FILE}"
    )"

    if [[ -n "${LINE}" ]]; then

        IFS=$'\t' read -r start end strand <<< "${LINE}"

        LOCUS_START["${locus}"]="${start}"
        LOCUS_END["${locus}"]="${end}"
        LOCUS_STRAND["${locus}"]="${strand}"

        ok "${locus} : ${start}-${end} (${strand})"

    else

        warn "${locus} non détecté."

    fi

done


# =============================================================================
# 8 — INDEXATION RÉFÉRENCE
# =============================================================================

section "6 — INDEXATION DE LA RÉFÉRENCE"

if [[ ! -f "${REFERENCE}.bwt" ]]; then

    info "Construction de l'index BWA..."

    bwa index "${REFERENCE}" \
        > "${LOG_DIR}/bwa_index.log" \
        2>&1

fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    samtools faidx "${REFERENCE}"
fi

ok "Référence indexée."


# =============================================================================
# 9 — DÉTECTION FASTQ
# =============================================================================

section "7 — DÉTECTION DES ÉCHANTILLONS"

mapfile -t R1_FILES < <(
    find "${RAW_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*_merged_R1.fastq.gz" \
        | sort
)

(( ${#R1_FILES[@]} > 0 )) \
    || die "Aucun FASTQ merged trouvé dans ${RAW_DIR}."

ok "${#R1_FILES[@]} échantillon(s) détecté(s)."


# =============================================================================
# 10 — FILTRE KRAKEN : B. CROCIDURAE
# =============================================================================

section "8 — FILTRE TAXONOMIQUE KRAKEN : B. CROCIDURAE"

KRAKEN_POSITIVE_LIST="${OUTPUT_DIR}/kraken_b_crocidurae_positive.txt"
KRAKEN_NEGATIVE_LIST="${OUTPUT_DIR}/kraken_b_crocidurae_negative.txt"
KRAKEN_SUMMARY="${OUTPUT_DIR}/kraken_taxonomic_filter.tsv"

# TaxIDs compatibles avec B. crocidurae
#
# 29520   = Borrelia crocidurae
# 1155096 = Borrelia crocidurae str. Achema
# 1293575 = Borrelia crocidurae DOU

B_CROCIDURAE_TAXIDS=(
    "29520"
    "1155096"
    "1293575"
)

: > "${KRAKEN_POSITIVE_LIST}"
: > "${KRAKEN_NEGATIVE_LIST}"

printf "sample\tB_crocidurae_reads\tstatus\n" \
    > "${KRAKEN_SUMMARY}"


for R1 in "${R1_FILES[@]}"; do

    SAMPLE="$(basename "${R1}" "_merged_R1.fastq.gz")"

    KRAKEN_REPORT="${KRAKEN_DIR}/${SAMPLE}/${SAMPLE}_report.txt"

    if [[ ! -r "${KRAKEN_REPORT}" ]]; then

        warn "${SAMPLE} : rapport Kraken absent."

        printf "%s\t0\tKRAKEN_REPORT_ABSENT\n" \
            "${SAMPLE}" \
            >> "${KRAKEN_SUMMARY}"

        continue

    fi


    # -------------------------------------------------------------------------
    # Recherche de B. crocidurae et de ses souches
    #
    # IMPORTANT :
    #
    # La colonne 2 du rapport Kraken correspond au nombre de reads
    # dans le clade.
    #
    # Pour le taxon parent 29520, ce nombre peut déjà inclure les
    # descendants Achema et DOU.
    #
    # On évite donc de sommer parent + descendants afin de ne pas
    # compter deux fois les mêmes reads.
    # -------------------------------------------------------------------------

    B_CROCIDURAE_READS="$(
        awk '
            $4 == "S" && ($5 == 29520 || $5 == 1155096 || $5 == 1293575) {

                # Priorité au taxon parent B. crocidurae
                if ($5 == 29520) {
                    parent = $2
                }

                # Conserver le maximum parmi les descendants
                if ($5 == 1155096 || $5 == 1293575) {
                    if ($2 > descendant_max)
                        descendant_max = $2
                }
            }

            END {

                # Si le parent existe, utiliser son nombre de reads.
                if (parent > 0) {
                    print parent
                }

                # Sinon utiliser le nombre du descendant détecté.
                else if (descendant_max > 0) {
                    print descendant_max
                }

                else {
                    print 0
                }
            }
        ' "${KRAKEN_REPORT}"
    )"

    B_CROCIDURAE_READS="${B_CROCIDURAE_READS:-0}"


    # -------------------------------------------------------------------------
    # Classification
    # -------------------------------------------------------------------------

    if [[ "${B_CROCIDURAE_READS}" =~ ^[0-9]+$ ]] &&
       (( B_CROCIDURAE_READS > 0 )); then

        echo "${SAMPLE}" >> "${KRAKEN_POSITIVE_LIST}"

        printf "%s\t%s\tPOSITIVE\n" \
            "${SAMPLE}" \
            "${B_CROCIDURAE_READS}" \
            >> "${KRAKEN_SUMMARY}"

        ok "${SAMPLE} : ${B_CROCIDURAE_READS} reads B. crocidurae → CANDIDAT"

    else

        echo "${SAMPLE}" >> "${KRAKEN_NEGATIVE_LIST}"

        printf "%s\t0\tNEGATIVE\n" \
            "${SAMPLE}" \
            >> "${KRAKEN_SUMMARY}"

        warn "${SAMPLE} : aucun read B. crocidurae → EXCLU"

    fi

done

# =============================================================================
# 11 — VÉRIFICATION DES CANDIDATS
# =============================================================================

section "9 — ÉCHANTILLONS RETENUS APRÈS KRAKEN"

echo
echo "Échantillons B. crocidurae candidats :"
cat "${KRAKEN_POSITIVE_LIST}" || true

echo

POSITIVE_COUNT=$(awk 'NF {n++} END {print n+0}' "${KRAKEN_POSITIVE_LIST}")
NEGATIVE_COUNT=$(awk 'NF {n++} END {print n+0}' "${KRAKEN_NEGATIVE_LIST}")

echo
echo "B. crocidurae candidats : ${POSITIVE_COUNT}"
echo "Échantillons exclus      : ${NEGATIVE_COUNT}"


(( POSITIVE_COUNT > 0 )) \
    || die "Aucun échantillon B. crocidurae détecté par Kraken."


# =============================================================================
# 12 — MAPPING UNIQUEMENT DES POSITIFS KRAKEN
# =============================================================================

section "10 — MAPPING DES CANDIDATS B. CROCIDURAE"

while read -r SAMPLE; do

    [[ -n "${SAMPLE}" ]] || continue

    R1="${RAW_DIR}/${SAMPLE}_merged_R1.fastq.gz"
    R2="${RAW_DIR}/${SAMPLE}_merged_R2.fastq.gz"

    [[ -r "${R1}" ]] || die "R1 absent : ${R1}"
    [[ -r "${R2}" ]] || die "R2 absent : ${R2}"

    BAM="${WORK_DIR}/${SAMPLE}.bam"

    info "Mapping : ${SAMPLE}"

    if [[ ! -f "${BAM}" ]]; then

        progress "${SAMPLE} → BWA MEM en cours..."

        bwa mem \
            -t "${THREADS}" \
            "${REFERENCE}" \
            "${R1}" \
            "${R2}" \
            2> "${LOG_DIR}/${SAMPLE}_bwa.log" \
        | samtools sort \
            -@ "${THREADS}" \
            -o "${BAM}" -

        ok "${SAMPLE} → mapping BWA terminé."

    else

        ok "${SAMPLE} → BAM déjà présent, mapping ignoré."

    fi

    samtools index "${BAM}"

    TOTAL_READS="$(
        samtools view -c "${BAM}"
    )"

    MAPPED_READS="$(
        samtools view -c -F 4 "${BAM}"
    )"

    if (( TOTAL_READS > 0 )); then

        MAPPED_PCT="$(
            awk -v m="${MAPPED_READS}" -v t="${TOTAL_READS}" \
                'BEGIN {printf "%.2f", (m/t)*100}'
        )"

    else

        MAPPED_PCT="0"

    fi

    ok "${SAMPLE} : ${MAPPED_READS}/${TOTAL_READS} reads mappés (${MAPPED_PCT}%)."

done < "${KRAKEN_POSITIVE_LIST}"


# =============================================================================
# 13 — ÉVALUATION DES LOCI
# =============================================================================

section "11 — ÉVALUATION AUTOMATIQUE DES LOCI"

SUMMARY="${OUTPUT_DIR}/locus_evaluation.tsv"

printf "locus\tlength\tsamples_total\tsamples_usable\tmean_coverage\tmin_depth\tmean_depth\tvariable_sites\tscore\n" \
    > "${SUMMARY}"

TOTAL_SAMPLES="${POSITIVE_COUNT}"


for locus in 16S 23S 5S flaB; do

    if [[ -z "${LOCUS_START[${locus}]:-}" ]]; then
        continue
    fi


    START="${LOCUS_START[${locus}]}"
    END="${LOCUS_END[${locus}]}"

    LENGTH=$((END - START + 1))

    info "Évaluation : ${locus} (${START}-${END})"
    progress "Locus ${locus} → évaluation de ${TOTAL_SAMPLES} échantillon(s)..."


    USABLE=0
    TOTAL_COVERAGE=0
    TOTAL_DEPTH=0


    while read -r SAMPLE; do

        [[ -n "${SAMPLE}" ]] || continue

        BAM="${WORK_DIR}/${SAMPLE}.bam"

        DEPTH_FILE="${WORK_DIR}/${SAMPLE}_${locus}_depth.tmp"

        samtools depth \
            -a \
            -r "CP003426.1:${START}-${END}" \
            "${BAM}" \
            > "${DEPTH_FILE}"


        COVERED="$(
            awk -v min="${MIN_DEPTH}" '
                $3 >= min {n++}
                END {print n+0}
            ' "${DEPTH_FILE}"
        )"


        MEAN_DEPTH="$(
            awk '
                {
                    sum += $3
                    n++
                }
                END {
                    if(n>0)
                        printf "%.2f", sum/n
                    else
                        print "0"
                }
            ' "${DEPTH_FILE}"
        )"


        COVERAGE="$(
            awk -v c="${COVERED}" -v l="${LENGTH}" \
                'BEGIN {
                    if(l>0)
                        printf "%.2f", (c/l)*100
                    else
                        print "0"
                }'
        )"


        if awk -v c="${COVERAGE}" -v min="${MIN_COVERAGE}" \
            'BEGIN {exit !(c >= min)}'; then

            USABLE=$((USABLE + 1))

        fi


        TOTAL_COVERAGE=$(
            awk -v a="${TOTAL_COVERAGE}" -v b="${COVERAGE}" \
                'BEGIN {printf "%.4f", a+b}'
        )


        TOTAL_DEPTH=$(
            awk -v a="${TOTAL_DEPTH}" -v b="${MEAN_DEPTH}" \
                'BEGIN {printf "%.4f", a+b}'
        )


        echo
        echo "    ${SAMPLE}"
        echo "       ${locus} couverture : ${COVERAGE}%"
        echo "       profondeur        : ${MEAN_DEPTH}x"

    done < "${KRAKEN_POSITIVE_LIST}"
    ok "Évaluation du locus ${locus} terminée."


    MEAN_COVERAGE=$(
        awk -v x="${TOTAL_COVERAGE}" -v n="${TOTAL_SAMPLES}" \
            'BEGIN {
                if(n>0)
                    printf "%.2f", x/n
                else
                    print "0"
            }'
    )


    MEAN_DEPTH_ALL=$(
        awk -v x="${TOTAL_DEPTH}" -v n="${TOTAL_SAMPLES}" \
            'BEGIN {
                if(n>0)
                    printf "%.2f", x/n
                else
                    print "0"
            }'
    )


    # -------------------------------------------------------------------------
    # SCORE
    #
    # priorité :
    #   - fraction d'échantillons exploitables
    #   - couverture
    #   - profondeur
    # -------------------------------------------------------------------------

    SCORE=$(
        awk \
            -v usable="${USABLE}" \
            -v total="${TOTAL_SAMPLES}" \
            -v cov="${MEAN_COVERAGE}" \
            -v depth="${MEAN_DEPTH_ALL}" \
            'BEGIN {

                if(total > 0)
                    fraction = usable / total
                else
                    fraction = 0

                score = (fraction * 50) + \
                        (cov * 0.40) + \
                        (depth * 1.0)

                printf "%.4f", score
            }'
    )


    printf "%s\t%d\t%d\t%d\t%s\t%d\t%s\t%d\t%s\n" \
        "${locus}" \
        "${LENGTH}" \
        "${TOTAL_SAMPLES}" \
        "${USABLE}" \
        "${MEAN_COVERAGE}" \
        "${MIN_DEPTH}" \
        "${MEAN_DEPTH_ALL}" \
        "0" \
        "${SCORE}" \
        >> "${SUMMARY}"


done


# =============================================================================
# 14 — AFFICHAGE CLASSEMENT
# =============================================================================

section "12 — CLASSEMENT DES LOCI"

sort \
    -t $'\t' \
    -k9,9nr \
    "${SUMMARY}" \
    | column -t -s $'\t'


# =============================================================================
# 15 — SÉLECTION DU MEILLEUR LOCUS
# =============================================================================

section "13 — SÉLECTION DU LOCUS"

BEST_LINE="$(
    awk -F '\t' '
        NR > 1 && $4 > 0 {
            print
        }
    ' "${SUMMARY}" |
    sort -t $'\t' -k9,9nr |
    head -n 1
)"

[[ -n "${BEST_LINE}" ]] \
    || die "Aucun locus ne possède une couverture suffisante."


IFS=$'\t' read -r \
    BEST_LOCUS \
    BEST_LENGTH \
    BEST_TOTAL \
    BEST_USABLE \
    BEST_COVERAGE \
    BEST_MIN_DEPTH \
    BEST_MEAN_DEPTH \
    BEST_VARIABLE \
    BEST_SCORE \
    <<< "${BEST_LINE}"


BEST_FRACTION=$(
    awk \
        -v u="${BEST_USABLE}" \
        -v t="${BEST_TOTAL}" \
        'BEGIN {
            if(t>0)
                printf "%.3f", u/t
            else
                print "0"
        }'
)


if ! awk \
    -v f="${BEST_FRACTION}" \
    -v min="${MIN_SAMPLE_FRACTION}" \
    'BEGIN {exit !(f >= min)}'; then

    die "Fraction insuffisante : ${BEST_USABLE}/${BEST_TOTAL}"

fi


START="${LOCUS_START[${BEST_LOCUS}]}"
END="${LOCUS_END[${BEST_LOCUS}]}"
STRAND="${LOCUS_STRAND[${BEST_LOCUS}]}"

progress "Comparaison des loci 16S / 23S / 5S / flaB..."
echo
echo "============================================================"
echo "  LOCUS RETENU : ${BEST_LOCUS}"
echo "============================================================"
echo
echo "  Longueur             : ${BEST_LENGTH} bp"
echo "  Échantillons         : ${BEST_USABLE}/${BEST_TOTAL}"
echo "  Couverture moyenne   : ${BEST_COVERAGE}%"
echo "  Profondeur moyenne   : ${BEST_MEAN_DEPTH}x"
echo "  Score                : ${BEST_SCORE}"
echo "  Coordonnées Achema   : CP003426.1:${START}-${END}"
echo "  Orientation          : ${STRAND}"
echo


# =============================================================================
# 16 — EXTRACTION LOCUS DES RÉFÉRENCES
# =============================================================================

section "14 — EXTRACTION DES LOCI DES RÉFÉRENCES"

REFERENCE_FASTA="${OUTPUT_DIR}/${BEST_LOCUS}_references.fasta"

: > "${REFERENCE_FASTA}"


QUERY="${WORK_DIR}/${BEST_LOCUS}_query.fasta"

samtools faidx \
    "${REFERENCE}" \
    "CP003426.1:${START}-${END}" \
    > "${QUERY}"


QUERY_LEN="$(
    awk '
        !/^>/ {
            sum += length($0)
        }
        END {
            print sum
        }
    ' "${QUERY}"
)"


ok "Séquence requête : ${QUERY_LEN} bp"


# =============================================================================
# 17 — NORMALISATION DES NOMS
# =============================================================================

SOUCHES_NORMALIZED="${WORK_DIR}/souches_normalized"

mkdir -p "${SOUCHES_NORMALIZED}"

mapfile -t SOUCHE_FILES < <(
    find "${SOUCHES_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.fasta" \
        | sort
)


(( ${#SOUCHE_FILES[@]} > 0 )) \
    || die "Aucun génome trouvé dans ${SOUCHES_DIR}."


NORMALIZED_FILES=()


for SOUCHE in "${SOUCHE_FILES[@]}"; do

    ORIGINAL_NAME="$(basename "${SOUCHE}")"

    NORMALIZED_NAME="${ORIGINAL_NAME// /_}"

    SYMLINK="${SOUCHES_NORMALIZED}/${NORMALIZED_NAME}"


    if [[ ! -e "${SYMLINK}" ]]; then

        ln -sf \
            "$(realpath "${SOUCHE}")" \
            "${SYMLINK}"

    fi


    NORMALIZED_FILES+=("${SYMLINK}")

done


ok "${#NORMALIZED_FILES[@]} génome(s) de référence détecté(s)."


# =============================================================================
# 18 — BLAST
# =============================================================================

for SOUCHE in "${NORMALIZED_FILES[@]}"; do

    SNAME="$(basename "${SOUCHE}" .fasta)"

    progress "BLAST du locus ${BEST_LOCUS} contre ${SNAME}..."

    BLAST_DB="${WORK_DIR}/${SNAME}_blastdb"

    if [[ ! -f "${BLAST_DB}.nhr" ]]; then

        makeblastdb \
            -in "${SOUCHE}" \
            -dbtype nucl \
            -out "${BLAST_DB}" \
            > "${LOG_DIR}/${SNAME}_makeblastdb.log" \
            2>&1

    fi

    BLAST_OUT="${WORK_DIR}/${SNAME}_${BEST_LOCUS}_blast.tsv"

    blastn \
        -query "${QUERY}" \
        -db "${BLAST_DB}" \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
        -max_target_seqs 1 \
        -evalue 1e-10 \
        > "${BLAST_OUT}" \
        2> "${LOG_DIR}/${SNAME}_${BEST_LOCUS}_blast.log"

    if [[ ! -s "${BLAST_OUT}" ]]; then

        warn "${SNAME} : aucun hit ${BEST_LOCUS}."

        continue

    fi

    BEST_HIT="$(
        sort -k10,10nr "${BLAST_OUT}" |
        head -n 1
    )"

    IFS=$'\t' read -r \
        QSEQ \
        SSEQ \
        PID \
        ALIGNLEN \
        QSTART \
        QEND \
        SSTART \
        SEND \
        EVALUE \
        BITSCORE \
        <<< "${BEST_HIT}"

    QUERY_COV=$(
        awk \
            -v a="${ALIGNLEN}" \
            -v l="${BEST_LENGTH}" \
            'BEGIN {
                if (l > 0)
                    printf "%.1f", (a/l)*100
                else
                    print "0"
            }'
    )

    if ! awk \
        -v c="${QUERY_COV}" \
        'BEGIN {exit !(c >= 70)}'; then

        warn "${SNAME} : couverture BLAST insuffisante (${QUERY_COV}%)."

        continue

    fi

    SOUCHE_EXTRACT="${WORK_DIR}/${SNAME}_${BEST_LOCUS}_raw.fasta"

    [[ -f "${SOUCHE}.fai" ]] ||
        samtools faidx "${SOUCHE}"

    if (( SSTART <= SEND )); then

        samtools faidx \
            "${SOUCHE}" \
            "${SSEQ}:${SSTART}-${SEND}" \
            > "${SOUCHE_EXTRACT}"

    else

        samtools faidx \
            "${SOUCHE}" \
            "${SSEQ}:${SEND}-${SSTART}" |
        seqkit seq -r -p \
            > "${SOUCHE_EXTRACT}"

    fi

    [[ -s "${SOUCHE_EXTRACT}" ]] ||
        die "${SNAME} : extraction de séquence échouée."

    awk \
        -v name="${SNAME}" \
        -v locus="${BEST_LOCUS}" '
        /^>/ {
            print ">" name "|" locus
            next
        }
        {
            print
        }
    ' "${SOUCHE_EXTRACT}" \
    >> "${REFERENCE_FASTA}"

    ok "${SNAME} : ${BEST_LOCUS} extrait (${PID}% identité, ${QUERY_COV}% couverture)."
    progress "BLAST du locus ${BEST_LOCUS} contre ${SNAME}..."

done


REF_COUNT="$(
    grep -c '^>' "${REFERENCE_FASTA}" || true
)"

(( REF_COUNT > 0 )) \
    || die "Aucune séquence de référence extraite."

ok "${REF_COUNT} séquence(s) de référence."


# =============================================================================
# 19 — EXTRACTION DES CONSENSUS
# =============================================================================

section "15 — EXTRACTION DES CONSENSUS DES ÉCHANTILLONS"

SAMPLE_FASTA="${OUTPUT_DIR}/${BEST_LOCUS}_samples.fasta"

: > "${SAMPLE_FASTA}"


while read -r SAMPLE; do

    [[ -n "${SAMPLE}" ]] || continue
    progress "${SAMPLE} → extraction du consensus ${BEST_LOCUS}..."

    BAM="${WORK_DIR}/${SAMPLE}.bam"

    DEPTH_FILE="${WORK_DIR}/${SAMPLE}_${BEST_LOCUS}_depth.txt"

    CONSENSUS="${WORK_DIR}/${SAMPLE}_${BEST_LOCUS}.fasta"


    samtools depth \
        -a \
        -r "CP003426.1:${START}-${END}" \
        "${BAM}" \
        > "${DEPTH_FILE}"


    COVERED="$(
        awk \
            -v min="${MIN_DEPTH}" '
            $3 >= min {n++}
            END {print n+0}
        ' "${DEPTH_FILE}"
    )"


    COVERAGE=$(
        awk \
            -v c="${COVERED}" \
            -v l="${BEST_LENGTH}" \
            'BEGIN {
                if(l>0)
                    printf "%.2f", (c/l)*100
                else
                    print "0"
            }'
    )


    if ! awk \
        -v c="${COVERAGE}" \
        -v min="${MIN_COVERAGE}" \
        'BEGIN {exit !(c >= min)}'; then

        warn "${SAMPLE} : couverture insuffisante sur ${BEST_LOCUS} (${COVERAGE}%) — EXCLU."

        continue

    fi


    samtools consensus \
        -r "CP003426.1:${START}-${END}" \
        --min-MQ 20 \
        --min-BQ 20 \
        "${BAM}" \
        > "${CONSENSUS}" \
        2> "${LOG_DIR}/${SAMPLE}_${BEST_LOCUS}_consensus.log"


    SEQ="$(
        awk '
            !/^>/ {
                printf "%s", $0
            }
        ' "${CONSENSUS}"
    )"


    [[ -n "${SEQ}" ]] ||
        die "${SAMPLE} : consensus vide."


    echo ">${SAMPLE}|${BEST_LOCUS}" \
        >> "${SAMPLE_FASTA}"

    echo "${SEQ}" \
        >> "${SAMPLE_FASTA}"


    ok "${SAMPLE} : ${BEST_LOCUS} extrait (${COVERAGE}%)."


done < "${KRAKEN_POSITIVE_LIST}"


SAMPLE_COUNT="$(
    grep -c '^>' "${SAMPLE_FASTA}" || true
)"


(( SAMPLE_COUNT >= 1 )) \
    || die "Aucun candidat B. crocidurae n'a une couverture suffisante."


ok "${SAMPLE_COUNT} séquence(s) échantillon."


# =============================================================================
# 20 — FASTA FINAL
# =============================================================================

section "16 — CONSTRUCTION DU FASTA FINAL"

COMBINED="${OUTPUT_DIR}/${BEST_LOCUS}_all_sequences.fasta"

cat \
    "${REFERENCE_FASTA}" \
    "${SAMPLE_FASTA}" \
    > "${COMBINED}"


FINAL_COUNT="$(
    grep -c '^>' "${COMBINED}"
)"


echo
echo "Séquences finales :"
grep '^>' "${COMBINED}"
echo


ok "${FINAL_COUNT} séquences finales."

(( FINAL_COUNT >= 3 )) \
    || die "Moins de 3 séquences."


# =============================================================================
# 21 — ALIGNEMENT MAFFT
# =============================================================================

section "17 — ALIGNEMENT MAFFT"

ALIGNMENT="${OUTPUT_DIR}/${BEST_LOCUS}_alignment.fasta"
progress "MAFFT → alignement de ${FINAL_COUNT} séquences..."

mafft \
    --auto \
    --thread "${THREADS}" \
    "${COMBINED}" \
    > "${ALIGNMENT}" \
    2> "${LOG_DIR}/mafft_${BEST_LOCUS}.log"


ALN_COUNT="$(
    grep -c '^>' "${ALIGNMENT}"
)"


ok "Alignement terminé : ${ALN_COUNT} séquences."


# =============================================================================
# 22 — IQ-TREE
# =============================================================================

section "18 — CONSTRUCTION DE L'ARBRE"

TREE_PREFIX="${TREE_DIR}/borrelia_crocidurae_${BEST_LOCUS}"


OUTGROUP_ID="Borrelia_miyamotoi_FR64b_Japon_Asia|${BEST_LOCUS}"


info "Outgroup : ${OUTGROUP_ID}"


OUTGROUP_COUNT="$(
    grep "^>" "${COMBINED}" |
    sed 's/^>//' |
    awk -v id="${OUTGROUP_ID}" '
        $0 == id {
            count++
        }
        END {
            print count+0
        }
    '
)"


if [[ "${OUTGROUP_COUNT}" -eq 0 ]]; then

    warn "Outgroup ${OUTGROUP_ID} absent du FASTA."

    warn "L'arbre sera construit sans enracinement explicite."

    progress "IQ-TREE → construction de l'arbre ${BEST_LOCUS}..."

    iqtree2 \
        -s "${ALIGNMENT}" \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T "${THREADS}" \
        --prefix "${TREE_PREFIX}" \
        --redo \
        > "${LOG_DIR}/iqtree_${BEST_LOCUS}.log" \
        2>&1

    ok "IQ-TREE → arbre phylogénétique terminé."

else

    ok "Outgroup trouvé exactement une fois."

    iqtree2 \
        -s "${ALIGNMENT}" \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T "${THREADS}" \
        -o "${OUTGROUP_ID}" \
        --prefix "${TREE_PREFIX}" \
        --redo \
        > "${LOG_DIR}/iqtree_${BEST_LOCUS}.log" \
        2>&1

fi


[[ -f "${TREE_PREFIX}.treefile" ]] \
    || die "IQ-TREE n'a pas produit le treefile."


ok "Arbre phylogénétique construit."


# =============================================================================
# 23 — RÉSUMÉ FINAL
# =============================================================================

section "19 — RÉSUMÉ FINAL"

echo
echo "============================================================"
echo "  FILTRE TAXONOMIQUE"
echo "============================================================"

echo
cat "${KRAKEN_SUMMARY}"


echo
echo "============================================================"
echo "  LOCI"
echo "============================================================"

cat "${SUMMARY}"


echo
echo "============================================================"
echo "  RÉSULTAT FINAL"
echo "============================================================"

echo
echo "Locus retenu          : ${BEST_LOCUS}"
echo "Candidats Kraken      : ${POSITIVE_COUNT}"
echo "Échantillons dans arbre : ${SAMPLE_COUNT}"
echo "Références            : ${REF_COUNT}"
echo "Séquences finales     : ${FINAL_COUNT}"
echo


echo "Fichiers :"

echo "  Filtre Kraken       : ${KRAKEN_SUMMARY}"

echo "  Liste positifs      : ${KRAKEN_POSITIVE_LIST}"

echo "  Liste négatifs      : ${KRAKEN_NEGATIVE_LIST}"

echo "  Évaluation loci     : ${SUMMARY}"

echo "  FASTA final         : ${COMBINED}"

echo "  Alignement          : ${ALIGNMENT}"

echo "  Arbre               : ${TREE_PREFIX}.treefile"

echo "  Bootstrap           : ${TREE_PREFIX}.contree"

echo "  Logs                : ${LOG_DIR}/"


echo
echo "============================================================"
echo "  09_locus_phylogeny.sh terminé"
echo "============================================================"

