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
# chargement références phylogeny_locus/
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
 
PHYLOGENY_LOCUS_DIR="data/reference/phylogeny_locus"
 
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
 
 
    B_CROCIDURAE_READS="$(
        awk '
            $4 == "S" && ($5 == 29520 || $5 == 1155096 || $5 == 1293575) {
 
                if ($5 == 29520) {
                    parent = $2
                }
 
                if ($5 == 1155096 || $5 == 1293575) {
                    if ($2 > descendant_max)
                        descendant_max = $2
                }
            }
 
            END {
 
                if (parent > 0) {
                    print parent
                }
 
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
# 16 — CHARGEMENT DES RÉFÉRENCES DEPUIS phylogeny_locus
# =============================================================================
#
# phylogeny_locus/ contient déjà les séquences du locus retenu (ex. 16S),
# une par souche. Pas besoin de BLAST ni d'extraction : on reformate juste
# les en-têtes et on concatène.
#
# Format attendu du header GenBank :
#   >ACC Genus species strain NOM 16S ribosomal RNA ...
#
# Header produit :
#   >ACC_Genus_species_strain_NOM|16S
#
# Collisions inter-fichiers : suffixe _2, _3, ...
# =============================================================================
 
section "14 — CHARGEMENT DES RÉFÉRENCES DEPUIS phylogeny_locus"
 
REFERENCE_FASTA="${OUTPUT_DIR}/${BEST_LOCUS}_references.fasta"
 
: > "${REFERENCE_FASTA}"
 
 
mapfile -t LOCUS_REF_FILES < <(
    find "${PHYLOGENY_LOCUS_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.fasta" \
        | sort
)
 
 
(( ${#LOCUS_REF_FILES[@]} > 0 )) \
    || die "Aucun FASTA trouvé dans ${PHYLOGENY_LOCUS_DIR}."
 
 
ok "${#LOCUS_REF_FILES[@]} fichier(s) de référence détecté(s)."
 
 
# Filtrer les fichiers non vides
VALID_FILES=()
for REF_FILE in "${LOCUS_REF_FILES[@]}"; do
    FNAME="$(basename "${REF_FILE}" .fasta)"
    if [[ -s "${REF_FILE}" ]]; then
        VALID_FILES+=("${REF_FILE}")
    else
        warn "${FNAME} : fichier vide, ignoré."
    fi
done
 
(( ${#VALID_FILES[@]} > 0 )) \
    || die "Tous les fichiers de référence sont vides."
 
 
# Un seul appel awk sur tous les fichiers valides :
# → le tableau count[] est global → déduplication inter-fichiers
#
# Deux formats de header supportés :
#   1. Pipe-délimité : >ACC | Genus species strain NOM 16S | Pays
#   2. GenBank classique : >ACC Genus species strain NOM 16S ribosomal RNA...
#
# Résultat : >ACC_Genus_species_strain_NOM|16S
# Collision inter-fichiers : suffixe _2, _3, ...
awk \
    -v locus="${BEST_LOCUS}" '
    /^>/ {
        sub(/^>/, "")
        line = $0

        # Détecter le format : pipe-délimité ou GenBank classique
        if (line ~ / \| /) {
            # Format : "ACC | description taxonomique | geo"
            split(line, p, / \| /)
            acc  = p[1]; gsub(/^ +| +$/, "", acc)
            taxo = p[2]; gsub(/^ +| +$/, "", taxo)
            geo  = p[3]; gsub(/^ +| +$/, "", geo)
        } else {
            # Format GenBank : "ACC description..."
            split(line, p, " ")
            acc  = p[1]
            taxo = substr(line, length(p[1]) + 2)
            geo  = ""
        }

        # Tronquer la description avant le premier mot de locus
        split(taxo, words, " ")
        desc = ""
        for (i = 1; i <= length(words); i++) {
            w = words[i]
            if (w ~ /^(16S|23S|5S|flaB|ribosomal|rRNA|RNA|gene|intergenic|partial|complete|sequence|genomic)$/)
                break
            desc = (desc == "") ? w : desc "_" w
        }

        gsub(/[^A-Za-z0-9_.\-]/, "_", acc)
        gsub(/[^A-Za-z0-9_.\-]/, "_", desc)
        gsub(/_+/, "_", desc)
        gsub(/_$/, "", desc)
        gsub(/[^A-Za-z0-9_.\-]/, "_", geo)
        gsub(/_+/, "_", geo)
        gsub(/_$/, "", geo)

        base_id = (desc != "") ? acc "_" desc : acc
        if (geo != "") base_id = base_id "_" geo
        full_id  = base_id "|" locus

        count[full_id]++
        if (count[full_id] > 1) {
            full_id = base_id "_" count[full_id] "|" locus
        }

        print ">" full_id
        next
    }
    { print }
' "${VALID_FILES[@]}" \
>> "${REFERENCE_FASTA}"
 
 
REF_COUNT="$(
    grep -c '^>' "${REFERENCE_FASTA}" || true
)"
 
(( REF_COUNT > 0 )) \
    || die "Aucune séquence de référence extraite."
 
ok "${REF_COUNT} séquence(s) de référence."
 
 
# =============================================================================
# 17 — EXTRACTION DES CONSENSUS DES ÉCHANTILLONS
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
# 18 — FASTA FINAL
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
# 19 — ALIGNEMENT MAFFT
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
# 20 — IQ-TREE
# =============================================================================
 
section "18 — CONSTRUCTION DE L'ARBRE"
 
TREE_PREFIX="${TREE_DIR}/borrelia_crocidurae_${BEST_LOCUS}"
 
 
# -------------------------------------------------------------------------
# Détection dynamique de l'outgroup :
# On cherche la première séquence "miyamotoi" dans le FASTA combiné.
# Cela évite tout hardcoding du nom exact de la souche.
# -------------------------------------------------------------------------
OUTGROUP_ID="$(
    grep '^>' "${COMBINED}" \
    | sed 's/^>//' \
    | grep -i 'miyamotoi' \
    | head -n 1 \
    || true
)"
 
 
if [[ -n "${OUTGROUP_ID}" ]]; then
 
    info "Outgroup détecté : ${OUTGROUP_ID}"
 
    progress "IQ-TREE → construction de l'arbre ${BEST_LOCUS} (outgroup : ${OUTGROUP_ID})..."
 
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
 
    ok "IQ-TREE → arbre enraciné sur ${OUTGROUP_ID}."
 
else
 
    warn "Aucune séquence miyamotoi dans le FASTA."
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
 
fi
 
 
[[ -f "${TREE_PREFIX}.treefile" ]] \
    || die "IQ-TREE n'a pas produit le treefile."
 
 
ok "Arbre phylogénétique construit."
 
 
# =============================================================================
# 21 — RÉSUMÉ FINAL
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
echo "Locus retenu            : ${BEST_LOCUS}"
echo "Candidats Kraken        : ${POSITIVE_COUNT}"
echo "Échantillons dans arbre : ${SAMPLE_COUNT}"
echo "Références              : ${REF_COUNT}"
echo "Séquences finales       : ${FINAL_COUNT}"
echo "Outgroup                : ${OUTGROUP_ID:-non défini}"
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
