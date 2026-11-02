#!/bin/bash

# ==============================================================================
#                               SLURM CONFIGURATION
# ==============================================================================

#SBATCH -p nadal-q
#SBATCH --chdir=/pool_sas/crh/RetroAtlas/Scripts/variantcalling
#SBATCH --job-name=Paco_test
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisco.iniestam1@um.es
#SBATCH --error=log/%j.err
#SBATCH --output=log/%j.out


# ==============================================================================
#                             ENVIRONMENT
# ==============================================================================

source ~/.bashrc
mamba activate sv_calling_env


# ==============================================================================
#                         USER CONFIGURATION
# ==============================================================================

# Reference genome
REF="/pool_sas/crh/ref/hg38.fa"

# Tandem repeats / VNTR file
TANDEM_REPEATS="SimpleRepeats.bed"

# Input CSV
CSV="/pool_sas/crh/RetroAtlas/bench.csv"

# Main output directory
OUTPATH="/pool_sas/crh/RetroAtlas/Scripts/variantcalling/bench"


# ------------------------------------------------------------------------------
# SUPPORT THRESHOLDS
# ------------------------------------------------------------------------------

# Minimum read support for *_full.vcf.gz
MINSUP=1

# Minimum read support for *_filtered.vcf.gz
FILTERSUP=3


# ------------------------------------------------------------------------------
# SV FILTERING
# ------------------------------------------------------------------------------

# Minimum insertion length
MIN_SVLEN=30


# ------------------------------------------------------------------------------
# OTHER PARAMETERS
# ------------------------------------------------------------------------------

THREADS="${SLURM_CPUS_PER_TASK:-24}"

# Parameters used by merge.R
MERGE_PARAM=64
MERGE_DISTANCE=500


# ==============================================================================
#                              CHECK PARAMETERS
# ==============================================================================

if (( FILTERSUP < MINSUP )); then
    echo "ERROR: FILTERSUP ($FILTERSUP) cannot be lower than MINSUP ($MINSUP)." >&2
    exit 1
fi

if (( MIN_SVLEN < 1 )); then
    echo "ERROR: MIN_SVLEN must be >= 1." >&2
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference not found: $REF" >&2
    exit 1
fi

if [[ ! -f "$CSV" ]]; then
    echo "ERROR: CSV not found: $CSV" >&2
    exit 1
fi

if [[ ! -f "$TANDEM_REPEATS" ]]; then
    echo "ERROR: Tandem repeat BED not found: $TANDEM_REPEATS" >&2
    exit 1
fi

mkdir -p "$OUTPATH"


# ==============================================================================
#                    FUNCTION: FILTER + COMPRESS VCF
# ==============================================================================
#
# All final VCFs:
#
#   - contain only insertions
#   - have |SVLEN| >= MIN_SVLEN
#   - are bgzip compressed
#   - are tabix indexed
#
# ==============================================================================

filter_and_compress_vcf() {

    local INPUT_VCF="$1"
    local OUTPUT_VCF_GZ="$2"

    bcftools filter \
        -i "INFO/SVTYPE=\"INS\" && abs(INFO/SVLEN)>=${MIN_SVLEN}" \
        "$INPUT_VCF" \
        -Oz \
        -o "$OUTPUT_VCF_GZ"

    tabix -f -p vcf "$OUTPUT_VCF_GZ"
}


# ==============================================================================
#                         FUNCTION: RUN SNIFFLES2
# ==============================================================================

run_sniffles() {

    local BAM="$1"
    local SUPPORT="$2"
    local OUTPUT_VCF_GZ="$3"

    local RAW_VCF="${OUTPUT_VCF_GZ%.vcf.gz}.raw.vcf"

    rm -f "$RAW_VCF"

    sniffles \
        --threads "$THREADS" \
        --reference "$REF" \
        --output-rnames \
        --minsupport "$SUPPORT" \
        --minsvlen "$MIN_SVLEN" \
        --phase \
        --phase-conflict-threshold 0.1 \
        --input "$BAM" \
        --vcf "$RAW_VCF" \
        --tandem-repeats "$TANDEM_REPEATS"

    filter_and_compress_vcf \
        "$RAW_VCF" \
        "$OUTPUT_VCF_GZ"

    rm -f "$RAW_VCF"
}


# ==============================================================================
#                           FUNCTION: RUN CUTESV
# ==============================================================================

run_cutesv() {

    local BAM="$1"
    local SUPPORT="$2"
    local OUTPUT_VCF_GZ="$3"
    local TMP_DIR="$4"

    local RAW_VCF="${OUTPUT_VCF_GZ%.vcf.gz}.raw.vcf"

    rm -rf "$TMP_DIR"
    mkdir -p "$TMP_DIR"

    rm -f "$RAW_VCF"

    cuteSV \
        -t "$THREADS" \
        -s "$SUPPORT" \
        -l "$MIN_SVLEN" \
        --report_readid \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --genotype \
        "$BAM" \
        "$REF" \
        "$RAW_VCF" \
        "$TMP_DIR"

    filter_and_compress_vcf \
        "$RAW_VCF" \
        "$OUTPUT_VCF_GZ"

    rm -f "$RAW_VCF"

    rm -rf "$TMP_DIR"
}


# ==============================================================================
#                          FUNCTION: RUN SEVERUS
# ==============================================================================

run_severus() {

    local BAM="$1"
    local SUPPORT="$2"
    local OUTPUT_VCF_GZ="$3"
    local TMP_DIR="$4"

    rm -rf "$TMP_DIR"
    mkdir -p "$TMP_DIR"

    severus \
        --target-bam "$BAM" \
        --out-dir "$TMP_DIR" \
        --threads "$THREADS" \
        --min-sv-size "$MIN_SVLEN" \
        --min-support "$SUPPORT" \
        --output-read-ids \
        --vntr-bed "$TANDEM_REPEATS"

    local RAW_VCF="${TMP_DIR}/all_SVs/severus_all.vcf"

    if [[ ! -f "$RAW_VCF" ]]; then
        echo "ERROR: Severus VCF not found: $RAW_VCF" >&2
        return 1
    fi

    filter_and_compress_vcf \
        "$RAW_VCF" \
        "$OUTPUT_VCF_GZ"
}


# ==============================================================================
#                             PRINT PARAMETERS
# ==============================================================================

echo
echo "======================================================================"
echo " SV CALLING PIPELINE"
echo "======================================================================"
echo "Reference:          $REF"
echo "CSV:                $CSV"
echo "Output:             $OUTPATH"
echo "Threads:            $THREADS"
echo "MINSUP:             $MINSUP"
echo "FILTERSUP:          $FILTERSUP"
echo "Minimum SVLEN:      $MIN_SVLEN"
echo "======================================================================"


# ==============================================================================
#                              MAIN LOOP
# ==============================================================================
#
# Expected CSV format:
#
# BAM ID,BAM
# sample1,/path/to/sample1.bam
# sample2,/path/to/sample2.bam
#
# tail -n +2 skips the header.
#
# ==============================================================================

while IFS=, read -r bam_id bam_path; do

    # Remove possible Windows carriage return
    bam_path="${bam_path%$'\r'}"

    # Ignore empty lines
    [[ -z "$bam_path" ]] && continue


    # ==========================================================================
    #                             SAMPLE SETUP
    # ==========================================================================

    BAM="$bam_path"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM/CRAM not found. Skipping: $BAM" >&2
        continue
    fi

    NAME="$(basename "$BAM")"
    NAME="${NAME%.bam}"
    NAME="${NAME%.cram}"

    SAMPLE_DIR="${OUTPATH}/${NAME}"

    mkdir -p "$SAMPLE_DIR"


    echo
    echo "======================================================================"
    echo " SAMPLE: $NAME"
    echo " BAM ID: $bam_id"
    echo " BAM:    $BAM"
    echo " OUT:    $SAMPLE_DIR"
    echo "======================================================================"


    # ==========================================================================
    #                          OUTPUT FILE NAMES
    # ==========================================================================

    # Sniffles
    SNIFFLES_FULL="${SAMPLE_DIR}/${NAME}_sniffles2_full.vcf.gz"
    SNIFFLES_FILTERED="${SAMPLE_DIR}/${NAME}_sniffles2_filtered.vcf.gz"

    # cuteSV
    CUTESV_FULL="${SAMPLE_DIR}/${NAME}_cutesv_full.vcf.gz"
    CUTESV_FILTERED="${SAMPLE_DIR}/${NAME}_cutesv_filtered.vcf.gz"

    # Severus
    SEVERUS_FULL="${SAMPLE_DIR}/${NAME}_severus_full.vcf.gz"
    SEVERUS_FILTERED="${SAMPLE_DIR}/${NAME}_severus_filtered.vcf.gz"


    # ==========================================================================
    #                              SNIFFLES2
    # ==========================================================================

    echo
    echo "----------------------------------------------------------------------"
    echo " SNIFFLES2 - FULL"
    echo " Support >= $MINSUP"
    echo "----------------------------------------------------------------------"

    run_sniffles \
        "$BAM" \
        "$MINSUP" \
        "$SNIFFLES_FULL"


    echo
    echo "----------------------------------------------------------------------"
    echo " SNIFFLES2 - FILTERED"
    echo " Support >= $FILTERSUP"
    echo "----------------------------------------------------------------------"

    run_sniffles \
        "$BAM" \
        "$FILTERSUP" \
        "$SNIFFLES_FILTERED"


    # ==========================================================================
    #                                CUTESV
    # ==========================================================================

    echo
    echo "----------------------------------------------------------------------"
    echo " CUTESV - FULL"
    echo " Support >= $MINSUP"
    echo "----------------------------------------------------------------------"

    run_cutesv \
        "$BAM" \
        "$MINSUP" \
        "$CUTESV_FULL" \
        "${SAMPLE_DIR}/tmp_cutesv_full"


    echo
    echo "----------------------------------------------------------------------"
    echo " CUTESV - FILTERED"
    echo " Support >= $FILTERSUP"
    echo "----------------------------------------------------------------------"

    run_cutesv \
        "$BAM" \
        "$FILTERSUP" \
        "$CUTESV_FILTERED" \
        "${SAMPLE_DIR}/tmp_cutesv_filtered"


    # ==========================================================================
    #                                SEVERUS
    # ==========================================================================

    echo
    echo "----------------------------------------------------------------------"
    echo " SEVERUS - FULL"
    echo " Support >= $MINSUP"
    echo "----------------------------------------------------------------------"

    run_severus \
        "$BAM" \
        "$MINSUP" \
        "$SEVERUS_FULL" \
        "${SAMPLE_DIR}/tmp_severus_full"


    echo
    echo "----------------------------------------------------------------------"
    echo " SEVERUS - FILTERED"
    echo " Support >= $FILTERSUP"
    echo "----------------------------------------------------------------------"

    run_severus \
        "$BAM" \
        "$FILTERSUP" \
        "$SEVERUS_FILTERED" \
        "${SAMPLE_DIR}/tmp_severus_filtered"


    # ==========================================================================
    #                                 MERGE
    # ==========================================================================

    echo
    echo "----------------------------------------------------------------------"
    echo " MERGING FILTERED VCF FILES"
    echo "----------------------------------------------------------------------"

    Rscript merge.R \
        "$MERGE_PARAM" \
        "$SAMPLE_DIR" \
        "$SNIFFLES_FILTERED" \
        "$CUTESV_FILTERED" \
        "$SEVERUS_FILTERED" \
        "$MERGE_DISTANCE"


    # ==========================================================================
    #                            SAMPLE COMPLETED
    # ==========================================================================

    echo
    echo "======================================================================"
    echo " COMPLETED: $NAME"
    echo "======================================================================"
    echo " Sniffles FULL:      $SNIFFLES_FULL"
    echo " Sniffles FILTERED:  $SNIFFLES_FILTERED"
    echo " cuteSV FULL:        $CUTESV_FULL"
    echo " cuteSV FILTERED:    $CUTESV_FILTERED"
    echo " Severus FULL:       $SEVERUS_FULL"
    echo " Severus FILTERED:   $SEVERUS_FILTERED"
    echo "======================================================================"

done < <(tail -n +2 "$CSV")


echo
echo "======================================================================"
echo " ALL SAMPLES COMPLETED"
echo "======================================================================"
