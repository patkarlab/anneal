#!/bin/bash
# =============================================================================
# stage2_variant_calling.sh -- Variant calling on consensus BAMs
#
# Runs samtools mpileup followed by the Rust variant caller on both
# SSCS and DCS BAMs for a given sample.
#
# Usage:
#   bash stage2_variant_calling.sh <sample_name> <output_dir>
#
# Expects Stage 1 outputs at:
#   <output_dir>/<sample>/consensus/<sample>.sscs.sc.sorted.bam
#   <output_dir>/<sample>/consensus/<sample>.dcs.sc.sorted.bam
#
# Produces:
#   <output_dir>/<sample>/variants/<sample>.sscs.vcf
#   <output_dir>/<sample>/variants/<sample>.dcs.vcf
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

if [ $# -ne 2 ]; then
    echo "Usage: $0 <sample_name> <output_dir>"
    exit 1
fi

SAMPLE="$1"
OUTPUT_DIR="$2"

CONSENSUS_DIR="${OUTPUT_DIR}/${SAMPLE}/consensus"
VARIANT_DIR="${OUTPUT_DIR}/${SAMPLE}/variants"
mkdir -p "${VARIANT_DIR}"

echo "================================================================"
echo "  Stage 2: Variant calling -- ${SAMPLE}"
echo "  $(date)"
echo "================================================================"

# ---- Function: mpileup + Rust variant caller for one BAM ----
call_variants() {
    local bam="$1"
    local label="$2"
    local mpileup="${VARIANT_DIR}/${SAMPLE}.${label}.mpileup"
    local vcf="${VARIANT_DIR}/${SAMPLE}.${label}.vcf"

    if [ ! -f "${bam}" ]; then
        echo "WARNING: BAM not found, skipping: ${bam}"
        return
    fi

    echo "[$(date '+%H:%M:%S')] ${SAMPLE} ${label}: generating mpileup..."
    samtools mpileup \
        -A \
        -l "${BEDFILE}" \
        -f "${REFERENCE}" \
        -d "${MAX_DEPTH}" \
        "${bam}" > "${mpileup}"

    echo "[$(date '+%H:%M:%S')] ${SAMPLE} ${label}: calling variants..."
    "${VARIANT_CALLER}" "${mpileup}" "${MIN_BASE_QUAL}" "${vcf}"

    local variant_count
    variant_count=$(grep -cv "^#" "${vcf}" 2>/dev/null || echo 0)
    echo "[$(date '+%H:%M:%S')] ${SAMPLE} ${label}: ${variant_count} variant records -> ${vcf}"

    # Clean up large intermediate mpileup
    rm -f "${mpileup}"
}

# ---- Run on both consensus types ----
SSCS_BAM="${CONSENSUS_DIR}/${SAMPLE}.sscs.sc.sorted.bam"
DCS_BAM="${CONSENSUS_DIR}/${SAMPLE}.dcs.sc.sorted.bam"

call_variants "${SSCS_BAM}" "sscs"
call_variants "${DCS_BAM}" "dcs"

echo ""
echo "[$(date '+%H:%M:%S')] Stage 2 complete for ${SAMPLE}"
ls -lh "${VARIANT_DIR}/"*.vcf 2>/dev/null
