#!/bin/bash
# =============================================================================
# stage1_consensus.sh -- Run Anneal consensus generation for a single sample
#
# Usage:
#   bash stage1_consensus.sh <sample_name> <fastq1> <fastq2> <output_dir>
#
# Produces:
#   <output_dir>/<sample>/consensus/<sample>.sscs.sc.sorted.bam
#   <output_dir>/<sample>/consensus/<sample>.dcs.sc.sorted.bam
#   <output_dir>/<sample>/consensus/<sample>.stats.txt
#   <output_dir>/<sample>/consensus/<sample>.family_sizes.tsv
#   <output_dir>/<sample>/consensus/<sample>.family_sizes.png
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

if [ $# -ne 4 ]; then
    echo "Usage: $0 <sample_name> <fastq1> <fastq2> <output_dir>"
    exit 1
fi

SAMPLE="$1"
FASTQ1="$2"
FASTQ2="$3"
OUTPUT_DIR="$4"

SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
CONSENSUS_DIR="${SAMPLE_DIR}/consensus"

echo "================================================================"
echo "  Stage 1: Anneal consensus -- ${SAMPLE}"
echo "  $(date)"
echo "================================================================"

# ---- Run Anneal ----
SC_FLAG=""
if [ "${SINGLETON_CORRECTION}" = true ]; then
    SC_FLAG="--singleton-correction"
fi

GPU_FLAG=""
if [ "${USE_GPU}" = false ]; then
    GPU_FLAG="--no-gpu"
fi

"${ANNEAL}" run \
    --fastq1 "${FASTQ1}" \
    --fastq2 "${FASTQ2}" \
    --reference "${REFERENCE}" \
    --output "${SAMPLE_DIR}" \
    --aligner "${ALIGNER}" \
    --bpattern "${BPATTERN}" \
    --cutoff "${CUTOFF}" \
    ${SC_FLAG} \
    ${GPU_FLAG} \
    -vv

# ---- Rename outputs with sample prefix ----
# Anneal outputs generic names; we rename for multi-sample clarity
cd "${CONSENSUS_DIR}"
for bam in sscs.sc.sorted.bam dcs.sc.sorted.bam; do
    if [ -f "${bam}" ] && [ ! -f "${SAMPLE}.${bam}" ]; then
        mv "${bam}" "${SAMPLE}.${bam}"
        [ -f "${bam}.bai" ] && mv "${bam}.bai" "${SAMPLE}.${bam}.bai"
    fi
done

# ---- Rename auxiliary outputs with sample prefix ----
for f in stats.txt family_sizes.tsv; do
    if [ -f "${f}" ] && [ ! -f "${SAMPLE}.${f}" ]; then
        mv "${f}" "${SAMPLE}.${f}"
    fi
done

# ---- Clean up intermediate BAMs (keep only the 2 recommended outputs) ----
for bam in sscs.sorted.bam singleton.sorted.bam sscs.rescue.sorted.bam \
           singleton.rescue.sorted.bam sscs.singleton.sorted.bam \
           rescue.remaining.sorted.bam; do
    if [ -f "${bam}" ]; then
        rm -f "${bam}" "${bam}.bai"
        echo "[cleanup] Removed intermediate: ${bam}"
    fi
done

# ---- Generate family size distribution plot ----
if [ -f "${SAMPLE}.family_sizes.tsv" ] && [ -f "${PLOT_SCRIPT}" ]; then
    echo "[$(date '+%H:%M:%S')] Generating family size distribution plot..."
    python3 "${PLOT_SCRIPT}" \
        --input "${SAMPLE}.family_sizes.tsv" \
        --output "${SAMPLE}.family_sizes.png" \
        --sample "${SAMPLE}" 2>/dev/null || echo "WARNING: plot generation failed (non-critical)"
fi

echo ""
echo "[$(date '+%H:%M:%S')] Stage 1 complete for ${SAMPLE}"
echo "  SSCS BAM:     ${CONSENSUS_DIR}/${SAMPLE}.sscs.sc.sorted.bam"
echo "  DCS BAM:      ${CONSENSUS_DIR}/${SAMPLE}.dcs.sc.sorted.bam"
echo "  Statistics:   ${CONSENSUS_DIR}/${SAMPLE}.stats.txt"
echo "  Family sizes: ${CONSENSUS_DIR}/${SAMPLE}.family_sizes.tsv"
echo "  Plot:         ${CONSENSUS_DIR}/${SAMPLE}.family_sizes.png"
