#!/bin/bash
# =============================================================================
# run_pipeline.sh -- Anneal pipeline for a single sample
#
# Stage 1: Duplex consensus generation (Anneal)
# Stage 2: Variant calling (samtools mpileup + Rust caller)
# Stage 3: Variant annotation (VEP + ANNOVAR, optional)
#
# Usage:
#   bash run_pipeline.sh <sample_name> <fastq1> <fastq2> <output_dir> [options]
#
# Options:
#   --stages 1,2     Run only these stages (comma-separated)
#   --annotate       Include stage 3 annotation (stages become 1,2,3)
#   --skip-vv        Skip VariantValidator in stage 3
#
# Profiles:
#   Core (default):    (no flag)             (stages 1,2: consensus + variants)
#   With annotation:   --annotate            (stages 1,2,3)
#   Resume annotation: --stages 3 --skip-vv  (annotate existing VCFs)
#
# Examples:
#   # Consensus + variant calling only (default)
#   bash run_pipeline.sh SAMPLE R1 R2 outdir
#
#   # Full run including annotation
#   bash run_pipeline.sh SAMPLE R1 R2 outdir --annotate --skip-vv
#
#   # Re-run annotation on existing VCFs
#   bash run_pipeline.sh SAMPLE R1 R2 outdir --stages 3 --skip-vv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
activate_conda

# ---- Parse arguments ----
STAGES=""
SKIP_VV=""
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --stages)
            STAGES="$2"
            shift 2
            ;;
        --annotate)
            if [ -z "${STAGES}" ]; then
                STAGES="1,2,3"
            else
                STAGES="${STAGES},3"
            fi
            shift
            ;;
        --skip-vv)
            SKIP_VV="--skip-vv"
            shift
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

set -- "${POSITIONAL[@]}"

# Default: stages 1,2 (no annotation unless requested)
if [ -z "${STAGES}" ]; then
    STAGES="1,2"
fi

if [ $# -ne 4 ]; then
    echo "Usage: $0 <sample_name> <fastq1> <fastq2> <output_dir> [options]"
    echo ""
    echo "Options:"
    echo "  --stages 1,2    Run only these stages (comma-separated)"
    echo "  --annotate      Include stage 3 annotation (stages become 1,2,3)"
    echo "  --skip-vv       Skip VariantValidator in annotation"
    echo ""
    echo "Profiles:"
    echo "  Core (default):    $0 SAMPLE R1 R2 outdir"
    echo "  With annotation:   $0 SAMPLE R1 R2 outdir --annotate --skip-vv"
    echo "  Resume annotation: $0 SAMPLE R1 R2 outdir --stages 3 --skip-vv"
    exit 1
fi

SAMPLE="$1"
FASTQ1="$2"
FASTQ2="$3"
OUTPUT_DIR="$4"

# Helper: check if a stage is in the comma-separated list
run_stage() {
    echo ",${STAGES}," | grep -q ",${1},"
}

PIPELINE_START=$(date +%s)

echo "################################################################"
echo "  Anneal 0.1.0 -- Pipeline"
echo "  Sample:  ${SAMPLE}"
echo "  Output:  ${OUTPUT_DIR}/${SAMPLE}/"
echo "  Stages:  ${STAGES}"
echo "  Started: $(date)"
echo "################################################################"
echo ""

# ---- Stage 1: Consensus ----
if run_stage 1; then
    bash "${SCRIPT_DIR}/stage1_consensus.sh" "${SAMPLE}" "${FASTQ1}" "${FASTQ2}" "${OUTPUT_DIR}"
    echo ""
else
    echo "[SKIP] Stage 1: consensus generation"
    echo ""
fi

# ---- Stage 2: Variant calling ----
if run_stage 2; then
    bash "${SCRIPT_DIR}/stage2_variant_calling.sh" "${SAMPLE}" "${OUTPUT_DIR}"
    echo ""
else
    echo "[SKIP] Stage 2: variant calling"
    echo ""
fi

# ---- Stage 3: Annotation ----
if run_stage 3; then
    bash "${SCRIPT_DIR}/stage3_annotate.sh" "${SAMPLE}" "${OUTPUT_DIR}" ${SKIP_VV}
    echo ""
else
    echo "[SKIP] Stage 3: variant annotation"
    echo ""
fi

PIPELINE_END=$(date +%s)
ELAPSED=$((PIPELINE_END - PIPELINE_START))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo "################################################################"
echo "  Pipeline complete: ${SAMPLE}"
echo "  Stages:  ${STAGES}"
echo "  Elapsed: ${MINUTES}m ${SECONDS}s"
echo "  Finished: $(date)"
echo "################################################################"
