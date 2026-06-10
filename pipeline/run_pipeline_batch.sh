#!/bin/bash
# =============================================================================
# run_pipeline_batch.sh -- Run the pipeline on all samples in a manifest.
#
# Usage:
#   bash run_pipeline_batch.sh <manifest.tsv> <output_dir> [options]
#
# Options:
#   --stages 1,2    Run only these stages (comma-separated)
#   --annotate      Run variant annotation (stage 3)
#   --skip-vv       Skip VariantValidator in annotation
#
# Profiles:
#   Core (default): bash run_pipeline_batch.sh manifest.tsv outdir
#   With annotation: bash run_pipeline_batch.sh manifest.tsv outdir --annotate --skip-vv
#
# Manifest format (tab-separated, header required):
#   sample_name    fastq1                                 fastq2
#   25NGS1071      /path/to/25NGS1071_R1.fastq.gz         /path/to/25NGS1071_R2.fastq.gz
#
# Generate a manifest with:
#   anneal manifest --dir /path/to/fastqs/ -o manifest.tsv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
activate_conda

EXTRA_ARGS=()
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotate|--skip-vv)
            EXTRA_ARGS+=("$1")
            shift
            ;;
        --stages)
            EXTRA_ARGS+=("$1" "$2")
            shift 2
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

set -- "${POSITIONAL[@]}"

if [ $# -ne 2 ]; then
    echo "Usage: $0 <manifest.tsv> <output_dir> [--stages 1,2] [--annotate] [--skip-vv]"
    exit 1
fi

MANIFEST="$1"
OUTPUT_DIR="$2"

if [ ! -f "${MANIFEST}" ]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Copy manifest for reproducibility
cp "${MANIFEST}" "${OUTPUT_DIR}/manifest.tsv"

TOTAL=0
SUCCESS=0
FAILED=0
FAILED_SAMPLES=""
BATCH_START=$(date +%s)

echo "################################################################"
echo "  Anneal 0.1.0 -- Batch Pipeline"
echo "  Manifest: ${MANIFEST}"
echo "  Output:   ${OUTPUT_DIR}/"
echo "  Options:  ${EXTRA_ARGS[*]:-none}"
echo "  Started:  $(date)"
echo "################################################################"
echo ""

# Count samples (skip header)
N_SAMPLES=$(tail -n +2 "${MANIFEST}" | grep -cv "^$" || echo 0)
echo "Samples to process: ${N_SAMPLES}"
echo ""

# ---- Process each sample ----
while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2; do
    # Skip header
    if [ "${SAMPLE}" = "sample_name" ]; then
        continue
    fi

    # Skip empty lines
    if [ -z "${SAMPLE}" ]; then
        continue
    fi

    TOTAL=$((TOTAL + 1))

    echo ""
    echo "============================================================"
    echo "  [${TOTAL}/${N_SAMPLES}] ${SAMPLE}"
    echo "============================================================"

    if bash "${SCRIPT_DIR}/run_pipeline.sh" "${SAMPLE}" "${FASTQ1}" "${FASTQ2}" "${OUTPUT_DIR}" "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"; then
        SUCCESS=$((SUCCESS + 1))
    else
        FAILED=$((FAILED + 1))
        FAILED_SAMPLES="${FAILED_SAMPLES}  ${SAMPLE}\n"
        echo "FAILED: ${SAMPLE}"
    fi

done < "${MANIFEST}"

BATCH_END=$(date +%s)
ELAPSED=$((BATCH_END - BATCH_START))
HOURS=$((ELAPSED / 3600))
MINUTES=$(( (ELAPSED % 3600) / 60 ))

# ---- Write batch summary ----
SUMMARY="${OUTPUT_DIR}/batch_summary.txt"
cat > "${SUMMARY}" << EOF
# Anneal 0.1.0 Batch Summary
Date: $(date '+%Y-%m-%d %H:%M:%S')
Manifest: ${MANIFEST}
Options: ${EXTRA_ARGS[*]:-none}
Total samples: ${TOTAL}
Completed: ${SUCCESS}
Failed: ${FAILED}
Total time: ${ELAPSED}s (${HOURS}h ${MINUTES}m)
EOF

if [ ${FAILED} -gt 0 ]; then
    echo "" >> "${SUMMARY}"
    echo "Failed samples:" >> "${SUMMARY}"
    echo -e "${FAILED_SAMPLES}" >> "${SUMMARY}"
fi

echo ""
echo "################################################################"
echo "  Batch complete"
echo "  Processed: ${SUCCESS}/${TOTAL}"
if [ ${FAILED} -gt 0 ]; then
    echo "  Failed: ${FAILED}"
    echo -e "  ${FAILED_SAMPLES}"
fi
echo "  Time: ${HOURS}h ${MINUTES}m"
echo "  Summary: ${SUMMARY}"
echo "################################################################"
