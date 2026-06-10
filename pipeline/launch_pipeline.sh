#!/bin/bash
# =============================================================================
# launch_pipeline.sh -- Start the batch pipeline in the background
#
# Usage:
#   bash launch_pipeline.sh <manifest.tsv> <output_dir>
#
# Runs the full pipeline detached via nohup so it survives SSH disconnects.
# Logs to <output_dir>/pipeline.log
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <manifest.tsv> <output_dir> [options]"
    echo ""
    echo "Options:  --stages 1,2, --annotate, --skip-vv"
    echo ""
    echo "Profiles:"
    echo "  Core (default):  $0 manifest.tsv outdir"
    echo "  With annotation: $0 manifest.tsv outdir --annotate --skip-vv"
    echo ""
    echo "Launches the full pipeline in the background (survives SSH disconnect)."
    echo "Monitor with: tail -f <output_dir>/pipeline.log"
    exit 1
fi

MANIFEST="$1"
OUTPUT_DIR="$2"
shift 2
EXTRA_ARGS="$*"

mkdir -p "${OUTPUT_DIR}"
LOG="${OUTPUT_DIR}/pipeline.log"

echo "================================================================"
echo "  Anneal 0.1.0 -- Launching pipeline"
echo "  Manifest: ${MANIFEST}"
echo "  Output:   ${OUTPUT_DIR}/"
echo "  Log:      ${LOG}"
echo "================================================================"

activate_conda

nohup bash "${SCRIPT_DIR}/run_pipeline_batch.sh" "${MANIFEST}" "${OUTPUT_DIR}" ${EXTRA_ARGS} \
    > "${LOG}" 2>&1 &

PID=$!
echo ""
echo "Pipeline started (PID: ${PID})"
echo "Monitor progress:"
echo "  tail -f ${LOG}"
echo ""
echo "Check if still running:"
echo "  ps -p ${PID}"
