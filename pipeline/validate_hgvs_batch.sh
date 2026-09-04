#!/bin/bash
# =============================================================================
# validate_hgvs_batch.sh -- HGVS validation through the public VariantValidator
#                           API, post hoc, on the login node.
#
# Usage:
#   bash validate_hgvs_batch.sh <outdir> [sample ...]
#
# Runs scripts/validate_hgvs.py on every <sample>.<track>.clinical.tsv under
# <outdir>/<sample>/annotated/ (all samples in <outdir> if none are named),
# sharing one cache so each unique HGVS query hits the API once across the
# batch. About one request per second; a cohort's unique clinical variants
# take minutes, not hours.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
if declare -F activate_conda >/dev/null 2>&1; then
    activate_conda anneal
fi
set -euo pipefail

OUTPUT_DIR="${1:?usage: validate_hgvs_batch.sh <outdir> [sample ...]}"
shift
if [ "$#" -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=()
    for d in "${OUTPUT_DIR}"/*/annotated; do
        [ -d "${d}" ] && SAMPLES+=("$(basename "$(dirname "${d}")")")
    done
fi

VV_URL="${VV_URL:-https://rest.variantvalidator.org}"
VV_CACHE="${VV_CACHE:-${ANNEAL_ROOT}/vv_cache.json}"
echo "VariantValidator: ${VV_URL}   cache: ${VV_CACHE}   samples: ${#SAMPLES[@]}"

for SAMPLE in "${SAMPLES[@]}"; do
    for track in dcs sscs; do
        tsv="${OUTPUT_DIR}/${SAMPLE}/annotated/${SAMPLE}.${track}.clinical.tsv"
        [ -f "${tsv}" ] || continue
        n=$(($(wc -l < "${tsv}") - 1))
        [ "${n}" -gt 0 ] || { echo "--- ${SAMPLE} ${track}: no clinical variants"; continue; }
        echo "--- ${SAMPLE} ${track}: ${n} variants"
        python3 "${ANNEAL_ROOT}/scripts/validate_hgvs.py" \
            -i "${tsv}" -o "${OUTPUT_DIR}/${SAMPLE}/annotated" \
            --vv-url "${VV_URL}" --cache "${VV_CACHE}"
    done
done
echo "done: $(date)"
