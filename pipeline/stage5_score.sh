#!/bin/bash
# =============================================================================
# stage5_score.sh -- Score each track's blind candidate list against the BNC
#                    background model (call_mrd_markers.py).
#
# Usage:
#   bash stage5_score.sh <sample_name> <output_dir>
#
# Inputs (per track in SCORE_TRACKS, default "dcs sscs"):
#   <output_dir>/<sample>/consensus/<sample>.<track>.sc.sorted.bam
#   <output_dir>/<sample>/variants/<sample>.<track>.vcf
#   <output_dir>/<sample>/variants/<sample>.<track>.indels.tsv
#   <output_dir>/<sample>/annotated/<sample>.<track>.filtered.tsv         (labels only, optional)
#   <output_dir>/<sample>/annotated/<sample>.<track>.indels.annotated.tsv (labels only, optional)
#
# Outputs:
#   <output_dir>/<sample>/scored/<sample>.<track>.candidates.tsv
#   <output_dir>/<sample>/scored/<sample>.<track>.calls.tsv
#
# The candidate list is the union of the Pisces calls and the CIGAR indel
# scan, derived from the sample's own reads only. No diagnosis variant list
# is read here; marker overlay is a separate post-hoc step.
#
# The artifact mask is deliberately not passed: it was retired 27 Aug 2026
# (1% floor against a measured 0.005% background) and the per-site beta
# matrix models those positions.
#
# Config (config.sh):
#   BETA_MATRIX_DCS / BETA_MATRIX_SSCS         per-track background model
#   INDEL_BLOCKLIST_DCS / INDEL_BLOCKLIST_SSCS per-track indel recurrence
#   INDEL_MIN_CONTROLS   (default 6)           BNCs an indel must recur in to be blocked
#   SCORE_TRACKS         (default "dcs sscs")
#   SCORE_EXTRA_ARGS     (default empty)       extra call_mrd_markers.py arguments
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
if declare -F activate_conda >/dev/null 2>&1; then
    activate_conda anneal
fi
set -euo pipefail

SAMPLE="${1:?usage: stage5_score.sh <sample_name> <output_dir>}"
OUTPUT_DIR="${2:?usage: stage5_score.sh <sample_name> <output_dir>}"

SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
SCORED_DIR="${SAMPLE_DIR}/scored"
ERROR_MODEL_DIR="${ANNEAL_ROOT}/scripts/error_model"

SCORE_TRACKS="${SCORE_TRACKS:-dcs sscs}"
INDEL_MIN_CONTROLS="${INDEL_MIN_CONTROLS:-6}"
SCORE_EXTRA_ARGS="${SCORE_EXTRA_ARGS:-}"

mkdir -p "${SCORED_DIR}"

echo "=========================================="
echo "  Stage 5: Background scoring -- ${SAMPLE}"
echo "=========================================="
echo "  Tracks:            ${SCORE_TRACKS}"
echo "  Indel min controls: ${INDEL_MIN_CONTROLS}"

for track in ${SCORE_TRACKS}; do
    upper="${track^^}"
    model_var="BETA_MATRIX_${upper}"
    model="${!model_var:-}"
    bl_var="INDEL_BLOCKLIST_${upper}"
    bl="${!bl_var:-}"

    bam="${SAMPLE_DIR}/consensus/${SAMPLE}.${track}.sc.sorted.bam"
    vcf="${SAMPLE_DIR}/variants/${SAMPLE}.${track}.vcf"
    indels="${SAMPLE_DIR}/variants/${SAMPLE}.${track}.indels.tsv"
    annot="${SAMPLE_DIR}/annotated/${SAMPLE}.${track}.filtered.tsv"
    iannot="${SAMPLE_DIR}/annotated/${SAMPLE}.${track}.indels.annotated.tsv"
    cand="${SCORED_DIR}/${SAMPLE}.${track}.candidates.tsv"
    calls="${SCORED_DIR}/${SAMPLE}.${track}.calls.tsv"

    echo ""
    echo "--- ${track} ---"
    echo "  model:     ${model}"
    echo "  blocklist: ${bl:-none}"

    for f in "${bam}" "${vcf}" "${indels}"; do
        if [ ! -f "${f}" ]; then
            echo "ERROR: missing input: ${f}"
            exit 1
        fi
    done
    if [ -z "${model}" ] || [ ! -f "${model}" ]; then
        echo "ERROR: ${model_var} not set or not found: ${model}"
        exit 1
    fi
    if [ -n "${bl}" ] && [ ! -f "${bl}" ]; then
        echo "ERROR: ${bl_var} not found: ${bl}"
        exit 1
    fi

    annot_args=()
    [ -f "${annot}" ]  && annot_args+=(--annotated "${annot}")
    [ -f "${iannot}" ] && annot_args+=(--indels-annotated "${iannot}")

    python "${ERROR_MODEL_DIR}/build_candidates.py" \
        --sample "${SAMPLE}" --track "${track}" \
        --vcf "${vcf}" --indels "${indels}" \
        ${annot_args[@]+"${annot_args[@]}"} \
        --out "${cand}"

    if [ ! -s "${cand}" ]; then
        echo "  no candidates for ${track}; nothing to score"
        : > "${calls}"
        continue
    fi

    python "${ERROR_MODEL_DIR}/call_mrd_markers.py" \
        --sample "${SAMPLE}" \
        --bam "${bam}" \
        --markers "${cand}" \
        --model "${model}" \
        ${bl:+--indel-blocklist "${bl}"} \
        --indel-min-controls "${INDEL_MIN_CONTROLS}" \
        ${SCORE_EXTRA_ARGS} \
        --out "${calls}"

    # Tier is presentational; the call is unchanged. From VAF and evidence:
    #   high_vaf   DETECTED at >= TIER_HIGH_VAF %  (germline or a major clone)
    #   mid_vaf    DETECTED at >= TIER_MID_VAF %
    #   mrd        DETECTED below that
    #   mrd_floor  DETECTED indel at <= 3 reads (minimum evidence, no background model)
    if ! head -1 "${calls}" | grep -q $'\ttier$'; then
        awk -F'\t' -v OFS='\t' -v hi="${TIER_HIGH_VAF:-20}" -v mid="${TIER_MID_VAF:-5}" '
            NR == 1 { print $0, "tier"; next }
            { t = ""
              if ($18 == "DETECTED") {
                  v = $10 + 0
                  if (v >= hi) t = "high_vaf"
                  else if (v >= mid) t = "mid_vaf"
                  else if ($7 != "snv" && $8 <= 3) t = "mrd_floor"
                  else t = "mrd" }
              print $0, t }' "${calls}" > "${calls}.tmp" && mv "${calls}.tmp" "${calls}"
    fi
    # Readable view: label split into columns, protein-altering consequences
    # only, annotation joined; the calls table stays the record.
    report="${SCORED_DIR}/${SAMPLE}.${track}.report.tsv"
    python "${ERROR_MODEL_DIR}/report_calls.py" --calls "${calls}" --out "${report}" \
        ${annot_args[@]+"${annot_args[@]}"}
    echo "  calls: ${calls}"
    echo "  report: ${report}"
    awk -F'\t' 'NR==1 { for (i = 1; i <= NF; i++) if ($i == "call") ci = i; next }
                ci { n[$ci]++ }
                END { for (k in n) printf "    %-16s %d\n", k, n[k] }' "${calls}"
done

echo ""
echo "Stage 5 complete: ${SCORED_DIR}"
