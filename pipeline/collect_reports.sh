#!/bin/bash
# =============================================================================
# collect_reports.sh -- Gather a batch's readable outputs into one zip for
#                       download.
#
# Usage:
#   bash pipeline/collect_reports.sh <outdir> [zip_name]
#
# Collects, for every sample under <outdir>, both tracks:
#   scored/*.report.tsv           readable calls (protein-altering, split, annotated)
#   scored/*.calls.tsv            full calls table with tier
#   annotated/*.filtered.tsv      stage 3 substitution table
#   annotated/*.clinical.tsv      stage 3 clinical table (+ VariantValidator columns)
#   annotated/*.indels.annotated.tsv
#   flt3/*.tsv                    getITD
#   consensus/*stats.txt          consensus statistics (QC)
#
# Writes ~/inbox/to_excel/<zip_name>.zip (default <batch>_reports_<date>).
# VCFs, BAMs and HTML are not included.
# =============================================================================
set -euo pipefail

OUTDIR="${1:?usage: collect_reports.sh <outdir> [zip_name]}"
BATCH="$(basename "${OUTDIR}")"
ZIP_NAME="${2:-${BATCH}_reports_$(date +%Y%m%d)}"
DEST="${HOME}/inbox/to_excel"
mkdir -p "${DEST}"
ZIP="${DEST}/${ZIP_NAME}.zip"

cd "${OUTDIR}"
n=$(ls -d */scored 2>/dev/null | wc -l)
[ "${n}" -gt 0 ] || { echo "ERROR: no scored/ directories under ${OUTDIR}"; exit 1; }

rm -f "${ZIP}"
zip -q "${ZIP}" \
    */scored/*.report.tsv \
    */scored/*.calls.tsv \
    */annotated/*.filtered.tsv \
    */annotated/*.clinical.tsv \
    */annotated/*.indels.annotated.tsv \
    */flt3/*.tsv \
    */consensus/*stats.txt 2>/dev/null || true

echo "samples: ${n}"
echo "files:   $(unzip -l "${ZIP}" | tail -1 | awk '{print $2}')"
echo "zip:     ${ZIP} ($(du -h "${ZIP}" | cut -f1))"
