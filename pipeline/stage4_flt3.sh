#!/bin/bash
# =============================================================================
# stage4_flt3.sh -- FLT3-ITD detection on consensus BAMs
#
# Extracts the FLT3 exon 14/15 window from the DCS consensus BAM and runs
# getITD on it. Emits getITD's output unfiltered; matching against a patient's
# known ITD happens downstream.
#
# Consensus reads are single-end (read pairs are merged during collapsing), so
# getITD is given one FASTQ. Tools requiring paired input, such as
# FLT3_ITD_ext, cannot read consensus BAMs at all.
#
# -min_read_copies 1: getITD deduplicates identical reads by default, which
# assumes raw data full of PCR duplicates. UMI families already did that, and
# leaving it on discards real signal.
#
# Primer filtering is left ON. It is nominally an amplicon-data filter, but on
# hybrid capture it removes reads with indels near read edges, which are the
# main source of spurious ITD calls. Turning it off quadruples usable reads and
# increases false positives.
#
# Note on discovery: untargeted getITD on DCS returns 1-5 spurious ITDs per
# sample with read counts indistinguishable from a true call. Filter the output
# to the patient's known ITD length and insertion site.
#
# getITD reports hg19 coordinates against its own FLT3 reference, while the BAM
# is queried in the build of BEDFILE. Both are correct; they are different
# coordinate spaces.
#
# Usage:
#   bash stage4_flt3.sh <sample_name> <output_dir>
#
# Expects Stage 1 output at:
#   <output_dir>/<sample>/consensus/<sample>.dcs.sc.sorted.bam
#
# Produces:
#   <output_dir>/<sample>/flt3/<sample>.flt3_itds.tsv
#   <output_dir>/<sample>/flt3/<sample>.flt3_insertions.tsv
#   <output_dir>/<sample>/flt3/<sample>.flt3_stats.txt
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
FLT3_DIR="${OUTPUT_DIR}/${SAMPLE}/flt3"
mkdir -p "${FLT3_DIR}"

echo "================================================================"
echo "  Stage 4: FLT3-ITD -- ${SAMPLE}"
echo "  $(date)"
echo "================================================================"

# ---- Locate the DCS BAM (stage 1 prefixes with the sample name; a bare
#      `anneal consensus` run does not) ----
DCS_BAM="${CONSENSUS_DIR}/${SAMPLE}.dcs.sc.sorted.bam"
if [ ! -f "${DCS_BAM}" ]; then
    DCS_BAM="${CONSENSUS_DIR}/dcs.sc.sorted.bam"
fi
if [ ! -f "${DCS_BAM}" ]; then
    echo "WARNING: no DCS BAM under ${CONSENSUS_DIR}, skipping stage 4"
    exit 0
fi

if [ ! -d "${GETITD_DIR}" ]; then
    echo "WARNING: getITD not found at ${GETITD_DIR}, skipping stage 4"
    echo "  install: git clone https://github.com/tjblaette/getitd.git"
    exit 0
fi

# ---- FLT3 region from BED columns 2-3. Column 4 is a legacy label and is
#      used only to select the rows. ----
REGION=$(awk -F'\t' -v pat="${FLT3_PROBE_PATTERN}" -v fl="${FLT3_FLANK}" '
    $4 ~ pat {
        if (min == "" || $2 < min) min = $2
        if ($3 > max) max = $3
        chrom = $1
    }
    END {
        if (chrom == "") exit 1
        lo = min - fl; if (lo < 0) lo = 0
        print chrom ":" lo "-" (max + fl)
    }' "${BEDFILE}")

if [ -z "${REGION}" ]; then
    echo "ERROR: no probe in ${BEDFILE} column 4 matching '${FLT3_PROBE_PATTERN}'"
    exit 1
fi
echo "[$(date '+%H:%M:%S')] FLT3 region: ${REGION}"

# ---- Extract region reads to single-end FASTQ ----
TMP=$(mktemp -d "${FLT3_DIR}/tmp.XXXXXX")
trap 'rm -rf "${TMP}"' EXIT

FQ="${TMP}/${SAMPLE}.flt3.fq"
samtools view -@ "${FLT3_THREADS}" -b "${DCS_BAM}" "${REGION}" -o "${TMP}/region.bam"
samtools sort -@ "${FLT3_THREADS}" -n "${TMP}/region.bam" -o "${TMP}/region.ns.bam"
samtools fastq -@ "${FLT3_THREADS}" -n \
    -0 /dev/null -s "${FQ}" -1 /dev/null -2 /dev/null \
    "${TMP}/region.ns.bam" > /dev/null 2>&1

NREADS=$(( $(wc -l < "${FQ}") / 4 ))
echo "[$(date '+%H:%M:%S')] ${NREADS} consensus reads over FLT3"

# getITD runs from its own directory, so the FASTQ path must be absolute
FQ_ABS="$(cd "$(dirname "${FQ}")" && pwd)/$(basename "${FQ}")"

if [ "${NREADS}" -lt "${FLT3_MIN_READS}" ]; then
    echo "WARNING: only ${NREADS} reads (< ${FLT3_MIN_READS}), skipping getITD"
    echo "not_evaluable: ${NREADS} reads" > "${FLT3_DIR}/${SAMPLE}.flt3_stats.txt"
    exit 0
fi

# ---- getITD ----
RUN="${SAMPLE}_flt3"
echo "[$(date '+%H:%M:%S')] running getITD..."
(
    cd "${GETITD_DIR}"
    rm -rf "${RUN}_getitd"
    python getitd.py \
        -reference anno/amplicon.txt \
        -anno anno/amplicon_kayser.tsv \
        -min_read_copies 1 \
        -nkern "${FLT3_THREADS}" \
        "${RUN}" "${FQ_ABS}"
) || { echo "ERROR: getITD failed for ${SAMPLE}"; exit 1; }

GOUT="${GETITD_DIR}/${RUN}_getitd"
HC="itds_collapsed-is-same_is-similar_is-close_is-same_trailing_hc.tsv"
INS="insertions_collapsed-is-same_is-similar_is-close_is-same_trailing_hc.tsv"

if [ -f "${GOUT}/${HC}" ]; then
    cp "${GOUT}/${HC}" "${FLT3_DIR}/${SAMPLE}.flt3_itds.tsv"
else
    # getITD writes no file when nothing survives filtering
    echo -e "sample\tlength\tstart\tvaf\tar\tcoverage\tcounts\ttrailing\tseq" \
        > "${FLT3_DIR}/${SAMPLE}.flt3_itds.tsv"
fi
[ -f "${GOUT}/${INS}" ] && cp "${GOUT}/${INS}" "${FLT3_DIR}/${SAMPLE}.flt3_insertions.tsv"
[ -f "${GOUT}/stats.txt" ] && cp "${GOUT}/stats.txt" "${FLT3_DIR}/${SAMPLE}.flt3_stats.txt"
rm -rf "${GOUT}"

NITD=$(( $(wc -l < "${FLT3_DIR}/${SAMPLE}.flt3_itds.tsv") - 1 ))
echo ""
echo "[$(date '+%H:%M:%S')] Stage 4 complete for ${SAMPLE}: ${NITD} ITD call(s)"
if [ "${NITD}" -gt 0 ]; then
    echo "  length  insertion_site  vaf%      reads"
    awk -F'\t' 'NR>1 {printf "  %-7s %-15s %-9s %s\n", $2, $19, $4, $7}' \
        "${FLT3_DIR}/${SAMPLE}.flt3_itds.tsv"
    echo ""
    echo "  Coordinates are hg19. Match against the patient's known ITD length"
    echo "  and insertion site; unmatched calls are not reliable."
fi
ls -lh "${FLT3_DIR}/"*.tsv 2>/dev/null
