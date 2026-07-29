#!/bin/bash
# =============================================================================
# stage2_variant_calling.sh -- Variant calling on consensus BAMs (Pisces)
#
# Runs Pisces on the SSCS and DCS consensus BAMs. Blind within the panel: the
# interval list comes from the panel BED, and no diagnosis loci are supplied.
#
# Two Pisces-specific steps:
#   1. Strip XV. anneal writes XV:Z:SSCS / XV:Z:DCS as a string; Pisces'
#      IsCollapsedRead() expects an integer and throws. Each BAM is
#      single-track, so XV carries nothing Pisces needs.
#   2. Panel-restricted masked genome. Pisces exhausts memory on the full
#      3,366-contig hg38 reference.
#
# Pisces is a .NET 2.0 app on a .NET 8 runtime, hence DOTNET_ROLL_FORWARD.
# Compute nodes have no dotnet on PATH, so it comes from the conda env.
#
# Pisces emits a gVCF with a record at every interval position; a variants-only
# VCF is written alongside for stage 3.
#
# Usage:  bash stage2_variant_calling.sh <sample_name> <output_dir>
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
echo "  Stage 2: Variant calling (Pisces) -- ${SAMPLE}"
echo "  $(date)"
echo "================================================================"

export DOTNET_ROLL_FORWARD=Major
export PATH="${DOTNET_DIR}:${PATH}"

command -v dotnet >/dev/null || {
    echo "FATAL: dotnet not found on $(hostname); expected under ${DOTNET_DIR}"; exit 3; }
[ -f "${PISCES_DLL}" ] || { echo "FATAL: Pisces not found: ${PISCES_DLL}"; exit 3; }
[ -f "${PISCES_GENOME}/GenomeSize.xml" ] || {
    echo "FATAL: GenomeSize.xml missing in ${PISCES_GENOME}"; exit 3; }

build_intervals() {
    local bam="$1" out="$2"
    samtools view -H "${bam}" | grep -E '^@HD|^@SQ' > "${out}"
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
            name = (NF>=4 ? $4 : $1":"$2"-"$3)
            print $1, $2+1, $3, "+", name
         }' "${BEDFILE}" >> "${out}"
}

call_variants() {
    local bam="$1" label="$2"

    if [ ! -f "${bam}" ]; then
        echo "WARNING: BAM not found, skipping: ${bam}"
        return
    fi

    local pbam="${VARIANT_DIR}/${SAMPLE}.${label}.pisces.bam"
    local intervals="${VARIANT_DIR}/${SAMPLE}.${label}.interval_list"
    local outdir="${VARIANT_DIR}/pisces_${label}"
    local gvcf="${outdir}/${SAMPLE}.${label}.pisces.genome.vcf"
    local final_gvcf="${VARIANT_DIR}/${SAMPLE}.${label}.pisces.genome.vcf"
    local vcf="${VARIANT_DIR}/${SAMPLE}.${label}.vcf"

    echo "[$(date '+%H:%M:%S')] ${SAMPLE} ${label}: stripping XV tag..."
    samtools view -h --remove-tag XV "${bam}" -b -o "${pbam}"
    samtools index "${pbam}"

    build_intervals "${pbam}" "${intervals}"
    echo "[$(date '+%H:%M:%S')] ${SAMPLE} ${label}: $(grep -vc '^@' "${intervals}") intervals"

    rm -rf "${outdir}"; mkdir -p "${outdir}"

    echo "[$(date '+%H:%M:%S')] ${SAMPLE} ${label}: running Pisces..."
    dotnet "${PISCES_DLL}" \
        -bam "${pbam}" -g "${PISCES_GENOME}" -i "${intervals}" \
        --minbq "${PISCES_MINBQ}" --minmq "${PISCES_MINMQ}" \
        --minvf "${PISCES_MINVF}" -c "${PISCES_MINCOV}" \
        --filterduplicates false -o "${outdir}"

    if [ ! -f "${gvcf}" ]; then
        echo "ERROR: expected gVCF not produced: ${gvcf}"; ls -la "${outdir}"; return 1
    fi
    mv "${gvcf}" "${final_gvcf}"

    awk -F'\t' '/^#/ || ($5 != "." && $5 != "<M>")' "${final_gvcf}" > "${vcf}"

    local n_all n_pass
    n_all=$(awk -F'\t' '!/^#/' "${vcf}" | wc -l)
    n_pass=$(awk -F'\t' '!/^#/ && $7=="PASS"' "${vcf}" | wc -l)
    echo "[$(date '+%H:%M:%S')] ${SAMPLE} ${label}: ${n_all} variants (${n_pass} PASS) -> ${vcf}"

    rm -rf "${outdir}" "${pbam}" "${pbam}.bai" "${intervals}"
}

for label in sscs dcs; do
    bam="${CONSENSUS_DIR}/${SAMPLE}.${label}.sc.sorted.bam"
    [ -f "${bam}" ] || bam="${CONSENSUS_DIR}/${label}.sc.sorted.bam"
    call_variants "${bam}" "${label}"
done

echo ""
echo "[$(date '+%H:%M:%S')] Stage 2 complete for ${SAMPLE}"
ls -lh "${VARIANT_DIR}/"*.vcf 2>/dev/null
