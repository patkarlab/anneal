#!/bin/bash
# =============================================================================
# make_pisces_reference.sh -- Build the panel-restricted reference Pisces needs.
#
# Usage:
#   bash scripts/make_pisces_reference.sh [output_dir]
#
# Pisces on the full 3,366-contig hg38 runs out of memory; on a reference
# holding only the chromosomes the panel touches it does not. This writes:
#   <output_dir>/genome.fa, genome.fa.fai, genome.dict, GenomeSize.xml
# Chromosomes come from column 1 of the panel BED; sequence from the masked
# reference (the same one used for alignment). GenomeSize.xml is produced by
# Pisces' CreateGenomeSizeFile on the .NET runtime, as stage 2 runs Pisces.
#
# Reads REFERENCE, BEDFILE, PISCES_DLL, DOTNET_DIR, PISCES_GENOME from
# pipeline/config.sh; the output directory defaults to PISCES_GENOME.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../pipeline/config.sh"
if declare -F activate_conda >/dev/null 2>&1; then
    activate_conda anneal
fi
set -euo pipefail

OUT="${1:-${PISCES_GENOME}}"
PISCES_BIN_DIR="$(dirname "${PISCES_DLL}")"
CREATE_DLL="${PISCES_BIN_DIR}/CreateGenomeSizeFile.dll"
export DOTNET_ROOT="${DOTNET_DIR}"
export DOTNET_ROLL_FORWARD=Major
DOTNET="${DOTNET_DIR}/dotnet"

for f in "${REFERENCE}" "${BEDFILE}" "${CREATE_DLL}" "${DOTNET}"; do
    [ -e "${f}" ] || { echo "ERROR: not found: ${f}"; exit 1; }
done
mkdir -p "${OUT}"

CHROMS=$(cut -f1 "${BEDFILE}" | grep -v '^#' | sort -u | sort -V | tr '\n' ' ')
echo "panel chromosomes: ${CHROMS}"
echo "reference:         ${REFERENCE}"
echo "output:            ${OUT}"

samtools faidx "${REFERENCE}" ${CHROMS} > "${OUT}/genome.fa"
samtools faidx "${OUT}/genome.fa"
samtools dict "${OUT}/genome.fa" > "${OUT}/genome.dict"
"${DOTNET}" "${CREATE_DLL}" -g "${OUT}" -s "Homo sapiens (hg38 masked, panel)" > "${OUT}/CreateGenomeSizeFile.log" 2>&1

echo "contigs in GenomeSize.xml: $(grep -c '<chromosome ' "${OUT}/GenomeSize.xml")"
ls -l "${OUT}"
