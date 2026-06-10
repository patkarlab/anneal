#!/bin/bash
# =============================================================================
# stage3_annotate.sh -- Variant annotation pipeline
#
# Runs VEP + ANNOVAR annotation, clinical filtering, and HGVS validation
# on VCFs from the Rust mpileup variant caller (stage 2).
#
# Usage:
#   bash stage3_annotate.sh <sample_name> <output_dir> [--skip-vv]
#
# Expects Stage 2 outputs at:
#   <output_dir>/<sample>/variants/<sample>.dcs.vcf
#   <output_dir>/<sample>/variants/<sample>.sscs.vcf
#
# Produces:
#   <output_dir>/<sample>/annotated/<sample>.dcs.annotated.tsv
#   <output_dir>/<sample>/annotated/<sample>.dcs.filtered.tsv
#   <output_dir>/<sample>/annotated/<sample>.dcs.clinical.tsv
#   <output_dir>/<sample>/annotated/<sample>.dcs.validated.tsv
#   (same for sscs)
#
# Requires:
#   - VEP in conda env 'vep'
#   - ANNOVAR with hg38 databases
#   - VariantValidator Docker (for HGVS step; skip with --skip-vv)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
activate_conda
SCRIPTS_DIR="${ANNEAL_ROOT}/scripts"

if [ $# -lt 2 ]; then
    echo "Usage: $0 <sample_name> <output_dir> [--skip-vv]"
    exit 1
fi

SAMPLE="$1"
OUTPUT_DIR="$2"
SKIP_VV=false
if [ "${3:-}" = "--skip-vv" ]; then
    SKIP_VV=true
fi

VARIANT_DIR="${OUTPUT_DIR}/${SAMPLE}/variants"
ANNOTATED_DIR="${OUTPUT_DIR}/${SAMPLE}/annotated"
mkdir -p "${ANNOTATED_DIR}"

echo "================================================================"
echo "  Stage 3: Variant Annotation -- ${SAMPLE}"
echo "  $(date)"
echo "================================================================"

# ---- Annotate each consensus type ----
for LABEL in dcs sscs; do
    VCF="${VARIANT_DIR}/${SAMPLE}.${LABEL}.vcf"

    if [ ! -f "${VCF}" ]; then
        echo "WARNING: VCF not found, skipping: ${VCF}"
        continue
    fi

    VARIANT_COUNT=$(grep -cv "^#" "${VCF}" 2>/dev/null || echo 0)
    echo ""
    echo "--- ${LABEL^^}: ${VARIANT_COUNT} variants ---"

    # Step 1: Annotate with VEP + ANNOVAR
    echo "[$(date '+%H:%M:%S')] Annotating ${LABEL}..."
    python3 "${SCRIPTS_DIR}/annotate_variants.py" \
        --vcf "${VCF}" \
        -s "${SAMPLE}" \
        --consensus "${LABEL}" \
        -o "${ANNOTATED_DIR}" 2>&1

    ANNOTATED_TSV="${ANNOTATED_DIR}/${SAMPLE}.${LABEL}.annotated.tsv"
    if [ ! -f "${ANNOTATED_TSV}" ]; then
        echo "WARNING: Annotation failed for ${LABEL}, skipping filter/validate"
        continue
    fi

    # Step 2: Filter
    echo "[$(date '+%H:%M:%S')] Filtering ${LABEL}..."
    python3 "${SCRIPTS_DIR}/filter_variants.py" \
        -i "${ANNOTATED_TSV}" \
        -o "${ANNOTATED_DIR}" 2>&1

    CLINICAL_TSV="${ANNOTATED_DIR}/${SAMPLE}.${LABEL}.clinical.tsv"

    # Step 3: HGVS validation (optional)
    if [ "${SKIP_VV}" = false ] && [ -f "${CLINICAL_TSV}" ]; then
        CLINICAL_COUNT=$(wc -l < "${CLINICAL_TSV}")
        if [ "${CLINICAL_COUNT}" -gt 1 ]; then
            echo "[$(date '+%H:%M:%S')] Validating HGVS for ${LABEL} (${CLINICAL_COUNT} variants)..."
            python3 "${SCRIPTS_DIR}/validate_hgvs.py" \
                -i "${CLINICAL_TSV}" \
                -o "${ANNOTATED_DIR}" 2>&1 || \
                echo "WARNING: HGVS validation failed (non-fatal)"
        else
            echo "[$(date '+%H:%M:%S')] No PASS variants for ${LABEL}, skipping HGVS"
        fi
    elif [ "${SKIP_VV}" = true ]; then
        echo "[$(date '+%H:%M:%S')] Skipping HGVS validation (--skip-vv)"
    fi

    echo "[$(date '+%H:%M:%S')] ${LABEL^^} complete"
done

# ---- Summary ----
echo ""
echo "================================================================"
echo "  Stage 3 complete for ${SAMPLE}"
echo "  $(date)"
echo "================================================================"
echo "Output directory: ${ANNOTATED_DIR}/"
ls -lh "${ANNOTATED_DIR}"/ 2>/dev/null
