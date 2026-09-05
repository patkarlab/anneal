#!/bin/bash
# run_consensus_G.sh — build the blind three-caller consensus + ANNOVAR annotation
# for the G rung. Runs ANNOVAR on the union of caller variants, then the matcher.
#
# BLIND: no diagnosis-variant list anywhere. Caller-count is computed with no
# knowledge of expected variants. Recovery scoring is a separate downstream step.

set -euo pipefail
eval "$(/home/patkarlab-clinical/miniconda3/bin/conda shell.bash hook)"
conda activate anneal

# ---- Paths ----
ANNOVAR=/home/patkarlab-clinical/programs/annovar
HUMANDB=/home/patkarlab-clinical/references/humandb
SAMPLE=DIL-A-G-Duplex
DIR=results_dilution_gpu/${SAMPLE}

VARDICT=$DIR/${SAMPLE}.sscs.vardict.vcf
LOFREQ=$DIR/${SAMPLE}.sscs.lofreq.baq_off.vcf
PISCES=$DIR/pisces_out/${SAMPLE}.sscs.panel.pisces.genome.vcf

OUTPREFIX=$DIR/${SAMPLE}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "== Consensus build: ${SAMPLE} =="
for f in "$VARDICT" "$LOFREQ" "$PISCES"; do
  [ -f "$f" ] || { echo "FATAL: missing caller VCF: $f" >&2; exit 1; }
done

# Step 1: matcher pass 1 -> produces the union VCF for ANNOVAR (no annotation yet).
python3 "$SCRIPT_DIR/consensus_callers.py" \
  --vardict "$VARDICT" --lofreq "$LOFREQ" --pisces "$PISCES" \
  --out-prefix "$OUTPREFIX" --indel-window 5

# Step 2: ANNOVAR on the union variant set.
ANNVCF=${OUTPREFIX}.for_annovar.vcf
perl $ANNOVAR/table_annovar.pl "$ANNVCF" "$HUMANDB" -buildver hg38 \
  -out ${OUTPREFIX}.annovar -remove -protocol refGene -operation g \
  -nastring . -vcfinput
MULTIANNO=${OUTPREFIX}.annovar.hg38_multianno.txt

# Step 3: matcher pass 2 -> merge annotation onto consensus.
python3 "$SCRIPT_DIR/consensus_callers.py" \
  --vardict "$VARDICT" --lofreq "$LOFREQ" --pisces "$PISCES" \
  --out-prefix "$OUTPREFIX" --indel-window 5 \
  --multianno "$MULTIANNO"

echo ""
echo "== Outputs =="
echo "  consensus:           ${OUTPREFIX}.consensus.tsv"
echo "  consensus+annotation: ${OUTPREFIX}.consensus.annotated.tsv"
echo ""
echo "== caller_count distribution =="
awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) print "  "k" caller(s): "c[k]}' ${OUTPREFIX}.consensus.annotated.tsv | sort -rn -k2
