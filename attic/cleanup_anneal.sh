#!/bin/bash
#
# cleanup_anneal.sh
#
# Cleans the anneal working directory ahead of pipeline lock.
#
#   Decisions applied:
#     - Rust mpileup variant caller: REMOVED (superseded by Pisces)
#     - CPU-only execution path: REMOVED (GPU-only pipeline)
#     - Rejected callers (Mutect2, SiNVICT, 3-caller matcher): moved to attic/
#
#   Nothing under results_*/, *_fastq*/, .git/, src/, or target/ is touched.
#
# Usage:
#     bash cleanup_anneal.sh              # dry run, prints actions only
#     APPLY=1 bash cleanup_anneal.sh      # actually performs them
#
set -uo pipefail

ROOT="${ANNEAL_ROOT:-$HOME/pipelines/anneal}"
APPLY="${APPLY:-0}"

cd "$ROOT" || { echo "FATAL: cannot cd to $ROOT"; exit 1; }

if [ ! -d .git ]; then
    echo "FATAL: $ROOT is not a git repository. Refusing to run."
    exit 1
fi

echo "Working directory : $ROOT"
if [ "$APPLY" = "1" ]; then
    echo "Mode              : APPLY (changes will be made)"
else
    echo "Mode              : DRY RUN (no changes; re-run with APPLY=1)"
fi
echo

N_DEL=0
N_MOV=0

do_rm() {
    # $1 = path
    [ -e "$1" ] || return 0
    echo "  DELETE  $1"
    N_DEL=$((N_DEL + 1))
    [ "$APPLY" = "1" ] && rm -rf -- "$1"
    return 0
}

do_mv() {
    # $1 = source, $2 = destination directory
    [ -e "$1" ] || return 0
    echo "  MOVE    $1  ->  $2/"
    N_MOV=$((N_MOV + 1))
    if [ "$APPLY" = "1" ]; then
        mkdir -p "$2"
        mv -- "$1" "$2"/
    fi
    return 0
}

# ---------------------------------------------------------------------------
# 1. PBS job logs and stray run logs
# ---------------------------------------------------------------------------
echo "[1] PBS job logs and stray logs"
while IFS= read -r f; do
    do_rm "$f"
done < <(find . -maxdepth 1 -type f \( \
            -name '*.o[0-9]*' -o \
            -name '*.e[0-9]*' -o \
            -name 'dilution_cpu.out' -o \
            -name 'dilution_gpu.out' -o \
            -name 'dilution_gpu.err' -o \
            -name 'script.log' -o \
            -name 'vardict.log' \) | sort)
echo

# ---------------------------------------------------------------------------
# 2. Timestamped backup files
# ---------------------------------------------------------------------------
echo "[2] Backup files (.bak / .bak_YYYYMMDD / .bak.YYYYMMDD)"
while IFS= read -r f; do
    do_rm "$f"
done < <(find . -path ./.git -prune -o -path './results_*' -prune -o \
            -type f \( -name '*.bak' -o -name '*.bak_*' -o -name '*.bak.*' \) \
            -print | sort)
echo

# ---------------------------------------------------------------------------
# 3. Zero-byte debris and Python caches
# ---------------------------------------------------------------------------
echo "[3] Zero-byte files and __pycache__"
while IFS= read -r f; do
    do_rm "$f"
done < <(find . -maxdepth 1 -type f -size 0 | sort)

while IFS= read -r d; do
    do_rm "$d"
done < <(find . -path ./.git -prune -o -type d -name '__pycache__' -print | sort)
echo

# ---------------------------------------------------------------------------
# 4. Rust mpileup variant caller  (DECISION: remove, superseded by Pisces)
# ---------------------------------------------------------------------------
echo "[4] Rust mpileup variant caller (superseded by Pisces)"
do_rm "mpileup_variant_caller"
echo

# ---------------------------------------------------------------------------
# 5. CPU-only execution path  (DECISION: GPU-only pipeline)
# ---------------------------------------------------------------------------
echo "[5] CPU-only execution artifacts"
do_rm "batch_cpu.pbs"
do_rm "target_cpu"
echo "  NOTE    removing the --no-gpu CLI flag is a source edit in src/;"
echo "          not handled here. Requires rebuild + version bump."
echo

# ---------------------------------------------------------------------------
# 6. Rejected / superseded callers -> attic (preserved, gitignored)
# ---------------------------------------------------------------------------
echo "[6] Rejected and superseded caller artifacts -> attic/"
do_mv "mutect2_G_sscs.pbs"   "attic/rejected_callers"
do_mv "sinvict_G_sscs.pbs"   "attic/rejected_callers"
do_mv "consensus_callers.py" "attic/rejected_callers"
do_mv "run_consensus_G.sh"   "attic/rejected_callers"
echo

# ---------------------------------------------------------------------------
# 7. One-off deployment tests -> attic
# ---------------------------------------------------------------------------
echo "[7] Deployment test scripts -> attic/"
do_mv "run_test.pbs"     "attic/deployment_tests"
do_mv "run_test_gpu.pbs" "attic/deployment_tests"
do_mv "run_pb_test.pbs"  "attic/deployment_tests"
do_mv "bwa_index.pbs"    "attic/deployment_tests"
echo

# ---------------------------------------------------------------------------
# 8. Stray analysis scripts -> scripts/
# ---------------------------------------------------------------------------
echo "[8] Stray Python scripts -> scripts/"
do_mv "build_background_pileup.py"    "scripts/error_model"
do_mv "score_calls_vs_background.py"  "scripts/error_model"
do_mv "npm1_cebpa_dilution.py"        "scripts/validation"
echo

# ---------------------------------------------------------------------------
# 9. PBS submission templates -> jobs/
# ---------------------------------------------------------------------------
echo "[9] PBS submission templates -> jobs/"
do_mv "batch_gpu.pbs"               "jobs"
do_mv "batch_gpu_bnc.pbs"           "jobs"
do_mv "annotate_dilution_array.pbs" "jobs"
do_mv "pisces_G_sscs.pbs"           "jobs"
do_mv "run_bnc_calibration.sh"      "jobs"
do_mv "run_dilution_rungs.sh"       "jobs"
echo

# ---------------------------------------------------------------------------
# 10. Sample / marker site lists -> assets/
# ---------------------------------------------------------------------------
echo "[10] Manifests and marker lists -> assets/"
do_mv "npm1_marker.sites.tsv"  "assets/marker_lists"
do_mv "dilution_manifest.tsv"  "assets/manifests"
do_mv "annot_manifest.tsv"     "assets/manifests"
do_mv "samples.tsv"            "assets/manifests"
echo

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "-----------------------------------------------------------"
printf "Deletions : %d\n" "$N_DEL"
printf "Moves     : %d\n" "$N_MOV"
echo "-----------------------------------------------------------"

if [ "$APPLY" != "1" ]; then
    echo
    echo "Dry run complete. Nothing was changed."
    echo "To apply:   APPLY=1 bash cleanup_anneal.sh"
else
    echo
    echo "Applied. Review before staging anything:"
    echo "    git status --short"
    echo
    echo "Reminder: stage files individually. Never 'git add -A' or 'git add .'"
fi
