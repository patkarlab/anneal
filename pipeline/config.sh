#!/bin/bash
# =============================================================================
# config.sh -- Shared configuration for Anneal 0.1.0 pipeline
#
# Directory layout:
#
#   /home/hemat/anneal/              <-- ANNEAL_ROOT (auto-resolved)
#     pipeline/                      <-- shell scripts + this config
#     src/                           <-- Rust source
#     target/release/anneal          <-- compiled binary
#     mpileup_variant_caller/        <-- Rust variant caller
#     scripts/                       <-- plotting + annotation scripts
#     bin/                           <-- bwa-mem2 wrapper
#
#   /goast/hemat_data/
#     duplex_fastqs/dilution/        <-- input FASTQs
#     duplex_results/                <-- pipeline outputs
#
# Edit the paths below if your layout differs.
# =============================================================================

# ---- Auto-resolve directories ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANNEAL_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"   # /home/hemat/anneal/

# ---- Workspace directories ----
SEQUENCES_DIR="/home/patkarlab-clinical/pipelines/anneal/sample_fastq"
RESULTS_DIR="/home/patkarlab-clinical/pipelines/anneal/results"

# ---- Reference genome (U2AF1-fixed hg38 from targeted-seq-pipeline) ----
REFERENCE="${REFERENCE:-/home/patkarlab-clinical/references/hg38_broad/Homo_sapiens_assembly38.masked.fasta}"

# ---- Target panel BED file ----
BEDFILE="${ANNEAL_ROOT}/AML_MRD_DUPLEX_probes_hg38_sortd.bed"

# ---- Annotation: VEP cache + ANNOVAR (install scripts vs database dir) ----
# ANNOVAR scripts live under programs/annovar; the CURRENT databases live under
# references/humandb (NOT programs/annovar/humandb, which holds an old clinvar).
VEP_CACHE="${VEP_CACHE:-/home/patkarlab-clinical/references/vep_cache}"
ANNOVAR_DIR="${ANNOVAR_DIR:-/home/patkarlab-clinical/programs/annovar}"
ANNOVAR_DB="${ANNOVAR_DB:-/home/patkarlab-clinical/references/humandb}"

# ---- Binaries ----
ANNEAL="${ANNEAL:-${ANNEAL_ROOT}/target/release/anneal}"

# ---- Family size plot script ----
PLOT_SCRIPT="${ANNEAL_ROOT}/scripts/plot_family_sizes.py"

# ---- Conda environment ----
CONDA_ENV="anneal"

# ---- Anneal consensus parameters ----
ALIGNER="${ALIGNER:-bwa-mem2}"
BPATTERN="NNNSS"
CUTOFF=0.6
SINGLETON_CORRECTION=true

# ---- Variant calling parameters ----
MIN_BASE_QUAL="5"    # ASCII character for Phred+33 quality threshold
MAX_DEPTH=100000

# ---- GPU ----
USE_GPU="${USE_GPU:-false}"

# ---- Activate conda environment ----
activate_conda() {
    # conda's own activation hooks are not written to survive `set -u`; the
    # dotnet hook installed for Pisces dereferences unset variables. Relax
    # nounset across activation only, then restore whatever the caller had.
    local _u_was_set=0
    case "$-" in *u*) _u_was_set=1; set +u ;; esac

    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
    conda activate "${CONDA_ENV}"

    [ "${_u_was_set}" -eq 1 ] && set -u
    export PATH="${ANNEAL_ROOT}/bin:$PATH"
}
PARABRICKS_FLAGS="${PARABRICKS_FLAGS:-}"

# ---- Stage 4: FLT3-ITD ----
# getITD checkout: git clone https://github.com/tjblaette/getitd.git
# Requires biopython < 1.84 (Bio.pairwise2 was removed in 1.84).
GETITD_DIR="${GETITD_DIR:-${HOME}/tools/getitd}"

# Substring of BED column 4 selecting the FLT3 probes to extract. Column 4
# carries legacy hg19 labels; only the gene/exon name is matched, coordinates
# are always read from columns 2-3.
FLT3_PROBE_PATTERN="${FLT3_PROBE_PATTERN:-FLT3_Ex_1[45]}"

# Padding either side of the probe span, to catch reads overhanging the target.
FLT3_FLANK="${FLT3_FLANK:-300}"

# Below this many consensus reads over FLT3 the locus is not evaluable.
FLT3_MIN_READS="${FLT3_MIN_READS:-500}"

FLT3_THREADS="${FLT3_THREADS:-8}"

# ---- Stage 2: Pisces ----
# Pisces 5.2.10.49. It is a .NET 2.0 application; the available runtime is
# .NET 8, so stage 2 exports DOTNET_ROLL_FORWARD=Major. Compute nodes carry no
# dotnet on PATH, hence the conda copy.
PISCES_DLL="${PISCES_DLL:-${HOME}/programs/pisces/Pisces_5.2.10.49/Pisces.dll}"
DOTNET_DIR="${DOTNET_DIR:-${HOME}/miniconda3/envs/anneal/lib/dotnet}"

# Panel-restricted masked genome plus GenomeSize.xml. Pisces exhausts memory on
# the full 3,366-contig hg38 reference.
PISCES_GENOME="${PISCES_GENOME:-${HOME}/references/pisces_hg38_panel}"

# --minbq 30  base-call quality floor. Q60 consensus bases clear it; this is
#             the parameter that sets the low-VAF detection limit.
# --minvf 1e-4  minimum variant frequency.
# --minmq 0   consensus reads are pre-validated, no MAPQ refiltering.
# -c 1        minimum coverage; depth is high everywhere on the panel.
PISCES_MINBQ="${PISCES_MINBQ:-30}"
PISCES_MINMQ="${PISCES_MINMQ:-0}"
PISCES_MINVF="${PISCES_MINVF:-0.0001}"
PISCES_MINCOV="${PISCES_MINCOV:-1}"

# ---- Background error model (built from biological negative controls) ----
# See scripts/background_model/ for how these are regenerated.
BETA_MATRIX_DCS="${BETA_MATRIX_DCS:-${ANNEAL_ROOT}/results_bnc/beta_matrix_DCS.txt}"
BETA_MATRIX_SSCS="${BETA_MATRIX_SSCS:-${ANNEAL_ROOT}/results_bnc/beta_matrix_SSCS.txt}"
ARTIFACT_MASK="${ARTIFACT_MASK:-${ANNEAL_ROOT}/results_bnc/artifact_mask.combined.bed}"
# Track-matched indel blocklists. SSCS artifacts that are strand-specific do
# not survive into DCS (C>A falls 29x between tracks), so applying an SSCS
# blocklist to DCS over-filters at TP53, RUNX1, GATA2 and STAG2 hot spots.
INDEL_BLOCKLIST_SSCS="${INDEL_BLOCKLIST_SSCS:-${ANNEAL_ROOT}/results_bnc/indel_blocklist.SSCS.patched.tsv}"
INDEL_BLOCKLIST_DCS="${INDEL_BLOCKLIST_DCS:-${ANNEAL_ROOT}/results_bnc/indel_blocklist.DCS.patched.tsv}"
# Fallback for anything still reading the single-variable form.
INDEL_BLOCKLIST="${INDEL_BLOCKLIST:-${ANNEAL_ROOT}/results_bnc/indel_blocklist.tsv}"

# ---- Extra arguments passed straight to the aligner ----
# Parabricks defaults to every GPU the driver reports. Apptainer's --nv honours
# the PBS ngpus allocation, so CUDA sees fewer devices than Parabricks asks for
# and pbrun aborts with "Number of GPUs requested (2) is more than number of
# GPUs (1) in the system". Docker's --gpus all exposed both, which is why this
# did not surface before the move to apptainer.
ALIGNER_ARGS="${ALIGNER_ARGS:---num-gpus 1}"

# ---- Indel scanning (stage 2, alongside Pisces) ----
# Anchor and deleted bases are read from real sequence, so this must be the
# UNMASKED reference, not the masked copy used for alignment.
REFERENCE_UNMASKED="${REFERENCE_UNMASKED:-/home/patkarlab-clinical/references/hg38_broad/Homo_sapiens_assembly38.fasta}"
INDEL_MIN_ALT="${INDEL_MIN_ALT:-2}"
