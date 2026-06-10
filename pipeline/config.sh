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
SEQUENCES_DIR="/goast/hemat_data/duplex_fastqs/dilution"
RESULTS_DIR="/goast/hemat_data/duplex_results"

# ---- Reference genome (U2AF1-fixed hg38 from targeted-seq-pipeline) ----
REFERENCE="/home/hemat/targeted-seq-pipeline/references/hg38_broad/Homo_sapiens_assembly38.masked.fasta"

# ---- Target panel BED file ----
BEDFILE="${ANNEAL_ROOT}/AML_MRD_DUPLEX_probes_hg38_sortd.bed"

# ---- Binaries ----
ANNEAL="${ANNEAL_ROOT}/target/release/anneal"
VARIANT_CALLER="${ANNEAL_ROOT}/mpileup_variant_caller/target/release/call_variants"

# ---- Family size plot script ----
PLOT_SCRIPT="${ANNEAL_ROOT}/scripts/plot_family_sizes.py"

# ---- Conda environment ----
CONDA_ENV="anneal"

# ---- Anneal consensus parameters ----
ALIGNER="bwa-mem2"
BPATTERN="NNNSS"
CUTOFF=0.6
SINGLETON_CORRECTION=true

# ---- Variant calling parameters ----
MIN_BASE_QUAL="5"    # ASCII character for Phred+33 quality threshold
MAX_DEPTH=100000

# ---- GPU ----
USE_GPU=false

# ---- Activate conda environment ----
activate_conda() {
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
    conda activate "${CONDA_ENV}"
    export PATH="${ANNEAL_ROOT}/bin:$PATH"
}
