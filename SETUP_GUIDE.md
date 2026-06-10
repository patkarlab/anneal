# Anneal 0.1.0 -- Setup Guide

Step-by-step instructions for deploying the Anneal duplex sequencing pipeline
(consensus + variant calling, with optional annotation) on a Linux server.

Paths shown under `/goast/hemat_data/...` are examples from one deployment.
Substitute your own server's paths throughout.

---

## Prerequisites

- Linux (Ubuntu/CentOS/Rocky)
- gcc/g++ (used as the Rust linker; `sudo apt install build-essential` or `sudo yum install gcc gcc-c++`)
- At least 128 GB RAM (for bwa-mem2 alignment of hg38)
- At least 500 GB free disk (reference genomes + sequencing data + outputs)
- Internet access (for initial tool downloads)

---

## Step 1: Install Miniconda (skip if already installed)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

---

## Step 2: Create the conda environment

Use `--override-channels` if default channels require TOS acceptance:

```bash
conda create -n anneal -y --override-channels -c conda-forge -c bioconda \
    python=3.11 \
    samtools \
    bwa-mem2 \
    matplotlib \
    pandas \
    numpy

conda activate anneal

# Verify
samtools --version | head -1
python3 --version
```

### bwa-mem2 wrapper (only if your CPU lacks AVX512)

The bwa-mem2 auto-launcher selects `avx512bw` by default. On CPUs without
AVX512 (e.g. the gandalf node) this fails, and you must create a wrapper that
forces a supported SIMD binary:

```bash
mkdir -p ~/anneal/bin

# Test which SIMD level works (try avx2 first, then sse42)
/path/to/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2 version

# Create wrapper pointing to the working binary
cat > ~/anneal/bin/bwa-mem2 << 'EOF'
#!/bin/bash
exec /path/to/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2 "$@"
EOF
chmod +x ~/anneal/bin/bwa-mem2

# Verify
~/anneal/bin/bwa-mem2 version
```

The pipeline's `activate_conda()` function in config.sh already prepends
`${ANNEAL_ROOT}/bin` to PATH, so the wrapper is picked up automatically.
If avx2 also fails, try `.sse42`. If you installed bwa-mem2 from conda
(Step 2) and your CPU supports it, no wrapper is needed.

---

## Step 3: Install Rust (skip if already installed)

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source $HOME/.cargo/env

# Verify
rustc --version
cargo --version
```

---

## Step 4: Set up directory structure

```bash
# Input FASTQs and output results directories.
# Example layout from one deployment:
#   /goast/hemat_data/duplex_fastqs/dilution/   (input FASTQs)
#   /goast/hemat_data/duplex_results/           (pipeline outputs)

# On a fresh server:
mkdir -p /path/to/sequences
mkdir -p /path/to/results
```

---

## Step 5: Clone and build Anneal

```bash
cd ~

# Clone from GitHub
git clone https://github.com/<your-org>/anneal.git
cd anneal

# Build Anneal (CPU-only, no GPU needed)
bash deploy.sh

# Build the Rust variant caller
cd mpileup_variant_caller
cargo build --release
cd ..

# Verify both binaries
./target/release/anneal --help
./mpileup_variant_caller/target/release/call_variants
```

---

## Step 6: Reference genome

You need a hg38 reference FASTA with bwa-mem2 indexes (`.0123`, `.amb`,
`.ann`, `.bwt.2bit.64`, `.pac`, `.fai`). If you already have one, point
`config.sh` to it. Otherwise:

```bash
# Download and index (one-time, ~1 hour)
mkdir -p ~/references/hg38_broad
cd ~/references/hg38_broad
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta
bwa-mem2 index Homo_sapiens_assembly38.fasta
```

---

## Step 7: Copy your BED file and sequences

The BED file is bundled in the repo at `AML_MRD_DUPLEX_probes_hg38_sortd.bed`.

Copy your FASTQs into the sequences directory:

```bash
cp /path/to/fastqs/*.fastq.gz /path/to/sequences/
```

---

## Step 8: Edit config.sh

Open `~/anneal/pipeline/config.sh` and verify the paths match your setup.
Most paths auto-resolve from the directory structure. Key lines to check:

```bash
nano ~/anneal/pipeline/config.sh
```

```bash
# -- These auto-resolve, usually no edit needed --
ANNEAL="${ANNEAL_ROOT}/target/release/anneal"
VARIANT_CALLER="${ANNEAL_ROOT}/mpileup_variant_caller/target/release/call_variants"
BEDFILE="${ANNEAL_ROOT}/AML_MRD_DUPLEX_probes_hg38_sortd.bed"

# -- Edit these to match your environment --
REFERENCE="/path/to/references/hg38_broad/Homo_sapiens_assembly38.fasta"
SEQUENCES_DIR="/path/to/sequences"
RESULTS_DIR="/path/to/results"
```

---

## Step 9: Create a manifest

```bash
cd ~/anneal

# Auto-generate from your sequences directory
./target/release/anneal manifest \
    --dir /path/to/sequences/ \
    -o manifest.tsv

# Verify
cat manifest.tsv
```

---

## Step 10: Test on a single sample

```bash
cd ~/anneal
conda activate anneal

bash pipeline/run_pipeline.sh \
    SAMPLE_NAME \
    /path/to/sequences/SAMPLE_R1_001.fastq.gz \
    /path/to/sequences/SAMPLE_R2_001.fastq.gz \
    /path/to/results/
```

Check the output:

```bash
ls -lh /path/to/results/SAMPLE_NAME/consensus/
ls -lh /path/to/results/SAMPLE_NAME/variants/
```

---

## Step 11: Run full batch

```bash
cd ~/anneal

# Foreground (blocks terminal, shows progress)
conda activate anneal
bash pipeline/run_pipeline_batch.sh manifest.tsv /path/to/results/

# OR background (survives SSH disconnect)
bash pipeline/launch_pipeline.sh manifest.tsv /path/to/results/

# Monitor background run
tail -f /path/to/results/pipeline.log
```

---

## Optional: Stage 3 annotation

Stage 3 (VEP + ANNOVAR + VariantValidator) is off by default and has heavier
dependencies. Enable it once those tools are installed:

```bash
# Single sample
bash pipeline/run_pipeline.sh SAMPLE R1 R2 /path/to/results/ --annotate --skip-vv

# Batch
bash pipeline/run_pipeline_batch.sh manifest.tsv /path/to/results/ --annotate --skip-vv
```

`--skip-vv` skips the VariantValidator HGVS step (which requires its Docker
container). Drop `--skip-vv` once VariantValidator is available.

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `samtools: command not found` | Run `conda activate anneal` first |
| `cargo: command not found` | Run `source $HOME/.cargo/env` |
| bwa-mem2 avx512bw error | Create `bin/bwa-mem2` wrapper forcing `.avx2` binary (see Step 2) |
| Conda TOS error | Add `--override-channels` or run `conda tos accept ...` |
| `bwa-mem2 index` killed by OOM | Need at least 64 GB RAM for hg38 indexing |
| Anneal fails at alignment | Check that bwa-mem2 index files exist alongside the FASTA |
| Plot not generated | `pip install matplotlib pandas numpy` in the anneal env |
| Permission denied on scripts | `chmod +x ~/anneal/pipeline/*.sh` |

---

## Expected directory layout after setup

```
~/anneal/                              # ANNEAL_ROOT (pipeline code)
  pipeline/config.sh                   # all paths configured here
  target/release/anneal                # compiled binary
  mpileup_variant_caller/target/release/call_variants
  scripts/plot_family_sizes.py
  bin/bwa-mem2                         # AVX2 wrapper (only if CPU lacks AVX512)
  manifest.tsv                         # sample sheet
  AML_MRD_DUPLEX_probes_hg38_sortd.bed # target panel (bundled)

/path/to/sequences/                    # input FASTQs
/path/to/results/                      # pipeline outputs
  SAMPLE_NAME/
    consensus/
    variants/
    annotated/                         # only if --annotate used

/path/to/references/
  hg38_broad/Homo_sapiens_assembly38.fasta  # hg38 + bwa-mem2 indexes
```
