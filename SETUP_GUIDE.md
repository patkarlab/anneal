# anneal: setup guide

Installing the pipeline on a PBS cluster with NVIDIA GPU nodes. Version
0.3.x. Running it once installed is in `SOP.md`.

The reference installation is `~/pipelines/anneal` on the BioInfinix
cluster (`patkarlab-clinical`, login node `ln1`): PBS Pro, A40 GPU nodes,
Apptainer as an environment module, no Docker daemon and no internet on the
compute nodes. Paths below are that installation's; every one of them is set
in `pipeline/config.sh` and can be changed there.

## Prerequisites

- A CUDA GPU node for alignment (Parabricks 4.3.1 runs on A40s; 16 CPUs and
  100 GB RAM per sample are what the jobs request).
- CPU nodes for everything else.
- Disk: about 150–200 GB of transient space per sample in flight, plus the
  references (~20 GB), the VEP cache (~20 GB) and ANNOVAR databases.
- gcc (Rust linker), Miniconda, git.
- Internet access on the login node (conda, VEP cache, VariantValidator API).

## 1. Conda environments

The exact environments of the reference installation are committed as
`envs/anneal.yml` and `envs/vep.yml` (conda exports, no build strings).
Create from those:

```bash
conda env create -n anneal -f envs/anneal.yml
conda env create -n vep -f envs/vep.yml
```

If a package in the export is unavailable on the new machine, the
environments can be built from names instead; the pinned versions in the
yml files are the ones to match:

```bash
conda create -n anneal -y --override-channels -c conda-forge -c bioconda \
    python=3.11 samtools pysam numpy pandas matplotlib requests dotnet-runtime=8
conda create -n vep -y --override-channels -c conda-forge -c bioconda ensembl-vep
```

The .NET 8 runtime in the `anneal` environment (`lib/dotnet` under the env)
is what runs Pisces; `config.sh` points `DOTNET_DIR` at it and stage 2 sets
`DOTNET_ROLL_FORWARD=Major` because Pisces was built against .NET 2.

VEP has its own environment so that its Perl does not collide with the
analysis one; `run_vep()` in `scripts/annotate_variants.py` strips the
analysis environment's Perl from PATH and clears `PERL5LIB` when it invokes
VEP. `config.sh` activates `anneal` through `activate_conda()`, which relaxes
`set -u` across activation (the .NET activation hook references an unset
variable) and prepends `bin/` to PATH for the Parabricks shim.

## 2. Build the consensus binary

```bash
git clone git@github.com:patkarlab/anneal.git ~/pipelines/anneal      # or https://github.com/patkarlab/anneal.git
cd ~/pipelines/anneal
source ~/.cargo/env 2>/dev/null || (curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && source ~/.cargo/env)
bash build.sh            # CPU build into target_cpu/release/anneal, runs the unit tests
```

The CPU build is the validated one; `config.sh` points `ANNEAL` at it. There
is a `gpu` cargo feature for a CUDA consensus path, but it is not used:
consensus runs on the CPU in every validated configuration, and stage 1
passes `--no-gpu`.

## 3. Parabricks through Apptainer

The binary calls Parabricks with a hard-coded `docker run`. On a cluster
without a Docker daemon, a shim named `docker` translates that into
`apptainer exec`.

```bash
# once: the Parabricks image as a SIF (from the Docker archive; no NGC login needed if you have the tar)
cd ~/pipelines
apptainer build parabricks_4.3.1.sif docker-archive://parabricks_4.3.1.tar

# the shim, on the PATH that activate_conda() prepends
cp ~/pipelines/anneal/pipeline/docker-apptainer-shim.sh ~/pipelines/anneal/bin/docker
chmod +x ~/pipelines/anneal/bin/docker
```

Two things every GPU job must do, and the reference jobs in `jobs/` do:
`module load apptainer/1.5.1` (the module is not on the default PATH of
compute nodes, so interactive tests on the login node pass where batch jobs
fail) and use `config.sh` defaults `ALIGNER=parabricks`,
`PARABRICKS_FLAGS=--parabricks-docker`, `ALIGNER_ARGS="--num-gpus 1"` (the
container sees one GPU; without the flag Parabricks asks for both).

Check from a GPU node:

```bash
module load apptainer/1.5.1
docker run --rm --gpus all nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1 pbrun version   # expect: pbrun: 4.3.1-1
```

Fallback without a GPU: `ALIGNER=bwa-mem2` with a bwa-mem2 index of the
masked reference and a CPU queue; alignment then takes an hour or more per
sample and the run is not the validated configuration.

## 4. References and annotation tools

| Purpose | Location (config variable) |
|---------|----------------------------|
| Alignment, background model | `~/references/hg38_broad_bwa/Homo_sapiens_assembly38.masked.fasta` (`REFERENCE`) |
| Indel scan, annotation | `~/references/hg38_broad/Homo_sapiens_assembly38.fasta` (`REFERENCE_UNMASKED`) |
| Pisces | `~/references/pisces_hg38_panel/` (`PISCES_GENOME`) |
| Pisces binaries | `~/programs/pisces/Pisces_5.2.10.49/Pisces.dll` (`PISCES_DLL`) |
| VEP cache | `~/references/vep_cache/` (`VEP_CACHE`) |
| ANNOVAR | `~/programs/annovar/` (`ANNOVAR_DIR`), databases in `~/references/humandb/` (`ANNOVAR_DB`) |
| getITD | `~/tools/getitd/` (`GETITD_DIR`) |
| Panel | `AML_MRD_DUPLEX_probes_hg38_sortd.bed`, in the repository (`BEDFILE`) |

**hg38.** The Broad hg38 (`Homo_sapiens_assembly38.fasta`) and the same
sequence with the U2AF1 duplicate region masked. Parabricks needs the classic
BWA index of the masked file; the unmasked one needs only `samtools faidx`:

```bash
mkdir -p ~/references/hg38_broad ~/references/hg38_broad_bwa
cd ~/references/hg38_broad
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta
# masked copy: obtain Homo_sapiens_assembly38.masked.fasta from the reference installation
cd ~/references/hg38_broad_bwa && samtools faidx Homo_sapiens_assembly38.masked.fasta && bwa index Homo_sapiens_assembly38.masked.fasta
```

**Pisces.** Release 5.2.10.49 from https://github.com/Illumina/Pisces/releases,
unpacked under `~/programs/pisces/`. Then the panel-restricted reference,
which the script builds from the masked FASTA and the BED (one FASTA of the
panel chromosomes, index, dictionary, and `GenomeSize.xml` from Pisces' own
`CreateGenomeSizeFile`):

```bash
bash scripts/make_pisces_reference.sh          # writes to PISCES_GENOME
```

**VEP cache**, GRCh38, for offline use, from the `vep` environment:

```bash
conda activate vep
vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/references/vep_cache --CONVERT
```

**ANNOVAR.** Registered download from https://annovar.openbioinformatics.org,
unpacked to `~/programs/annovar/`. `scripts/annotate_variants.py` uses the
hg38 tables `refGene`, `cosmic103`, `gnomad211_exome`, `clinvar_20250721`,
`avsnp151`, and skips any that are missing. Three come from ANNOVAR's server:

```bash
cd ~/programs/annovar
for db in refGene gnomad211_exome avsnp151; do
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar $db ~/references/humandb/
done
```

`cosmic103` and `clinvar_20250721` are user-prepared tables (ANNOVAR does
not distribute COSMIC; ClinVar dates on its server lag): build them with
ANNOVAR's `prepare_annovar_user.pl` from the COSMIC and ClinVar downloads, or
copy `hg38_cosmic103.txt` and `hg38_clinvar_20250721.txt` (with their
`.idx` files) from the reference installation.

**getITD** and its Python dependencies (per its README: numpy, pandas,
biopython), into the `anneal` environment:

```bash
git clone https://github.com/tjblaette/getitd.git ~/tools/getitd
conda install -n anneal -y -c conda-forge biopython
```

## 5. Configuration

`pipeline/config.sh` holds every path, tool and parameter, each as
`VAR="${VAR:-default}"` so it can be overridden from the environment. After
installation, check the block that names the machine:

```bash
grep -n "^REFERENCE=\|^REFERENCE_UNMASKED=\|^ANNEAL=\|^VEP_CACHE=\|^ANNOVAR_DIR=\|^ANNOVAR_DB=\|^GETITD_DIR=\|^PISCES_DLL=\|^DOTNET_DIR=\|^PISCES_GENOME=" pipeline/config.sh
```

Do not change the analysis parameters (cutoff, quality, singleton
correction, alpha, blocklist threshold, estimator settings): they are the
validated configuration, and `SOP.md` section 3 lists them.

## 6. Background model assets

Stage 5 needs the per-track background matrices and indel blocklists, built
once from eight biological negative controls:

```
results_bnc/beta_matrix_DCS.txt   results_bnc/beta_matrix_SSCS.txt      (+ .report.tsv)
results_bnc/indel_blocklist.DCS.patched.tsv   results_bnc/indel_blocklist.SSCS.patched.tsv
```

They are data, not code, and are not in the repository. On a new
installation they are copied from the reference installation, or rebuilt:
run the eight control FASTQs through stage 1 (`jobs/anneal_e2e.pbs` with
`--stages 1`), place the consensus BAMs under `results_bnc_patched/<control>/`,
then `qsub jobs/rebuild_background.pbs` (matrices) and
`scripts/error_model/build_indel_blocklist_v2.py` (blocklists). A rebuild
must reproduce the validation in `CHANGELOG.md` before it is used.

## 7. Verify the installation

Run one sample end to end and compare with the expectations:

```bash
cd ~/pipelines/anneal
mkdir -p logs
qsub -v FULL=<sample>_S<n>,OUT=/scratch/<user>/verify,FQ=/scratch/<user>/<fastq_dir> jobs/anneal_e2e.pbs
```

On completion the log ends with `run_pipeline exit: 0`, `<sample>/scored/`
holds four tables, and `consensus/<sample>.stats.txt` shows a singleton rate
of 55–60% and DCS recovery of 15–22% on a typical library. A sample run
twice from FASTQ gives identical consensus statistics and identical calls;
the pipeline is deterministic on identical input and code.

## Other schedulers, or a single server

The stage scripts and `run_pipeline.sh` do not depend on the scheduler; only
the files in `jobs/` do. On Slurm, replace the `#PBS` headers with `#SBATCH`
equivalents (`--array` for the cohort job, `--gres=gpu:1` for stage 1) and
replace `PBS_ARRAY_INDEX` with `SLURM_ARRAY_TASK_ID`; the quota gate uses
`lfs quota`, which is Lustre-specific, and becomes `df` elsewhere. On a
single server without a scheduler, `pipeline/run_pipeline_batch.sh` runs a
manifest of samples one after another, which matches the throughput the
cluster achieves anyway. The hard dependency is Parabricks: an NVIDIA GPU
with at least 16 GB of memory and either Apptainer or a Docker daemon; with
Docker available, the shim is not needed and `PARABRICKS_FLAGS` can be
empty. Without a GPU the bwa-mem2 fallback in section 3 works but is not the
validated configuration.

## Troubleshooting installation

| Symptom | Cause | Fix |
|---------|-------|-----|
| `apptainer: not found` inside a job | module not loaded | `module load apptainer/1.5.1` in the job |
| `Aligner 'parabricks' not found at 'pbrun'` | shim not on PATH or flag missing | `bin/docker` present and executable; `PARABRICKS_FLAGS=--parabricks-docker` |
| Parabricks: "requires 2 GPUs" | container sees one | `ALIGNER_ARGS="--num-gpus 1"` |
| Pisces: framework version error | roll-forward | stage 2 sets `DOTNET_ROLL_FORWARD=Major`; check `DOTNET_DIR` |
| Pisces killed, out of memory | full reference | use the panel reference in `PISCES_GENOME` |
| VEP: `Compilation failed in require ... base.pm` | Perl from another environment | VEP is invoked only through `run_vep()`, which strips the analysis env's Perl from PATH and clears `PERL5LIB` |
| `cargo` not found | Rust not on PATH | `source ~/.cargo/env` |
| Linker errors building the binary | no gcc | `build.sh` sets `RUSTFLAGS="-C linker=gcc"`; install `build-essential` or `gcc gcc-c++` |
