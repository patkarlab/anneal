# CLAUDE.md -- Anneal Architecture & Design Notes

Guidance for working in this repository. Anneal is a duplex consensus pipeline
for UMI-based error suppression and ultra-sensitive variant detection.

## Project Overview

Anneal takes UMI-tagged paired-end FASTQs and produces duplex consensus BAMs,
then calls SNVs and indels from those high-confidence consensus reads. An
optional annotation stage adds functional annotation and HGVS validation.

The core engine is written in Rust (CPU, with optional CUDA acceleration). The
variant caller is a separate Rust crate. Orchestration is a set of bash stage
scripts driven by a single shared `config.sh`.

## Pipeline Stages

| Stage | Name | Tool | Output |
|-------|------|------|--------|
| 1 | Consensus generation | `anneal` (Rust) | `{sample}.sscs.sc.sorted.bam`, `{sample}.dcs.sc.sorted.bam`, stats, family sizes |
| 2 | Variant calling | samtools mpileup + `call_variants` (Rust) | `{sample}.sscs.vcf`, `{sample}.dcs.vcf` |
| 3 | Annotation (optional) | VEP + ANNOVAR + VariantValidator | `{sample}.{dcs,sscs}.annotated/filtered/clinical.tsv` |

Stage 1 internals: barcode extraction -> BWA-MEM2 alignment -> family grouping
-> SSCS (single-strand consensus) -> singleton correction -> DCS (duplex
consensus). VCF FORMAT is `GT:ALT:TOT:FRAC` (alt reads, total depth, VAF).

Default run is stages 1,2. Stage 3 is opt-in via `--annotate` because it
depends on VEP, ANNOVAR, and the VariantValidator Docker container.

## File Structure

```
src/                          # Rust consensus engine
  barcode/                    # UMI extraction
  grouping/                   # family grouping
  consensus/                  # SSCS, DCS, pipeline, config
  singleton/                  # singleton correction
  cuda/                       # optional GPU kernels (PTX + .cu)
  manifest.rs                 # manifest subcommand
  main.rs
mpileup_variant_caller/       # separate Rust crate (call_variants)
pipeline/
  config.sh                   # shared config -- edit paths here
  stage1_consensus.sh
  stage2_variant_calling.sh
  stage3_annotate.sh          # optional
  run_pipeline.sh             # single sample
  run_pipeline_batch.sh       # manifest-driven batch
  launch_pipeline.sh          # background launcher (nohup)
scripts/                      # python: plotting + annotation helpers
  plot_family_sizes.py
  annotate_variants.py
  filter_variants.py
  validate_hgvs.py
deploy.sh                     # build script (cpu | gpu)
AML_MRD_DUPLEX_probes_hg38_sortd.bed
```

## Key Design Decisions

### Why only 2 BAM outputs?
SSCS (with singleton correction) and DCS are the two consensus types that
matter clinically. SSCS maximizes sensitivity; DCS maximizes specificity
(both strands must agree). Intermediate BAMs are cleaned up automatically.

### Consensus cutoff 0.6
The base-call agreement fraction within a family. 0.6 balances retaining
real low-frequency signal against over-calling noise; tunable via `--cutoff`.

### Rust variant caller
Replaces an older Perl mpileup parser, roughly an order of magnitude faster,
and emits clean VCF with explicit alt-read/depth/VAF fields for MRD work.

### UMI-aware callers are unnecessary downstream
Because consensus is already built per UMI family upstream, the Stage 2 caller
operates on consensus reads and does not need to be UMI-aware itself.

## Configuration (config.sh)

`ANNEAL_ROOT` auto-resolves from the script location, and most binary/reference
paths derive from it. The variables you actually edit per server:

- `REFERENCE` -- hg38 FASTA with bwa-mem2 indexes
- `BEDFILE` -- target panel (bundled)
- `SEQUENCES_DIR`, `RESULTS_DIR` -- input/output
- `USE_GPU` -- false (CPU) by default

`activate_conda()` sources conda, activates the `anneal` env, and prepends
`${ANNEAL_ROOT}/bin` to PATH (for the optional bwa-mem2 SIMD wrapper).

## Build & Run

```bash
bash deploy.sh                                   # build anneal (cpu); use `gpu` for CUDA
cd mpileup_variant_caller && cargo build --release && cd ..
bash pipeline/run_pipeline.sh SAMPLE R1 R2 out/  # stages 1,2
bash pipeline/run_pipeline.sh SAMPLE R1 R2 out/ --annotate --skip-vv   # + stage 3
```

## Changelog

### 0.1.0 -- first public release
- Removed FLT3-ITD detection (former Stage 3: getITD + FiLT3R + concordance).
  FLT3 exon regions remain in the panel BED for standard SNV/indel calling.
- Annotation renumbered from Stage 4 to Stage 3.
- Server-agnostic paths; no PBS dependency.

## References

- Wang TT, Abelson S, et al. (2019) Nucleic Acids Research, 47(15), e87
- Kennedy SR, et al. (2014) Nature Protocols, 9(11), 2586-2606
- ConsensusCruncher: https://github.com/pughlab/ConsensusCruncher
