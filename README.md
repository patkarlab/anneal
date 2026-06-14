# Anneal 0.1.0

GPU-accelerated duplex consensus pipeline for UMI-based error suppression in NGS data.
Designed for ultra-sensitive variant detection in liquid biopsy and MRD monitoring.

## Overview

Anneal builds duplex consensus sequences from UMI-tagged reads to suppress sequencing
error, then calls SNVs and indels from the resulting high-confidence consensus BAMs.
An optional annotation stage adds VEP + ANNOVAR functional annotation and HGVS validation.

- **Duplex consensus engine** (Rust): barcode extraction, family grouping, SSCS,
  singleton correction, and DCS, with optional CUDA acceleration.
- **Rust mpileup variant caller**: fast SNV + indel detection, VCF output with
  GT:ALT:TOT:FRAC.
- **Optional annotation**: VEP + ANNOVAR + VariantValidator integration (Stage 3).
- **Server-agnostic deployment**: auto-resolved paths, no PBS dependency.
- **Automated cleanup**: intermediate BAMs and mpileup files removed automatically.

## Quick Start

```bash
# Build
bash deploy.sh
cd mpileup_variant_caller && cargo build --release && cd ..

# Edit paths
vi pipeline/config.sh

# Single sample -- consensus + variant calling (default)
bash pipeline/run_pipeline.sh SAMPLE R1.fastq.gz R2.fastq.gz output/

# Single sample -- include annotation (Stage 3)
bash pipeline/run_pipeline.sh SAMPLE R1.fastq.gz R2.fastq.gz output/ --annotate --skip-vv

# Resume annotation only on existing VCFs
bash pipeline/run_pipeline.sh SAMPLE R1.fastq.gz R2.fastq.gz output/ --stages 3 --skip-vv

# Batch processing
anneal manifest --dir /path/to/fastqs/ -o samples.tsv
bash pipeline/run_pipeline_batch.sh samples.tsv output/
bash pipeline/run_pipeline_batch.sh samples.tsv output/ --annotate --skip-vv
```

## Pipeline Stages

### Stage 1: Consensus Generation
Barcode extraction -> BWA-MEM2 alignment -> Family grouping -> SSCS -> Singleton correction -> DCS

```
Output: {sample}.sscs.sc.sorted.bam, {sample}.dcs.sc.sorted.bam
        {sample}.stats.txt, {sample}.family_sizes.tsv, {sample}.family_sizes.png
```

### Stage 2: Variant Calling
samtools mpileup -> Rust variant caller (quality-filtered SNV + indel detection)

```
Output: {sample}.sscs.vcf, {sample}.dcs.vcf
VCF FORMAT: GT:ALT:TOT:FRAC (ALT=alt reads, TOT=depth, FRAC=VAF)
```

### Stage 3: Variant Annotation (optional)
VEP + ANNOVAR -> Clinical filtering -> HGVS validation (VariantValidator)

```
Output: {sample}.{dcs|sscs}.annotated.tsv -> .filtered.tsv -> .clinical.tsv
Requires: VEP conda env, ANNOVAR, VariantValidator Docker
Skip HGVS step with: --skip-vv
```

## Output Structure

```
output/
  {sample}/
    consensus/
      {sample}.sscs.sc.sorted.bam      # SSCS (with singleton correction)
      {sample}.dcs.sc.sorted.bam       # DCS (highest confidence)
      {sample}.stats.txt               # Consensus statistics
      {sample}.family_sizes.tsv        # Family size distribution
      {sample}.family_sizes.png        # Distribution plot
    variants/
      {sample}.sscs.vcf                # SSCS variant calls
      {sample}.dcs.vcf                 # DCS variant calls
    annotated/                         # Stage 3 (optional)
      {sample}.dcs.annotated.tsv
      {sample}.dcs.filtered.tsv
      {sample}.dcs.clinical.tsv
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| --cutoff | 0.6 | Consensus base call fraction (0.5-1.0) |
| --min-qual | 30 | Minimum base quality for consensus |
| --min-family | 2 | Minimum reads per family for SSCS |
| --bpattern | NNNSS | UMI pattern: 3bp random + 2bp spacer |
| --aligner | bwa-mem2 | bwa-mem2, parabricks, minimap2, bwa |
| --singleton-correction | on | Singleton rescue strategies |

## Dependencies

**Rust** (1.93+): Anneal core + variant caller
**Conda** (anneal env): samtools, bwa-mem2, python, matplotlib, pandas, numpy
**Annotation** (optional): VEP, ANNOVAR, VariantValidator Docker

## Configuration

Edit `pipeline/config.sh` with your local paths before running. Key variables:

- `REFERENCE`: hg38 masked FASTA (U2AF1-fixed Broad assembly)
- `BEDFILE`: Target panel BED file
- `SEQUENCES_DIR`: Input FASTQ directory
- `RESULTS_DIR`: Output directory
- `USE_GPU`: false for CPU-only (default), true for CUDA

All other paths auto-resolve from the installation directory.

## Documentation

- `SETUP_GUIDE.md`: Step-by-step fresh installation guide
- `ARCHITECTURE.md`: Architecture notes and design decisions

## References

- Wang TT, Abelson S, et al. (2019) Nucleic Acids Research, 47(15), e87
- Kennedy SR, et al. (2014) Nature Protocols, 9(11), 2586-2606
- ConsensusCruncher: https://github.com/pughlab/ConsensusCruncher
