# Changelog

## 0.1.0 -- first public release

First public release of Anneal: a duplex consensus pipeline for UMI-based
error suppression and ultra-sensitive SNV/indel detection.

### Changed
- Removed FLT3-ITD detection (formerly Stage 3: getITD + FiLT3R with a
  concordance report). FLT3 exon regions remain in the target panel BED, so
  FLT3 point/TKD mutations are still called as standard SNVs/indels in Stage 2.
- Renumbered variant annotation from Stage 4 to Stage 3. It remains optional
  (`--annotate`) and unchanged in behavior.
- Default single/batch runs now execute Stages 1-2 (consensus + variant
  calling). Removed the `--no-flt3` flag.
- Server-agnostic paths and documentation; legacy private install paths removed
  from the build script.

### Pipeline
- Stage 1: Duplex consensus generation (Rust; optional CUDA)
- Stage 2: Variant calling (samtools mpileup + Rust caller)
- Stage 3: Variant annotation (VEP + ANNOVAR + VariantValidator; optional)
