## v0.2.1 - validation corrected against measured counts

No code change. Corrects the dilution validation recorded in v0.2.0.

### Corrected
- The 25NGS1601 series is a two-fold serial dilution, not five-fold. The
  previously reported "PTPN11 steps 4.95x against a nominal 5x" was a
  coincidence read as a result.
- Expected VAFs on the dilution worksheet are nominal, derived from dilution
  factors rather than measured. They run 1.5x to 4x from observed at the lower
  rungs. Validation is now quoted against read counts taken from the consensus
  BAMs.
- IDH2 at rung I was briefly recorded as a false negative above the stated LoD.
  It is not: measured content there is 0.024% against a 0.060% LoD, so the
  non-detection is correct. The nominal 0.097% was the source of the error.
- PTPN11 at rung G was briefly recorded as a probable false positive. It is
  real: background at that position is one read in 18,000, flat across the
  lower three rungs, and the 15 SSCS reads at G are a 15-fold excess.

### Result
Every DCS call and non-call across four markers and four rungs is correct when
assessed against measured counts. No false positives, no false negatives.

Deepest true detections: NPM1 0.056% (3 reads in 5,402), IDH2 0.056% (2 reads
in 3,594).

### Note on SSCS
SSCS at the lowest rungs sits at its own noise floor rather than tracking
dilution - IDH2 reads 10 at rung J against 5 at rung I on equal depth, when J
should be half of I. DCS correctly reported non-detection at both. This is the
practical case for DCS as the primary substitution track.
## v0.2.0 - MRD marker tracking

Locks the MRD calling configuration as a marker-tracking assay. Untargeted
discovery at MRD sensitivity is explicitly out of scope.

### Added
- `scripts/error_model/call_mrd_markers.py` - scores a patient's diagnosis
  variants against the BNC beta matrix, reporting VAF, p-value, alt reads,
  depth and per-locus limit of detection. Counts from the BAM directly rather
  than from a VCF, so markers below the caller's Q-score floor are not lost.
- `scripts/error_model/build_indel_blocklist.py`
- `scripts/error_model/derive_artifact_mask.py`
- `Cargo.lock` now tracked, for reproducible builds.

### Removed
- `mpileup_variant_caller/` - superseded by Pisces.
- CPU-only pipeline scripts. The pipeline is GPU-only; `--no-gpu` remains in the
  binary but is unsupported and unvalidated.

### Changed
- Background model is `beta_matrix_DCS.txt` (Waalkes per-position alpha/beta,
  8 BNCs, Pisces-derived). DCS is the primary track.
- Sensitivity is reported per arm rather than as a single figure. Previous
  documentation quoted ~0.02% without qualification; validated performance is
  0.056% for the indel arm and 0.138% for the substitution arm at current depth,
  with per-locus LoD reported on every call. See ARCHITECTURE.md.
- `.gitignore` excludes sample manifests and marker lists. These contain lab
  accessions and diagnosis variant coordinates and must not be published.

### Validation
Sample 25NGS1601, 1/5 serial dilution across four rungs. Markers detected fall
3 to 2 to 1 to 0 across the series while untargeted calls at the same threshold
stay flat at 195, 190, 160, 190 - establishing that untargeted survivors are
per-sample artifact and that marker tracking is the supported mode. Every marker
dropout occurs where the expected VAF falls below that locus's reported LoD.

### Known issues
- CEBPA probe delivers 133-222x against 2,800-7,000x elsewhere. Non-evaluable at
  every rung; panel issue.
- `beta_distribution.py` computes the default error rate as
  `default_error_rate / No_of_samples`, placing the floor 8x below the cited
  literature value (1/120,000 rather than 1/15,000 for SSCS; 1/1,600,000 rather
  than 1/200,000 for DCS). Measured effect at DCS depth is negligible - two alt
  reads clear significance under either floor below ~20,000x - but the value
  does not match the documentation. Fix pending upstream.
- Scripts that build the background matrix are not yet in this repository.
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
