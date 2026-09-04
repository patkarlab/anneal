## 0.3.0 - 2026-09-04

Validated release: corrected duplex engine, rebuilt background model, and a
scoring stage. Every number below was reproduced through `run_pipeline.sh` on
the DIL-A two-fold dilution series (G, H, I, J).

### Consensus engine (commit 691a92f)
- Four defects fixed: CIGAR restored to the family key; duplex consensus
  requires both strands to call the same non-N base at >= Q30 (no gap-fill);
  rescued singletons report family_size 1; per-strand depth emitted as
  `XA`/`XB`. DCS error rate on DIL-A-G 2.903e-04 -> 7.635e-05 at -1.1%
  depth; C>A (8-oxoG) 2.939e-05 -> 1.000e-06. Eight BNCs: DCS 6.22e-05,
  CV 12.4%.
- Tiering by min(XA, XB) is flat: `ANNEAL_MIN_READS_PER_STRAND` stays 1 and
  singleton correction stays on. Settled.
- Binary is `target_cpu/release/anneal`; `config.sh` defaults to it, to the
  BWA-indexed masked reference, and to the rebuilt background files.

### Background model
- `build_background.py` replaces the Pisces-VCF fit. BNC consensus BAMs are
  pileuped directly, so no site is censored by a caller; zero-count sites are
  anchored by a Jeffreys prior on pooled depth.
- Estimator corrections in this release:
  - binomial sampling variance is subtracted before between-control
    dispersion is derived, so sampling noise on a few reads no longer
    registers as dispersion;
  - method of moments only when estimable (>= 20 pooled alt reads in >= 3
    non-zero controls), otherwise the `--dispersion 3.0` fallback;
  - a control whose count is a Poisson outlier (p < 1e-3) at >= 10x the
    pooled rate of the other controls, with >= 3 reads, is dropped for that
    substitution and recorded in the report (`outlier_dropped:k/d`). The
    ratio keeps systematic artifact sites, where every control is high and
    one is merely highest, with the moment estimator.
- DCS matrix: 1,310 substitutions with a control dropped, 118 with measured
  dispersion; small-count MoM sites fell from 7,451 to 82. SSCS matrix:
  5,327 dropped, 6,555 measured.
- Effect at IDH2 R140 C>T (chr15:90088702): one control carried 5 reads
  (0.13%), the other seven none. Before: alpha 0.113, beta 738, 12 reads
  needed at 3,205x (0.37%). After: that control excluded, alpha 0.167, beta
  10,695, 3 reads (0.094%).
- Indel blocklists are per track (`build_indel_blocklist_v2.py`, from the
  `scan_indels.py` tables, consumed with `--indel-min-controls 6`). The
  artifact mask is retired: a 1% floor against a 0.005% background, and the
  per-site matrix now does its job.
- Files: `results_bnc/beta_matrix_{DCS,SSCS}.txt` (+ `.report.tsv`, now with
  a `note` column and `n_samples` = controls used). The previous
  Pisces-derived files are archived as `*.pisces_20260731.txt`; the
  intermediate `*.patched.txt` files are the comparator for this change.

### Scoring
- Stage 5 (`pipeline/stage5_score.sh`; `--stages 5` or `--score`). Per
  track, candidates are the union of the Pisces VCF and the CIGAR indel scan,
  VCF-anchored, deduplicated, N-containing alleles excluded, labelled from
  the annotation tables (`scripts/error_model/build_candidates.py`), then
  scored by `call_mrd_markers.py` against the per-track matrix and blocklist
  into `scored/{sample}.{track}.calls.tsv`. Previously the indel scan never
  entered scoring, so a marker below Pisces' Q20 gate had no route to a
  call. The Pisces records reproduce the 27 Aug by-hand verdicts exactly.
- Strand filter is relative. With alt >= `--strand-min-alt` and the alt
  fraction >= `--strand-thresh` one-sided, a call is rejected only if alt
  orientation differs from reference orientation at the same site (Fisher
  exact, two-sided, p < `--strand-p` 1e-3). Duplex reads are both-strand by
  construction, so orientation alone had censored sites covered by one mate
  (IDH1, TP53, PHF6 at 10-13 reads); the PTPN11 A>C cluster against balanced
  reference reads still fails. Rows that are one-sided but not biased carry
  the note "one-sided coverage: ref F/R".
- Marker overlay stays a post-hoc `call_mrd_markers.py` run with the
  patient's diagnosis list. Nothing in stages 1-5 reads that list.

### Validation: DIL-A, two-fold series, DCS (reported track)

| marker        | G                    | H                            | I                        | J                        |
|---------------|----------------------|------------------------------|--------------------------|--------------------------|
| NPM1 type A   | 22/7017 0.314% DET   | 14/8855 0.158% DET           | 6/6772 0.089% DET        | 0/6841 ND (LoD 0.044%)   |
| IDH2 R140Q    | 11/3205 0.343% DET   | 1/3492 0.029% ND (LoD 0.086%)| 0/3262 ND (LoD 0.092%)   | 0/2796 ND (LoD 0.072%)   |
| CEBPA H24Afs  | 4/219 1.83% DET      | 0/132 ND (LoD 2.3%)          | 0/222 ND (LoD 1.4%)      | 0/221 ND (LoD 1.4%)      |

G's candidate list scored on each rung: DETECTED 111 / 45 / 41 / 45. The
flat 45 are germline and diluent constants (SF3B1 intronic 7%, the FLT3
intronic SNPs, the RAD21 homopolymer ladder).

SSCS confirmatory track, same markers:

| marker        | G                     | H                     | I                     | J                     |
|---------------|-----------------------|-----------------------|-----------------------|-----------------------|
| NPM1 type A   | 147/51919 0.283% DET  | 49/52030 0.094% DET   | 83/49465 0.168% DET   | 9/47218 0.019% DET    |
| IDH2 R140Q    | 134/30397 0.441% DET  | 24/25312 0.095% DET   | 6/29377 0.020% ND     | 12/27277 0.044% DET   |
| CEBPA H24Afs  | not evaluable at any rung: SSCS blocklist (homopolymer slippage in 8/8 BNCs) |

IDH2 at H: SSCS shows the molecules at 0.095%; DCS has 1 of an expected 6
(p ~ 0.02); the pre-fix engine had 2 at the same site. G converts normally
(0.44% SSCS -> 0.34% DCS) and NPM1 in the same H library converts normally,
so this is molecule sampling at ~3,500x duplex depth, not an engine defect.

### Reporting policy (settled)
DCS is the reported track. SSCS is confirmatory only: a known marker
positive in SSCS and negative in DCS is reported as "detected below duplex
sensitivity" with both counts. SSCS is not quantitative (NPM1 SSCS 0.28 ->
0.094 -> 0.17 -> 0.019 while DCS halves cleanly) and is damage-prone for
C>T/G>A near its limit (IDH2 SSCS: I 0.020% ND, J 0.044% DET).

### Front-door reproduction
DX-1 run from FASTQ through `run_pipeline.sh --stages 1,2,3,4,5` on
`v0.3.0` with no overrides (`jobs/anneal_e2e.pbs`, A40, 100 min): consensus
statistics identical to the 27 Aug by-hand run (246,401,504 reads, 295,090
rescued, 788,886 DCS), Stage 5 identical (155 candidates, 54 DETECTED, no
differing call, NPM1 4/6526).

### Known issues
- CEBPA probe delivers 130-220x DCS: non-evaluable at every rung. Panel.
- TP53 chr17:7676341 T>C: controls disagree (5-15%), fitted concentration
  15; effectively non-evaluable.
- `CUTOFF=0.6` and the `read_len` tail behaviour (a position covered by one
  read in a multi-read family is called at proportion 1.0) are unchanged
  from the validated runs. Changing either means redoing the BNCs and the
  dilution.
- `scan_indels.py` reports insertions whose inserted base is N (strand
  disagreement inside the insertion) as separate records with the same
  support as the resolved allele. Excluded at candidate build, not yet at
  source. `strand_frac` in the indel table is read orientation, not strand
  of origin.
- DX-3 stage 1 outstanding (/scratch quota during alignment).
- Depth: ~2,800x DCS gives a 3-molecule limit near 0.1%; 0.05% needs
  ~6,000x. Run the saturation test (1601-G subsampled to 50% and 25%)
  before the next sequencing batch.
## v0.3.0 - Pisces in stage 2, FLT3-ITD as stage 4

### Fixed
- **Stage 2 could not run.** It invoked `mpileup_variant_caller`, the Rust
  caller removed in v0.2.0, so `config.sh` pointed at a binary that no longer
  existed. Stage 2 now runs Pisces, matching what the documentation already
  described and what the dilution validation was performed with.
- `activate_conda()` relaxes `set -u` across conda activation. The dotnet hooks
  installed for Pisces dereference unset variables and abort any script using
  `set -euo pipefail`, which silently killed several batch jobs.

### Added
- **Stage 4: FLT3-ITD** (`--flt3`, or `--stages 4`). Extracts the FLT3 exon
  14/15 window from the DCS consensus BAM and runs getITD. Output is getITD's
  own, unfiltered; matching against a patient's known ITD happens downstream.
  - `-min_read_copies 1`: getITD deduplicates identical reads by default,
    assuming raw data with PCR duplicates. UMI families already did that.
  - Primer filtering left on. Nominally an amplicon filter, but on hybrid
    capture it removes reads with indels near read edges, the main source of
    spurious ITD calls. Disabling it quadruples usable reads and increases
    false positives.
  - Consensus reads are single-end, so getITD receives one FASTQ. Tools
    requiring paired input, such as FLT3_ITD_ext, cannot read consensus BAMs.
- `bin/docker` - shim translating anneal's hardcoded `docker run` for
  Parabricks into `apptainer exec`. Docker is unusable from batch jobs on this
  cluster (permission denied on the daemon socket); apptainer is not. `bin/` is
  already prepended to PATH by `activate_conda()`.
- `scripts/qc_panel.py` - per-gene and per-exon coverage with limit of
  detection, alongside the existing `coverage_plot.py`.

### Stage 2 notes
Pisces needs two accommodations, both handled inside the stage:
- The XV tag is stripped. anneal writes `XV:Z:SSCS` / `XV:Z:DCS` as a string;
  Pisces expects an integer and throws. Each BAM is single-track, so XV carries
  nothing Pisces needs.
- A panel-restricted masked genome (`PISCES_GENOME`) is used. Pisces exhausts
  memory on the full 3,366-contig hg38 reference.

The interval list is regenerated per sample from the panel BED, so no derived
file needs to be tracked. Pisces runs on the full consensus BAM - the interval
list does the restriction, measured at 1.0 GB peak and 25 s per track.

Both a gVCF and a variants-only VCF are written; stage 3 consumes the latter at
its existing filename.

### Known limitations
- **FLT3-ITD discovery is not supported.** Untargeted getITD on DCS returns 1-5
  spurious ITDs per sample with read counts indistinguishable from a true call.
  On DIL-B the confirmed 63 bp ITD was found at rung G with 9 reads, while
  artifacts at other rungs carried 5-10. Filter to the patient's known ITD
  length and insertion site; every artifact observed sat 44-60 bp away.
- **Long ITDs are invisible to consensus.** Family grouping needs an alignment
  position, so unmapped reads are dropped: a raw BAM carried 944,067, SSCS and
  DCS exactly zero. Soft-clipped reads survive at 0.77 the rate of all reads, so
  ITDs short enough to still align are retained. 63 bp is confirmed; the ceiling
  above that is untested.
- getITD reports hg19 coordinates against its own FLT3 reference, while the BAM
  is queried in the build of `BEDFILE`. Both are correct and they are different
  coordinate spaces.
## v0.2.1 - validation corrected against measured counts

No code change. Corrects the dilution validation recorded in v0.2.0.

### Corrected
- The DIL-A series is a two-fold serial dilution, not five-fold. The
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
Sample DIL-A, 1/5 serial dilution across four rungs. Markers detected fall
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
