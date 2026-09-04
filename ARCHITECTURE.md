# Anneal -- Architecture & Design Notes

Guidance for working in this repository. Anneal is a duplex consensus pipeline
for UMI-based error suppression and ultra-sensitive variant detection.

## Project Overview

Anneal takes UMI-tagged paired-end FASTQs and produces duplex consensus BAMs,
then calls SNVs and indels from those high-confidence consensus reads. An
optional annotation stage adds functional annotation and HGVS validation.

The core engine is written in Rust with CUDA acceleration; GPU is the only
supported and validated configuration. Variant calling uses Pisces. MRD calling
is marker tracking against a per-site background model -- see the MRD section
below. Orchestration is a set of bash stage scripts driven by a single shared
`config.sh`.

## Pipeline Stages

| Stage | Name | Tool | Output |
|-------|------|------|--------|
| 1 | Consensus generation | `anneal` (Rust) | `{sample}.sscs.sc.sorted.bam`, `{sample}.dcs.sc.sorted.bam`, stats, family sizes |
| 2 | Variant calling | Pisces 5.2.10.49 | `{sample}.{sscs,dcs}.pisces.genome.vcf` |
| 3 | Annotation (optional) | VEP + ANNOVAR + VariantValidator | `{sample}.{dcs,sscs}.annotated/filtered/clinical.tsv` |
| 4 | FLT3-ITD (optional) | getITD | `{sample}.flt3_itds.tsv` |

Stage 1 internals: barcode extraction -> BWA-MEM2 alignment -> family grouping
-> SSCS (single-strand consensus) -> singleton correction -> DCS (duplex
consensus).

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
pipeline/
  config.sh                   # shared config -- edit paths here
  stage1_consensus.sh
  stage2_variant_calling.sh
  stage3_annotate.sh          # optional
  stage4_flt3.sh              # optional, FLT3-ITD via getITD
  docker-apptainer-shim.sh    # copy to ~/bin/docker; see SETUP_GUIDE
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

### Pisces as the variant caller
Chosen over LoFreq, VarDict, Mutect2 and SiNVICT on the dilution series. Two
accommodations, both inside stage 2: the XV tag is stripped, since anneal writes
it as a string and Pisces expects an integer; and a panel-restricted masked
genome is used, since Pisces exhausts memory on the full 3,366-contig hg38
reference.

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
bash deploy.sh                                   # build anneal (gpu)
bash pipeline/run_pipeline.sh SAMPLE R1 R2 out/  # stages 1,2
bash pipeline/run_pipeline.sh SAMPLE R1 R2 out/ --annotate --skip-vv   # + stage 3
bash pipeline/run_pipeline.sh SAMPLE R1 R2 out/ --flt3                  # + stage 4
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
# MRD calling: marker-tracking configuration

## Scope

anneal calls MRD by tracking a patient's known diagnosis variants into follow-up
samples. It does not perform untargeted variant discovery at MRD sensitivity.

This is a deliberate restriction, established by dilution experiment rather than
assumed. Scoring a DCS sample against the background model without restricting
to a marker list yields roughly 2,000 calls at p < 0.005. That count does not
change across a 125-fold dilution series:

| rung | untargeted calls (p < 1e-9) | markers detected |
|------|------------------------------|------------------|
| G    | 195                          | 3                |
| H    | 190                          | 2                |
| I    | 160                          | 1                |
| J    | 190                          | 0                |

Tumour content falls 5x per rung. The marker count tracks it; the untargeted
count does not. The surviving untargeted calls are therefore per-sample
systematic artifact, not signal, and no p-value threshold separates them —
tightening from 0.005 to 1e-9 reduces the count but leaves it equally flat.

Restricting to the diagnosis-variant list removes the problem rather than
solving it: with 3-5 positions under test, the artifacts at other positions are
never examined.

## Configuration

| component | setting |
|-----------|---------|
| consensus track | DCS (`dcs.sc.sorted.bam`) |
| variant caller | Pisces 5.2.10.49 (discovery); marker scoring reads the BAM directly |
| background model | `beta_matrix_DCS.txt`, Waalkes per-position alpha/beta from 8 BNCs |
| substitution test | beta-binomial upper tail |
| indel test | recurrence blocklist (`indel_blocklist.tsv`) |
| artifact mask | `artifact_mask.combined.bed`, 624 loci |
| strand filter | reject if alt >= 10 reads and >= 90% on one strand |
| alpha | 0.005 |
| minimum alt reads | 2 |
| scope | patient diagnosis variants only |

### Why the marker scorer does not read a VCF

Pisces applies `MinimumVariantQScore` (default 20) beneath any requested
`--minvf`. On the BNC panel this censors everything below roughly 7-12 alt
reads, and because the Q-score is depth-dependent the censoring point moves
between samples. A marker present at 3 reads in a remission sample would be
absent from the caller output entirely.

`call_mrd_markers.py` therefore counts from the BAM at each marker position
regardless of what any caller emitted. Absence in the report means absence in
the reads.

Reads are counted unfiltered, matching how the beta matrix was built (Pisces on
the unfiltered `*_pisces.bam`), so both sides of the comparison see the same
read population.

### Why alpha is 0.005 and not 1e-9

1e-9 was derived to suppress untargeted calls across ~3,600 tested positions. It
is a multiple-testing correction, and it does not apply when 3-5 positions are
under test with per-site background already modelled.

Carried into tracking mode it drops real markers. On the G rung:

| marker | p |
|--------|---|
| IDH2 R140Q | 8.75e-10 |
| PTPN11 F285L | 1.74e-09 |

At 1e-9 the first survives by a 14% margin and the second is lost. 0.005 is the
cutoff in the source method and is appropriate here. Sensitivity is bounded by
the reported per-locus LoD, not by alpha.

## Limit of detection

Every marker is reported with the smallest alt count that would have reached
significance at that locus, given the observed depth and that site's background,
and the VAF that count corresponds to.

A negative marker is only interpretable alongside that number. "Not detected,
sensitivity 0.03% at this locus" is a result; "not detected" alone is not.

LoD is depth-limited:

| DCS depth | LoD (2 reads) |
|-----------|---------------|
| 3,000     | 0.067%        |
| 5,000     | 0.040%        |
| 10,000    | 0.020%        |
| 20,000    | 0.010%        |

## Validated sensitivity

Two arms with different floors. The pipeline should not be described as having a
single sensitivity figure.

| arm | mechanism | demonstrated |
|-----|-----------|--------------|
| indel (NPM1 type A) | blocklist | 0.056% at 3 reads / 5,402 |
| substitution (IDH2 R140Q) | beta-binomial | 0.056% at 2 reads / 3,594 |

Reported per-locus LoD in that series ranged 0.043% to 0.075% at DCS depths of
2,800 to 7,000.

Reaching 0.02% requires roughly 10,000x DCS depth, i.e. 2-3x more sequencing
than the current configuration delivers.

## Dilution validation

Sample 25NGS1601, two-fold serial dilution across four rungs (G to J).

The expected VAFs recorded on the dilution worksheet are nominal, calculated
from dilution factors rather than measured. They run 1.5x to 4x away from
observed at the lower rungs and should not be used as ground truth. The table
below reports read counts taken directly from the consensus BAMs at each marker
position, base quality 20.

### IDH2 R140Q, chr15 C>T

| rung | nominal | SSCS | DCS | DCS call |
|------|---------|------|-----|----------|
| G | 0.387% | 0.417% (90/21,564) | 0.336% (11/3,274) | detected |
| H | 0.193% | 0.113% (24/21,255) | 0.056% (2/3,594) | detected |
| I | 0.097% | 0.024% (5/21,114) | 0 (0/3,324) | not detected |
| J | 0.048% | 0.049% (10/20,520) | 0 (0/2,829) | not detected |

True content at rung I is 0.024%, below the 0.060% LoD that rung supports, so
the DCS non-detection is correct. The nominal 0.097% would have implied a false
negative; the counts show otherwise.

SSCS at rungs I and J is at its own noise floor, not tracking dilution: 10 reads
at J exceeds 5 at I on equal depth, when J should be half of I. SSCS would have
returned two quantitative-looking values where DCS correctly reported nothing.

### PTPN11 F285L, chr12 T>C

| rung | nominal | SSCS | DCS | DCS call |
|------|---------|------|-----|----------|
| G | 0.032% | 0.083% (15/18,021) | 0.141% (4/2,832) | detected |
| H | 0.016% | 0.005% (1/20,067) | 0 (0/3,524) | not detected |
| I | 0.008% | 0.006% (1/17,360) | 0 (0/2,619) | not detected |
| J | 0.004% | 0.006% (1/17,970) | 0 (0/2,796) | not detected |

Background at this position is one read in roughly 18,000, flat across H, I and
J. The 15 reads at G are a 15-fold excess over that floor and both tracks agree,
so the detection is real despite the nominal value sitting below the LoD.

### NPM1 type A insertion, chr5

| rung | nominal | DCS | DCS call |
|------|---------|-----|----------|
| G | 0.284% | 0.214% (12/5,620) | detected |
| H | 0.142% | 0.156% (11/7,039) | detected |
| I | 0.071% | 0.056% (3/5,402) | detected |
| J | 0.035% | 0 (0/5,335) | not detected |

Scored through the indel blocklist, not the beta model. Quantitative across
three rungs; the dropout at J is below the 0.056% LoD for that rung.

### CEBPA H24Afs

Not evaluable at any rung. Depth is 133x to 222x against 2,800x to 7,000x
elsewhere on the panel, roughly a twenty-fold shortfall. Capture problem, not a
caller problem. CEBPA should be excluded as a trackable marker until the panel
is addressed.

### Summary

Every DCS call and non-call across the four markers and four rungs is correct
when assessed against measured counts. No false positives, no false negatives.

Deepest true detections: NPM1 at 0.056% (3 reads in 5,402), IDH2 at 0.056%
(2 reads in 3,594).

### Untargeted scoring, same samples

Scored without restricting to the marker list, at p < 1e-9:

| rung | calls |
|------|-------|
| G | 195 |
| H | 190 |
| I | 160 |
| J | 190 |

Flat across the series while marker detections fall 3, 2, 1, 0. The untargeted
survivors are per-sample systematic artifact. This is the basis for restricting
the assay to marker tracking.

## Known limitations

**CEBPA is not trackable.** The CEBPA probe returns 133-222x against 2,800-7,000x
elsewhere on the panel, a roughly 20-fold shortfall. It is non-evaluable at every
dilution rung and will remain so until the capture is addressed. This is a panel
issue, not a caller issue.

**Background model floors on a literature constant.** 83% of positions in the DCS
matrix carry the default rather than a fitted value, because 8 BNCs at ~3,000x
DCS provide roughly 26,000 molecules per site while confirming a rate of
1/200,000 requires an order of magnitude more. The default is taken from Wang et
al. NAR 2019 supplementary table 5 (SmallDeep panel). It is not measurable from
the current control panel.

**Untargeted discovery is not supported.** See Scope.

## Stage 5: background scoring (0.3.0)

Stage 5 is the point at which variant calls and the error model meet. It
runs per sample and per track (`SCORE_TRACKS`, default `dcs sscs`) and is
blind: the candidate list is derived from the sample's own reads and no
diagnosis variant list is read anywhere in stages 1-5.

**Candidates.** `scripts/error_model/build_candidates.py` takes the union of
the Pisces variants-only VCF and the stage 2b indel table
(`scan_indels.py`), both VCF-anchored, deduplicated on chrom/pos/ref/alt.
Multi-allelic VCF records are split. Alleles containing N are excluded: in
the corrected engine an N inside a consensus read means the two strands
disagreed at that base, so `G>GN` is an unresolved event, not an allele.
Each candidate carries a label `Gene|Consequence|HGVSp;src=P|S|PS` from the
stage 3 tables (P = Pisces, S = indel scan). Output is the header-less
`chrom pos ref alt label` format that `call_mrd_markers.py` reads.

**Scoring.** `call_mrd_markers.py` counts reads at every candidate directly
from the consensus BAM (Pisces' Q20 gate would otherwise remove a marker at
a few reads), scores substitutions against the per-track beta matrix with a
beta-binomial tail at `--alpha-level 0.005`, and indels against the
per-track blocklist (`--indel-min-controls 6`, `--min-indel-alt 3`). The
artifact mask is not passed. Output:
`scored/{sample}.{track}.calls.tsv`, one row per candidate with alt count,
depth, VAF, per-orientation counts, background, p-value, minimum callable
count, per-locus LoD, call and note.

**Strand test.** Duplex reads are both-strand by construction, so
`alt_fwd`/`alt_rev` on DCS is the alignment orientation of the consensus
read, not strand of origin. A site covered by one mate only is one-sided
for reference reads as well. The test is therefore relative: with alt >=
`--strand-min-alt` (10) and the alt fraction >= `--strand-thresh` (0.90)
one-sided, a call is rejected only if the alt orientation differs from the
reference orientation at the same site (Fisher exact, two-sided, p <
`--strand-p` 1e-3). With no reference reads the absolute rule applies.
Rows that are one-sided but not biased are annotated
`one-sided coverage: ref F/R`.

**Marker overlay.** Tracking a patient's diagnosis variants is a separate,
post-hoc run of `call_mrd_markers.py` with `--markers <patient.tsv>` on the
stage 1 BAM, using the same matrix and blocklist. For a dilution or
longitudinal series the baseline sample's candidate list is scored on every
later sample; scoring each sample's own list produces spurious zeros for
markers that fell below the caller's gate.

## Background estimator (0.3.0)

`build_background.py` fits one Beta per (position, alt base) from the eight
BNC consensus BAMs per track. Pileups are taken directly (`-Q 30`), so no
site is censored by a caller. The mean is the pooled rate under a Jeffreys
prior, `(k + 0.5) / (N + 1)`, so a site with no observed alt is anchored at
about `0.5/N` rather than at an external constant.

The concentration (`alpha + beta`) decides how much between-control
variation the model tolerates:

1. Outlier control. The control with the most alt reads is tested against
   the pooled rate of the others. If it carries at least `--outlier-min-alt`
   (3) reads, is a Poisson outlier at `--outlier-p` (1e-3), and its rate is
   at least `--outlier-min-ratio` (10) times the others', it is dropped for
   that substitution and the exclusion is written to the report. This keeps
   one control's clone or library-specific event from setting the site's
   limit for every patient (IDH2 R140: 5 reads in one control, none in
   seven). The ratio keeps systematic artifact sites, where every control is
   high and one is merely highest, with the moment estimator.
2. Method of moments, only when estimable: `--mom-min-alt` (20) pooled alt
   reads in at least three non-zero controls. Binomial sampling variance is
   subtracted from the between-control variance first; the estimate is used
   when it implies more dispersion than pure binomial sampling.
3. Otherwise the fallback: binomial concentration divided by `--dispersion`
   (3.0), a conservative widening.

The report (`beta_matrix_{TRACK}.report.tsv`) records, per substitution,
controls used, pooled depth and count, alpha, beta, model rate, and a note
(`outlier_dropped:k/d`, `mom`, or empty).

## Reporting policy

DCS is the reported track. SSCS is confirmatory only. For each known
marker the report carries a DCS line and an SSCS line; a marker positive in
SSCS and negative in DCS is reported as "detected below duplex sensitivity"
with both counts and both limits. SSCS is not quantitative and is
damage-prone for C>T/G>A near its limit; SSCS indels in homopolymers are
blocklisted (CEBPA H24Afs is not evaluable on SSCS). See the 0.3.0
CHANGELOG for the dilution series that fixed this policy.
