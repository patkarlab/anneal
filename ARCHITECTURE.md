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
| 4 | MRD marker tracking | `call_mrd_markers.py` | `{sample}.mrd_report.tsv` |

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
| substitution | beta-binomial | 0.138% at 5 reads / 3,619 |

Reported per-locus LoD in that series ranged 0.028% to 0.056% at DCS depths of
2,800 to 7,000.

Reaching 0.02% requires roughly 10,000x DCS depth, i.e. 2-3x more sequencing
than the current configuration delivers.

### Dilution behaviour

Sample 25NGS1601, 1/5 serial dilution:

| marker | G | H | I | J |
|--------|---|---|---|---|
| IDH2 R140Q | 0.363% | 0.138% | not detected | not detected |
| NPM1 type A | 0.214% | 0.156% | 0.056% | not detected |
| PTPN11 F285L | 0.140% | not detected | not detected | not detected |
| CEBPA H24Afs | not evaluable at any rung | | | |

Every dropout occurs where the expected VAF falls below that locus's reported
LoD. Nothing disappears above its stated detection limit and nothing appears
below it. PTPN11 steps 4.95x from G to H against a nominal 5x dilution.

At 5-12 alt reads counting noise is large; treat these as detection rather than
quantification.

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
