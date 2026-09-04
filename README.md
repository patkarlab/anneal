# anneal

Duplex-sequencing pipeline for measurable residual disease (MRD) in acute
myeloid leukaemia. UMI-tagged paired-end reads are aligned, collapsed into
single-strand (SSCS) and duplex (DCS) consensus reads, called for
substitutions and indels, annotated, and scored against a per-site
background model built from biological negative controls. Version 0.3.0.

## Pipeline

| Stage | Script | What it does | Output |
|-------|--------|--------------|--------|
| 1 | `stage1_consensus.sh` | Parabricks alignment on the GPU; UMI family grouping and SSCS/DCS consensus in the Rust `anneal` binary, with singleton correction | `consensus/<sample>.{sscs,dcs}.sc.sorted.bam`, statistics, family-size plot |
| 2 | `stage2_variant_calling.sh` | Pisces on each track (variants-only VCF plus gVCF); CIGAR-based indel scan, since Pisces does not report indels at MRD frequencies | `variants/<sample>.<track>.vcf`, `<sample>.<track>.indels.tsv` |
| 3 | `stage3_annotate.sh` | VEP and ANNOVAR for substitutions and indels | `annotated/<sample>.<track>.{annotated,filtered,clinical}.tsv`, `<sample>.<track>.indels.annotated.tsv` |
| 4 | `stage4_flt3.sh` | getITD on the DCS reads over FLT3 exons 14 and 15 | `flt3/` |
| 5 | `stage5_score.sh` | Blind candidate list (Pisces calls plus indel scan) scored against the per-track background model and indel blocklist | `scored/<sample>.<track>.calls.tsv` |

Stage 5 is where calls and the error model meet. Each row of the calls table
carries alt count, depth, VAF, per-orientation counts, the site's modelled
background, a beta-binomial p-value, the minimum alt count that would be
callable at that depth, the per-locus limit of detection, and a call
(`DETECTED`, `not_detected`, `not_evaluable`) with a note. No diagnosis
variant list is read anywhere in stages 1-5; marker tracking is a separate,
post-hoc step (see Reporting).

## Running

```bash
# one sample, all stages
bash pipeline/run_pipeline.sh SAMPLE R1.fastq.gz R2.fastq.gz OUTDIR --stages 1,2,3,4,5

# selected stages on an existing sample directory
bash pipeline/run_pipeline.sh SAMPLE x x OUTDIR --stages 5

# a manifest of samples (sample_name, fastq1, fastq2)
bash pipeline/run_pipeline_batch.sh manifest.tsv OUTDIR
```

`--annotate`, `--flt3` and `--score` add stages 3, 4 and 5 to the default
`1,2`. All paths, tool locations and parameters live in `pipeline/config.sh`;
the defaults are the validated configuration (Parabricks through the
Apptainer shim, the CPU consensus build, cutoff 0.6, singleton correction on,
the per-track background files). Every value can be overridden through the
environment. `jobs/anneal_e2e.pbs` runs one sample end to end on a GPU node
and is the reference invocation. `SETUP_GUIDE.md` covers installation.

Requirements: CUDA GPU node for stage 1 (Parabricks via Apptainer), Rust
toolchain for the consensus binary, conda environment `anneal` (pysam,
numpy, pandas), Pisces 5.2 with the .NET 8 runtime, a separate `vep`
environment, ANNOVAR, getITD. No scipy anywhere.

## Background model

`build_background.py` fits a Beta distribution per position and alternate
base from the consensus BAMs of eight biological negative controls, per
track. Pileups are taken directly, so no site is censored by a caller.
Zero-count sites are anchored by a Jeffreys prior on pooled depth. The
concentration comes from the controls: a control carrying a sample-specific
event at a site (a clone, or an artifact confined to one library) is dropped
for that substitution and recorded; between-control dispersion is estimated
by method of moments only where the counts support it, with binomial
sampling variance subtracted; otherwise a conservative fallback applies.
`jobs/rebuild_background.pbs` rebuilds both matrices. Indel artifacts are
handled by per-track recurrence blocklists (`build_indel_blocklist_v2.py`).

## Reporting

DCS is the reported track. SSCS is confirmatory: for a known marker the
report carries a DCS line and an SSCS line, and a marker positive in SSCS
and negative in DCS is reported as "detected below duplex sensitivity" with
both counts and both limits. SSCS is not quantitative and is damage-prone
for C>T/G>A near its limit. Marker tracking is run post hoc with
`scripts/error_model/call_mrd_markers.py --markers <patient markers>` on the
stage 1 BAMs; for a longitudinal series the baseline candidate list is
scored on every later sample.

## Validation

A diagnostic sample carrying NPM1 type A and IDH2 R140Q was diluted two-fold
across four rungs (G to J) and run end to end. On DCS, NPM1 halves per rung
(0.31%, 0.16%, 0.09%) and drops out at J against a 0.04% limit; IDH2 is
detected at G (0.34%) and lost from H where the library carried too few
duplex molecules, while SSCS still shows it (0.095%). Untargeted calls from
the top rung fall from 111 to a flat 45 that are germline and diluent
constants. A second diagnostic sample run from FASTQ through the front door
reproduced its earlier by-hand run exactly: identical consensus statistics
and identical calls. Details and the full tables are in `CHANGELOG.md`;
design notes are in `ARCHITECTURE.md`.

The corrected duplex engine (0.3.0) took the DCS error rate on the top rung
from 2.9e-04 to 7.6e-05 at unchanged depth, with C>A (8-oxoG) falling 29x,
which is the signature of duplex agreement being enforced. Per-locus limits
at the current ~2,800x duplex depth are near 0.1% for three molecules; 0.05%
needs roughly 6,000x.

## Known limits

- One low-coverage probe (CEBPA) is non-evaluable at every dilution rung.
- A TP53 intronic position where controls disagree widely is effectively
  non-evaluable.
- Insertion alleles containing N (strand disagreement inside the insertion)
  are excluded at candidate build, not yet at the scanner.
- Consensus cutoff 0.6 and the tail-position behaviour of the consensus
  caller are unchanged from the validated runs; changing either requires
  rebuilding the background model and revalidating.

## Layout

```
pipeline/            config.sh, run_pipeline.sh, stage1-5 scripts, batch runner
src/                 Rust consensus engine (anneal)
scripts/             annotation, indel scan, error_model/ (scoring, blocklists)
build_background.py  background matrix builder
jobs/                reference PBS jobs (end-to-end run, background rebuild)
ARCHITECTURE.md      design notes and rationale
CHANGELOG.md         release history with validation numbers
SETUP_GUIDE.md       installation
```

## Method background

Consensus construction follows ConsensusCruncher (Wang et al., NAR 2019),
including singleton correction. The per-site Beta background follows the
approach of Waalkes et al. (Haematologica 2017), rebuilt from pileups with
the estimator changes described above.
