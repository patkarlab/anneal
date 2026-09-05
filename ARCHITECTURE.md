# anneal: architecture and design notes

Version 0.3.x. What the pipeline does and why it does it that way. Running
it is in `SOP.md`; installing it in `SETUP_GUIDE.md`; the validation record
and history in `CHANGELOG.md`.

## 1. Overview

anneal turns UMI-tagged paired-end reads from the AML MRD hybrid-capture
panel (182 probes, hg38) into a scored table of variant calls per sample.
The chain is: alignment; grouping of reads into per-molecule, per-strand
families by UMI, position and alignment; single-strand consensus (SSCS);
duplex consensus (DCS) requiring the two strands to agree; substitution
calling and a CIGAR-based indel scan on both tracks; annotation; and
scoring of every candidate against a per-site background model built from
eight biological negative controls.

Two properties are fixed by design. The pipeline is blind: stages 1–5 never
read a patient's known mutations, and tracking diagnosis markers is a
separate command on the finished consensus BAMs. And DCS is the reported
track: SSCS exists to confirm, never to quantify.

## 2. Data flow

| Stage | Input | Tool | Output |
|-------|-------|------|--------|
| 1 | FASTQ (UMI in the read name, pattern `NNNSS`) | Parabricks fq2bam (GPU); `anneal` (Rust, CPU) | `consensus/<sample>.{sscs,dcs}.sc.sorted.bam`, `stats.txt`, family-size plot |
| 2 | consensus BAMs | Pisces 5.2 (`--minvf 1e-4`, XV tag stripped, panel reference); `scan_indels.py` | `variants/<sample>.<track>.vcf`, `.pisces.genome.vcf`, `.indels.tsv` |
| 3 | VCF, indel table | VEP (`--flag_pick --everything --hgvs`), ANNOVAR | `annotated/<sample>.<track>.{annotated,filtered,clinical}.tsv`, `.indels.annotated.tsv` |
| 4 | DCS BAM | getITD (`-min_read_copies 1`) on FLT3 exons 14–15 | `flt3/` |
| 5 | VCF, indel table, annotation, background model | `build_candidates.py`, `call_mrd_markers.py`, `report_calls.py` | `scored/<sample>.<track>.{candidates,calls,report}.tsv` |

Stage 1 is the only GPU step (alignment, ~25 min); consensus and stages 2–5
are CPU. `run_pipeline.sh` chains the stages; `jobs/anneal_e2e.pbs` runs one
sample end to end on a GPU node and `jobs/anneal_cohort.pbs` does the same
as an array over a sample sheet.

## 3. Consensus engine

The engine follows ConsensusCruncher (Wang et al. 2019) with singleton
correction, implemented in Rust.

**Family key.** A read belongs to a family by UMI pair, chromosome, start
positions of read and mate, strand, read number, and the ordered CIGAR
pair of read and mate. CIGAR in the key is what makes index-wise consensus
valid: reads with different soft-clipping or indels vote at shifted
positions otherwise, and the emitted consensus would carry the alignment
of an arbitrary member.

**SSCS.** Per family, a base is called where the majority fraction reaches
`CUTOFF` (0.6) among bases of quality ≥ 30; otherwise N. Quality is the sum
of the members' qualities, capped at 60.

**Singleton correction.** A molecule represented by one read on one strand
is not discarded. The lone read is paired with the complementary strand's
SSCS (strategy 1) or its lone read (strategy 2) and must agree with it base
by base; disagreement is N. This is where anneal departs from plain duplex
calling, and it is why families of size 1+N enter the DCS pool. Per-strand
read counts are emitted as `XA` and `XB`; `ANNEAL_MIN_READS_PER_STRAND=1`
keeps 1+N molecules, a choice the tier analysis in the CHANGELOG supports
(the error rate is flat in min(XA, XB)).

**DCS.** Complementary SSCS reads (UMI halves swapped, strand flipped, same
coordinates and CIGAR) are combined; a base is kept only where both call
the same non-N base at ≥ Q30, otherwise N. No gap-fill: a position where
one strand is N stays N, since passing the other strand's base through
would make single-strand damage indistinguishable from duplex-confirmed
sequence.

**Why cutoff 0.6.** A 2:1 split in a three-read family calls the majority
base; at 0.7 it would be N. The choice was made for yield on this panel's
family-size distribution and is part of the validated configuration.
Changing it means rebuilding the background model and revalidating.

**Why only two BAMs are kept.** `sscs.sc` and `dcs.sc` (singleton-corrected)
are the tracks everything downstream consumes; the intermediate SSCS,
singleton and rescue BAMs are analysis artefacts and stage 1 removes them.

**0.3.0.** Four defects in the earlier engine (CIGAR missing from the key,
gap-fill in the duplex step, rescued singletons paired with their own
corrector, double-counted family sizes) made DCS no better than SSCS. Their
correction took the DCS error rate on the dilution top rung from 2.9e-04 to
7.6e-05 at unchanged depth, with C>A falling 29×. Details in the CHANGELOG.

## 4. Candidate generation

**Substitutions: Pisces.** Chosen after a comparison with an in-house
mpileup caller and with SSCS-aware alternatives: it is deterministic,
reports every site at `--minvf 1e-4`, and its calls on the DCS track matched
manual counts. Two requirements: the consensus `XV` tag must be stripped
(Pisces expects it to be an integer), and the reference must be restricted
to the panel chromosomes (the full contig set exhausts memory). UMI-aware
callers add nothing downstream of consensus reads: each consensus read is
already one molecule.

**Indels: CIGAR scan.** Pisces applies a variant-quality gate beneath any
requested frequency, and at MRD frequencies an NPM1 insertion at 20 reads
in 50,000 scores VQ 0 and is dropped. `scan_indels.py` counts insertions
and deletions directly from the CIGAR strings of the consensus BAM, with
strand counts, and reports recurrence in the controls. Insertions whose
inserted sequence contains N are strand disagreements inside the insertion,
not alleles, and are excluded when candidates are built.

**Stage 5 candidates.** The union of the Pisces VCF and the indel table,
VCF-anchored, deduplicated, labelled from the annotation. Nothing else
enters; in particular no diagnosis list.

## 5. Background model

A Beta distribution per position and alternate base, per track, from the
consensus BAMs of eight biological negative controls (`build_background.py`).

- Counts come from pileups of the control BAMs at every panel position
  (`-Q 30`), so no site is censored by a caller. An earlier version fitted
  Pisces VCFs and had left 84% of scored calls against a default constant.
- The mean is the pooled rate under a Jeffreys prior, `(k + 0.5)/(N + 1)`:
  a site never seen mutated in 35,000 control reads is anchored near
  1.4e-5 rather than at an external constant.
- Outlier control. The control with the most alt reads at a site is
  dropped for that substitution if it carries ≥ 3 reads, is a Poisson
  outlier (p < 1e-3) against the pooled rate of the other controls, and
  exceeds that rate ≥ 10×. This removes one control's clone or
  library-specific event from the site's limit (IDH2 R140Q: five reads in
  one control, none in seven) while keeping systematic artefact sites, where
  every control is high, with the moment estimator.
- Dispersion. Between-control dispersion is estimated by method of moments,
  with binomial sampling variance subtracted, only when ≥ 20 pooled alt
  reads sit in ≥ 3 controls; otherwise the binomial concentration divided by
  3 (a conservative widening). Small-count MoM had been reading a few reads
  as dispersion at three quarters of the non-zero sites.
- Indels have no per-site model. Per-track recurrence blocklists
  (`build_indel_blocklist_v2.py`) list alleles seen in the controls; an
  allele present in ≥ 6 of 8 (`INDEL_MIN_CONTROLS`) makes a candidate
  `not_evaluable`. A common germline homopolymer polymorphism the eight
  controls happen not to carry is therefore called, and the tier column
  makes that visible.
- The artefact mask of earlier versions is retired: a 1% floor against a
  0.005% measured background, and the per-site model does its job.

The report beside each matrix lists, per substitution, controls used,
pooled depth and count, alpha, beta, the model rate, and a note
(`outlier_dropped:k/d`, `mom`, or empty).

## 6. Scoring

`call_mrd_markers.py` counts reads at each candidate directly from the
consensus BAM (min base quality 30, the same filter the model was built
with) and tests the alt count against the site's Beta with a beta-binomial
tail. `min_alt_callable` is the smallest count that would pass at that
depth, and `lod_pct` that count as a VAF: the per-locus limit of detection.

**Alpha 0.005.** Marker tracking tests a handful of pre-specified sites
per patient, so a per-site alpha of 0.005 is the appropriate error rate
and matches Waalkes et al. Earlier untargeted work used 1e-9 as a crude
multiple-testing correction across the panel; that is not the clinical
mode and is not used.

**Calls.** `DETECTED` when the p-value is below alpha and the count is ≥
`--min-alt`; `not_evaluable` when the site has no model (germline in the
controls), is blocklisted, or is too shallow; `not_detected` otherwise, with
the reason in `note`.

**Strand test.** Duplex reads are both-strand by construction, so read
orientation is alignment orientation, and a site covered by one mate is
one-sided for reference reads too. With ≥ 10 alt reads ≥ 90% one-sided, a
call is rejected only if alt orientation differs from reference orientation
at the same site (Fisher exact, p < 1e-3); rows that are one-sided but not
biased carry the note `one-sided coverage: ref F/R`.

**Tier.** Presentational, beside the call: `high_vaf` (≥ 20%), `mid_vaf`
(5–20%), `mrd` (< 5%), `mrd_floor` (an indel at ≤ 3 reads, the weakest
evidence reported). `report_calls.py` writes the readable view: protein-
altering consequences only, label split into gene, consequence, protein
change and source, annotation joined, DETECTED first.

## 7. Marker tracking and reporting

The clinical question is whether a patient's diagnosis variants are present
in a follow-up sample. That is answered post hoc with
`call_mrd_markers.py --markers`, counting each marker in the stage 1 BAMs
and scoring it against the same model. Counting from the BAM rather than
reading a VCF is deliberate: a marker at a few reads would be absent from
any caller's output. For a series, the baseline sample's candidate list is
scored on every later sample; scoring each sample's own list produces
spurious zeros for markers that fell below the caller's gate.

DCS is reported. SSCS is confirmatory: it can show a known marker that DCS
has too few molecules for, reported as "detected below duplex sensitivity"
with both counts and both limits. SSCS is not quantitative (its dilution
series do not step cleanly) and is damage-prone for C>T/G>A near its limit.

## 8. Depth and limit of detection

The error rate is no longer the constraint; molecules are. At ~2,800×
duplex depth, three molecules is ~0.1%; 0.05% needs ~6,000×, 0.03% ~10,000×.
About 15–22% of input molecules yield a duplex read, and a marker can be
lost from a rung not because of noise but because too few of its molecules
formed duplexes (IDH2 at rung H in the dilution series: one duplex read
against six expected, with SSCS showing the molecules at 0.095%). Whether
more sequencing or more input DNA buys depth is what the saturation test
(`jobs/anneal_saturation.pbs`) measures; note that the binary's
`consensus --bedfile` mode double-emits families across overlapping probes
and must not be used for yield comparisons, which is why the job runs the
subsamples genome-wide.

## 9. Configuration

Every path, tool and parameter is in `pipeline/config.sh` as
`VAR="${VAR:-default}"`. The defaults are the validated configuration:
`ALIGNER=parabricks` through the Apptainer shim, `target_cpu/release/anneal`,
`CUTOFF=0.6`, min quality 30, singleton correction on,
`ANNEAL_MIN_READS_PER_STRAND=1`, the rebuilt matrices and per-track
blocklists, `INDEL_MIN_CONTROLS=6`, alpha 0.005, tiers 20/5. Changing any
engine or model parameter invalidates the matrices and the validation.

## 10. Known limitations

- One low-coverage probe (CEBPA, 130–220× DCS) is non-evaluable at MRD
  frequencies; a panel issue.
- A TP53 intronic position where controls disagree widely (5–15%) is
  effectively non-evaluable.
- Indels have no per-site statistical model; a three-read indel absent from
  the controls is reported, tiered `mrd_floor`.
- N-containing insertion alleles are excluded at candidate build, not yet
  at the scanner, which also reports them with the same support as the
  resolved allele.
- The consensus caller calls tail positions covered by a single read in a
  multi-read family with full confidence; unchanged from the validated
  configuration, and a change requires revalidation.
- `strand_frac` in the indel table is read orientation, not strand of origin.

## 11. Validation

The 0.3.0 engine and model were validated on a two-fold dilution series of
a diagnostic sample and by reproducing a by-hand run from FASTQ through the
front door with identical consensus statistics and identical calls. The
tables are in `CHANGELOG.md` and are the numbers any rebuild must reproduce.

## 12. References

- Wang et al. ConsensusCruncher: a tool for consensus sequence generation of
  UMI-tagged reads. Nucleic Acids Research 2019.
- Waalkes et al. Ultrasensitive detection of acute myeloid leukemia minimal
  residual disease using single molecule molecular inversion probes.
  Haematologica 2017.
- Kennedy et al. Detecting ultralow-frequency mutations by Duplex Sequencing.
  Nature Protocols 2014.
