# anneal: standard operating procedure

How to run the pipeline, from FASTQ to the report tables. Version 0.3.x.
Design and rationale are in `ARCHITECTURE.md`; validation numbers and
release history in `CHANGELOG.md`.

## 1. What the pipeline does

anneal takes UMI-tagged paired-end duplex sequencing reads from the AML MRD
panel (182 probes, hg38) and produces, per sample, a table of variant calls
scored against a background model built from eight biological negative
controls. It runs in five stages. A sample takes about 100 minutes
(about 200 for a 500M-read sample).

| Stage | Step | Tool | Runs on | Time |
|-------|------|------|---------|------|
| 1a | Alignment | Parabricks (fq2bam) through Apptainer | **GPU** (A40) | ~25 min |
| 1b | UMI grouping, SSCS and DCS consensus, singleton correction | `anneal` (Rust, CPU build) | CPU, 16 threads | ~30 min |
| 2 | Substitution calling on each track; CIGAR indel scan | Pisces 5.2; `scan_indels.py` | CPU | ~10 min |
| 3 | Annotation | VEP, ANNOVAR | CPU | ~15 min |
| 4 | FLT3-ITD | getITD on the DCS reads | CPU | ~2 min |
| 5 | Scoring against the background model; tier; report | `call_mrd_markers.py`, `report_calls.py` | CPU | ~10 min (SSCS is the slow track) |

Only alignment uses the GPU. The consensus step and stages 2–5 are CPU
work, and in the standard job they run on the same GPU node right after
alignment so that one sample is one self-contained job. The GPU is idle for
about three quarters of each job; throughput comes from running two or three
samples at once, and the limit on that is disk (section 4), not GPUs.

Two tracks are produced: **DCS** (duplex consensus, both strands agree) is
the reported track; **SSCS** (single-strand consensus) is confirmatory. The
pipeline is blind: it never reads a patient's diagnosis variants. Tracking
known markers is a separate command on the finished BAMs (section 8).

## 2. Dependencies and locations

Account `patkarlab-clinical` on the BioInfinix cluster, login node `ln1`
(`10.100.95.23`), institute network only.

| Item | Location |
|------|----------|
| Repository | `~/pipelines/anneal` |
| Configuration (all paths, tools, parameters) | `pipeline/config.sh` |
| Pipeline driver | `pipeline/run_pipeline.sh` (stages 1–5); stage scripts `pipeline/stage{1..5}_*.sh` |
| Consensus binary | `target_cpu/release/anneal` (the validated CPU build) |
| Conda environments | `anneal` (pysam, numpy, pandas, samtools, requests, .NET runtime for Pisces); `vep` |
| Parabricks | `~/pipelines/parabricks_4.3.1.sif` via `module load apptainer/1.5.1`; the shim `bin/docker` makes the binary's `docker run` call use Apptainer |
| Pisces | 5.2.10.49, panel-restricted reference `~/references/pisces_hg38_panel/` |
| VEP cache, ANNOVAR, getITD | `~/references/vep_cache/` and the paths in `config.sh` |
| References | `~/references/hg38_broad_bwa/Homo_sapiens_assembly38.masked.fasta` (alignment, background model); `~/references/hg38_broad/Homo_sapiens_assembly38.fasta` (unmasked, indel scan and annotation) |
| Panel BED | `~/pipelines/anneal/AML_MRD_DUPLEX_probes_hg38_sortd.bed` |
| Background model | `results_bnc/beta_matrix_DCS.txt`, `beta_matrix_SSCS.txt`, `indel_blocklist.DCS.patched.tsv`, `indel_blocklist.SSCS.patched.tsv` |
| Control BAMs (for rebuilding the model) | `results_bnc_patched/AMLMRD_DUPLEX_BNC{1..8}/` |
| Control FASTQs | `/scratch/patkarlab-clinical/Duplex_bncs/` — the only copy; never delete |
| Input FASTQs | `/scratch/patkarlab-clinical/<batch>/` |
| Outputs | `/scratch/patkarlab-clinical/<outdir>/<sample>/` |
| Jobs | `~/pipelines/anneal/jobs/*.pbs` |
| Job logs | `~/pipelines/anneal/logs/<jobid>.hn1.OU` |
| Download folder | `~/inbox/to_excel/` |

Queues: `A40b` (one node, one GPU, usually free), `a40` (two nodes, two GPUs
each, often queued), `short` (CPU, 6 h), `medium` (CPU, 72 h). The login node
kills anything that runs more than a few minutes, so all pipeline work is
submitted as jobs; only the quick checks and the HGVS validation (which
needs internet, and compute nodes have none) run on the login node.

## 3. What to change for a run, and what never to change

**Per run**, on the `qsub` line, never inside the scripts:

| Variable | Meaning |
|----------|---------|
| `FQ` | directory holding the batch's FASTQs |
| `OUT` | output directory for the batch |
| `MANIFEST` | sample sheet (section 4) |
| `FULL` | one sample's FASTQ name, for a single-sample run |

FASTQ names are `<sample><R1_SUFFIX>` and `<sample><R2_SUFFIX>`; the defaults are
the Illumina `_R1_001.fastq.gz` / `_R2_001.fastq.gz`, and other naming is passed on the
`qsub` line (`R1_SUFFIX=_R1.fastq.gz,R2_SUFFIX=_R2.fastq.gz` for the cohort batches).
The sample name is everything before `_S<n>` if present, otherwise the whole prefix.

**Never change** without rebuilding the background model and revalidating
the dilution series (see `CHANGELOG.md` for the numbers that must be
reproduced): the consensus binary, `CUTOFF=0.6`, minimum base quality 30,
`SINGLETON_CORRECTION=true`, `ANNEAL_MIN_READS_PER_STRAND=1`, the reference
files, the background matrices and blocklists, alpha 0.005,
`INDEL_MIN_CONTROLS=6`, and the estimator settings of `build_background.py`.
`config.sh` holds the validated values as defaults; a run should need no
overrides.

**Safe to change**: the tier thresholds (`TIER_HIGH_VAF=20`,
`TIER_MID_VAF=5`), which only label rows; queue and resource lines in the
jobs; `SCORE_TRACKS` if a run needs DCS only.

## 4. Preparing a batch

**Disk.** Each stage 1 needs 100–184 GB of transient space on `/scratch`
(soft quota 800 GB, hard limit 1 TB). The jobs will not start with less than
220 GB free, and the cohort job waits for space rather than fail. Before a
batch:

```bash
lfs quota -h -u $USER /scratch | tail -1      # keep "used" under 600 GB
```

**Code.** The working tree must be clean and on a release tag; the jobs
record the tag in their logs.

```bash
cd ~/pipelines/anneal && git status --short && git describe --tags
```

**Sample sheet.** One FASTQ name per line, derived from the batch directory:

```bash
FQ=/scratch/patkarlab-clinical/<batch>
ls $FQ/*_R1*.fastq.gz | sed 's#.*/##; s#_R1.*##' > ~/cohort_<batch>.manifest
wc -l ~/cohort_<batch>.manifest; head -3 ~/cohort_<batch>.manifest
```

Sample sheets carry accession numbers and stay outside the repository.

## 5. Executing the pipeline

**One sample:**

```bash
cd ~/pipelines/anneal
qsub -v FULL=<sample>_S<n>,OUT=/scratch/patkarlab-clinical/<outdir>,FQ=/scratch/patkarlab-clinical/<batch> jobs/anneal_e2e.pbs
```

**A batch** (array job; each subjob is one sample):

```bash
cd ~/pipelines/anneal
N=$(wc -l < ~/cohort_<batch>.manifest)
qsub -J 0-$((N-1)) -v MANIFEST=$HOME/cohort_<batch>.manifest,OUT=/scratch/patkarlab-clinical/<outdir>,FQ=/scratch/patkarlab-clinical/<batch> jobs/anneal_cohort.pbs
```

The array runs on `A40b` by default; to use the `a40` nodes as well, submit
a second copy with `-q a40` and the same arguments, and let the subjobs
skip samples the other copy has finished (a subjob exits immediately when
the sample's `scored/` table already exists). Expect two samples in
progress at a time and about 110 samples in three to four days.

**Monitoring:**

```bash
qstat -u $USER                                          # running / queued subjobs
ls /scratch/patkarlab-clinical/<outdir>/                # one directory per sample
ls /scratch/patkarlab-clinical/<outdir>/<sample>/       # alignment -> consensus -> variants -> annotated -> flt3 -> scored
tail -12 ~/pipelines/anneal/logs/<jobid>.hn1.OU         # "run_pipeline exit: 0" and "alignment/ removed" on success
```

**Failures.** A failed sample keeps its `alignment/` directory for
inspection and its log ends with a non-zero exit. Fix the cause (section
9), remove the sample's output directory, and resubmit only that index:
`qsub -J i-i ...` with the same arguments.

## 6. After the batch

**HGVS validation** through the public VariantValidator API (RefSeq
nomenclature for the clinical tables). This is part of every batch and runs
on the login node because compute nodes have no internet. Coding and
splice-site variants only; a shared cache means each unique variant is
queried once, so a batch takes minutes:

```bash
cd ~/pipelines/anneal && conda activate anneal
bash pipeline/validate_hgvs_batch.sh /scratch/patkarlab-clinical/<outdir>
```

**Reports for download.** One zip per batch in `~/inbox/to_excel/`, holding
for every sample and both tracks the report, calls, filtered, clinical and
indel tables, the FLT3 result and the consensus statistics:

```bash
bash pipeline/collect_reports.sh /scratch/patkarlab-clinical/<outdir>
```

Download the zip from `~/inbox/to_excel/` (VS Code: right-click, Download).

## 7. Results and how to read them

Layout per sample:

```
<outdir>/<sample>/
  consensus/   <sample>.dcs.sc.sorted.bam, <sample>.sscs.sc.sorted.bam, <sample>.stats.txt, family-size plot
  variants/    <sample>.<track>.vcf (Pisces), <sample>.<track>.indels.tsv (indel scan)
  annotated/   <sample>.<track>.filtered.tsv, .clinical.tsv, .indels.annotated.tsv (+ VariantValidator columns)
  flt3/        <sample>.flt3_itds.tsv
  scored/      <sample>.<track>.calls.tsv, .report.tsv, .candidates.tsv
```

**QC, first.** `consensus/<sample>.stats.txt`: expect singleton rate
55–60%, DCS recovery 15–22%, DCS reads 0.8–2.3M per 250–540M input reads,
DCS depth 3,000–17,000× depending on input. Values far outside these mean a
library or run problem, and the calls should not be read until it is
understood.

**`scored/<sample>.dcs.report.tsv`** is what to read: protein-altering
calls only (missense, nonsense, frameshift, in-frame; intronic, synonymous,
splice, UTR removed), one row per candidate, DETECTED first by VAF, with:

| Column | Meaning |
|--------|---------|
| `gene`, `consequence`, `protein_change`, `protein_id` | the annotation, split |
| `source` | `P` Pisces call, `S` indel scan, `PS` both |
| `alt_count`, `depth`, `vaf_pct` | counts from the consensus BAM |
| `bg_pct`, `p_value`, `min_alt_callable`, `lod_pct` | the site's modelled background, the beta-binomial p-value, the smallest alt count callable at this depth, and that as a VAF (the per-locus limit of detection) |
| `call` | `DETECTED`, `not_detected`, or `not_evaluable` (blocklisted, masked, no background model, or germline in the controls) |
| `tier` | qualifies DETECTED: `high_vaf` ≥ 20% (germline or a major clone), `mid_vaf` 5–20%, `mrd` < 5%, `mrd_floor` an indel at ≤ 3 reads (no background model; weakest evidence reported) |
| `note` | why, when not DETECTED; `one-sided coverage` when a site is covered by one mate |
| `HGVSc`, `HGVSg`, `COSMIC`, `ClinVar`, `gnomAD_AF`, `rsID` | joined from the annotation tables |

The full record with every candidate, including the non-coding ones, is
`calls.tsv`. Common germline homopolymer polymorphisms that the eight
controls do not carry appear as `mid_vaf` indels (the RAD21 splice-acceptor
ladder is the usual one) and are not MRD.

**SSCS** files have the same layout. SSCS is not quantitative and is
damage-prone for C>T/G>A near its limit; use it to confirm a DCS call or
to see a known marker that DCS has too few molecules for (section 8), not
to discover.

## 8. Marker tracking

For a patient with known diagnosis variants, one markers file: hg38,
1-based, VCF-anchored, no header, `chrom pos ref alt [label]`. Then, on the
login node (seconds per sample):

```bash
cd ~/pipelines/anneal && source pipeline/config.sh && activate_conda
S=<sample>; O=/scratch/patkarlab-clinical/<outdir>
python scripts/error_model/call_mrd_markers.py --sample $S --bam $O/$S/consensus/$S.dcs.sc.sorted.bam \
    --markers <patient>.markers.tsv --model "$BETA_MATRIX_DCS" --indel-blocklist "$INDEL_BLOCKLIST_DCS" \
    --indel-min-controls 6 --out $O/$S/scored/$S.dcs.markers.tsv
python scripts/error_model/call_mrd_markers.py --sample $S --bam $O/$S/consensus/$S.sscs.sc.sorted.bam \
    --markers <patient>.markers.tsv --model "$BETA_MATRIX_SSCS" --indel-blocklist "$INDEL_BLOCKLIST_SSCS" \
    --indel-min-controls 6 --out $O/$S/scored/$S.sscs.markers.tsv
```

Report the DCS line. Where DCS is negative and SSCS positive, report
"detected below duplex sensitivity" with both counts and both limits. For a
longitudinal series, score the baseline sample's candidate list on every
follow-up sample rather than each sample's own list. Marker files are
patient data and stay outside the repository.

## 9. Troubleshooting

| Symptom | Cause | Action |
|---------|-------|--------|
| Job exits 75 at the start | less than 220 GB free on `/scratch` | free space; nothing was written |
| `Aligner 'parabricks' not found at 'pbrun'` | Parabricks flag missing | `PARABRICKS_FLAGS` must be `--parabricks-docker` in `config.sh` |
| Stage 1 dies with a Lustre quota error | disk filled while aligning | free space, remove the sample's output, resubmit |
| Pisces: .NET version error | roll-forward not set | stage 2 sets `DOTNET_ROLL_FORWARD=Major`; check the `anneal` env still has `lib/dotnet` |
| VEP: `Compilation failed in require ... base.pm` | another env's Perl on PATH | VEP only runs through `run_vep()` in `annotate_variants.py`, which sanitises the environment |
| `ModuleNotFoundError: pysam` | shell in `base` | `conda activate anneal` |
| A login-node process disappears | login node time limit | run it as a job |
| HGVS validation: `Cannot reach VariantValidator` | started on a compute node, or the public service is down | run on the login node; if `curl -s https://rest.variantvalidator.org` fails, retry later, the cache resumes |
| Front-door output differs from an earlier run of the same sample | different code or config | compare `git describe` in the two logs; the pipeline is deterministic on identical input and code |

## 10. Data handling

Outputs, logs, sample sheets, marker files and `vv_cache.json` contain
patient data and are never committed. Documents use neutral sample labels.
The control FASTQs on `/scratch` are the only copy and must be preserved.
