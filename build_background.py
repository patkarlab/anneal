#!/usr/bin/env python3
"""
build_background.py

Builds the per-site, per-substitution background model directly from the
biological negative control consensus BAMs.

WHY THIS REPLACES THE EXISTING CHAIN
------------------------------------
The current model is produced by: Pisces -> remove_variants_gtr_20.pl ->
print_multiple_variants_at_same_location.pl -> fill_empty_mips.pl ->
beta_distribution.py. That chain has a defect its own README documents:

    "Pisces applies MinimumVariantQScore (default 20) beneath the requested
     --minvf, censoring roughly everything below 7-12 alt reads. Positions
     below that are recorded as zeros and take the default."

So at nearly every site the BNCs contribute nothing and the model falls back
to default_error_rate, a literature constant of 1/200000 for DCS. The measured
DCS floor on this panel is 6.6e-05, which is 13x higher. At 3000x depth the
old model expects 0.015 background reads where reality is ~0.2, so a single
alt read scores p ~ 0.015 (apparently significant) when it is in fact routine.
That is a false-positive generator at every site where the default applies.

This script counts from the pileup instead, so no site is censored. Zero-count
sites are anchored by a Jeffreys prior on the pooled depth rather than by a
constant from another panel: with eight controls at ~3000x, a site with no
observed alt gets ~0.5/24000 = 2.1e-05 per base, which is the right order of
magnitude for the measured floor rather than 13x below it.

OUTPUT
------
Tab-separated, one row per panel position, matching the format that
call_mrd_markers.py::load_beta_matrix expects:

    chr  pos  alpha_A  beta_A  alpha_T  beta_T  alpha_G  beta_G  alpha_C  beta_C

The reference base at each position is written with alpha=beta=0 so that
load_beta_matrix skips it (it requires both > 0).

A companion report with the raw per-sample counts is written alongside, for
inspection and for deriving the artifact mask.

USAGE
-----
    python3 build_background.py \\
        --bams /scratch/.../AMLMRD_DUPLEX_BNC*/consensus/dcs.sc.sorted.bam \\
        --ref ~/references/hg38_broad/Homo_sapiens_assembly38.masked.fasta \\
        --bed AML_MRD_DUPLEX_probes_hg38_sortd.bed \\
        --out beta_matrix_DCS.txt \\
        --report beta_matrix_DCS.report.tsv
"""

import math
import argparse
import os
import subprocess
import sys
from collections import defaultdict

BASES = ("A", "T", "G", "C")          # column order required by the consumer

# A per-sample allele fraction at or above this is treated as germline (or a
# real variant in that control) and that sample is dropped for that position
# and alt base only. Panel data shows a clean gap between background (<0.5%)
# and germline (>40%), so the exact value is not load-bearing.
GERMLINE_AF = 0.20

# Samples contributing less than this depth at a position are ignored there.
MIN_SAMPLE_DEPTH = 100

# A position must have this many usable control samples to be modelled.
# Below it, the position is emitted with zeros and becomes not_evaluable
# downstream, which is safer than modelling it from one or two controls.
MIN_SAMPLES = 4


def parse_pileup_column(pileup):
    """Return (ref_match_count, {alt_base: count}) for one mpileup column.

    Indel and read-start markers must be stepped over so inserted or deleted
    sequence is not miscounted as substitutions.
    """
    alts = defaultdict(int)
    ref_matches = 0
    i = 0
    n = len(pileup)
    while i < n:
        c = pileup[i]
        if c == "^":
            i += 2                      # read start: skip the mapping quality
        elif c == "$":
            i += 1
        elif c in ".,":
            ref_matches += 1
            i += 1
        elif c in "+-":
            j = i + 1
            digits = ""
            while j < n and pileup[j].isdigit():
                digits += pileup[j]
                j += 1
            i = j + (int(digits) if digits else 0)
        elif c in "ACGTacgt":
            alts[c.upper()] += 1
            i += 1
        else:
            i += 1
    return ref_matches, alts


def pileup_sample(bam, ref, bed):
    """Yield (chrom, pos, ref_base, depth, {alt: count}) for one BAM."""
    cmd = ["samtools", "mpileup", "-A", "-Q", "30", "-d", "1000000",
           "-l", bed, "-f", ref, bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        ref_base = f[2].upper()
        if ref_base not in BASES:
            continue
        ref_matches, alts = parse_pileup_column(f[4])
        depth = ref_matches + sum(alts.values())
        yield f[0], int(f[1]), ref_base, depth, alts
    proc.wait()


def poisson_sf(k, lam):
    """P(X >= k) for X ~ Poisson(lam), without scipy."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    log_term = -lam
    cdf = math.exp(log_term)
    for i in range(1, k):
        log_term += math.log(lam) - math.log(i)
        cdf += math.exp(log_term)
    return max(0.0, 1.0 - cdf)


def fit_beta(counts, depths, dispersion, mom_min_alt=20, outlier_p=1e-3,
             outlier_min_alt=3, outlier_min_ratio=10.0):
    """Fit Beta parameters for one (position, alt base) across controls.

    counts, depths: per-sample alt counts and depths, germline already removed.

    The mean is the pooled rate under a Jeffreys prior, so a site with no
    observed alt is anchored at ~0.5/N rather than at an external constant.

    The concentration sets how much between-sample variability the model
    tolerates. If the controls are more dispersed than binomial sampling
    alone would explain, that dispersion is used (method of moments).
    Otherwise the binomial concentration is divided by --dispersion, which
    widens the distribution and makes the test conservative.
    """
    pairs = [(k, d) for k, d in zip(counts, depths) if d > 0]
    if not pairs:
        return None
    note = ""

    # Outlier control. The control with the most alt reads is tested against
    # the pooled rate of the others; a Poisson outlier is dropped for this
    # substitution so that one control's clone does not set the site's limit.
    if len(pairs) >= 3:
        k_max, d_max = max(pairs, key=lambda kd: kd[0])
        if k_max >= outlier_min_alt:
            k_rest = sum(k for k, _ in pairs) - k_max
            d_rest = sum(d for _, d in pairs) - d_max
            rate_rest = (k_rest + 0.5) / (d_rest + 1.0)
            # Both conditions: Poisson-improbable given the others, and a
            # large excess over them. The ratio keeps systematic artifact
            # sites, where every control is high and one is merely highest,
            # with the moment estimator instead of excising a control.
            if (poisson_sf(k_max, rate_rest * d_max) < outlier_p
                    and k_max / d_max >= outlier_min_ratio * rate_rest):
                pairs.remove((k_max, d_max))
                note = "outlier_dropped:%d/%d" % (k_max, d_max)

    n_tot = sum(d for _, d in pairs)
    k_tot = sum(k for k, _ in pairs)
    mean = (k_tot + 0.5) / (n_tot + 1.0)
    m_binomial = n_tot + 1.0
    concentration = m_binomial / dispersion

    # Method of moments, only when estimable. Binomial sampling variance is
    # subtracted from the between-control variance first.
    nonzero = sum(1 for k, _ in pairs if k > 0)
    if len(pairs) >= 3 and k_tot >= mom_min_alt and nonzero >= 3:
        rates = [k / d for k, d in pairs]
        r_mean = sum(rates) / len(rates)
        r_var = sum((r - r_mean) ** 2 for r in rates) / (len(rates) - 1)
        mean_depth = n_tot / len(pairs)
        var_binomial = r_mean * (1.0 - r_mean) / mean_depth
        excess = r_var - var_binomial
        if excess > 0 and 0.0 < r_mean < 1.0:
            candidate = r_mean * (1.0 - r_mean) / excess - 1.0
            if 1.0 < candidate < m_binomial:
                concentration = candidate
                note = (note + ";" if note else "") + "mom"

    alpha = mean * concentration
    beta = (1.0 - mean) * concentration
    if alpha <= 0 or beta <= 0:
        return None
    return alpha, beta, len(pairs), n_tot, k_tot, note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", nargs="+", required=True,
                    help="Control consensus BAMs (all DCS, or all SSCS)")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True,
                    help="Output beta matrix, e.g. beta_matrix_DCS.txt")
    ap.add_argument("--report",
                    help="Optional per-site raw count report")
    ap.add_argument("--mom-min-alt", type=int, default=20,
                    help="Pooled alt reads required before between-control "
                         "dispersion is estimated by method of moments; below "
                         "this the --dispersion fallback is used (default 20)")
    ap.add_argument("--outlier-p", type=float, default=1e-3,
                    help="Poisson tail probability below which the highest-"
                         "count control is dropped for a substitution "
                         "(default 1e-3)")
    ap.add_argument("--outlier-min-alt", type=int, default=3,
                    help="Minimum alt reads in a control before it can be "
                         "dropped as an outlier (default 3)")
    ap.add_argument("--outlier-min-ratio", type=float, default=10.0,
                    help="The dropped control's rate must exceed the pooled "
                         "rate of the others by this factor; keeps systematic "
                         "artifact sites intact for the moment estimator "
                         "(default 10)")
    ap.add_argument("--dispersion", type=float, default=3.0,
                    help="Variance inflation over binomial when the controls "
                         "give no dispersion estimate (default 3.0). Higher "
                         "is more conservative.")
    ap.add_argument("--germline-af", type=float, default=GERMLINE_AF)
    ap.add_argument("--min-sample-depth", type=int, default=MIN_SAMPLE_DEPTH)
    ap.add_argument("--min-samples", type=int, default=MIN_SAMPLES)
    args = ap.parse_args()

    for path in [args.ref, args.bed] + args.bams:
        if not os.path.isfile(path):
            sys.exit("ERROR: not found: %s" % path)
    if len(args.bams) < args.min_samples:
        sys.exit("ERROR: %d BAMs given but --min-samples is %d"
                 % (len(args.bams), args.min_samples))

    # site -> ref base
    ref_at = {}
    # (site, alt) -> list of (count, depth) across samples
    obs = defaultdict(list)
    # site -> list of depths, for the report
    depth_at = defaultdict(list)

    for idx, bam in enumerate(args.bams, 1):
        print("[%d/%d] %s" % (idx, len(args.bams), os.path.basename(bam)),
              file=sys.stderr)
        seen = 0
        for chrom, pos, ref_base, depth, alts in pileup_sample(
                bam, args.ref, args.bed):
            if depth < args.min_sample_depth:
                continue
            site = (chrom, pos)
            ref_at[site] = ref_base
            depth_at[site].append(depth)
            for alt in BASES:
                if alt == ref_base:
                    continue
                k = alts.get(alt, 0)
                if k / depth >= args.germline_af:
                    continue            # germline or real variant in control
                obs[(site, alt)].append((k, depth))
            seen += 1
        print("        %d positions" % seen, file=sys.stderr)

    if not ref_at:
        sys.exit("ERROR: no positions passed the depth filter")

    out = open(args.out, "w")
    print("chr\tpos\talpha A\tbeta A\talpha T\tbeta T"
          "\talpha G\tbeta G\talpha C\tbeta C", file=out)

    rep = None
    if args.report:
        rep = open(args.report, "w")
        print("chr\tpos\tref\tn_samples\tpooled_depth\t"
              "alt\talt_count\tpooled_rate\talpha\tbeta\tmodel_rate\tnote",
              file=rep)

    n_modelled = n_skipped = 0
    for site in sorted(ref_at, key=lambda s: (s[0], s[1])):
        chrom, pos = site
        ref_base = ref_at[site]
        fields = []
        for alt in BASES:
            if alt == ref_base:
                fields.extend(["0", "0"])   # consumer skips non-positive
                continue
            records = obs.get((site, alt), [])
            if len(records) < args.min_samples:
                fields.extend(["0", "0"])
                continue
            counts = [k for k, _ in records]
            depths = [d for _, d in records]
            fit = fit_beta(counts, depths, args.dispersion,
                           args.mom_min_alt, args.outlier_p, args.outlier_min_alt,
                           args.outlier_min_ratio)
            if fit is None:
                fields.extend(["0", "0"])
                continue
            alpha, beta, n_used, n_tot, k_tot, note = fit
            fields.extend(["%.6g" % alpha, "%.6g" % beta])
            if rep:
                print("%s\t%d\t%s\t%d\t%d\t%s\t%d\t%.3e\t%.6g\t%.6g\t%.3e\t%s"
                      % (chrom, pos, ref_base, n_used, n_tot,
                         alt, k_tot, k_tot / n_tot if n_tot else 0.0,
                         alpha, beta, alpha / (alpha + beta), note), file=rep)
        if any(f != "0" for f in fields):
            n_modelled += 1
        else:
            n_skipped += 1
        print("%s\t%d\t%s" % (chrom, pos, "\t".join(fields)), file=out)

    out.close()
    if rep:
        rep.close()

    print("", file=sys.stderr)
    print("positions modelled: %d" % n_modelled, file=sys.stderr)
    print("positions skipped:  %d (fewer than %d usable controls)"
          % (n_skipped, args.min_samples), file=sys.stderr)
    print("written: %s" % args.out, file=sys.stderr)
    if args.report:
        print("report:  %s" % args.report, file=sys.stderr)


if __name__ == "__main__":
    main()
