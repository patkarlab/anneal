#!/usr/bin/env python3
"""
snv_specificity.py  (Step 2 of 3)

Measure the false-positive rate of the substitution arm on the eight
marker-negative normals, using the per-site base-count tables already produced
by sscs_site_counts.py. For every panel site and substitution in each normal,
the normal's own alt/depth is scored against the site's error Beta with both
tests; any positive call in a marker-negative sample is a false positive.

This reuses existing counts rather than re-fetching from the BAMs, so it is
fast. Two caveats it does not cover, both of which make the count it reports an
UPPER bound: the strand filter is not applied here (the count tables carry no
strand), so any call found should be strand-checked with inspect_locus_reads;
and the model was built from these same normals, so a clean leave-one-out test
would rebuild the model without each normal before scoring it. Both refinements
only reduce the false-positive count, so a clean result here is conclusive.

Output: a summary to stderr (calls per normal) and a TSV of every called site
for inspection. Calls are most often real germline/CHIP variants rather than
sequencing artefacts; the TSV lets you tell them apart by gene and recurrence.
"""

import argparse
import math
import os
import sys

BASES = ("A", "C", "G", "T")


# --- statistics (identical to the apply step; validated against scipy) ---

def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 300, 1e-16, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        if abs(d * c - 1.0) < EPS:
            break
    return h


def beta_cdf(x, a, b):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    bt = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _logB(x, y):
    return math.lgamma(x) + math.lgamma(y) - math.lgamma(x + y)


def betabinom_sf(k, n, a, b):
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    logBab = _logB(a, b)
    lower = 0.0
    for j in range(0, k):
        logpmf = (math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
                  + _logB(j + a, n - j + b) - logBab)
        lower += math.exp(logpmf)
    return max(0.0, 1.0 - lower)


def load_model(path):
    sites = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            sites[(f[0], int(f[1]), f[2], f[3])] = (float(f[10]), float(f[11]))
    return sites


def load_mask(path):
    masked = set()
    if not path:
        return masked
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            masked.add((f[0], int(f[1]) + 1))
    return masked


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True)
    ap.add_argument("--counts", nargs="+", required=True,
                    help="Per-normal count TSVs from sscs_site_counts.py")
    ap.add_argument("--out", required=True, help="TSV of called sites")
    ap.add_argument("--mask")
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--alpha-level", type=float, default=0.005)
    ap.add_argument("--test", choices=["pbeta", "betabinom", "both"],
                    default="betabinom")
    args = ap.parse_args()

    sites = load_model(args.model)
    masked = load_mask(args.mask)

    per_sample = {}
    with open(args.out, "w") as out:
        out.write("#sample\tchrom\tpos\tref\talt\talt_count\tdepth\tvaf\t"
                  "pbeta_p\tbetabinom_p\n")
        for path in args.counts:
            sample = os.path.basename(path).split(".")[0]
            calls = 0
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    chrom, pos, ref, depth, a, c, g, t = line.rstrip("\n").split("\t")
                    pos = int(pos)
                    depth = int(depth)
                    if depth < args.min_depth or ref not in BASES:
                        continue
                    if (chrom, pos) in masked:
                        continue
                    cnt = {"A": int(a), "C": int(c), "G": int(g), "T": int(t)}
                    for alt in BASES:
                        if alt == ref:
                            continue
                        key = (chrom, pos, ref, alt)
                        if key not in sites:
                            continue
                        alpha, beta = sites[key]
                        ac = cnt[alt]
                        if ac == 0:
                            continue
                        vaf = ac / depth
                        pb = 1.0 - beta_cdf(vaf, alpha, beta)
                        bb = betabinom_sf(ac, depth, alpha, beta)
                        if args.test == "pbeta":
                            passed = pb <= args.alpha_level
                        elif args.test == "betabinom":
                            passed = bb <= args.alpha_level
                        else:
                            passed = pb <= args.alpha_level and bb <= args.alpha_level
                        if passed:
                            calls += 1
                            out.write(f"{sample}\t{chrom}\t{pos}\t{ref}\t{alt}\t"
                                      f"{ac}\t{depth}\t{vaf:.3e}\t{pb:.3e}\t"
                                      f"{bb:.3e}\n")
            per_sample[sample] = calls

    total = sum(per_sample.values())
    sys.stderr.write(f"Substitution-arm false-positive calls per normal "
                     f"(test={args.test}, alpha={args.alpha_level}):\n")
    for s in sorted(per_sample):
        sys.stderr.write(f"  {s}: {per_sample[s]}\n")
    sys.stderr.write(f"  TOTAL: {total} (listed in {args.out} for inspection)\n")


if __name__ == "__main__":
    main()
