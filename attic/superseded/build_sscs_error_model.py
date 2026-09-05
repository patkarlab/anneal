#!/usr/bin/env python3
"""
build_sscs_error_model.py

Stage 3 of the SSCS error model. Aggregates the per-BNC base-count tables
produced by sscs_site_counts.py and fits a per-site, per-substitution
background error model across the biological negative controls (BNCs).

Model
-----
For every panel position and every substitution (ref -> each other base):

  1. Per BNC, the observed error rate is alt_reads / depth. Any BNC whose
     rate exceeds --vaf-exclude (default 0.20) is dropped at that site as a
     germline / real variant rather than background error. This mirrors the
     reference model and protects the pool from a single germline carrier.

  2. Surviving BNCs are pooled (summed alt and summed depth).

  3. The per-site rate is shrunk toward a substitution-class prior by
     empirical Bayes:

        posterior  Beta(alpha, beta)
        alpha = alpha0 + pooled_alt
        beta  = beta0  + (pooled_depth - pooled_alt)

     where (alpha0, beta0) is a Beta fit (method of moments) to the
     distribution of per-site rates within that directed substitution class.
     This replaces the reference model's fixed 1/15000 floor: with only eight
     controls most sites carry zero observed alt reads, and the class prior
     supplies a data-driven, substitution-specific floor of appropriate
     strength. Sites with real systematic noise are dominated by their counts.

  4. The reported background is (posterior_mean + 3 * posterior_sd) * 100,
     keeping the reference model's reporting convention.

The posterior (alpha, beta) is what the apply step uses for the Beta-Binomial
test on patient SSCS calls. Per-class priors are written as ## comment lines
at the top of the output for reproducibility.

Diagnostic columns (n_bnc_used, max_bnc_vaf) are included so sites driven by a
single outlier control (e.g. intermediate-VAF clonal haematopoiesis that slips
under the 0.20 cutoff) are visible for review rather than silently fit.
"""

import argparse
import math
import sys
from collections import defaultdict

BASES = ("A", "C", "G", "T")


def read_counts(paths):
    """Read all per-BNC count tables produced by sscs_site_counts.py.

    Returns a dict keyed by (chrom, pos) -> {"ref": ref, "bnc": [entry, ...]}
    where each entry (one per input file, None if the position was absent) is
    {"A":.., "C":.., "G":.., "T":.., "depth":..}.
    """
    sites = {}
    n_files = len(paths)
    for fi, path in enumerate(paths):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                chrom, pos, ref, depth, a, c, g, t = line.rstrip("\n").split("\t")
                key = (chrom, int(pos))
                rec = sites.get(key)
                if rec is None:
                    rec = {"ref": ref, "bnc": [None] * n_files}
                    sites[key] = rec
                elif rec["ref"] != ref:
                    sys.stderr.write(
                        f"WARN: reference mismatch at {chrom}:{pos} "
                        f"({rec['ref']} vs {ref}); keeping first\n")
                rec["bnc"][fi] = {"A": int(a), "C": int(c), "G": int(g),
                                  "T": int(t), "depth": int(depth)}
    return sites


def moment_beta(rates):
    """Method-of-moments Beta(alpha, beta) from a list of rates in (0, 1).
    Returns None if a valid strictly-positive fit is not possible."""
    n = len(rates)
    if n < 2:
        return None
    m = sum(rates) / n
    if m <= 0 or m >= 1:
        return None
    v = sum((r - m) ** 2 for r in rates) / (n - 1)
    if v <= 0:
        return None
    s = m * (1.0 - m) / v - 1.0
    if s <= 0:
        return None
    return m * s, (1.0 - m) * s


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", nargs="+", required=True,
                    help="Per-BNC count TSVs from sscs_site_counts.py")
    ap.add_argument("--out", required=True, help="Output model TSV")
    ap.add_argument("--vaf-exclude", type=float, default=0.20,
                    help="Drop a BNC at a site if its rate exceeds this value "
                         "(germline/real, not error). Default 0.20")
    ap.add_argument("--prior-min-depth", type=int, default=1000,
                    help="Minimum pooled depth for a site to contribute to its "
                         "class prior fit. Default 1000")
    ap.add_argument("--prior-max-rate", type=float, default=0.01,
                    help="Exclude a site from its class prior fit if its pooled "
                         "rate exceeds this. Keeps systematic-artifact loci "
                         "(paralog/mapping) from inflating the genuine error "
                         "floor. Default 0.01")
    ap.add_argument("--fallback-kappa", type=float, default=1000.0,
                    help="Prior strength in pseudo-reads when a class prior "
                         "cannot be fit by moments. Default 1000")
    args = ap.parse_args()

    sites = read_counts(args.counts)
    sys.stderr.write(f"Loaded {len(sites)} positions from "
                     f"{len(args.counts)} BNC tables\n")

    # ---- Pass 1: per-site pooled counts, with per-BNC VAF exclusion ----
    records = []
    class_rates = defaultdict(list)          # sub -> [per-site rate, ...]
    class_totals = defaultdict(lambda: [0, 0])  # sub -> [alt_total, depth_total]

    for (chrom, pos), rec in sites.items():
        ref = rec["ref"]
        if ref not in BASES:
            continue
        for alt in BASES:
            if alt == ref:
                continue
            pooled_alt = 0
            pooled_depth = 0
            n_used = 0
            max_vaf = 0.0
            for entry in rec["bnc"]:
                if entry is None or entry["depth"] == 0:
                    continue
                a = entry[alt]
                d = entry["depth"]
                vaf = a / d
                if vaf > args.vaf_exclude:
                    continue  # germline/real in this control; drop it here
                pooled_alt += a
                pooled_depth += d
                n_used += 1
                if vaf > max_vaf:
                    max_vaf = vaf
            if pooled_depth == 0:
                continue
            sub = f"{ref}>{alt}"
            site_rate = pooled_alt / pooled_depth
            records.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                            "sub": sub, "n_used": n_used,
                            "pooled_alt": pooled_alt,
                            "pooled_depth": pooled_depth,
                            "site_rate": site_rate, "max_vaf": max_vaf})
            class_totals[sub][0] += pooled_alt
            class_totals[sub][1] += pooled_depth
            if (pooled_depth >= args.prior_min_depth
                    and site_rate <= args.prior_max_rate):
                class_rates[sub].append(site_rate)

    # ---- Fit a prior per substitution class ----
    priors = {}
    for sub in sorted(class_totals):
        fit = moment_beta(class_rates.get(sub, []))
        if fit is not None:
            a0, b0 = fit
            source = "moments"
        else:
            alt_tot, dep_tot = class_totals[sub]
            mean = (alt_tot + 0.5) / (dep_tot + 1.0)   # Jeffreys class rate
            k = args.fallback_kappa
            a0, b0 = mean * k, (1.0 - mean) * k
            source = "fallback"
        priors[sub] = (a0, b0, source)

    # ---- Pass 2: posterior per site; write output ----
    cols = ["#chrom", "pos", "ref", "alt", "sub", "n_bnc_used",
            "pooled_alt", "pooled_depth", "site_rate", "max_bnc_vaf",
            "alpha", "beta", "post_mean", "post_sd", "background_pct"]
    with open(args.out, "w") as out:
        for sub in sorted(priors):
            a0, b0, source = priors[sub]
            out.write(f"## prior {sub} alpha0={a0:.6g} beta0={b0:.6g} "
                      f"mean={a0 / (a0 + b0):.3e} source={source}\n")
        out.write("\t".join(cols) + "\n")
        for r in records:
            a0, b0, _ = priors[r["sub"]]
            alpha = a0 + r["pooled_alt"]
            beta = b0 + (r["pooled_depth"] - r["pooled_alt"])
            tot = alpha + beta
            post_mean = alpha / tot
            post_var = (alpha * beta) / (tot * tot * (tot + 1.0))
            post_sd = math.sqrt(post_var)
            background_pct = (post_mean + 3.0 * post_sd) * 100.0
            out.write("\t".join(str(x) for x in [
                r["chrom"], r["pos"], r["ref"], r["alt"], r["sub"],
                r["n_used"], r["pooled_alt"], r["pooled_depth"],
                f"{r['site_rate']:.3e}", f"{r['max_vaf']:.3e}",
                f"{alpha:.6g}", f"{beta:.6g}", f"{post_mean:.3e}",
                f"{post_sd:.3e}", f"{background_pct:.4f}"]) + "\n")

    sys.stderr.write(f"Wrote model for {len(records)} site-substitutions "
                     f"to {args.out}\n")


if __name__ == "__main__":
    main()
