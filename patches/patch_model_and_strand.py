#!/usr/bin/env python3
"""
patch_model_and_strand.py

Two exact-anchor patches, applied only inside the function each one targets.
Every anchor must match exactly once inside its span or the file is left
untouched. Backups go to <file>.bak.<timestamp>. Re-running is a no-op.

1. build_background.py :: fit_beta
   - Outlier control: if the control with the most alt reads carries at
     least --outlier-min-alt reads and that count is a Poisson outlier
     (p < --outlier-p) against the rate of the remaining controls, it is
     dropped for that substitution. One control with a low-level clone at a
     hotspot (IDH2 R140Q: 5 reads in one BNC, 0 in six) must not set the
     panel's detection limit at that site.
   - Method of moments only when estimable: pooled alt >= --mom-min-alt and
     at least three non-zero controls. Otherwise the --dispersion fallback.
   - Binomial sampling variance is subtracted from the between-control
     variance before the concentration is derived, so sampling noise on a
     handful of reads no longer registers as dispersion.
   - The report gains a `note` column ("outlier_dropped:k/d", "mom", or
     empty) and n_samples reflects the controls actually used.

2. call_mrd_markers.py :: score_snv
   - The strand test becomes relative: with alt >= --strand-min-alt and the
     alt fraction >= --strand-thresh one-sided, a call is rejected only if the
     alt orientation differs from the reference orientation at the same site
     (Fisher exact, two-sided, p < --strand-p, default 1e-3). Duplex reads
     are both-strand by construction, so fwd/rev is read orientation; a site
     covered by one mate only used to fail regardless of truth. If there are
     no reference reads to compare against, the absolute rule still applies.
   - Rows that are one-sided but not biased carry the note
     "one-sided coverage: ref F/R" when no other note applies.

Usage:
    python patch_model_and_strand.py --root ~/pipelines/anneal [--dry-run]
"""

import argparse
import os
import re
import shutil
import sys
import time


# --------------------------------------------------------------------------
# build_background.py
# --------------------------------------------------------------------------

BB_FIT_OLD = '''    n_tot = sum(depths)
    k_tot = sum(counts)
    if n_tot <= 0:
        return None

    mean = (k_tot + 0.5) / (n_tot + 1.0)
    m_binomial = n_tot + 1.0

    m_moment = None
    if len(counts) >= 3 and k_tot > 0:
        rates = [k / d for k, d in zip(counts, depths) if d > 0]
        if len(rates) >= 3:
            r_mean = sum(rates) / len(rates)
            r_var = sum((r - r_mean) ** 2 for r in rates) / (len(rates) - 1)
            mean_depth = n_tot / len(depths)
            var_binomial = r_mean * (1.0 - r_mean) / mean_depth
            if r_var > var_binomial > 0 and 0.0 < r_mean < 1.0:
                candidate = r_mean * (1.0 - r_mean) / r_var - 1.0
                if candidate > 1.0:
                    m_moment = candidate

    if m_moment is not None and m_moment < m_binomial:
        concentration = m_moment
    else:
        concentration = m_binomial / dispersion

    alpha = mean * concentration
    beta = (1.0 - mean) * concentration
    if alpha <= 0 or beta <= 0:
        return None
    return alpha, beta
'''

BB_FIT_NEW = '''    pairs = [(k, d) for k, d in zip(counts, depths) if d > 0]
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
            if poisson_sf(k_max, rate_rest * d_max) < outlier_p:
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
'''

BB_SIG_OLD = "def fit_beta(counts, depths, dispersion):\n"
BB_SIG_NEW = '''def poisson_sf(k, lam):
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
             outlier_min_alt=3):
'''

BB_LOOP_OLD = '''            fit = fit_beta(counts, depths, args.dispersion)
            if fit is None:
                fields.extend(["0", "0"])
                continue
            alpha, beta = fit
            fields.extend(["%.6g" % alpha, "%.6g" % beta])
            if rep:
                n_tot = sum(depths)
                k_tot = sum(counts)
                print("%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%.3e\\t%.6g\\t%.6g\\t%.3e"
                      % (chrom, pos, ref_base, len(records), n_tot,
                         alt, k_tot, k_tot / n_tot if n_tot else 0.0,
                         alpha, beta, alpha / (alpha + beta)), file=rep)
'''

BB_LOOP_NEW = '''            fit = fit_beta(counts, depths, args.dispersion,
                           args.mom_min_alt, args.outlier_p, args.outlier_min_alt)
            if fit is None:
                fields.extend(["0", "0"])
                continue
            alpha, beta, n_used, n_tot, k_tot, note = fit
            fields.extend(["%.6g" % alpha, "%.6g" % beta])
            if rep:
                print("%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%.3e\\t%.6g\\t%.6g\\t%.3e\\t%s"
                      % (chrom, pos, ref_base, n_used, n_tot,
                         alt, k_tot, k_tot / n_tot if n_tot else 0.0,
                         alpha, beta, alpha / (alpha + beta), note), file=rep)
'''

BB_ARG_ANCHOR = '    ap.add_argument("--dispersion", type=float, default=3.0,\n'
BB_ARG_NEW = '''    ap.add_argument("--mom-min-alt", type=int, default=20,
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
'''


# --------------------------------------------------------------------------
# call_mrd_markers.py
# --------------------------------------------------------------------------

CM_FUNCS = '''def fisher_two_sided(a, b, c, d):
    """Two-sided Fisher exact p for the 2x2 table [[a, b], [c, d]]."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    r1, c1 = a + b, a + c
    lg = math.lgamma

    def logp(x):
        return (lg(r1 + 1) - lg(x + 1) - lg(r1 - x + 1)
                + lg(n - r1 + 1) - lg(c1 - x + 1) - lg(n - r1 - c1 + x + 1)
                - (lg(n + 1) - lg(c1 + 1) - lg(n - c1 + 1)))

    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    p_obs = logp(a)
    total = 0.0
    for x in range(lo, hi + 1):
        lp = logp(x)
        if lp <= p_obs + 1e-9:
            total += math.exp(lp)
    return min(1.0, total)


def orientation_bias(alt_fwd, alt_rev, ref_fwd, ref_rev, p_thresh):
    """True if alt read orientation differs from reference orientation.

    Duplex consensus reads are both-strand by construction, so fwd/rev is
    the alignment orientation of the consensus read, not strand of origin.
    A site covered by one mate only is one-sided for reference reads too;
    only a departure from the reference orientation indicates a mapping
    artifact. With no reference reads the caller's absolute rule applies.
    """
    if ref_fwd + ref_rev == 0:
        return True
    return fisher_two_sided(alt_fwd, alt_rev, ref_fwd, ref_rev) < p_thresh


'''

CM_HEAD_OLD = '''    counts, strand = count_at(bam, chrom, pos - 1, args.min_bq)
    depth = sum(counts.values())
    ac = counts.get(alt, 0)
    fwd, rev = strand.get(alt, [0, 0])
    vaf = ac / depth if depth else 0.0
    sfrac = (max(fwd, rev) / ac) if ac else 0.0
'''

CM_HEAD_NEW = '''    counts, strand = count_at(bam, chrom, pos - 1, args.min_bq)
    depth = sum(counts.values())
    ac = counts.get(alt, 0)
    fwd, rev = strand.get(alt, [0, 0])
    ref_fwd, ref_rev = strand.get(m.get("ref", ""), [0, 0])
    vaf = ac / depth if depth else 0.0
    sfrac = (max(fwd, rev) / ac) if ac else 0.0
    one_sided = ac >= args.strand_min_alt and sfrac >= args.strand_thresh
    biased = one_sided and orientation_bias(fwd, rev, ref_fwd, ref_rev,
                                            args.strand_p)
'''

CM_ELIF_OLD = '''    elif ac >= args.strand_min_alt and sfrac >= args.strand_thresh:
        row["call"], row["note"] = "not_detected", "strand bias"
'''

CM_ELIF_NEW = '''    elif biased:
        row["call"], row["note"] = "not_detected", "strand bias"
'''

CM_TAIL_OLD = '''    else:
        row["call"], row["note"] = "not_detected", "below background"
    return row
'''

CM_TAIL_NEW = '''    else:
        row["call"], row["note"] = "not_detected", "below background"
    if one_sided and not biased and not row.get("note"):
        row["note"] = "one-sided coverage: ref %d/%d" % (ref_fwd, ref_rev)
    return row
'''

CM_ARG_ANCHOR = '    ap.add_argument("--strand-min-alt", type=int, default=10)\n'
CM_ARG_NEW = ('    ap.add_argument("--strand-min-alt", type=int, default=10)\n'
              '    ap.add_argument("--strand-p", type=float, default=1e-3,\n'
              '                    help="Fisher p below which alt orientation is "\n'
              '                         "called biased against reference orientation "\n'
              '                         "(default 1e-3)")\n')


# --------------------------------------------------------------------------
# machinery
# --------------------------------------------------------------------------

def span(text, start_pat, end_pat):
    """Return (i, j) delimiting the function that starts at start_pat."""
    i = text.find(start_pat)
    if i < 0:
        sys.exit("ERROR: anchor not found: %r" % start_pat[:50])
    j = text.find(end_pat, i + len(start_pat)) if end_pat else len(text)
    if j < 0:
        j = len(text)
    return i, j


def replace_in_span(text, i, j, old, new, label):
    seg = text[i:j]
    n = seg.count(old)
    if n != 1:
        sys.exit("ERROR: %s: anchor matched %d times inside its span, expected 1"
                 % (label, n))
    return text[:i] + seg.replace(old, new) + text[j:]


def ensure_import_math(text):
    if re.search(r"^import math$", text, re.M):
        return text
    m = re.search(r"^import \w+", text, re.M)
    if not m:
        sys.exit("ERROR: no import block found")
    return text[:m.start()] + "import math\n" + text[m.start():]


def patch_build_background(text):
    if "def poisson_sf(" in text:
        print("  build_background.py already patched, skipping")
        return text
    text = ensure_import_math(text)
    i, j = span(text, "def fit_beta(", "\ndef main(")
    text = replace_in_span(text, i, j, BB_FIT_OLD, BB_FIT_NEW, "fit_beta body")
    text = replace_in_span(text, i, j + len(BB_FIT_NEW) - len(BB_FIT_OLD),
                           BB_SIG_OLD, BB_SIG_NEW, "fit_beta signature")
    i, j = span(text, "\ndef main(", None)
    text = replace_in_span(text, i, j, BB_LOOP_OLD, BB_LOOP_NEW, "site loop")
    i, j = span(text, "\ndef main(", None)
    text = replace_in_span(text, i, j, BB_ARG_ANCHOR, BB_ARG_NEW + BB_ARG_ANCHOR,
                           "argparse")
    # report header: tolerate either a literal string or a list of names
    if 'model_rate\\tnote' not in text and '"model_rate", "note"' not in text:
        if 'model_rate"' in text:
            text = text.replace('model_rate"', 'model_rate\\tnote"', 1)
        elif '"model_rate"]' in text:
            text = text.replace('"model_rate"]', '"model_rate", "note"]', 1)
        else:
            print("  WARNING: report header anchor not found; add a 'note' "
                  "column to the header by hand")
    return text


def patch_call_mrd_markers(text):
    if "def orientation_bias(" in text:
        print("  call_mrd_markers.py already patched, skipping")
        return text
    text = ensure_import_math(text)
    k = text.find("def score_snv(")
    if k < 0:
        sys.exit("ERROR: def score_snv not found")
    text = text[:k] + CM_FUNCS + text[k:]
    i, j = span(text, "def score_snv(", "\ndef score_indel(")
    text = replace_in_span(text, i, j, CM_HEAD_OLD, CM_HEAD_NEW, "score_snv head")
    i, j = span(text, "def score_snv(", "\ndef score_indel(")
    text = replace_in_span(text, i, j, CM_ELIF_OLD, CM_ELIF_NEW, "strand elif")
    i, j = span(text, "def score_snv(", "\ndef score_indel(")
    text = replace_in_span(text, i, j, CM_TAIL_OLD, CM_TAIL_NEW, "score_snv tail")
    i, j = span(text, "\ndef main(", None)
    text = replace_in_span(text, i, j, CM_ARG_ANCHOR, CM_ARG_NEW, "argparse")
    return text


def apply(path, fn, dry_run):
    print("== %s" % path)
    if not os.path.isfile(path):
        sys.exit("ERROR: not found: %s" % path)
    original = open(path).read()
    patched = fn(original)
    if patched == original:
        print("  no change")
        return
    compile(patched, path, "exec")
    if dry_run:
        print("  dry run: would write %d -> %d bytes" % (len(original), len(patched)))
        return
    backup = "%s.bak.%s" % (path, time.strftime("%Y%m%d-%H%M%S"))
    shutil.copy2(path, backup)
    with open(path, "w") as fh:
        fh.write(patched)
    print("  written; backup at %s" % backup)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="anneal repo root")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    apply(os.path.join(args.root, "build_background.py"),
          patch_build_background, args.dry_run)
    apply(os.path.join(args.root, "scripts", "error_model", "call_mrd_markers.py"),
          patch_call_mrd_markers, args.dry_run)


if __name__ == "__main__":
    main()
