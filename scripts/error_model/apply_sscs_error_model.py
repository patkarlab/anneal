#!/usr/bin/env python3
"""
apply_sscs_error_model.py

Score candidate variants in a patient SSCS sample against the per-site error
model, applying all four layers consistently with how the model was built:

  1. NM <= 10 filter on the patient BAM (calmd + samtools view), so paralog and
     noisy reads are removed before counting, exactly as for the BNC panel.
  2. Per-candidate counts recomputed from the FILTERED BAM (the VCF's ALT/TOT is
     not trusted, since it was produced before filtering).
  3. Strand filter: a candidate is rejected if its alt reads are >= --strand-thresh
     one strand (paralog signature), applied only when there are enough alt reads
     to judge strand (--strand-min-alt).
  4. Mask: candidates at masked positions are reported as no-call.

Two statistical tests are computed for every candidate so they can be compared:
  - pbeta (Waalkes): call if 1 - BetaCDF(VAF, alpha, beta) <= --alpha-level.
    Tests the observed VAF against the site error-rate Beta; ignores depth.
  - betabinom: call if P(X >= alt | depth, alpha, beta) <= --alpha-level.
    Depth-aware; the appropriate test at very low VAF.

Inputs: patient SSCS BAM (unfiltered), the model TSV from
build_sscs_error_model.py, the mask BED from derive_artifact_mask.py, and a
--sites TSV of candidates (chrom, pos, ref, alt; 1-based). The candidate list is
typically the patient's diagnosis variants for MRD tracking.
"""

import argparse
import math
import os
import subprocess
import sys
import pysam

BASES = ("A", "C", "G", "T")


# ---------- statistics (validated against scipy to ~1e-11) ----------

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
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS:
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
    """P(X >= k) for X ~ BetaBinomial(n, a, b)."""
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


# ---------- IO ----------

def filter_patient_bam(in_bam, ref, out_bam, max_nm):
    """calmd to populate NM against ref, then keep reads with NM <= max_nm.
    Mirrors the BNC filtering exactly."""
    calmd_tmp = out_bam + ".calmd.tmp.bam"
    with open(calmd_tmp, "wb") as fh:
        subprocess.run(["samtools", "calmd", "-b", in_bam, ref],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)
    subprocess.run(["samtools", "view", "-b", "-e", f"[NM]<={max_nm}",
                    "-o", out_bam, calmd_tmp], check=True)
    subprocess.run(["samtools", "index", out_bam], check=True)
    os.remove(calmd_tmp)


def load_model(path):
    """Return (sites, priors). sites[(chrom,pos,ref,alt)] = (alpha,beta,bg_pct);
    priors[sub] = (alpha0,beta0) from the ## prior header lines."""
    sites, priors = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("## prior"):
                f = line.split()
                sub = f[2]
                a0 = float(f[3].split("=")[1])
                b0 = float(f[4].split("=")[1])
                priors[sub] = (a0, b0)
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            key = (f[0], int(f[1]), f[2], f[3])
            sites[key] = (float(f[10]), float(f[11]), float(f[14]))
    return sites, priors


def load_mask(path):
    masked = set()
    if not path:
        return masked
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            masked.add((f[0], int(f[1]) + 1))   # BED 0-based start -> 1-based pos
    return masked


def read_sites(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            yield f[0], int(f[1]), f[2].upper(), f[3].upper()


def count_at(bam, chrom, pos0, min_bq):
    """Per-base counts and per-base [fwd,rev] strand at a single position."""
    counts = {b: 0 for b in BASES}
    strand = {b: [0, 0] for b in BASES}
    for read in bam.fetch(chrom, pos0, pos0 + 1):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_qcfail or read.is_duplicate):
            continue
        if read.query_sequence is None:
            continue
        qpos = None
        for qp, rp in read.get_aligned_pairs(matches_only=True):
            if rp == pos0:
                qpos = qp
                break
        if qpos is None:
            continue
        quals = read.query_qualities
        if quals is not None and quals[qpos] < min_bq:
            continue
        base = read.query_sequence[qpos].upper()
        if base in counts:
            counts[base] += 1
            strand[base][1 if read.is_reverse else 0] += 1
    return counts, strand


def count_indel(bam, chrom, anchor_pos0, indel_type, indel_len, shimmer):
    """Count reads supporting an indel of the given type and length near the
    anchor, plus the spanning depth. The event sits just after the anchor base.
    Matching is by type (INS/DEL) and length within +/- shimmer bp, tolerant of
    alignment placement; the exact inserted/deleted sequence is not required."""
    event_pos = anchor_pos0 + 1
    support = 0
    depth = 0
    for read in bam.fetch(chrom, anchor_pos0, anchor_pos0 + 1):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_qcfail or read.is_duplicate):
            continue
        if read.cigartuples is None:
            continue
        if not (read.reference_start <= anchor_pos0 < read.reference_end):
            continue
        depth += 1
        rp = read.reference_start
        found = False
        for op, ln in read.cigartuples:
            if op in (0, 7, 8):          # M, =, X: consume reference and query
                rp += ln
            elif op == 2:                # D: deletion, consumes reference
                if (indel_type == "DEL" and ln == indel_len
                        and abs(rp - event_pos) <= shimmer):
                    found = True
                rp += ln
            elif op == 1:                # I: insertion, consumes query only
                if (indel_type == "INS" and ln == indel_len
                        and abs(rp - event_pos) <= shimmer):
                    found = True
            elif op == 3:                # N: ref skip
                rp += ln
            # S(4), H(5), P(6): no reference advance
        if found:
            support += 1
    return support, depth


def load_indel_blocklist(path):
    """Return {(chrom,pos,ref,alt): n_controls} from the recurrence table."""
    bl = {}
    if not path:
        return bl
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            bl[(f[0], int(f[1]), f[2].upper(), f[3].upper())] = int(f[4])
    return bl


def score_indel(bam, chrom, pos, ref, alt, blocklist, masked, args):
    """Score one indel candidate and return its 18-column output row.
    Strand and the Beta tests do not apply to indels; those columns are NA."""
    if "N" in ref or "N" in alt:
        support = depth = 0
        vaf = 0.0
        final = "unsupported_n_allele"
    elif len(alt) > len(ref):
        support, depth = count_indel(bam, chrom, pos - 1, "INS",
                                     len(alt) - len(ref), args.indel_shimmer)
        vaf = support / depth if depth else 0.0
        final = None
    elif len(ref) > len(alt):
        support, depth = count_indel(bam, chrom, pos - 1, "DEL",
                                     len(ref) - len(alt), args.indel_shimmer)
        vaf = support / depth if depth else 0.0
        final = None
    else:
        support = depth = 0
        vaf = 0.0
        final = "unsupported_mnv"

    if final is None:
        n_ctrl = blocklist.get((chrom, pos, ref, alt), 0)
        if (chrom, pos) in masked:
            final = "no_call_masked"
        elif n_ctrl >= args.indel_min_controls:
            final = f"reject_blocklist(n={n_ctrl})"
        elif depth < args.min_depth:
            final = "no_call_lowdepth"
        elif support < args.min_indel_alt or vaf < args.min_indel_vaf:
            final = "background"
        else:
            final = "call"

    return "\t".join([chrom, str(pos), ref, alt, str(support), str(depth),
                      f"{vaf:.3e}"] + ["NA"] * 10 + [final])


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="Patient SSCS BAM (unfiltered)")
    ap.add_argument("--ref", required=True, help="Reference FASTA (unmasked)")
    ap.add_argument("--model", required=True, help="Error model TSV")
    ap.add_argument("--sites", required=True,
                    help="Candidate TSV: chrom, pos, ref, alt (1-based)")
    ap.add_argument("--out", required=True, help="Output calls TSV")
    ap.add_argument("--mask", help="Artifact mask BED (no-call positions)")
    ap.add_argument("--filtered-bam",
                    help="Where to write the NM-filtered patient BAM "
                         "(default: alongside --out)")
    ap.add_argument("--skip-filter", action="store_true",
                    help="Treat --bam as already NM-filtered; skip calmd step")
    ap.add_argument("--max-nm", type=int, default=10,
                    help="NM ceiling for the patient filter. Default 10")
    ap.add_argument("--min-bq", type=int, default=20,
                    help="Minimum base quality. Default 20")
    ap.add_argument("--min-depth", type=int, default=100,
                    help="Minimum recomputed depth to attempt a call. Default 100")
    ap.add_argument("--strand-thresh", type=float, default=0.90,
                    help="Reject if alt reads are >= this fraction one strand. "
                         "Default 0.90")
    ap.add_argument("--strand-min-alt", type=int, default=10,
                    help="Apply the strand filter only when alt count is at least "
                         "this, to protect ultra-low-VAF calls. Default 10")
    ap.add_argument("--alpha-level", type=float, default=0.005,
                    help="P-value threshold for a positive call. Default 0.005")
    ap.add_argument("--test", choices=["pbeta", "betabinom", "both"],
                    default="betabinom",
                    help="Which test drives the final call. 'both' requires both "
                         "to pass. Default betabinom")
    ap.add_argument("--indel-blocklist",
                    help="Indel recurrence TSV from build_indel_blocklist.py")
    ap.add_argument("--indel-min-controls", type=int, default=2,
                    help="Block an indel allele present in at least this many "
                         "controls. Default 2")
    ap.add_argument("--min-indel-alt", type=int, default=3,
                    help="Minimum supporting reads for an indel call. Default 3")
    ap.add_argument("--min-indel-vaf", type=float, default=1e-4,
                    help="Minimum VAF for an indel call. Default 1e-4")
    ap.add_argument("--indel-shimmer", type=int, default=3,
                    help="Allowed bp offset when matching an indel near the "
                         "anchor. Default 3")
    args = ap.parse_args()

    # 1. filter the patient BAM
    if args.skip_filter:
        filt_bam = args.bam
    else:
        filt_bam = args.filtered_bam or (os.path.splitext(args.out)[0]
                                         + ".patient.nm{}.bam".format(args.max_nm))
        sys.stderr.write(f"Filtering patient BAM (NM<={args.max_nm}) -> {filt_bam}\n")
        filter_patient_bam(args.bam, args.ref, filt_bam, args.max_nm)

    sites, priors = load_model(args.model)
    masked = load_mask(args.mask)
    indel_blocklist = load_indel_blocklist(args.indel_blocklist)
    bam = pysam.AlignmentFile(filt_bam, "rb")

    cols = ["#chrom", "pos", "ref", "alt", "alt_count", "depth", "vaf",
            "alt_fwd", "alt_rev", "strand_frac", "alpha", "beta",
            "site_bg_pct", "pbeta_p", "pbeta_call", "betabinom_p",
            "betabinom_call", "final_call"]

    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for chrom, pos, ref, alt in read_sites(args.sites):
            # substitutions go through the Beta model; everything else (indels,
            # and any non-substitution allele) goes through the indel path.
            if not (len(ref) == 1 and len(alt) == 1
                    and ref in BASES and alt in BASES):
                out.write(score_indel(bam, chrom, pos, ref, alt,
                                      indel_blocklist, masked, args) + "\n")
                continue
            counts, strand = count_at(bam, chrom, pos - 1, args.min_bq)
            depth = sum(counts.values())
            alt_count = counts.get(alt, 0)
            fwd, rev = strand.get(alt, [0, 0])
            vaf = alt_count / depth if depth else 0.0
            strand_frac = (max(fwd, rev) / alt_count) if alt_count else 0.0

            # site Beta: per-site if present, else the substitution-class prior
            sub = f"{ref}>{alt}"
            if (chrom, pos, ref, alt) in sites:
                alpha, beta, bg = sites[(chrom, pos, ref, alt)]
            elif sub in priors:
                alpha, beta = priors[sub]
                bg = float("nan")
            else:
                alpha = beta = bg = float("nan")

            # the two tests
            if depth and not math.isnan(alpha):
                pbeta_p = 1.0 - beta_cdf(vaf, alpha, beta)
                bb_p = betabinom_sf(alt_count, depth, alpha, beta)
            else:
                pbeta_p = bb_p = float("nan")
            pbeta_call = (not math.isnan(pbeta_p)) and pbeta_p <= args.alpha_level
            bb_call = (not math.isnan(bb_p)) and bb_p <= args.alpha_level

            # decision, with mask and strand overriding the statistics
            if (chrom, pos) in masked:
                final = "no_call_masked"
            elif depth < args.min_depth:
                final = "no_call_lowdepth"
            elif alt_count >= args.strand_min_alt and strand_frac >= args.strand_thresh:
                final = "reject_strand"
            elif math.isnan(alpha):
                final = "no_model"
            else:
                if args.test == "pbeta":
                    passed = pbeta_call
                elif args.test == "betabinom":
                    passed = bb_call
                else:
                    passed = pbeta_call and bb_call
                final = "call" if passed else "background"

            def fmt(x):
                return "NA" if isinstance(x, float) and math.isnan(x) else (
                    f"{x:.3e}" if isinstance(x, float) else str(x))

            out.write("\t".join([
                chrom, str(pos), ref, alt, str(alt_count), str(depth),
                f"{vaf:.3e}", str(fwd), str(rev), f"{strand_frac:.3f}",
                fmt(alpha), fmt(beta), fmt(bg), fmt(pbeta_p), str(pbeta_call),
                fmt(bb_p), str(bb_call), final]) + "\n")

    sys.stderr.write(f"Wrote calls to {args.out}\n")


if __name__ == "__main__":
    main()
