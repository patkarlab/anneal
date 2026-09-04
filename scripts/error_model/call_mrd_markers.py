#!/usr/bin/env python3
"""
Marker-tracking MRD caller (DCS track).

Counts reads directly from the consensus BAM at the patient's diagnosis variant
positions and scores them against the BNC beta matrix. Does not read a VCF:
Pisces applies MinimumVariantQScore beneath any requested --minvf, so a marker
present at a few reads would be absent from the caller output.

Reads are counted unfiltered, matching how the beta matrix was built.

--markers is a 1-based TSV, no header:  chrom  pos  ref  alt  [label]

    call_mrd_markers.py --sample S --bam S.dcs.sc.sorted.bam \\
        --markers patient.tsv --model beta_matrix_DCS.txt \\
        --mask artifact_mask.combined.bed \\
        --indel-blocklist indel_blocklist.tsv --out S.mrd_report.tsv
"""

import math
import argparse
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apply_sscs_error_model import (          # noqa: E402
    betabinom_sf, count_at, count_indel, load_mask, load_indel_blocklist)

BASES = ("A", "C", "G", "T")
MATRIX_ORDER = ("A", "T", "G", "C")

COLS = ["sample", "label", "chrom", "pos", "ref", "alt", "type",
        "alt_count", "depth", "vaf_pct", "alt_fwd", "alt_rev",
        "strand_frac", "bg_pct", "p_value", "min_alt_callable",
        "lod_pct", "call", "note"]


def load_beta_matrix(path):
    """{(chrom, pos, alt): (alpha, beta)}. The matrix has no ref column, so
    entries key on the alt base alone."""
    table = {}
    with open(path) as fh:
        if not fh.readline().lower().startswith("chr"):
            fh.seek(0)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            try:
                chrom, pos = f[0], int(f[1])
            except ValueError:
                continue
            for i, base in enumerate(MATRIX_ORDER):
                try:
                    a, b = float(f[2 + 2 * i]), float(f[3 + 2 * i])
                except (ValueError, IndexError):
                    continue
                if a > 0 and b > 0:
                    table[(chrom, pos, base)] = (a, b)
    return table


def load_markers(path):
    out = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 4:
                sys.stderr.write(f"markers line {n}: need 4 columns, skipped\n")
                continue
            try:
                pos = int(f[1])
            except ValueError:
                if n > 1:
                    sys.stderr.write(f"markers line {n}: bad position, skipped\n")
                continue
            out.append({"chrom": f[0], "pos": pos, "ref": f[2].upper().strip(),
                        "alt": f[3].upper().strip(),
                        "label": f[4] if len(f) > 4 else ""})
    return out


def min_alt_callable(depth, alpha, beta, level, cap=500):
    """Smallest alt count reaching p < level. None if the tail stops
    decreasing first, i.e. the numerical floor sits above level."""
    if not depth or alpha is None:
        return None
    prev = None
    for k in range(1, min(cap, depth) + 1):
        p = betabinom_sf(k, depth, alpha, beta)
        if p < level:
            return k
        if prev is not None and p >= prev:
            return None
        prev = p
    return None


def fmt(x, spec=".3e"):
    return "NA" if x is None else format(x, spec)


def fisher_two_sided(a, b, c, d):
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


def score_snv(bam, row, m, matrix, masked, args):
    chrom, pos, alt = m["chrom"], m["pos"], m["alt"]
    counts, strand = count_at(bam, chrom, pos - 1, args.min_bq)
    depth = sum(counts.values())
    ac = counts.get(alt, 0)
    fwd, rev = strand.get(alt, [0, 0])
    ref_fwd, ref_rev = strand.get(m.get("ref", ""), [0, 0])
    vaf = ac / depth if depth else 0.0
    sfrac = (max(fwd, rev) / ac) if ac else 0.0
    one_sided = ac >= args.strand_min_alt and sfrac >= args.strand_thresh
    biased = one_sided and orientation_bias(fwd, rev, ref_fwd, ref_rev,
                                            args.strand_p)

    ab = matrix.get((chrom, pos, alt))
    alpha, beta = ab if ab else (None, None)
    bg = alpha / (alpha + beta) if ab else None
    p = betabinom_sf(ac, depth, alpha, beta) if (depth and ab and ac) else None
    kmin = min_alt_callable(depth, alpha, beta, args.alpha_level)

    row.update({
        "alt_count": ac, "depth": depth, "vaf_pct": fmt(100 * vaf, ".4f"),
        "alt_fwd": fwd, "alt_rev": rev,
        "strand_frac": fmt(sfrac, ".3f") if ac else "NA",
        "bg_pct": fmt(100 * bg, ".5f") if bg is not None else "NA",
        "p_value": fmt(p), "min_alt_callable": kmin or "NA",
        "lod_pct": fmt(100 * kmin / depth, ".4f") if (kmin and depth) else "NA",
    })

    if (chrom, pos) in masked:
        row["call"], row["note"] = "not_evaluable", "artifact mask"
    elif depth < args.min_depth:
        row["call"], row["note"] = "not_evaluable", f"depth < {args.min_depth}"
    elif ab is None:
        row["call"], row["note"] = "not_evaluable", "no background model"
    elif ac == 0:
        row["call"] = "not_detected"
    elif biased:
        row["call"], row["note"] = "not_detected", "strand bias"
    elif ac < args.min_alt:
        row["call"], row["note"] = "not_detected", f"alt < {args.min_alt}"
    elif p is not None and p < args.alpha_level:
        row["call"] = "DETECTED"
    else:
        row["call"], row["note"] = "not_detected", "below background"
    if one_sided and not biased and not row.get("note"):
        row["note"] = "one-sided coverage: ref %d/%d" % (ref_fwd, ref_rev)
    return row


def score_indel(bam, row, m, masked, blocklist, args):
    chrom, pos, ref, alt = m["chrom"], m["pos"], m["ref"], m["alt"]
    kind = "INS" if len(alt) > len(ref) else "DEL"
    support, depth = count_indel(bam, chrom, pos - 1, kind,
                                 abs(len(alt) - len(ref)), args.indel_shimmer)
    vaf = support / depth if depth else 0.0
    n_ctrl = blocklist.get((chrom, pos, ref, alt), 0)

    row.update({
        "alt_count": support, "depth": depth,
        "vaf_pct": fmt(100 * vaf, ".4f"), "bg_pct": f"blocklist_n={n_ctrl}",
        "min_alt_callable": args.min_indel_alt,
        "lod_pct": fmt(100 * args.min_indel_alt / depth, ".4f") if depth else "NA",
    })

    if (chrom, pos) in masked:
        row["call"], row["note"] = "not_evaluable", "artifact mask"
    elif n_ctrl >= args.indel_min_controls:
        row["call"] = "not_evaluable"
        row["note"] = f"blocklisted in {n_ctrl} controls"
    elif depth < args.min_depth:
        row["call"], row["note"] = "not_evaluable", f"depth < {args.min_depth}"
    elif support >= args.min_indel_alt and vaf >= args.min_indel_vaf:
        row["call"] = "DETECTED"
    else:
        row["call"] = "not_detected"
    return row


def score_marker(bam, m, matrix, masked, blocklist, args):
    row = {c: "NA" for c in COLS}
    row.update({"sample": args.sample, "label": m["label"], "chrom": m["chrom"],
                "pos": m["pos"], "ref": m["ref"], "alt": m["alt"], "note": ""})
    ref, alt = m["ref"], m["alt"]

    if len(ref) == 1 and len(alt) == 1 and ref in BASES and alt in BASES:
        row["type"] = "snv"
        return score_snv(bam, row, m, matrix, masked, args)
    if len(ref) != len(alt):
        row["type"] = "indel"
        return score_indel(bam, row, m, masked, blocklist, args)
    row["type"] = "mnv"
    row["call"], row["note"] = "not_evaluable", "MNV not supported"
    return row


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True, help="DCS consensus BAM")
    ap.add_argument("--markers", required=True)
    ap.add_argument("--model", required=True, help="beta_matrix_DCS.txt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mask")
    ap.add_argument("--indel-blocklist")
    # 0.005 is the source-method cutoff, appropriate at 3-5 tests. The 1e-9 used
    # for untargeted scoring is a multiple-testing correction and drops real
    # markers here (PTPN11 scored 1.74e-09 on the G rung). See ARCHITECTURE.md.
    ap.add_argument("--alpha-level", type=float, default=5e-3)
    ap.add_argument("--min-alt", type=int, default=2)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--strand-thresh", type=float, default=0.90)
    ap.add_argument("--strand-min-alt", type=int, default=10)
    ap.add_argument("--strand-p", type=float, default=1e-3,
                    help="Fisher p below which alt orientation is "
                         "called biased against reference orientation "
                         "(default 1e-3)")
    ap.add_argument("--indel-min-controls", type=int, default=2)
    ap.add_argument("--min-indel-alt", type=int, default=3)
    ap.add_argument("--min-indel-vaf", type=float, default=1e-4)
    ap.add_argument("--indel-shimmer", type=int, default=3)
    args = ap.parse_args()

    markers = load_markers(args.markers)
    if not markers:
        sys.exit("FATAL: no usable markers")
    matrix = load_beta_matrix(args.model)
    masked = load_mask(args.mask) if args.mask else set()
    blocklist = load_indel_blocklist(args.indel_blocklist) if args.indel_blocklist else {}

    bam = pysam.AlignmentFile(args.bam, "rb")
    rows = [score_marker(bam, m, matrix, masked, blocklist, args) for m in markers]
    bam.close()

    with open(args.out, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in COLS) + "\n")

    det = sum(1 for r in rows if r["call"] == "DETECTED")
    nev = sum(1 for r in rows if r["call"] == "not_evaluable")
    sys.stderr.write(f"\n{args.sample}: {det} detected, "
                     f"{len(rows) - det - nev} not detected, {nev} not evaluable\n")
    for r in rows:
        tag = r["label"] or f"{r['chrom']}:{r['pos']}"
        sys.stderr.write(f"  {tag:24} {r['ref']}>{r['alt']:<6} {r['call']:<14} "
                         f"VAF {r['vaf_pct']:>9}%  alt {r['alt_count']:>5}/"
                         f"{r['depth']:<7} LoD {r['lod_pct']:>8}%  p {r['p_value']}\n")
    sys.stderr.write(f"\nwrote {args.out}\n")


if __name__ == "__main__":
    main()
