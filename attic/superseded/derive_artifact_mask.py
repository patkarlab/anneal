#!/usr/bin/env python3
"""
derive_artifact_mask.py

Classify the inflated loci in an SSCS error model into two groups, using the
strand distribution of the alt reads in a representative filtered (NM<=10) BNC
BAM:

  - strand-clearable: the alt reads are strongly one-strand (>= --strand-thresh),
    the signature of a paralog / mis-mapped copy that aligns in a fixed
    orientation. These are handled at calling time by the strand-bias filter,
    because a real variant would inherit the locus strand balance.

  - mask: the alt reads are NOT strand-skewed, so the strand filter cannot tell
    a paralog read from a real variant read. These positions are genuinely
    irreducible by read features and are written to a mask BED for no-call.

Input is the model TSV from build_sscs_error_model.py and one filtered BNC BAM
(the paralog signal is systematic across controls, so a single representative
BAM suffices). Output is a BED of masked positions plus a per-locus report.

The strand criterion here is absolute (one-strand fraction), which is the right
test for these known-artifact loci. The calling-time strand filter uses a
relative test (alt strand vs the locus ref strand) so it is robust to ordinary
capture-strand bias on real candidates; that lives in the apply step.
"""

import argparse
import sys
import pysam

BASES = ("A", "C", "G", "T")


def alt_strand_counts(bam, chrom, pos0, alt, min_bq):
    """Count alt-supporting reads on forward and reverse strands at a position."""
    fwd = rev = 0
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
        if read.query_sequence[qpos].upper() == alt:
            if read.is_reverse:
                rev += 1
            else:
                fwd += 1
    return fwd, rev


def read_inflated_loci(model_path, min_bg):
    """Yield (chrom, pos, alt, background_pct) for rows above the background cut.
    Column order matches build_sscs_error_model.py output."""
    with open(model_path) as fh:
        for line in fh:
            if line.startswith(("#", "##")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            try:
                bg = float(f[14])
            except ValueError:
                continue
            if bg > min_bg:
                yield f[0], int(f[1]), f[3], bg


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="Error model TSV")
    ap.add_argument("--bam", required=True,
                    help="One representative filtered (NM<=10) BNC BAM")
    ap.add_argument("--out", required=True, help="Output mask BED")
    ap.add_argument("--report", help="Optional per-locus classification TSV")
    ap.add_argument("--min-bg", type=float, default=1.0,
                    help="Only classify positions with background_pct above this "
                         "(percent). Default 1.0")
    ap.add_argument("--strand-thresh", type=float, default=0.90,
                    help="Alt one-strand fraction at or above which a locus is "
                         "treated as strand-clearable paralog (not masked). "
                         "Default 0.90")
    ap.add_argument("--min-bq", type=int, default=20,
                    help="Minimum base quality at the locus. Default 20")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    masked_positions = {}   # (chrom, pos) -> reason
    rows = []

    for chrom, pos, alt, bg in read_inflated_loci(args.model, args.min_bg):
        fwd, rev = alt_strand_counts(bam, chrom, pos - 1, alt, args.min_bq)
        total = fwd + rev
        if total == 0:
            frac = 0.0
            verdict = "MASK"      # inflated but no surviving alt reads: uninterpretable
        else:
            frac = max(fwd, rev) / total
            verdict = "clearable" if frac >= args.strand_thresh else "MASK"
        rows.append((chrom, pos, alt, bg, fwd, rev, frac, verdict))
        if verdict == "MASK":
            masked_positions[(chrom, pos)] = f"{alt}:strand{frac:.2f}"

    # write mask BED (unique positions, 0-based half-open)
    with open(args.out, "w") as out:
        for (chrom, pos) in sorted(masked_positions):
            out.write(f"{chrom}\t{pos - 1}\t{pos}\t{masked_positions[(chrom, pos)]}\n")

    if args.report:
        with open(args.report, "w") as rep:
            rep.write("#chrom\tpos\talt\tbackground_pct\talt_fwd\talt_rev\t"
                      "one_strand_frac\tverdict\n")
            for chrom, pos, alt, bg, fwd, rev, frac, verdict in rows:
                rep.write(f"{chrom}\t{pos}\t{alt}\t{bg}\t{fwd}\t{rev}\t"
                          f"{frac:.3f}\t{verdict}\n")

    n_mask = len(masked_positions)
    n_clear = sum(1 for r in rows if r[7] == "clearable")
    sys.stderr.write(f"Classified {len(rows)} inflated loci: {n_clear} strand-"
                     f"clearable, {n_mask} masked positions written to {args.out}\n")


if __name__ == "__main__":
    main()
