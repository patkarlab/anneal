#!/usr/bin/env python3
"""
inspect_locus_reads.py

Diagnostic for a single suspect locus. Partitions the consensus reads covering
a position into those carrying the suspect ALT base and those carrying the REF
base, then reports, per group: read count, strand balance, and how many OTHER
mismatches against the reference each read carries.

Purpose: separate a genuine per-site signal from a paralog / misalignment
artifact. Reads originating from a divergent paralogous copy carry several
LINKED mismatches at consistent nearby positions; genuine variant- or
error-bearing reads differ from the reference only at the locus itself. A large
excess of co-occurring mismatches among the ALT reads, clustered at a few
recurrent positions, is the paralog signature.

No MD tag is required; mismatches are computed directly against the FASTA.
"""

import argparse
from collections import Counter
import pysam

BASES = ("A", "C", "G", "T")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="Consensus BAM (indexed)")
    ap.add_argument("--ref", required=True, help="Reference FASTA (indexed)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--pos", type=int, required=True,
                    help="1-based position of the suspect locus")
    ap.add_argument("--alt", required=True, choices=BASES,
                    help="Suspect ALT base to partition on")
    ap.add_argument("--min-bq", type=int, default=20,
                    help="Minimum base quality at the locus to use a read "
                         "(default 20)")
    ap.add_argument("--flank", type=int, default=150,
                    help="Report co-mismatch positions within this many bp of "
                         "the locus (default 150)")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.ref)
    pos0 = args.pos - 1
    ref_base = fasta.fetch(args.chrom, pos0, pos0 + 1).upper()

    groups = {"ALT": [], "REF": []}   # each entry: (n_other_mismatches, is_reverse)
    co_positions = Counter()          # co-mismatch ref positions among ALT reads

    for read in bam.fetch(args.chrom, pos0, pos0 + 1):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_qcfail or read.is_duplicate):
            continue
        if read.query_sequence is None:
            continue
        seq = read.query_sequence.upper()
        quals = read.query_qualities

        pairs = read.get_aligned_pairs(matches_only=True)  # (qpos, rpos), no indels

        # base carried at the target position
        target_qpos = None
        for qp, rp in pairs:
            if rp == pos0:
                target_qpos = qp
                break
        if target_qpos is None:
            continue  # locus deleted/skipped in this read
        if quals is not None and quals[target_qpos] < args.min_bq:
            continue

        base_here = seq[target_qpos]
        if base_here == args.alt:
            grp = "ALT"
        elif base_here == ref_base:
            grp = "REF"
        else:
            continue  # third allele; ignore for this two-way contrast

        # count OTHER substitution mismatches against the reference
        ref_span = fasta.fetch(args.chrom, read.reference_start,
                               read.reference_end).upper()
        n_other = 0
        for qp, rp in pairs:
            if rp == pos0:
                continue
            if seq[qp] != ref_span[rp - read.reference_start]:
                n_other += 1
                if abs(rp - pos0) <= args.flank:
                    co_positions[rp + 1] += 1  # 1-based for reporting
        groups[grp].append((n_other, read.is_reverse))

    def summarize(name):
        g = groups[name]
        n = len(g)
        if n == 0:
            return f"{name}: 0 reads"
        rev = sum(1 for _, r in g if r)
        mean_mm = sum(mm for mm, _ in g) / n
        return (f"{name}: {n} reads (fwd {n - rev}, rev {rev}); "
                f"mean other mismatches/read = {mean_mm:.2f}")

    print(f"Locus {args.chrom}:{args.pos}  ref={ref_base}  alt={args.alt}")
    print(summarize("REF"))
    print(summarize("ALT"))
    print(f"Top co-occurring mismatch positions among ALT reads "
          f"(within {args.flank} bp):")
    for p, c in co_positions.most_common(8):
        print(f"  {args.chrom}:{p}  in {c} ALT reads")


if __name__ == "__main__":
    main()
