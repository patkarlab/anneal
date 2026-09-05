#!/usr/bin/env python3
"""
make_homopolymer_sites.py  (Step 1 of 3)

Generate indel candidates at every homopolymer run in the panel, for an indel
specificity sweep. Homopolymer and short-repeat contexts are where polymerase
slippage produces artefactual indels, so they are the worst case for the
alt-count/VAF gate. Running the apply step on the eight normals with these
candidates measures the false-positive rate of the indel arm where it is most
vulnerable.

For each run of >= --min-run identical bases, two candidates are emitted in
left-aligned anchored form: a 1 bp deletion and a 1 bp insertion of the run
base. The apply step's shimmer-tolerant matching catches slippage events placed
anywhere in the run.

Panel BED coordinates are read from columns 2-3 only (column 4 is a legacy
label).

Output: a sites TSV (chrom, pos, ref, alt) plus a context note column.
"""

import argparse
import sys
import pysam


def read_bed_intervals(path):
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            yield f[0], int(f[1]), int(f[2])


def homopolymer_runs(seq, min_run):
    """Yield (start_index, base, length) for runs of >= min_run identical bases."""
    n = len(seq)
    i = 0
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        run_len = j - i
        if seq[i] in "ACGT" and run_len >= min_run:
            yield i, seq[i], run_len
        i = j


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="Panel BED (coords from cols 2-3)")
    ap.add_argument("--ref", required=True, help="Reference FASTA (indexed)")
    ap.add_argument("--out", required=True, help="Output sites TSV")
    ap.add_argument("--min-run", type=int, default=5,
                    help="Minimum homopolymer length to test. Default 5")
    args = ap.parse_args()

    fasta = pysam.FastaFile(args.ref)
    contigs = set(fasta.references)
    seen = set()
    n_runs = 0

    with open(args.out, "w") as out:
        out.write("#chrom\tpos\tref\talt\tcontext\n")
        for chrom, start, stop in read_bed_intervals(args.bed):
            if chrom not in contigs:
                continue
            seq = fasta.fetch(chrom, start, stop).upper()
            for idx, base, run_len in homopolymer_runs(seq, args.min_run):
                s0 = start + idx          # genomic 0-based start of the run
                if s0 == 0:
                    continue              # need the base before the run
                before = fasta.fetch(chrom, s0 - 1, s0).upper()
                if before not in "ACGT":
                    continue
                pos = s0                  # 1-based anchor = (s0-1)+1
                note = f"{base}x{run_len}"
                # 1 bp deletion: REF = before+base, ALT = before
                key_d = (chrom, pos, before + base, before)
                if key_d not in seen:
                    seen.add(key_d)
                    out.write(f"{chrom}\t{pos}\t{before + base}\t{before}\tdel_{note}\n")
                # 1 bp insertion: REF = before, ALT = before+base
                key_i = (chrom, pos, before, before + base)
                if key_i not in seen:
                    seen.add(key_i)
                    out.write(f"{chrom}\t{pos}\t{before}\t{before + base}\tins_{note}\n")
                n_runs += 1

    sys.stderr.write(f"Found {n_runs} homopolymer runs (>= {args.min_run} bp); "
                     f"wrote {len(seen)} indel candidates to {args.out}\n")


if __name__ == "__main__":
    main()
