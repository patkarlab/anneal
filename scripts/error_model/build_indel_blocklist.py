#!/usr/bin/env python3
"""
build_indel_blocklist.py

Build an indel recurrence table from the BNC SSCS VCFs. Each indel allele is
keyed by its anchored representation (chrom, pos, ref, alt) exactly as anneal
writes it, and the table records how many distinct control samples carry it.

The apply step treats an allele as a systematic artifact when its control count
reaches a threshold k, applied at load time so k stays tunable without
rebuilding. Recurrence across many controls is reproducible chemistry/alignment
behaviour; a real patient marker such as NPM1 type A is absent from the controls
and so never enters the table.

Indels whose inserted or deleted sequence contains N are excluded: those are
consensus base-masking no-calls, not real events.

Input: the BNC SSCS VCFs. VCF FORMAT is GT:ALT:TOT:FRAC; an indel is any record
where REF and ALT differ in length. Output: TSV of chrom, pos, ref, alt,
n_controls, sorted by descending recurrence.
"""

import argparse
import sys
from collections import defaultdict


def indel_alleles_in_vcf(path):
    """Return the set of (chrom, pos, ref, alt) indel alleles present in one VCF,
    excluding any allele whose REF or ALT contains N."""
    alleles = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[3].upper(), f[4].upper()
            if len(ref) == len(alt):
                continue                      # substitution, not an indel
            if "N" in ref or "N" in alt:
                continue                      # base-masking no-call, not real
            alleles.add((chrom, int(pos), ref, alt))
    return alleles


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcfs", nargs="+", required=True,
                    help="BNC SSCS VCFs")
    ap.add_argument("--out", required=True, help="Output recurrence TSV")
    args = ap.parse_args()

    control_counts = defaultdict(int)
    for path in args.vcfs:
        for allele in indel_alleles_in_vcf(path):
            control_counts[allele] += 1

    rows = sorted(control_counts.items(),
                  key=lambda kv: (-kv[1], kv[0][0], kv[0][1]))
    with open(args.out, "w") as out:
        out.write("#chrom\tpos\tref\talt\tn_controls\n")
        for (chrom, pos, ref, alt), n in rows:
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{n}\n")

    n_total = len(rows)
    n_vcfs = len(args.vcfs)
    by_k = defaultdict(int)
    for _, n in rows:
        by_k[n] += 1
    spectrum = ", ".join(f"{k}:{by_k[k]}" for k in range(1, n_vcfs + 1) if by_k[k])
    sys.stderr.write(f"{n_total} distinct non-N indel alleles across {n_vcfs} "
                     f"controls written to {args.out}\n")
    sys.stderr.write(f"recurrence spectrum (n_controls:count): {spectrum}\n")


if __name__ == "__main__":
    main()
