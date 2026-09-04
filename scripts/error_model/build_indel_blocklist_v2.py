#!/usr/bin/env python3
"""
build_indel_blocklist_v2.py

Builds the indel recurrence table from scan_indels.py TSV output.

WHY THIS REPLACES build_indel_blocklist.py
------------------------------------------
The original reads VCF and indexes f[0], f[1], f[3], f[4] as chrom, pos, ref,
alt. scan_indels.py has emitted a TSV since commit 6812a01 (31 Jul 2026), so
those indices now land on sample, track, pos, ref and the parse dies on
int('SSCS').

More importantly, the blocklist in results_bnc/ is dated 23 Jun 2026 -- five
weeks before scan_indels.py existed -- so it was built from Pisces VCFs. That
same commit message records why that is fatal:

    "Pisces cannot report indels at MRD frequencies. Lowest indel Pisces
     reports is 0.128% against a substitution floor near 0.01%."

So every allele in the old table is a high-frequency indel, and the systematic
low-frequency artifacts the blocklist exists to catch could never have entered
it. Rebuilding from BAM-derived indels is the point of this script.

WHAT IT ADDS
------------
The original recorded recurrence only, so an allele seen at 3 reads in all
eight controls and one seen at 5000 reads in all eight were indistinguishable.
Panel data shows both extremes: chrX:134393658 and chr17:76737017 run to
thousands of reads, while most recurrent alleles sit near the 2-read floor.

Columns 1-5 keep the original meaning and order, so
apply_sscs_error_model.py::load_indel_blocklist reads it unchanged (it takes
f[0..4] positionally and ignores the rest). Columns 6+ carry magnitude for
apply steps that want to threshold on more than recurrence.

Output TSV:
  #chrom pos ref alt n_controls median_vaf_pct max_vaf_pct median_count
  max_count median_strand_frac type length

USAGE
-----
    python3 build_indel_blocklist_v2.py \\
        --tsvs results_bnc/indel_vcfs/*.sscs.indels.tsv \\
        --out results_bnc/indel_blocklist.patched.tsv
"""

import argparse
import os
import sys
from collections import defaultdict

REQUIRED = ("chrom", "pos", "ref", "alt")


def median(values):
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2:
        return float(s[mid])
    return (s[mid - 1] + s[mid]) / 2.0


def read_indel_tsv(path):
    """Yield dicts for indel rows in one scan_indels.py TSV.

    Column positions are resolved from the header rather than hardcoded, so a
    future column reorder does not silently shift the parse the way the VCF
    indices did.
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in REQUIRED if c not in idx]
        if missing:
            sys.exit("ERROR: %s lacks column(s) %s. Header was:\n  %s"
                     % (path, ", ".join(missing), "\t".join(header)))

        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue                      # truncated row
            try:
                ref = f[idx["ref"]].upper()
                alt = f[idx["alt"]].upper()
                if len(ref) == len(alt):
                    continue                  # substitution, not an indel
                if "N" in ref or "N" in alt:
                    continue                  # consensus no-call, not an event
                rec = {
                    "chrom": f[idx["chrom"]],
                    "pos": int(f[idx["pos"]]),
                    "ref": ref,
                    "alt": alt,
                }
            except (ValueError, IndexError):
                continue

            for name, cast in (("count", int), ("depth", int),
                               ("vaf_pct", float), ("strand_frac", float),
                               ("length", int)):
                try:
                    rec[name] = cast(f[idx[name]])
                except (ValueError, IndexError, KeyError):
                    rec[name] = 0
            rec["type"] = f[idx["type"]] if "type" in idx else ""
            yield rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsvs", nargs="+", required=True,
                    help="scan_indels.py output, one per control sample")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-controls", type=int, default=1,
                    help="Only emit alleles seen in at least this many "
                         "controls (default 1; the apply step does the real "
                         "thresholding, so keep this low)")
    args = ap.parse_args()

    for path in args.tsvs:
        if not os.path.isfile(path):
            sys.exit("ERROR: not found: %s" % path)

    # allele -> {sample: record}. Keyed by sample so a duplicated input file
    # cannot inflate n_controls.
    seen = defaultdict(dict)

    for path in args.tsvs:
        sample = os.path.basename(path).split(".")[0]
        n = 0
        for rec in read_indel_tsv(path):
            key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
            seen[key][sample] = rec
            n += 1
        print("  %-40s %6d indels" % (os.path.basename(path), n),
              file=sys.stderr)

    n_samples = len({os.path.basename(p).split(".")[0] for p in args.tsvs})
    print("\n  %d control samples, %d distinct alleles"
          % (n_samples, len(seen)), file=sys.stderr)

    rows = []
    for key, per_sample in seen.items():
        recs = list(per_sample.values())
        if len(recs) < args.min_controls:
            continue
        rows.append((
            key[0], key[1], key[2], key[3], len(recs),
            median([r["vaf_pct"] for r in recs]),
            max(r["vaf_pct"] for r in recs),
            median([r["count"] for r in recs]),
            max(r["count"] for r in recs),
            median([r["strand_frac"] for r in recs]),
            recs[0]["type"],
            recs[0]["length"],
        ))

    # Most recurrent first, then largest, so the head of the file is the
    # set of alleles that matter most.
    rows.sort(key=lambda r: (-r[4], -r[8]))

    with open(args.out, "w") as out:
        print("#chrom\tpos\tref\talt\tn_controls\tmedian_vaf_pct\t"
              "max_vaf_pct\tmedian_count\tmax_count\tmedian_strand_frac\t"
              "type\tlength", file=out)
        for r in rows:
            print("%s\t%d\t%s\t%s\t%d\t%.4f\t%.4f\t%.1f\t%d\t%.3f\t%s\t%d"
                  % r, file=out)

    print("\n  recurrence distribution:", file=sys.stderr)
    dist = defaultdict(int)
    for r in rows:
        dist[r[4]] += 1
    for k in sorted(dist):
        print("    n_controls=%-2d  %6d" % (k, dist[k]), file=sys.stderr)
    print("\n  written: %s (%d alleles)" % (args.out, len(rows)),
          file=sys.stderr)


if __name__ == "__main__":
    main()
