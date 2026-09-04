#!/usr/bin/env python3
"""
build_candidates.py

Build the blind per-sample candidate list that stage 5 scores with
call_mrd_markers.py.

Candidates are the union of
  1. every record in the Pisces variants-only VCF for the track, and
  2. every record in the stage 2b indel table (scan_indels.py output),
both already VCF-anchored (chrom, 1-based pos, ref, alt), deduplicated on
(chrom, pos, ref, alt). No diagnosis list is read: the list is derived from
the sample's own reads only.

The optional annotation tables contribute a label (Gene|Consequence|HGVSp)
so it travels through call_mrd_markers.py into the calls table. A candidate
with no annotation gets label NA. Every label ends with the source of the
candidate: ;src=P (Pisces VCF), ;src=S (indel scan) or ;src=PS (both).

Output is the --markers format of call_mrd_markers.py, 1-based, no header:

    chrom  pos  ref  alt  label

Usage:
    build_candidates.py --sample S --track dcs \
        --vcf S.dcs.vcf --indels S.dcs.indels.tsv \
        [--annotated S.dcs.filtered.tsv] \
        [--indels-annotated S.dcs.indels.annotated.tsv] \
        --out S.dcs.candidates.tsv
"""

import argparse
import csv
import os
import sys

MISSING = {"", ".", "NA", "-"}

# Alleles containing N are not candidates. In the patched engine an N inside a
# consensus read means the two strands disagreed at that base, so an insertion
# such as G>GN is an unresolved event, not an allele. The indel scan reports
# these as separate records with the same support as the resolved allele.
SKIPPED_N = {"count": 0}


def has_n(ref, alt):
    if "N" in ref.upper() or "N" in alt.upper():
        SKIPPED_N["count"] += 1
        return True
    return False


def natural_chrom_key(chrom):
    """Sort chromosomes chr1..chr22, chrX, chrY, chrM, then anything else."""
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name.upper() in order:
        return (0, order[name.upper()], "")
    return (1, 0, chrom)


def read_vcf(path):
    """Yield (chrom, pos, ref, alt) for every ALT allele in a VCF."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alts = f[0], f[1], f[3], f[4]
            try:
                pos = int(pos)
            except ValueError:
                continue
            for alt in alts.split(","):
                if alt in MISSING or alt == "*" or "<" in alt or alt == ref:
                    continue
                if has_n(ref, alt):
                    continue
                yield chrom, pos, ref, alt


def read_indels(path):
    """Yield (chrom, pos, ref, alt) from the scan_indels.py TSV (header-driven)."""
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        needed = {"chrom", "pos", "ref", "alt"}
        if not needed.issubset(reader.fieldnames or []):
            sys.exit("ERROR: %s lacks columns %s" % (path, sorted(needed)))
        for row in reader:
            try:
                pos = int(row["pos"])
            except ValueError:
                continue
            if row["ref"] in MISSING or row["alt"] in MISSING:
                continue
            if has_n(row["ref"], row["alt"]):
                continue
            yield row["chrom"], pos, row["ref"], row["alt"]


def make_label(row):
    parts = []
    for col in ("Gene", "Consequence", "HGVSp"):
        val = (row.get(col) or "").strip()
        if val not in MISSING:
            parts.append(val)
    return "|".join(parts) if parts else "NA"


def read_labels(path, chrom_col, pos_cols, ref_col, alt_col):
    """Map (chrom, pos, ref, alt) -> label. Multiple pos columns are all keyed."""
    labels = {}
    if not path:
        return labels
    if not os.path.isfile(path):
        sys.stderr.write("WARNING: annotation not found, labels skipped: %s\n" % path)
        return labels
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        for col in (chrom_col, ref_col, alt_col):
            if col not in cols:
                sys.stderr.write("WARNING: %s lacks column %s, labels skipped\n" % (path, col))
                return labels
        for row in reader:
            label = make_label(row)
            if label == "NA":
                continue
            for pos_col in pos_cols:
                if pos_col not in cols:
                    continue
                try:
                    pos = int(row[pos_col])
                except (TypeError, ValueError):
                    continue
                key = (row[chrom_col], pos, row[ref_col], row[alt_col])
                labels.setdefault(key, label)
    return labels


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--track", required=True)
    ap.add_argument("--vcf", required=True, help="Pisces variants-only VCF")
    ap.add_argument("--indels", required=True, help="scan_indels.py TSV")
    ap.add_argument("--annotated", help="stage 3 substitution table (Chr/Start/End/Ref/Alt/Gene/...)")
    ap.add_argument("--indels-annotated", help="stage 3 indel table (chrom/pos/ref/alt/Gene/...)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sources = {}
    n_vcf = n_scan = 0
    for key in read_vcf(args.vcf):
        n_vcf += 1
        sources.setdefault(key, set()).add("P")
    for key in read_indels(args.indels):
        n_scan += 1
        sources.setdefault(key, set()).add("S")

    labels = read_labels(args.annotated, "Chr", ("End", "Start"), "Ref", "Alt")
    labels.update({k: v for k, v in
                   read_labels(args.indels_annotated, "chrom", ("pos",), "ref", "alt").items()
                   if k not in labels})

    n_labelled = 0
    ordered = sorted(sources, key=lambda k: (natural_chrom_key(k[0]), k[1], k[2], k[3]))
    with open(args.out, "w") as out:
        for key in ordered:
            label = labels.get(key, "NA")
            if label != "NA":
                n_labelled += 1
            src = "".join(s for s in "PS" if s in sources[key])
            out.write("%s\t%d\t%s\t%s\t%s;src=%s\n" % (key[0], key[1], key[2], key[3], label, src))

    both = sum(1 for s in sources.values() if len(s) == 2)
    sys.stderr.write(
        "%s %s: vcf=%d scan=%d both=%d skipped_N=%d candidates=%d labelled=%d -> %s\n"
        % (args.sample, args.track, n_vcf, n_scan, both, SKIPPED_N["count"],
           len(ordered), n_labelled, args.out))


if __name__ == "__main__":
    main()
