#!/usr/bin/env python3
"""
measure_duplex.py

Reports the two numbers that decide whether the duplex step is working.

  1. Per-base background error rate at Q30, SSCS vs DCS, plus the
     12 substitution classes. Same method and AF cap as the earlier
     run, so results are directly comparable to the pre-patch baseline:

         SSCS  depth 540,854,706  alt 151,435  rate 2.800e-04
         DCS   depth  61,572,753  alt  18,378  rate 2.985e-04

     The DCS rate is the target. Published duplex sits at 1e-7 to 1e-8.
     Anything above ~1e-6 means duplex is still not suppressing error.

  2. Per-strand support in the DCS pool, read from the XA/XB tags added
     by the patch. min(XA, XB) == 1 means one strand contributed a single
     read: that pairing is not duplex evidence.

Usage:
    python3 measure_duplex.py <sample_dir> [--ref PATH] [--bed PATH]

    e.g. python3 measure_duplex.py results_dilution_gpu/25NGS1601-G-Duplex
"""

import argparse
import glob
import os
import subprocess
import sys
from collections import defaultdict

DEFAULT_REF = os.path.expanduser(
    "~/references/hg38_broad/Homo_sapiens_assembly38.masked.fasta")
DEFAULT_BED = "AML_MRD_DUPLEX_probes_hg38_sortd.bed"

# Anything at or above this allele fraction is treated as real signal, not
# background. Matches the pre-patch baseline run.
AF_CAP = 0.05

# Depth below this is too shallow to judge a position.
MIN_DEPTH = 100


def find_bam(sample_dir, kind):
    """Locate the sscs or dcs BAM inside a sample directory."""
    patterns = [
        os.path.join(sample_dir, "consensus", "*.%s.sc.sorted.bam" % kind),
        os.path.join(sample_dir, "consensus", "%s.sc.sorted.bam" % kind),
        os.path.join(sample_dir, "consensus", "*.%s.sorted.bam" % kind),
        os.path.join(sample_dir, "consensus", "%s.sorted.bam" % kind),
        os.path.join(sample_dir, "*.%s.sc.sorted.bam" % kind),
    ]
    for pattern in patterns:
        hits = sorted(glob.glob(pattern))
        if hits:
            return hits[0]
    return None


def parse_pileup_column(pileup):
    """Return (ref_match_count, {alt_base: count}) for one mpileup column.

    The pileup string uses markers that must be stepped over so that
    inserted and deleted sequence is not miscounted as substitutions.
    """
    alts = defaultdict(int)
    ref_matches = 0
    i = 0
    n = len(pileup)
    while i < n:
        c = pileup[i]
        if c == "^":
            # Read start: the next character is the mapping quality.
            i += 2
        elif c == "$":
            i += 1
        elif c in ".,":
            ref_matches += 1
            i += 1
        elif c in "+-":
            # Indel: a length prefix followed by that many sequence bases.
            j = i + 1
            digits = ""
            while j < n and pileup[j].isdigit():
                digits += pileup[j]
                j += 1
            i = j + (int(digits) if digits else 0)
        elif c in "ACGTacgt":
            alts[c.upper()] += 1
            i += 1
        else:
            i += 1
    return ref_matches, alts


def background_rate(bam, ref, bed):
    """Aggregate background depth, alt count, and per-substitution counts."""
    total_depth = 0
    total_alt = 0
    subs = defaultdict(int)

    cmd = ["samtools", "mpileup", "-A", "-Q", "30", "-d", "1000000",
           "-l", bed, "-f", ref, bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)

    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        ref_base = fields[2].upper()
        if ref_base not in "ACGT":
            continue

        ref_matches, alts = parse_pileup_column(fields[4])
        depth = ref_matches + sum(alts.values())
        if depth < MIN_DEPTH:
            continue

        total_depth += depth
        for alt_base, count in alts.items():
            if alt_base == ref_base:
                continue
            if count / depth >= AF_CAP:
                continue          # real signal, not background
            total_alt += count
            subs["%s>%s" % (ref_base, alt_base)] += count

    proc.wait()
    return total_depth, total_alt, subs


def strand_support(bam):
    """Distribution of min(XA, XB) across the DCS pool."""
    hist = defaultdict(int)
    missing_tags = 0

    proc = subprocess.Popen(["samtools", "view", bam],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        xa = xb = None
        for field in fields[11:]:
            if field.startswith("XA:i:"):
                xa = int(field[5:])
            elif field.startswith("XB:i:"):
                xb = int(field[5:])
        if xa is None or xb is None:
            missing_tags += 1
            continue
        hist[min(xa, xb)] += 1
    proc.wait()
    return hist, missing_tags


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_dir")
    parser.add_argument("--ref", default=DEFAULT_REF)
    parser.add_argument("--bed", default=DEFAULT_BED)
    args = parser.parse_args()

    for path, label in ((args.ref, "reference"), (args.bed, "BED")):
        if not os.path.isfile(path):
            sys.exit("ERROR: %s not found at %s" % (label, path))

    print("# sample=%s  Q30  AF_cap=%.2f  min_depth=%d"
          % (args.sample_dir, AF_CAP, MIN_DEPTH))
    print("")

    dcs_bam = None
    for kind in ("sscs", "dcs"):
        bam = find_bam(args.sample_dir, kind)
        if bam is None:
            print("%-5s NO_BAM" % kind)
            continue
        if kind == "dcs":
            dcs_bam = bam

        depth, alt, subs = background_rate(bam, args.ref, args.bed)
        if depth == 0:
            print("%-5s no positions passed the depth filter" % kind)
            continue

        print("%-5s depth=%-14d alt=%-10d rate=%.3e   %s"
              % (kind, depth, alt, alt / depth, os.path.basename(bam)))
        for name, count in sorted(subs.items(), key=lambda kv: -kv[1]):
            print("        %-5s %-10d %.3e" % (name, count, count / depth))
        print("")

    if dcs_bam is None:
        return

    print("# DCS per-strand support: min(XA, XB)")
    hist, missing = strand_support(dcs_bam)
    if missing and not hist:
        print("  XA/XB tags absent. This BAM predates the patch; regenerate")
        print("  stage 1 to get per-strand depths.")
        return

    total = sum(hist.values())
    cumulative = 0
    for support in sorted(hist):
        count = hist[support]
        cumulative += count
        label = "%d" % support if support < 10 else "10+"
        print("  min_strand=%-4s %-10d %5.1f%%   (cum %5.1f%%)"
              % (label, count, 100 * count / total, 100 * cumulative / total))
    if missing:
        print("  untagged reads: %d" % missing)

    single = hist.get(1, 0)
    if single:
        print("")
        print("  %.1f%% of the DCS pool has only one read on its weaker strand."
              % (100 * single / total))
        print("  Set ANNEAL_MIN_READS_PER_STRAND=2 to exclude these.")


if __name__ == "__main__":
    main()
