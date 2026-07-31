#!/usr/bin/env python3
"""
Panel-wide indel counts from a consensus BAM.

Pisces cannot report indels at MRD frequencies. Its variant Q-score model treats
indels as higher-error than substitutions, so the NPM1 type A insertion in
25NGS1071 (20 reads in 50,462, 0.040%) scores VQ 0 and is dropped at the default
--minvq 20. Lowering the threshold does not help incrementally: minvq 5, 10 and
15 all miss it, and only minvq 0 recovers it, which inflates the call set 22x
(2,640 to 57,521) and turns the caller into a pileup dump. The lowest indel
Pisces does report is 0.128%, against a substitution floor near 0.01%.

Pisces also places insertions ambiguously in tandem repeats. The validated NPM1
result was obtained by counting from the BAM with +/-3 bp positional tolerance
for that reason.

So indels are counted here directly off the CIGARs. No caller, no Q-score gate.

Coordinates are VCF-anchored: POS is the base before the event and REF and ALT
both carry it, so output joins to a diagnostic VCF and to the indel blocklist.

Panel coordinates come from BED columns 2-3. Column 4 is a legacy hg19 label and
is not read.

Output:
    {sample}.{track}.indels.tsv

    scan_indels.py --sample S --bam S.dcs.sc.sorted.bam --track dcs \\
        --bed panel.bed --ref genome.fa \\
        --indel-blocklist indel_blocklist.tsv --mask artifact_mask.bed \\
        --out S.dcs.indels.tsv
"""

import argparse
import os
import sys
from collections import defaultdict

import pysam

COLS = ["sample", "track", "chrom", "pos", "ref", "alt", "type", "length",
        "count", "depth", "vaf_pct", "fwd", "rev", "strand_frac",
        "blocklist_n", "masked"]


def load_regions(bed):
    out = []
    for line in open(bed):
        if line.startswith(("#", "track")):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def load_blocklist(path):
    d = {}
    if not path:
        return d
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5:
            try:
                d[(f[0], int(f[1]), f[2].upper(), f[3].upper())] = int(f[4])
            except ValueError:
                pass
    return d


def load_mask(path):
    s = set()
    if not path:
        return s
    for line in open(path):
        f = line.split()
        if len(f) >= 3:
            for p in range(int(f[1]) + 1, int(f[2]) + 1):
                s.add((f[0], p))
    return s


def scan(bam, regions, panel_pos, min_count):
    """Tally indels from CIGARs.

    pysam reference_start is 0-based. After consuming R reference bases the
    cursor sits at 0-based R, so an I or D there is anchored on 0-based R-1,
    which is 1-based R. VCF POS is therefore the cursor value itself.
    """
    tally = defaultdict(lambda: [0, 0, 0])          # count, fwd, rev
    seen = set()
    for chrom, start, end in regions:
        try:
            it = bam.fetch(chrom, start, end)
        except ValueError:
            continue
        for r in it:
            key = (r.query_name, r.reference_start, r.flag)
            if key in seen:
                continue
            seen.add(key)
            if r.is_unmapped or r.cigartuples is None:
                continue
            rev = r.is_reverse
            refpos = r.reference_start
            qpos = 0
            for op, ln in r.cigartuples:
                if op in (0, 7, 8):                       # M / = / X
                    refpos += ln
                    qpos += ln
                elif op == 1:                             # I
                    if (chrom, refpos) in panel_pos:
                        seq = (r.query_sequence[qpos:qpos + ln].upper()
                               if r.query_sequence else "")
                        t = tally[(chrom, refpos, "INS", ln, seq)]
                        t[0] += 1
                        t[2 if rev else 1] += 1
                    qpos += ln
                elif op == 2:                             # D
                    if (chrom, refpos) in panel_pos:
                        t = tally[(chrom, refpos, "DEL", ln, "")]
                        t[0] += 1
                        t[2 if rev else 1] += 1
                    refpos += ln
                elif op == 3:                             # N
                    refpos += ln
                elif op == 4:                             # S
                    qpos += ln
    return {k: v for k, v in tally.items() if v[0] >= min_count}


def depth_at(bam, chrom, pos, min_bq):
    """Total base coverage at a 1-based position."""
    try:
        cov = bam.count_coverage(chrom, pos - 1, pos,
                                 quality_threshold=min_bq)
    except ValueError:
        return 0
    return int(sum(c[0] for c in cov))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True, help="consensus BAM")
    ap.add_argument("--track", required=True, choices=["sscs", "dcs"])
    ap.add_argument("--bed", required=True, help="panel BED")
    ap.add_argument("--ref", required=True, help="reference FASTA")
    ap.add_argument("--out", required=True)
    ap.add_argument("--indel-blocklist", help="annotation only")
    ap.add_argument("--mask", help="artifact mask BED, annotation only")
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-alt", type=int, default=2,
                    help="minimum supporting reads to emit a row. Default 2")
    args = ap.parse_args()

    regions = load_regions(args.bed)
    if not regions:
        sys.exit("FATAL: no regions in BED")
    blocklist = load_blocklist(args.indel_blocklist)
    masked = load_mask(args.mask)

    fa = pysam.FastaFile(args.ref)
    bam = pysam.AlignmentFile(args.bam, "rb")

    panel_pos = set()
    for chrom, start, end in regions:
        for p in range(start, end + 1):
            panel_pos.add((chrom, p))

    sys.stderr.write(f"{len(regions)} probes, {len(panel_pos):,} positions | "
                     f"blocklist {len(blocklist)} | mask {len(masked)}\n")

    tally = scan(bam, regions, panel_pos, args.min_alt)

    rows = []
    for (chrom, cur, kind, ln, seq), (n, fwd, rev) in sorted(tally.items()):
        pos = cur                                        # VCF POS, 1-based
        try:
            anchor = fa.fetch(chrom, pos - 1, pos).upper()
        except (KeyError, ValueError):
            anchor = "N"
        if kind == "INS":
            ref, alt = anchor, anchor + seq
        else:
            try:
                deleted = fa.fetch(chrom, pos, pos + ln).upper()
            except (KeyError, ValueError):
                deleted = "N" * ln
            ref, alt = anchor + deleted, anchor
        depth = depth_at(bam, chrom, pos, args.min_bq)
        rows.append({
            "sample": args.sample, "track": args.track.upper(),
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "type": kind, "length": ln, "count": n, "depth": depth,
            "vaf_pct": f"{100*n/depth:.4f}" if depth else "NA",
            "fwd": fwd, "rev": rev,
            "strand_frac": f"{max(fwd,rev)/n:.3f}" if n else "NA",
            "blocklist_n": blocklist.get((chrom, pos, ref, alt), 0),
            "masked": "yes" if (chrom, pos) in masked else "no"})

    with open(args.out, "w") as out:
        out.write("\t".join(COLS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in COLS) + "\n")

    bam.close()
    fa.close()

    ins = sum(1 for r in rows if r["type"] == "INS")
    sys.stderr.write(f"{args.sample} {args.track}: {len(rows)} indels "
                     f"({ins} INS, {len(rows)-ins} DEL) with >= {args.min_alt} reads\n")
    sys.stderr.write(f"wrote {args.out}\n")


if __name__ == "__main__":
    main()
