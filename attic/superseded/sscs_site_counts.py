#!/usr/bin/env python3
"""
sscs_site_counts.py

Stage 1 of the SSCS error model: produce a per-site base-count table from a
single SSCS consensus BAM over the panel BED. One table per BNC; the fitter
(Stage 3) aggregates these across all biological negative controls.

Output (TSV, one row per panel position):
    chrom  pos  ref  depth  A  C  G  T

  - pos is 1-based.
  - depth is the count of A+C+G+T bases passing the base-quality and read
    filters (deletions, N, and indels are not counted here; indels are handled
    separately by an indel blocklist).
  - Counts come from pysam.count_coverage, which does NOT cap at 8000, so the
    ~30k SSCS depth at these loci is reported in full.

Panel BED note: genomic coordinates are read from columns 2-3 only. Column 4
of AML_MRD_DUPLEX_probes_hg38_sortd.bed carries legacy/hg19 labels and must
never be used as a coordinate.
"""

import argparse
import sys
import pysam


def make_read_filter(min_mapq):
    """Read-level filter applied by count_coverage (in addition to base quality).

    Excludes unmapped, secondary, supplementary, QC-fail and duplicate reads,
    then applies a minimum mapping-quality cutoff. Consensus reads should pass
    these cleanly; the filter is here so the same script is safe on any BAM.
    """
    def keep(read):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_qcfail or read.is_duplicate):
            return False
        if read.mapping_quality < min_mapq:
            return False
        return True
    return keep


def read_bed_intervals(bed_path):
    """Yield (chrom, start, stop) from a BED file using columns 1-3 (0-based,
    half-open). Skips blank lines and track/browser/comment headers."""
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            chrom = f[0]
            start = int(f[1])   # 0-based start
            stop = int(f[2])    # exclusive end  -- never use column 4
            if stop > start:
                yield chrom, start, stop


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="SSCS consensus BAM (indexed)")
    ap.add_argument("--bed", required=True, help="Panel BED (coords from cols 2-3)")
    ap.add_argument("--ref", required=True, help="Reference FASTA (hg38, indexed)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--min-bq", type=int, default=20,
                    help="Minimum consensus base quality to count (default 20)")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="Minimum mapping quality to count (default 0)")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    fasta = pysam.FastaFile(args.ref)
    keep = make_read_filter(args.min_mapq)

    bam_contigs = set(bam.references)
    fasta_contigs = set(fasta.references)

    n_pos = 0
    warned = set()
    with open(args.out, "w") as out:
        out.write("#chrom\tpos\tref\tdepth\tA\tC\tG\tT\n")
        for chrom, start, stop in read_bed_intervals(args.bed):
            if chrom not in bam_contigs or chrom not in fasta_contigs:
                if chrom not in warned:
                    sys.stderr.write(f"WARN: {chrom} not in BAM and/or FASTA; "
                                     f"skipping its intervals\n")
                    warned.add(chrom)
                continue

            # count_coverage returns four arrays (A, C, G, T), one count per
            # position in [start, stop). No 8000 depth cap.
            cov = bam.count_coverage(chrom, start, stop,
                                     quality_threshold=args.min_bq,
                                     read_callback=keep)
            ref_seq = fasta.fetch(chrom, start, stop).upper()

            for i in range(stop - start):
                a, c, g, t = cov[0][i], cov[1][i], cov[2][i], cov[3][i]
                depth = a + c + g + t
                ref_base = ref_seq[i] if i < len(ref_seq) else "N"
                out.write(f"{chrom}\t{start + i + 1}\t{ref_base}\t{depth}\t"
                          f"{a}\t{c}\t{g}\t{t}\n")
                n_pos += 1

    sys.stderr.write(f"Wrote {n_pos} panel positions to {args.out}\n")


if __name__ == "__main__":
    main()
