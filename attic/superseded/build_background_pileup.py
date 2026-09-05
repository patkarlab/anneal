#!/usr/bin/env python3
"""
build_background_pileup.py - caller-independent per-site-per-alt SNV background.

For every panel position, for every non-reference base, collect the alt fraction
across N normal SSCS BAMs -> background mean, SD, and a POOLED error rate for the
binomial test. Uses pysam count_coverage with a base-quality floor matching the
Pisces --minbq setting, so the background sees reads the same way the caller does.

HARD RULES honored:
  - No hardcoded coordinates: positions come from the panel BED at runtime.
  - BED columns 2-3 ONLY (col 4 is a legacy hg19 label) for genomic coords.
  - Blind: background is computed at EVERY panel site/alt; no diagnosis list.
  - No scipy. pysam + stdlib only.

Output TSV (one row per panel-position x non-ref alt):
  chrom  pos  ref  alt  n_normals  mean_vaf  sd_vaf  max_vaf
  pooled_alt  pooled_depth  bg_error_rate  germline_flag
"""
import argparse, sys, math
import pysam

BASES = "ACGT"  # pysam count_coverage array order is A,C,G,T

def read_bed_regions(path):
    """(chrom, start0, end) from cols 1-3 only. Col 4 ignored on purpose."""
    regions = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            regions.append((f[0], int(f[1]), int(f[2])))
    return regions

def sample_sd(vals, mean):
    n = len(vals)
    if n < 2:
        return 0.0
    return math.sqrt(sum((v - mean) ** 2 for v in vals) / (n - 1))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", nargs="+", required=True,
                    help="Normal SSCS BAM paths (the biological negative controls)")
    ap.add_argument("--bed", required=True, help="Panel BED (cols 2-3 = hg38 coords)")
    ap.add_argument("--ref", required=True, help="Reference FASTA (indexed)")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--min-bq", type=int, default=30,
                    help="Base-quality floor; match Pisces --minbq (default 30)")
    ap.add_argument("--germline-vaf", type=float, default=0.30,
                    help="Flag site/alt as germline if mean normal VAF exceeds this (default 0.30)")
    args = ap.parse_args()

    fasta = pysam.FastaFile(args.ref)
    bams = [pysam.AlignmentFile(b) for b in args.bams]
    regions = read_bed_regions(args.bed)
    sys.stderr.write(f"normals: {len(bams)} | panel regions: {len(regions)} | min_bq={args.min_bq}\n")

    cols = ["chrom","pos","ref","alt","n_normals","mean_vaf","sd_vaf","max_vaf",
            "pooled_alt","pooled_depth","bg_error_rate","germline_flag"]
    n_rows = 0
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for chrom, start, end in regions:
            width = end - start
            if width <= 0:
                continue
            ref_seq = fasta.fetch(chrom, start, end).upper()
            # one count_coverage per region per BAM (fast); arrays are length=width
            per_bam = [b.count_coverage(chrom, start, end, quality_threshold=args.min_bq)
                       for b in bams]
            for j in range(width):
                ref_base = ref_seq[j]
                if ref_base not in BASES:
                    continue
                pos1 = start + j + 1  # 1-based VCF coordinate
                for alt in BASES:
                    if alt == ref_base:
                        continue
                    ai = BASES.index(alt)
                    fracs, pooled_alt, pooled_dp = [], 0, 0
                    for counts in per_bam:
                        a, c, g, t = (counts[0][j], counts[1][j], counts[2][j], counts[3][j])
                        dp = a + c + g + t
                        if dp == 0:
                            continue
                        altc = counts[ai][j]
                        fracs.append(altc / dp)
                        pooled_alt += altc
                        pooled_dp += dp
                    if not fracs:
                        continue
                    mean = sum(fracs) / len(fracs)
                    sd = sample_sd(fracs, mean)
                    bg_rate = pooled_alt / pooled_dp if pooled_dp else 0.0
                    germ = 1 if mean > args.germline_vaf else 0
                    out.write("\t".join(str(x) for x in [
                        chrom, pos1, ref_base, alt, len(fracs),
                        f"{mean:.8f}", f"{sd:.8f}", f"{max(fracs):.8f}",
                        pooled_alt, pooled_dp, f"{bg_rate:.3e}", germ]) + "\n")
                    n_rows += 1
    sys.stderr.write(f"wrote {n_rows} site/alt background rows -> {args.out}\n")

if __name__ == "__main__":
    main()
