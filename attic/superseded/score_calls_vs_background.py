#!/usr/bin/env python3
"""
score_calls_vs_background.py - blind per-site-per-alt MRD scorer.

Tumor evidence comes from a PILEUP of the tumor BAM (pysam count_coverage,
base-quality-filtered to match the background), NOT from a caller's VCF. This is
deliberate: at the detection floor a sensitive caller (Pisces) collapses a real
low-VAF driver onto a gVCF reference line, so its per-alt count is not in the VCF.
Scoring the pileup directly tests every panel position x 3 alts and lets the
per-alt background - not a caller threshold - decide significance.

For each panel position, for each non-ref alt:
  tumor_alt, tumor_depth  := count_coverage(quality>=min_bq) on the tumor BAM
  bg_rate                 := matched (chrom,pos,ref,alt) row in the background table
  p_value                 := P(X >= tumor_alt | n=tumor_depth, p=bg_rate)   [exact binomial]
Significance: -log10 p >= --sig AND not germline-flagged.

A Pisces gVCF may be supplied with --pisces purely to ANNOTATE each row with the
caller's verdict (PASS / filtered / reference-line) for cross-checking. It never
gates what is scored and never supplies the evidence.

HARD RULES: no hardcoded coordinates (panel BED drives positions; col 2-3 only,
col 4 ignored); blind (no diagnosis list); no scipy (math only).
"""
import argparse, sys, math
import pysam

BASES = "ACGT"

# ---------- scipy-free exact binomial upper tail: P(X>=k) = I_p(k, n-k+1) ----------
def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 300, 1.0e-14, 1.0e-300
    qab, qap, qam = a+b, a+1.0, a-1.0
    c = 1.0; d = 1.0 - qab*x/qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0/d; h = d
    for m in range(1, MAXIT+1):
        m2 = 2*m
        aa = m*(b-m)*x/((qam+m2)*(a+m2))
        d = 1.0+aa*d; d = FPMIN if abs(d)<FPMIN else d
        c = 1.0+aa/c; c = FPMIN if abs(c)<FPMIN else c
        d = 1.0/d; h *= d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d = 1.0+aa*d; d = FPMIN if abs(d)<FPMIN else d
        c = 1.0+aa/c; c = FPMIN if abs(c)<FPMIN else c
        d = 1.0/d; de = d*c; h *= de
        if abs(de-1.0) < EPS: break
    return h

def _betai(a, b, x):
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    lb = math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)
    bt = math.exp(lb + a*math.log(x) + b*math.log1p(-x))
    return bt*_betacf(a,b,x)/a if x < (a+1.0)/(a+b+2.0) else 1.0 - bt*_betacf(b,a,1.0-x)/b

def binom_sf_ge(k, n, p):
    if k <= 0: return 1.0
    if k > n:  return 0.0
    return _betai(k, n-k+1, p)

def strand_bias_p(fwd, rev):
    """Two-sided binomial p that alt-strand split deviates from 50/50.
       Low p = strand-biased (artifact-like); high p = balanced (real duplex signal).
       Depth-aware: a 353:56 split (p~1e-53) is caught where a raw ratio would not."""
    n = fwd + rev
    if n == 0:
        return 1.0
    k = min(fwd, rev)
    return min(1.0, 2.0 * binom_sf_ge(n - k, n, 0.5))

def neg_log10(k, n, p):
    pv = binom_sf_ge(k, n, p)
    if pv > 0.0:
        return -math.log10(pv), pv
    logpmf = (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
              + k*math.log(p) + (n-k)*math.log1p(-p))
    return -logpmf/math.log(10.0), 0.0

# ---------- inputs ----------
def read_bed_regions(path):
    regions = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#","track","browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3: continue
            regions.append((f[0], int(f[1]), int(f[2])))   # cols 1-3 only
    return regions

def load_background(path):
    bg = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c:k for k,c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[i["chrom"]], f[i["pos"]], f[i["ref"]], f[i["alt"]])
            bg[key] = dict(n=int(f[i["n_normals"]]),
                           mean=float(f[i["mean_vaf"]]), sd=float(f[i["sd_vaf"]]),
                           pooled_depth=int(f[i["pooled_depth"]]),
                           rate=float(f[i["bg_error_rate"]]),
                           germ=int(f[i["germline_flag"]]))
    return bg

def load_pisces_annotation(path):
    """Optional cross-check only: map (chrom,pos,alt)->FILTER, and ref-line VF presence."""
    ann = {}
    if not path:
        return ann
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8: continue
            chrom, pos, _id, ref, alt, qual, filt = f[:7]
            if alt == ".":
                ann[(chrom, pos, "REFLINE")] = filt          # position was emitted as reference
            else:
                ann[(chrom, pos, alt)] = filt                # explicit variant call
    return ann

def sample_sd(vals, mean):
    n = len(vals)
    return math.sqrt(sum((v-mean)**2 for v in vals)/(n-1)) if n >= 2 else 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tumor-bam", required=True, help="Tumor SSCS consensus BAM (evidence source)")
    ap.add_argument("--background", required=True, help="per_site_per_alt_background.tsv")
    ap.add_argument("--bed", required=True, help="Panel BED (cols 2-3 = hg38 coords)")
    ap.add_argument("--ref", required=True, help="Reference FASTA matching the BAM alignment")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pisces", default=None, help="Optional Pisces gVCF for cross-check annotation only")
    ap.add_argument("--min-bq", type=int, default=30, help="Base-quality floor; match background build")
    ap.add_argument("--sig", type=float, default=6.0, help="Significance cutoff as -log10 p (default 6)")
    ap.add_argument("--bg-floor", type=float, default=1.0e-6,
                    help="Error-rate floor for a site/alt ABSENT from the background table")
    ap.add_argument("--strand-bias-p", type=float, default=1.0e-3,
                    help="Min two-sided strand-split binomial p to accept a call; "
                         "below this the alt is strand-biased (artifact) and rejected (default 1e-3). "
                         "Set to 0 to disable the strand filter.")
    ap.add_argument("--emit", choices=["significant","all"], default="significant",
                    help="Write only significant rows (default) or every tested site/alt")
    args = ap.parse_args()

    fasta = pysam.FastaFile(args.ref)
    tbam  = pysam.AlignmentFile(args.tumor_bam)
    bg    = load_background(args.background)
    ann   = load_pisces_annotation(args.pisces)
    regions = read_bed_regions(args.bed)
    sys.stderr.write(f"tumor={args.tumor_bam}\nbackground site/alts: {len(bg)} | panel regions: {len(regions)} "
                     f"| pisces_annot: {len(ann)} | emit={args.emit}\n")

    cols = ["chrom","pos","ref","alt","tumor_alt","tumor_depth","tumor_vaf_pct",
            "alt_fwd","alt_rev","strand_bias_p",
            "bg_n_normals","bg_mean_vaf_pct","bg_sd_vaf_pct","bg_error_rate","rate_used",
            "germline_flag","pisces_filter","p_value","neg_log10_p","significant"]
    n_tested = n_sig = 0
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for chrom, start, end in regions:
            width = end - start
            if width <= 0: continue
            ref_seq = fasta.fetch(chrom, start, end).upper()
            # strand-resolved coverage: forward and reverse separately (same BQ filter).
            # fwd[base]+rev[base] reconstructs the total count_coverage exactly.
            cov_f = tbam.count_coverage(chrom, start, end, quality_threshold=args.min_bq,
                                        read_callback=lambda r: not r.is_reverse)
            cov_r = tbam.count_coverage(chrom, start, end, quality_threshold=args.min_bq,
                                        read_callback=lambda r: r.is_reverse)
            for j in range(width):
                ref_base = ref_seq[j]
                if ref_base not in BASES: continue
                pos1 = start + j + 1
                fwd = {b: cov_f[k][j] for k, b in enumerate(BASES)}
                rev = {b: cov_r[k][j] for k, b in enumerate(BASES)}
                depth = sum(fwd[b] + rev[b] for b in BASES)
                if depth == 0: continue
                for alt in BASES:
                    if alt == ref_base: continue
                    af, ar = fwd[alt], rev[alt]
                    talt = af + ar
                    n_tested += 1
                    key = (chrom, str(pos1), ref_base, alt)
                    b = bg.get(key)
                    if b is not None:
                        floor = 0.5/(b["pooled_depth"]+1)
                        rate = max(b["rate"], floor)
                        bn, bmean, bsd, brate, germ = b["n"], b["mean"], b["sd"], b["rate"], b["germ"]
                    else:
                        rate = args.bg_floor
                        bn, bmean, bsd, brate, germ = -1, -1.0, -1.0, -1.0, 0
                    nlp, pv = neg_log10(talt, depth, rate)
                    sb_p = strand_bias_p(af, ar)
                    # significance: above background AND strand-balanced AND not germline
                    above_bg     = nlp >= args.sig
                    strand_ok    = (args.strand_bias_p <= 0.0) or (sb_p >= args.strand_bias_p)
                    sig = 1 if (above_bg and strand_ok and germ != 1) else 0
                    # pisces cross-check annotation (does NOT affect scoring)
                    pf = ann.get((chrom, str(pos1), alt))
                    if pf is None and (chrom, str(pos1), "REFLINE") in ann:
                        pf = "refline"
                    pf = pf if pf is not None else "NA"
                    if sig: n_sig += 1
                    if args.emit == "all" or sig:
                        out.write("\t".join(str(x) for x in [
                            chrom, pos1, ref_base, alt, talt, depth,
                            f"{100*talt/depth:.4f}", af, ar, f"{sb_p:.3e}",
                            bn, f"{bmean*100:.4f}" if bmean>=0 else -1,
                            f"{bsd*100:.4f}" if bsd>=0 else -1,
                            f"{brate:.3e}" if brate>=0 else -1, f"{rate:.3e}",
                            germ, pf, f"{pv:.3e}", f"{nlp:.2f}", sig]) + "\n")
    sys.stderr.write(f"tested {n_tested} site/alts | {n_sig} significant (-log10p>={args.sig}) -> {args.out}\n")

if __name__ == "__main__":
    main()
