#!/usr/bin/env python3
"""
score_multianno.py

Post-annotation scoring stage. Reads the ANNOVAR .hg38_multianno.txt that feeds
filter_variants.py, scores each called variant against the SSCS error model
(substitutions) or the indel blocklist (indels), and writes the same table with
scoring columns appended, as a .scored.txt the filter then thresholds on.

Allele representation: ANNOVAR rewrites indels into '-' form in its Ref/Alt
columns, which does not match the matrix. The original anchored VCF columns are
preserved in the Otherinfo block, so this stage locates the VCF FORMAT field
(GT:ALT:TOT:FRAC) and reads CHROM/POS/REF/ALT by offset from it, ignoring
ANNOVAR's columns entirely. Rows whose anchored REF or ALT contain N, or whose
ALT is ANNOVAR's '0', are skipped (base-masking no-calls, not real alleles).

Counts are recomputed from the patient SSCS BAM, NM-filtered the same way as the
model was built, so the p-value is measured consistently. Substitutions get both
the pbeta and beta-binomial p-values and a strand-bias flag; indels get the
blocklist verdict and the alt-count/VAF gate. The pass column reflects the
chosen test and the cutoff, which should be calibrated on the negative controls
rather than left at the 0.005 default.

This stage reuses the functions in apply_sscs_error_model.py, so that script must
sit in the same directory.
"""

import argparse
import os
import sys
import pysam

from apply_sscs_error_model import (
    beta_cdf, betabinom_sf, count_at, count_indel,
    filter_patient_bam, load_model, load_mask, load_indel_blocklist)

def load_beta_matrix(path):
    """Read a Waalkes-format beta matrix into {(chrom, pos, alt): (alpha, beta)}.

    Header: chr, pos, then alpha/beta pairs for A, T, G, C in that order.
    Rows whose alpha/beta cannot be parsed are skipped rather than aborting the
    run, so a partially written matrix degrades to no_model instead of crashing.
    """
    order = ("A", "T", "G", "C")
    table = {}
    with open(path) as fh:
        header = fh.readline()
        if not header.lower().startswith("chr"):
            fh.seek(0)
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 10:
                continue
            chrom = f[0]
            try:
                pos = int(f[1])
            except ValueError:
                continue
            for i, base in enumerate(order):
                try:
                    alpha = float(f[2 + 2 * i])
                    beta = float(f[3 + 2 * i])
                except (ValueError, IndexError):
                    continue
                if alpha > 0 and beta > 0:
                    table[(chrom, pos, base)] = (alpha, beta)
    return table

BASES = ("A", "C", "G", "T")

NEW_COLS = ["anneal_ref", "anneal_alt", "anneal_alt_count", "anneal_depth",
            "anneal_vaf", "anneal_strand_frac", "anneal_strand_flag",
            "anneal_pbeta_p", "anneal_betabinom_p", "anneal_call", "anneal_pass"]


def fmt(x):
    return "NA" if x is None else f"{x:.3e}"


def find_format_index(fields):
    """Locate the VCF FORMAT field (GT:ALT:TOT:FRAC) in a multianno row, scanning
    from the right. CHROM/POS/REF/ALT then sit at known offsets before it."""
    for i in range(len(fields) - 1, -1, -1):
        if fields[i].startswith("GT:") and "FRAC" in fields[i]:
            return i
    return None


def score_call(bam, chrom, pos, ref, alt, sites, priors, masked, blocklist, args):
    """Return the list of NEW_COLS values for one anchored allele."""
    if ("N" in ref or "N" in alt or ref in ("", "0", ".")
            or alt in ("", "0", ".")):
        return [ref, alt, "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "skipped_n_allele", "False"]

    is_snv = (len(ref) == 1 and len(alt) == 1
              and ref in BASES and alt in BASES)

    if is_snv:
        counts, strand = count_at(bam, chrom, pos - 1, args.min_bq)
        depth = sum(counts.values())
        ac = counts.get(alt, 0)
        fwd, rev = strand.get(alt, [0, 0])
        vaf = ac / depth if depth else 0.0
        sfrac = (max(fwd, rev) / ac) if ac else 0.0
        strand_flag = ac >= args.strand_min_alt and sfrac >= args.strand_thresh

        sub = f"{ref}>{alt}"
        if (chrom, pos, alt) in sites:
            alpha, beta = sites[(chrom, pos, alt)]
        elif (chrom, pos, ref, alt) in sites:
            alpha, beta = sites[(chrom, pos, ref, alt)][:2]
        elif sub in priors:
            alpha, beta = priors[sub]
        else:
            alpha = beta = None

        if depth and alpha is not None:
            pb = 1.0 - beta_cdf(vaf, alpha, beta)
            bb = betabinom_sf(ac, depth, alpha, beta)
        else:
            pb = bb = None

        if (chrom, pos) in masked:
            call = "no_call_masked"
        elif depth < args.min_depth:
            call = "no_call_lowdepth"
        elif strand_flag:
            call = "reject_strand"
        elif alpha is None:
            call = "no_model"
        else:
            if args.test == "pbeta":
                passed = pb <= args.alpha_level
            elif args.test == "betabinom":
                passed = bb <= args.alpha_level
            else:
                passed = pb <= args.alpha_level and bb <= args.alpha_level
            call = "call" if passed else "background"

        return [ref, alt, str(ac), str(depth), f"{vaf:.3e}", f"{sfrac:.3f}",
                str(strand_flag), fmt(pb), fmt(bb), call, str(call == "call")]

    # indel
    if len(alt) > len(ref):
        support, depth = count_indel(bam, chrom, pos - 1, "INS",
                                     len(alt) - len(ref), args.indel_shimmer)
    elif len(ref) > len(alt):
        support, depth = count_indel(bam, chrom, pos - 1, "DEL",
                                     len(ref) - len(alt), args.indel_shimmer)
    else:
        return [ref, alt, "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "unsupported_mnv", "False"]

    vaf = support / depth if depth else 0.0
    n_ctrl = blocklist.get((chrom, pos, ref, alt), 0)
    if (chrom, pos) in masked:
        call = "no_call_masked"
    elif n_ctrl >= args.indel_min_controls:
        call = f"reject_blocklist(n={n_ctrl})"
    elif depth < args.min_depth:
        call = "no_call_lowdepth"
    elif support < args.min_indel_alt or vaf < args.min_indel_vaf:
        call = "background"
    else:
        call = "call"

    return [ref, alt, str(support), str(depth), f"{vaf:.3e}", "NA", "NA",
            "NA", "NA", call, str(call == "call")]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--multianno", required=True,
                    help="ANNOVAR .hg38_multianno.txt feeding filter_variants.py")
    ap.add_argument("--bam", required=True, help="Patient SSCS BAM (unfiltered)")
    ap.add_argument("--ref", required=True, help="Reference FASTA (unmasked)")
    ap.add_argument("--model", required=True,
                    help="beta_matrix_SSCS.txt / beta_matrix_DCS.txt, "
                         "or a long-format error model TSV")
    ap.add_argument("--model-format", choices=["waalkes", "long"],
                    default="waalkes",
                    help="waalkes = per-position alpha/beta matrix "
                         "(default). long = legacy per-site-per-alt TSV")
    ap.add_argument("--out", required=True, help="Output scored table (.scored.txt)")
    ap.add_argument("--mask", help="Artifact mask BED")
    ap.add_argument("--indel-blocklist", help="Indel recurrence TSV")
    ap.add_argument("--filtered-bam", help="Where to write the NM-filtered BAM")
    ap.add_argument("--skip-filter", action="store_true",
                    help="Treat --bam as already NM-filtered")
    ap.add_argument("--max-nm", type=int, default=10)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--strand-thresh", type=float, default=0.90)
    ap.add_argument("--strand-min-alt", type=int, default=10)
    ap.add_argument("--alpha-level", type=float, default=0.005,
                    help="P-value cutoff for the pass column; calibrate on the "
                         "negative controls. Default 0.005")
    ap.add_argument("--test", choices=["pbeta", "betabinom", "both"],
                    default="betabinom")
    ap.add_argument("--indel-min-controls", type=int, default=2)
    ap.add_argument("--min-indel-alt", type=int, default=3)
    ap.add_argument("--min-indel-vaf", type=float, default=1e-4)
    ap.add_argument("--indel-shimmer", type=int, default=3)
    args = ap.parse_args()

    if args.skip_filter:
        filt_bam = args.bam
    else:
        filt_bam = args.filtered_bam or (os.path.splitext(args.out)[0]
                                         + f".patient.nm{args.max_nm}.bam")
        sys.stderr.write(f"Filtering patient BAM (NM<={args.max_nm}) -> {filt_bam}\n")
        filter_patient_bam(args.bam, args.ref, filt_bam, args.max_nm)

    if args.model_format == "waalkes":
        sites = load_beta_matrix(args.model)
        priors = {}
        sys.stderr.write(
            f"Loaded {len(sites)} site-substitutions from "
            f"{args.model}\n")
    else:
        sites, priors = load_model(args.model)
    masked = load_mask(args.mask)
    blocklist = load_indel_blocklist(args.indel_blocklist)
    bam = pysam.AlignmentFile(filt_bam, "rb")

    n_rows = n_scored = n_skipped = 0
    with open(args.multianno) as inp, open(args.out, "w") as out:
        header = inp.readline().rstrip("\n")
        out.write(header + "\t" + "\t".join(NEW_COLS) + "\n")
        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue
            n_rows += 1
            fields = line.split("\t")
            f_idx = find_format_index(fields)
            if f_idx is None or f_idx < 8:
                out.write(line + "\t" + "\t".join(["NA"] * len(NEW_COLS)) + "\n")
                n_skipped += 1
                continue
            chrom = fields[f_idx - 8]
            ref = fields[f_idx - 5].upper()
            alt = fields[f_idx - 4].upper()
            try:
                pos = int(fields[f_idx - 7])
            except ValueError:
                out.write(line + "\t" + "\t".join(["NA"] * len(NEW_COLS)) + "\n")
                n_skipped += 1
                continue
            sys.stderr.write(f"row {n_scored} {chrom}:{pos} {ref}>{alt}\n"); sys.stderr.flush()
            result = score_call(bam, chrom, pos, ref, alt, sites, priors,
                                masked, blocklist, args)
            out.write(line + "\t" + "\t".join(result) + "\n")
            n_scored += 1

    sys.stderr.write(f"Scored {n_scored} of {n_rows} rows "
                     f"({n_skipped} could not be located); wrote {args.out}\n")


if __name__ == "__main__":
    main()
