#!/usr/bin/env python3
"""
call_mrd_eln.py

Apply the ELN 2021 NGS-MRD positivity definition to a scored sample.

The ELN consensus (Heuser et al., Blood 2021, supplementary "bioinformatics
analysis") defines NGS-MRD positivity with background error correction as:

    VAF > mean_background + 3 * SD_background

where the background is measured per site from negative controls. This script
takes the per-site background thresholds computed from the BNC panel
(per_site_background_threshold.tsv) and applies them to a patient/dilution
.scored.txt, calling each substitution positive if its VAF exceeds its site's
mean+3SD threshold.

Substitutions are called against the per-site background. Indels (e.g. the NPM1
insertion) are not substitution sites and carry NA VAF/threshold here; they are
reported with their existing call from the indel blocklist + gate logic.

No gene-level filtering is applied. The model reports, per variant, whether the
VAF exceeds the per-site background (with the global floor and the strand verdict
from score_multianno). Clinical interpretation -- germline, CHIP, gene relevance --
is left entirely downstream.

Column positions are addressed explicitly (not by NF arithmetic) to avoid the
off-by-one errors that variable column counts can introduce. VAF is column 60.
"""

import argparse
import sys

# 1-based -> 0-based explicit indices in the .scored.txt
C_CHR = 0
C_POS = 1            # ANNOVAR Start
C_FUNC = 5
C_GENE = 6
C_EXONICFUNC = 8
C_VAF = 59           # anneal_vaf  (column 60)
C_ALT_COUNT = 57     # column 58
C_DEPTH = 58         # column 59
C_CALL = 64          # column 65
A_REF = -11          # anneal_ref (anchored)
A_ALT = -10          # anneal_alt (anchored)

KEEP_EXONICFUNC = {
    "nonsynonymous SNV", "stopgain", "stoploss",
    "frameshift insertion", "frameshift deletion",
    "nonframeshift insertion", "nonframeshift deletion", "startloss",
}

OUT_HEADER = [
    "gene", "chrom", "pos", "ref", "alt", "exonic_func",
    "alt_count", "depth", "vaf",
    "bg_mean", "bg_sd", "threshold_3sd",
    "mrd_call", "site_evaluable",
]


def load_thresholds(path):
    """site 'chrom:pos' -> (mean, sd, threshold)."""
    thr = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                thr[f[0]] = (float(f[2]), float(f[3]), float(f[4]))
            except ValueError:
                continue
    return thr


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scored", required=True, help="Sample .scored.txt")
    ap.add_argument("--thresholds", required=True,
                    help="per_site_background_threshold.tsv from the BNC panel")
    ap.add_argument("--out", required=True, help="Output MRD call table")
    ap.add_argument("--max-evaluable-bg", type=float, default=0.05,
                    help="Sites whose mean+3SD threshold exceeds this are "
                         "flagged not evaluable (background too high). Default 0.05")
    ap.add_argument("--include-noncoding", action="store_true")
    ap.add_argument("--positive-only", action="store_true",
                    help="Output only variants called MRD-positive")
    args = ap.parse_args()

    thr = load_thresholds(args.thresholds)

    n_in = n_pos = n_noneval = 0
    with open(args.scored) as inp, open(args.out, "w") as out:
        inp.readline()  # original header
        out.write("\t".join(OUT_HEADER) + "\n")
        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < 66:
                continue

            chrom = f[C_CHR]
            ref = f[A_REF]
            alt = f[A_ALT]
            if ref in ("", "0", ".") or alt in ("", "0", ".") \
                    or "N" in ref or "N" in alt:
                continue

            try:
                pos = int(f[C_POS])
            except ValueError:
                continue

            func = f[C_FUNC]
            exonic_func = f[C_EXONICFUNC]
            if not args.include_noncoding:
                if func != "exonic" and "splicing" not in func:
                    continue
                if func == "exonic" and exonic_func not in KEEP_EXONICFUNC:
                    continue

            n_in += 1
            gene = f[C_GENE]
            is_indel = (len(ref) != len(alt))

            site = f"{chrom}:{pos}"
            if is_indel:
                # indel: no substitution background; carry its existing call
                bg_mean = bg_sd = threshold = "NA"
                vaf = f[C_VAF]
                call = f[C_CALL]  # from blocklist + gate logic
                evaluable = "indel_arm"
                is_pos = (call == "call")
            else:
                vaf_str = f[C_VAF]
                scored_call = f[C_CALL]  # anneal_call from score_multianno (strand/mask/background aware)
                if site in thr and vaf_str != "NA":
                    m, sd, t = thr[site]
                    bg_mean = f"{m:.4e}"
                    bg_sd = f"{sd:.4e}"
                    threshold = f"{t:.4e}"
                    vaf = vaf_str
                    evaluable = "yes" if t <= args.max_evaluable_bg else "no_high_bg"
                    # Positive requires BOTH: VAF clears the ELN per-site threshold,
                    # AND score_multianno did not reject it (strand/mask/background).
                    clears = float(vaf_str) > t
                    if scored_call == "reject_strand":
                        is_pos = False
                        call = "reject_strand"
                    elif not clears:
                        is_pos = False
                        call = "negative"
                    else:
                        is_pos = True
                        call = "positive"
                else:
                    bg_mean = bg_sd = threshold = "NA"
                    vaf = vaf_str
                    evaluable = "no_threshold"
                    is_pos = False
                    call = "no_threshold"

            if evaluable in ("no_high_bg",):
                n_noneval += 1
            if is_pos and evaluable in ("yes", "indel_arm"):
                n_pos += 1

            if args.positive_only and not (
                    is_pos and evaluable in ("yes", "indel_arm")):
                continue

            out.write("\t".join([
                gene, chrom, str(pos), ref, alt, exonic_func,
                f[C_ALT_COUNT], f[C_DEPTH], vaf,
                bg_mean, bg_sd, threshold,
                call, evaluable,
            ]) + "\n")

    sys.stderr.write(
        f"evaluated {n_in} coding variants; "
        f"above per-site background (evaluable): {n_pos}; "
        f"non-evaluable high-bg sites: {n_noneval}\n")


if __name__ == "__main__":
    main()
