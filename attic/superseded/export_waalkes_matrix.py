#!/usr/bin/env python3
"""
export_waalkes_matrix.py  (Step 3 of 3)

Reshape the per-site error model into the Waalkes MRD_alpha_beta_matrix layout,
so the methods comparison can run at the matrix level. The Waalkes format is one
row per position with four (alpha, beta) pairs, ordered ALT = A, T, G, C (the
order their calculate_beta_P_values_vcf.pl indexes), and is consumed by their
1 - pbeta(VAF, alpha, beta) <= 0.005 test.

Each non-reference substitution's column is filled from this model's per-site
Beta. The column whose base equals the reference is never queried by the test (a
called variant's ALT is never the reference), so it is filled with the
substitution-class transition prior purely to keep the layout a valid numeric
matrix that their Perl can read.

Input: the model TSV from build_sscs_error_model.py. Output: a tab-separated
matrix with the Waalkes header.
"""

import argparse
import sys
from collections import defaultdict

BASES_ORDER = ("A", "T", "G", "C")          # Waalkes column order
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="Model TSV")
    ap.add_argument("--out", required=True, help="Output Waalkes-format matrix")
    args = ap.parse_args()

    priors = {}
    # position -> {"ref": ref, alt: (alpha, beta)}
    by_pos = defaultdict(dict)

    with open(args.model) as fh:
        for line in fh:
            if line.startswith("## prior"):
                f = line.split()
                priors[f[2]] = (float(f[3].split("=")[1]),
                                float(f[4].split("=")[1]))
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[2], f[3]
            rec = by_pos[(chrom, pos)]
            rec["ref"] = ref
            rec[alt] = (float(f[10]), float(f[11]))

    n = 0
    with open(args.out, "w") as out:
        out.write("chr\tpos\talpha\tbeta\talpha\tbeta\talpha\tbeta\talpha\tbeta\n")
        for (chrom, pos) in sorted(by_pos):
            rec = by_pos[(chrom, pos)]
            ref = rec.get("ref")
            cells = []
            for base in BASES_ORDER:
                if base == ref:
                    # never queried; fill with the transition prior for this ref
                    ab = priors.get(f"{ref}>{TRANSITION.get(ref, 'A')}", (0.25, 18745.0))
                elif base in rec:
                    ab = rec[base]
                else:
                    ab = priors.get(f"{ref}>{base}", (0.25, 18745.0))
                cells.extend([f"{ab[0]:.8g}", f"{ab[1]:.8g}"])
            out.write(f"{chrom}\t{pos}\t" + "\t".join(cells) + "\n")
            n += 1

    sys.stderr.write(f"Wrote {n} positions in Waalkes matrix format to {args.out}\n")


if __name__ == "__main__":
    main()
