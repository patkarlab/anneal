#!/usr/bin/env python3
"""
patch_score_multianno.py

Teaches scripts/error_model/score_multianno.py to read the Waalkes-format beta
matrices (results_bnc/beta_matrix_SSCS.txt, beta_matrix_DCS.txt) instead of the
long-format error model.

The matrix layout is:

    chr  pos  alpha A  beta A  alpha T  beta T  alpha G  beta G  alpha C  beta C

There is no reference-base column, so entries are keyed on (chrom, pos, alt)
rather than (chrom, pos, ref, alt). Alpha/beta are per alt base regardless of
reference, so nothing is lost. The 3-tuple key cannot collide with the
long-format 4-tuple key, so both model types coexist.

Substitutions only. Indels are unaffected and continue through the blocklist.

Usage:
    python patch_score_multianno.py --file scripts/error_model/score_multianno.py
    python patch_score_multianno.py --file ... --dry-run
"""

import argparse
import os
import shutil
import sys
import time

LOADER = '''

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
            line = line.rstrip("\\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\\t")
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

'''

EDITS = [
    # 1. import: keep load_model available for the long format
    (
        "from apply_sscs_error_model import (\n"
        "    beta_cdf, betabinom_sf, count_at, count_indel,\n"
        "    filter_patient_bam, load_model, load_mask, load_indel_blocklist)",

        "from apply_sscs_error_model import (\n"
        "    beta_cdf, betabinom_sf, count_at, count_indel,\n"
        "    filter_patient_bam, load_model, load_mask, load_indel_blocklist)"
        + LOADER.rstrip("\n"),
    ),

    # 2. lookup: try the matrix key first, then fall through unchanged
    (
        '        sub = f"{ref}>{alt}"\n'
        "        if (chrom, pos, ref, alt) in sites:\n"
        "            alpha, beta = sites[(chrom, pos, ref, alt)][:2]\n"
        "        elif sub in priors:\n"
        "            alpha, beta = priors[sub]\n"
        "        else:\n"
        "            alpha = beta = None",

        '        sub = f"{ref}>{alt}"\n'
        "        if (chrom, pos, alt) in sites:\n"
        "            alpha, beta = sites[(chrom, pos, alt)]\n"
        "        elif (chrom, pos, ref, alt) in sites:\n"
        "            alpha, beta = sites[(chrom, pos, ref, alt)][:2]\n"
        "        elif sub in priors:\n"
        "            alpha, beta = priors[sub]\n"
        "        else:\n"
        "            alpha = beta = None",
    ),

    # 3. CLI flag
    (
        '    ap.add_argument("--model", required=True, help="Error model TSV")',

        '    ap.add_argument("--model", required=True,\n'
        '                    help="beta_matrix_SSCS.txt / beta_matrix_DCS.txt, "\n'
        '                         "or a long-format error model TSV")\n'
        '    ap.add_argument("--model-format", choices=["waalkes", "long"],\n'
        '                    default="waalkes",\n'
        '                    help="waalkes = per-position alpha/beta matrix "\n'
        '                         "(default). long = legacy per-site-per-alt TSV")',
    ),

    # 4. dispatch the loader
    (
        "    sites, priors = load_model(args.model)",

        "    if args.model_format == \"waalkes\":\n"
        "        sites = load_beta_matrix(args.model)\n"
        "        priors = {}\n"
        "        sys.stderr.write(\n"
        "            f\"Loaded {len(sites)} site-substitutions from \"\n"
        "            f\"{args.model}\\n\")\n"
        "    else:\n"
        "        sites, priors = load_model(args.model)",
    ),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--file", required=True, help="path to score_multianno.py")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.file):
        sys.exit(f"FATAL: not found: {args.file}")

    with open(args.file) as fh:
        src = fh.read()

    if "load_beta_matrix" in src:
        sys.exit("Already patched (load_beta_matrix present). Nothing to do.")

    # verify every anchor is present and unique before touching anything
    for i, (old, _) in enumerate(EDITS, 1):
        n = src.count(old)
        if n != 1:
            sys.stderr.write(f"\nFATAL: edit {i} anchor found {n} times, expected 1.\n")
            sys.stderr.write("Anchor was:\n" + old[:300] + "\n")
            sys.exit("Refusing to patch. The file differs from what was expected.")
        print(f"edit {i}: anchor found")

    for old, new in EDITS:
        src = src.replace(old, new, 1)

    if args.dry_run:
        print("\nDry run. All 4 anchors matched. Nothing written.")
        return

    backup = f"{args.file}.bak_{time.strftime('%Y%m%d_%H%M%S')}"
    shutil.copy2(args.file, backup)
    with open(args.file, "w") as fh:
        fh.write(src)

    print(f"\nPatched  : {args.file}")
    print(f"Backup   : {backup}")
    print("\nVerify:  python " + args.file + " --help | head -30")


if __name__ == "__main__":
    main()
