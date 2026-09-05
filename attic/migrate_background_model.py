#!/usr/bin/env python3
"""
migrate_background_model.py

Two jobs:

1. Copy the background-model build chain out of results_bnc/ErrorModel/, which
   is gitignored, into scripts/background_model/ so the published pipeline can
   regenerate its own beta matrix.

2. Correct the default error rate in beta_distribution.py.

The correction
--------------
    mean_default_error = default_error_rate / No_of_samples

The mean of N identical rates is the rate, not the rate over N. The floor
therefore sits 8x below the cited literature value: 1/120,000 instead of
1/15,000 for SSCS, 1/1,600,000 instead of 1/200,000 for DCS.

The trigger carries the same units error -- sum_x_error_list is a sum over the
eight controls and is compared against a single-sample rate. Both are fixed
together; fixing only the floor introduces a discontinuity at the threshold.

Rebuild required
----------------
Any existing matrix was produced by the unfixed code. After patching, rebuild
from the existing Pisces gVCFs and re-run validation.

Usage:
    python migrate_background_model.py --root ~/pipelines/anneal --dry-run
    python migrate_background_model.py --root ~/pipelines/anneal
"""

import argparse
import os
import re
import shutil
import sys
import time

FILES = [
    "beta_distribution.py",
    "fill_empty_mips.pl",
    "print_multiple_variants_at_same_location.pl",
    "remove_variants_gtr_20.pl",
    "non-overlapping_bed_regions.py",
    "error_model.md",
]

README = """# Background error model

Site- and substitution-specific error model built from biological negative
controls, after Waalkes et al. Haematologica 2017;102(9):1549-1557.

See `error_model.md` for the full procedure. Summary:

1. Pisces on the 8 BNC consensus BAMs, `--minbq 30 --minvf 0.0001 -c 1`
2. `remove_variants_gtr_20.pl`  -- drop VAF > 0.2 (germline)
3. `print_multiple_variants_at_same_location.pl` -- split multi-allelic
4. `fill_empty_mips.pl`  -- write zero counts at uncalled positions
5. `beta_distribution.py` -- fit per-position Beta, emit the matrix

Run once per consensus track. `default_error_rate` in `beta_distribution.py`
must be set per track: 1/15000 for SSCS, 1/200000 for DCS, both from Wang et al.
NAR 2019 supplementary table 5 (SmallDeep panel).

## Known issues

- `default_error_rate` is a module-level constant, so the two tracks currently
  require two copies of the script. It should become a command-line argument.

- `fill_empty_mips.pl` iterates `$start <= $i <= $end` against a 0-based BED
  start, so each probe beginning a contiguous block emits one extra position at
  its 5' end. 97 of 182 probes are block starts, giving 21,937 rows where the
  panel spans 21,840 positions. The surplus rows carry zero counts and take the
  default, so the model is unaffected, but the matrix cannot be joined to other
  tables on position without dropping them.

- Depth at filled positions is carried forward from the last matched variant
  line rather than being the true depth at that position. Error rates are
  unaffected (0/anything = 0) but the depth column is not usable, which rules
  out any depth-aware test built from this file.

- Pisces applies `MinimumVariantQScore` (default 20) beneath the requested
  `--minvf`, censoring roughly everything below 7-12 alt reads. Positions below
  that are recorded as zeros and take the default. The default is therefore
  doing the work at most positions, and it is a literature constant rather than
  a measurement from this panel.
"""


def patch_beta(src):
    """Return (patched_text, list_of_changes) or (None, reason)."""
    changes = []

    if "mean_default_error = default_error_rate\n" in src and \
       "/ No_of_samples" not in src.split("mean_default_error")[1][:40]:
        return None, "already patched"

    old_floor = "mean_default_error = default_error_rate / No_of_samples"
    if old_floor not in src:
        return None, "floor assignment not found"
    src = src.replace(old_floor, "mean_default_error = default_error_rate", 1)
    changes.append("floor: removed / No_of_samples")

    pat = re.compile(r"sum_([atgc])_error_list(\s*)<(\s*)default_error_rate")
    src, n = pat.subn(
        lambda m: f"sum_{m.group(1)}_error_list / No_of_samples"
                  f"{m.group(2)}<{m.group(3)}default_error_rate", src)
    if n != 4:
        return None, f"expected 4 trigger comparisons, found {n}"
    changes.append(f"trigger: divided by No_of_samples in {n} comparisons")

    return src, changes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.path.expanduser("~/pipelines/anneal"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = os.path.abspath(os.path.expanduser(args.root))
    srcdir = os.path.join(root, "results_bnc", "ErrorModel")
    dstdir = os.path.join(root, "scripts", "background_model")
    stamp = time.strftime("%Y%m%d_%H%M%S")

    if not os.path.isdir(srcdir):
        sys.exit(f"FATAL: not found: {srcdir}")

    print(f"source : {srcdir}")
    print(f"target : {dstdir}")
    print(f"mode   : {'DRY RUN' if args.dry_run else 'APPLY'}\n")

    print("[1] copy build chain into the repo")
    for f in FILES:
        s = os.path.join(srcdir, f)
        if not os.path.exists(s):
            print(f"  MISSING  {f}")
            continue
        print(f"  COPY     {f}")
        if not args.dry_run:
            os.makedirs(dstdir, exist_ok=True)
            shutil.copy2(s, os.path.join(dstdir, f))

    print(f"\n  WRITE    README.md")
    if not args.dry_run:
        os.makedirs(dstdir, exist_ok=True)
        with open(os.path.join(dstdir, "README.md"), "w") as fh:
            fh.write(README)

    print("\n[2] correct default error rate")
    targets = [os.path.join(dstdir, "beta_distribution.py"),
               os.path.join(srcdir, "beta_distribution.py"),
               os.path.join(root, "results_bnc", "ErrorModelDCS",
                            "beta_distribution.py")]

    for t in targets:
        label = os.path.relpath(t, root)
        if not os.path.exists(t):
            if args.dry_run and t.startswith(dstdir):
                print(f"  {label}: will exist after copy, patch deferred")
            else:
                print(f"  {label}: MISSING")
            continue
        patched, info = patch_beta(open(t).read())
        if patched is None:
            print(f"  {label}: SKIP ({info})")
            continue
        print(f"  {label}:")
        for c in info:
            print(f"      {c}")
        if not args.dry_run:
            shutil.copy2(t, f"{t}.bak_{stamp}")
            with open(t, "w") as fh:
                fh.write(patched)

    if args.dry_run:
        print("\nDry run. Nothing written. Re-run without --dry-run to apply.")
        return

    print("\nApplied. The existing matrices were built with the unfixed code "
          "and must be rebuilt before use.")


if __name__ == "__main__":
    main()
