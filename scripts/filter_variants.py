#!/usr/bin/env python3
"""
filter_variants.py -- Filter annotated duplex variants for clinical reporting.

Adapted from targeted-seq-pipeline/14_variant_filter.py for the Anneal duplex
pipeline. Removes SomaticSeq-specific logic (multi-caller count, verdict,
U2AF1 rescue) since the Rust mpileup caller is the sole variant caller.

Input:  results/{sample}/annotated/{sample}.dcs.annotated.tsv
Output: results/{sample}/annotated/{sample}.dcs.filtered.tsv   (all + Filter)
        results/{sample}/annotated/{sample}.dcs.clinical.tsv   (PASS only)

Workflow:
  1. Investigate ANNOVAR orphan rows (VariantCaller_Count = -1 equivalent)
  2. Deduplicate: same position + same gene keeps best-annotated row
  3. Apply clinical filters
  4. Report HIGH impact variants
"""

import argparse
import logging
import os
import sys

import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

MISSING = "-1"


def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter annotated duplex variants for clinical reporting.",
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Annotated TSV from annotate_variants.py")
    ap.add_argument("-o", "--outdir", default=None,
                    help="Output directory (default: same as input)")
    ap.add_argument("--min-alt-count", type=int, default=3,
                    help="Minimum ALT_COUNT for PASS (default: 3)")
    ap.add_argument("--min-depth", type=int, default=50,
                    help="Minimum Total_Depth for PASS (default: 50)")
    ap.add_argument("--max-pop-af", type=float, default=0.01,
                    help="Max population AF for common polymorphism (default: 0.01)")
    return ap.parse_args()


def safe_float(val, default=np.nan):
    if val is None or str(val).strip() in ("", "-1", ".", "nan", "NA"):
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val, default=-1):
    if val is None or str(val).strip() in ("", "-1", ".", "nan", "NA"):
        return default
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return default


def is_missing(val):
    return val is None or str(val).strip() in ("", "-1", ".", "nan", "NA")


def count_non_missing(row):
    return sum(1 for v in row if not is_missing(v))


# ---------------------------------------------------------------------------
# Step 1: Investigate orphan rows
# ---------------------------------------------------------------------------

def investigate_orphans(df):
    """Report on rows with missing caller info (ANNOVAR-only orphans)."""
    mask = df["ALT_COUNT"].astype(str).str.strip() == "-1"
    n_orphan = mask.sum()
    total = len(df)

    log.info("=" * 70)
    log.info("ORPHAN ROW INVESTIGATION")
    log.info("=" * 70)
    log.info("Total rows: %d", total)
    log.info("Rows with ALT_COUNT = -1: %d (%.1f%%)",
             n_orphan, 100 * n_orphan / total if total > 0 else 0)

    if n_orphan == 0:
        log.info("No orphan rows found.")
        return

    orphan_df = df[mask].copy()
    has_gene = orphan_df["Gene"].astype(str).str.strip() != "-1"
    has_hgvsc = orphan_df["HGVSc"].astype(str).str.strip() != "-1"
    log.info("  Has Gene:  %d / %d", has_gene.sum(), n_orphan)
    log.info("  Has HGVSc: %d / %d", has_hgvsc.sum(), n_orphan)

    ref_is_dash = orphan_df["Ref"].astype(str).str.strip() == "-"
    alt_is_dash = orphan_df["Alt"].astype(str).str.strip() == "-"
    log.info("  Ref = '-' (insertion):  %d", ref_is_dash.sum())
    log.info("  Alt = '-' (deletion):   %d", alt_is_dash.sum())

    log.info("These are ANNOVAR-only orphan indel representations.")
    log.info("They will be resolved during deduplication.")


# ---------------------------------------------------------------------------
# Step 2: Deduplication
# ---------------------------------------------------------------------------

def _is_orphan_indel(row):
    ref = str(row["Ref"]).strip()
    alt = str(row["Alt"]).strip()
    return ref == "-" or alt == "-"


def deduplicate(df):
    """Remove duplicate variants at same position + gene.

    Strategy:
      - ANNOVAR orphan indels (Ref='-' or Alt='-') with no VEP annotation
        are matched to VEP-annotated rows at same Chr:Start:Gene and dropped.
      - Remaining rows deduped by exact allele, keeping best-annotated.
    """
    log.info("")
    log.info("=" * 70)
    log.info("DEDUPLICATION")
    log.info("=" * 70)

    n_before = len(df)
    df = df.copy()

    df["_has_vep"] = df["HGVSc"].apply(lambda v: 0 if is_missing(v) else 1)
    df["_has_depth"] = df["ALT_COUNT"].apply(
        lambda v: 0 if str(v).strip() == "-1" else 1)
    df["_n_filled"] = df.apply(lambda row: count_non_missing(row), axis=1)
    df["_score"] = df["_has_vep"] * 1000 + df["_has_depth"] * 100 + df["_n_filled"]
    df["_is_orphan"] = df.apply(_is_orphan_indel, axis=1)

    # Phase 1: Remove ANNOVAR orphan indels with a VEP partner
    pos_key = (df["Chr"].astype(str) + ":" + df["Start"].astype(str) + ":"
               + df["Gene"].astype(str))
    df["_pos_key"] = pos_key
    vep_pos_keys = set(df.loc[df["_has_vep"] == 1, "_pos_key"].unique())
    orphan_mask = (df["_is_orphan"] & df["_pos_key"].isin(vep_pos_keys)
                   & (df["_has_vep"] == 0))
    n_orphans_removed = orphan_mask.sum()
    df = df[~orphan_mask].copy()
    log.info("Orphan indels with VEP partner removed: %d", n_orphans_removed)

    # Phase 2: Exact allele dedup
    df["_dedup_key"] = (df["Chr"].astype(str) + ":" + df["Start"].astype(str) + ":"
                        + df["Ref"].astype(str) + ":" + df["Alt"].astype(str) + ":"
                        + df["Gene"].astype(str))

    keep_idx = []
    n_groups_deduped = 0
    for key, group in df.groupby("_dedup_key"):
        if len(group) == 1:
            keep_idx.append(group.index[0])
        else:
            n_groups_deduped += 1
            best = group.sort_values("_score", ascending=False).iloc[0]
            keep_idx.append(best.name)

    df_deduped = df.loc[keep_idx].copy()
    helper_cols = ["_has_vep", "_has_depth", "_n_filled", "_score",
                   "_is_orphan", "_pos_key", "_dedup_key"]
    df_deduped.drop(columns=helper_cols, inplace=True)
    df_deduped.reset_index(drop=True, inplace=True)

    n_after = len(df_deduped)
    n_removed = n_before - n_after
    log.info("Exact-allele duplicate groups: %d", n_groups_deduped)
    log.info("Input rows:    %d", n_before)
    log.info("Rows removed:  %d (orphan: %d, exact-dup: %d)",
             n_removed, n_orphans_removed, n_removed - n_orphans_removed)
    log.info("Output rows:   %d", n_after)

    return df_deduped


# ---------------------------------------------------------------------------
# Step 3: Flag overlapping variants
# ---------------------------------------------------------------------------

def flag_overlapping_variants(df):
    """Flag variants at similar positions that may represent the same event."""
    log.info("")
    log.info("=" * 70)
    log.info("OVERLAPPING VARIANT DETECTION")
    log.info("=" * 70)

    df = df.copy()
    df["Dedup_Note"] = ""

    starts = df["Start"].apply(lambda v: safe_int(v, default=-1))
    vafs = df["VAF_pct"].apply(safe_float)
    alt_counts = df["ALT_COUNT"].apply(lambda v: safe_int(v, default=-1))
    genes = df["Gene"].astype(str)
    refs = df["Ref"].astype(str)
    alts = df["Alt"].astype(str)

    def _is_indel(ref, alt):
        return len(ref) != len(alt) or ref == "-" or alt == "-"

    n_flagged = 0
    for i in range(len(df)):
        if starts.iloc[i] < 0:
            continue
        for j in range(i + 1, len(df)):
            if starts.iloc[j] < 0:
                continue
            if genes.iloc[i] != genes.iloc[j]:
                continue
            if abs(starts.iloc[i] - starts.iloc[j]) > 10:
                continue
            if not _is_indel(refs.iloc[i], alts.iloc[i]) and \
               not _is_indel(refs.iloc[j], alts.iloc[j]):
                continue
            vaf_i, vaf_j = vafs.iloc[i], vafs.iloc[j]
            if np.isnan(vaf_i) or np.isnan(vaf_j):
                continue
            if abs(vaf_i - vaf_j) > 2.0:
                continue
            ac_i, ac_j = alt_counts.iloc[i], alt_counts.iloc[j]
            if ac_i > 0 and ac_j > 0:
                ac_mean = (ac_i + ac_j) / 2
                if abs(ac_i - ac_j) / ac_mean > 0.20:
                    continue

            gene = genes.iloc[i]
            note = (f"MANUAL_REVIEW: overlapping variants in {gene}, "
                    f"same VAF, likely single event -- verify in IGV")
            df.at[df.index[i], "Dedup_Note"] = note
            df.at[df.index[j], "Dedup_Note"] = note
            n_flagged += 1

    if n_flagged > 0:
        log.info("Flagged %d variant pair(s) for manual review", n_flagged)
    else:
        log.info("No overlapping variant pairs found.")

    return df


# ---------------------------------------------------------------------------
# Step 4: Filtering
# ---------------------------------------------------------------------------

def apply_filters(df, min_alt_count, min_depth, max_pop_af):
    """Apply clinical filters.

    Filters (priority order):
      COMMON_POLYMORPHISM  Max_AF > max_pop_af
      LOW_IMPACT           IMPACT == "MODIFIER" and not splice region
      LOW_DEPTH            ALT_COUNT < min_alt_count
      LOW_TOTAL_DEPTH      Total_Depth < min_depth
      NO_ANNOTATION        No gene annotation from VEP or ANNOVAR
      PASS                 Everything else
    """
    log.info("")
    log.info("=" * 70)
    log.info("FILTERING")
    log.info("=" * 70)
    log.info("Parameters: min_alt_count=%d  min_depth=%d  max_pop_af=%.4f",
             min_alt_count, min_depth, max_pop_af)

    df = df.copy()

    df["_max_af"] = df["Max_AF"].apply(safe_float)
    df["_alt_count"] = df["ALT_COUNT"].apply(safe_int)
    df["_total_depth"] = df["Total_Depth"].apply(safe_int)
    df["_impact"] = df["IMPACT"].astype(str).str.strip()
    df["_consequence"] = df["Consequence"].astype(str).str.strip()
    df["_gene"] = df["Gene"].astype(str).str.strip()

    filters = []
    for _, row in df.iterrows():
        # Priority 1: Common polymorphism
        if not np.isnan(row["_max_af"]) and row["_max_af"] > max_pop_af:
            filters.append("COMMON_POLYMORPHISM")
            continue

        # Priority 2: Low impact (MODIFIER, not splice)
        if row["_impact"] == "MODIFIER":
            csq = row["_consequence"].lower()
            if "splice" not in csq:
                filters.append("LOW_IMPACT")
                continue

        # Priority 3: Low alt depth
        alt_count = row["_alt_count"]
        if 0 <= alt_count < min_alt_count:
            filters.append("LOW_DEPTH")
            continue

        # Priority 4: Low total depth
        total_depth = row["_total_depth"]
        if 0 <= total_depth < min_depth:
            filters.append("LOW_TOTAL_DEPTH")
            continue

        # Priority 5: No annotation at all
        if row["_gene"] in ("-1", "", "nan") and row["_alt_count"] < 0:
            filters.append("NO_ANNOTATION")
            continue

        filters.append("PASS")

    df["Filter"] = filters

    df.drop(columns=["_max_af", "_alt_count", "_total_depth",
                      "_impact", "_consequence", "_gene"],
            inplace=True)

    counts = df["Filter"].value_counts()
    log.info("Filter results:")
    for filt in ["PASS", "COMMON_POLYMORPHISM", "LOW_IMPACT",
                 "LOW_DEPTH", "LOW_TOTAL_DEPTH", "NO_ANNOTATION"]:
        n = counts.get(filt, 0)
        if n > 0:
            log.info("  %-25s %d", filt, n)
    log.info("  %-25s %d", "TOTAL", len(df))

    return df


# ---------------------------------------------------------------------------
# Step 5: Report
# ---------------------------------------------------------------------------

def report_high_impact(df):
    """Report HIGH impact variants."""
    log.info("")
    log.info("=" * 70)
    log.info("HIGH IMPACT VARIANTS")
    log.info("=" * 70)

    high = df[df["IMPACT"].astype(str).str.strip() == "HIGH"].copy()
    if high.empty:
        log.info("  (none)")
    else:
        for _, row in high.iterrows():
            log.info(
                "  %-8s %-12s %-6s>%-6s  %-8s %-30s %-35s  VAF=%-8s  Filter=%s",
                row["Chr"], row["Start"],
                str(row["Ref"])[:6], str(row["Alt"])[:6],
                row["Gene"],
                str(row["Consequence"])[:30],
                str(row["HGVSp"])[:35],
                row["VAF_pct"],
                row["Filter"],
            )
    log.info("-" * 100)
    return high


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    input_path = args.input

    if not os.path.isfile(input_path):
        log.error("Input file not found: %s", input_path)
        sys.exit(1)

    outdir = args.outdir or os.path.dirname(input_path)
    os.makedirs(outdir, exist_ok=True)

    basename = os.path.splitext(os.path.basename(input_path))[0]
    # Remove .annotated suffix if present
    basename = basename.replace(".annotated", "")

    log.info("=" * 70)
    log.info("filter_variants.py -- Duplex Variant Filtering")
    log.info("=" * 70)
    log.info("Input:   %s", input_path)
    log.info("Output:  %s", outdir)

    df = pd.read_csv(input_path, sep="\t", dtype=str)
    log.info("Loaded %d variants", len(df))

    # Step 1: Investigate orphans
    investigate_orphans(df)

    # Step 2: Deduplicate
    df = deduplicate(df)

    # Step 2b: Flag overlapping variants
    df = flag_overlapping_variants(df)

    # Step 3: Apply filters
    df = apply_filters(df, args.min_alt_count, args.min_depth, args.max_pop_af)

    # Step 4: Report
    report_high_impact(df)

    # Write filtered output (all variants)
    filtered_path = os.path.join(outdir, f"{basename}.filtered.tsv")
    df.to_csv(filtered_path, sep="\t", index=False)
    log.info("")
    log.info("Written: %s (%d variants)", filtered_path, len(df))

    # Write clinical output (PASS only)
    clinical = df[df["Filter"] == "PASS"].copy()
    clinical.sort_values(["Gene", "Chr", "Start"], inplace=True)
    clinical.reset_index(drop=True, inplace=True)
    clinical_path = os.path.join(outdir, f"{basename}.clinical.tsv")
    clinical.to_csv(clinical_path, sep="\t", index=False)
    log.info("Written: %s (%d PASS variants)", clinical_path, len(clinical))

    # Summary
    log.info("")
    log.info("=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)
    log.info("Total input:     %d", pd.read_csv(input_path, sep="\t").shape[0])
    log.info("After dedup:     %d", len(df))
    for filt in ["PASS", "COMMON_POLYMORPHISM", "LOW_IMPACT",
                 "LOW_DEPTH", "LOW_TOTAL_DEPTH", "NO_ANNOTATION"]:
        n = (df["Filter"] == filt).sum()
        if n > 0:
            log.info("  %-25s %d", filt, n)
    log.info("Clinical (PASS): %d", len(clinical))


if __name__ == "__main__":
    main()
