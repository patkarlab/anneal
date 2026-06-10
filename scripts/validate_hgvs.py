#!/usr/bin/env python3
"""
validate_hgvs.py -- Validate HGVS nomenclature via local VariantValidator.

Adapted from targeted-seq-pipeline/17_variant_validator.py for the Anneal
duplex pipeline. Operates on the clinical or filtered TSV from
filter_variants.py.

For each variant with an HGVSc value, queries the local VariantValidator
REST API to validate and correct HGVS nomenclature. Replaces original
HGVSc/HGVSp/HGVSg with validated versions and adds:
  VV_Transcript     Reference transcript used
  VV_Valid          True/False
  VV_Warnings       Validation warnings or correction notes

Usage:
    python3 validate_hgvs.py \\
        -i results/SAMPLE/annotated/SAMPLE.dcs.clinical.tsv \\
        -o results/SAMPLE/annotated/

Requires:
  - Local VariantValidator Docker container running on localhost:5001
"""

import argparse
import logging
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULT_VV_URL = "http://localhost:5001"
GENOME_BUILD = "GRCh38"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate HGVS nomenclature via VariantValidator",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input TSV (clinical or filtered)")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Output directory (default: same as input)")
    parser.add_argument("--vv-url", default=DEFAULT_VV_URL,
                        help=f"VariantValidator base URL (default: {DEFAULT_VV_URL})")
    parser.add_argument("--threads", type=int, default=1,
                        help="Parallel query threads (default: 1)")
    parser.add_argument("--timeout", type=int, default=120,
                        help="Per-query timeout in seconds (default: 120)")
    return parser.parse_args()


def check_vv_connection(base_url):
    """Verify VariantValidator is reachable."""
    try:
        test_url = (f"{base_url}/VariantValidator/variantvalidator/"
                    f"GRCh38/NM_000088.4:c.589G>T/all?content-type=application/json")
        resp = requests.get(test_url, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            if data.get("flag") in ("gene_variant", "warning"):
                log.info("VariantValidator reachable at %s", base_url)
                return True
    except (requests.ConnectionError, requests.Timeout):
        pass
    except Exception:
        pass

    log.error(
        "Cannot connect to VariantValidator at %s\n"
        "  Is the Docker container running?\n"
        "  Start it with:\n"
        "    cd ~/targeted-seq-pipeline/software/rest_variantValidator\n"
        "    docker compose up -d\n"
        "    docker exec -d rest_variantvalidator-rest-variantvalidator-1 "
        "bash -c 'cd /app/rest_VariantValidator && "
        "gunicorn -b 0.0.0.0:5000 --timeout 600 app --threads=5'",
        base_url,
    )
    return False


def build_query_hgvs(hgvsc, mane_select="", hgvsg=""):
    """Build HGVS query string for VariantValidator.

    Strategy:
    1. If MANE_SELECT available, combine RefSeq transcript with c. notation
    2. If no MANE_SELECT but HGVSg exists, use genomic notation
    3. Otherwise try original HGVSc as-is
    """
    hgvsc = str(hgvsc).strip()
    mane = str(mane_select).strip()
    hgvsg = str(hgvsg).strip()

    if hgvsc in ("-1", "", "nan"):
        return None

    match = re.search(r':(c\..+)$', hgvsc)
    if match:
        c_change = match.group(1)
        if mane not in ("-1", "", "nan"):
            return f"{mane}:{c_change}"

    if hgvsg not in ("-1", "", "nan"):
        return hgvsg

    return hgvsc


def query_variant(hgvsc, base_url, timeout):
    """Query VariantValidator for a single variant."""
    result = {
        "VV_HGVSc": "",
        "VV_HGVSp": "",
        "VV_HGVSg": "",
        "VV_Transcript": "",
        "VV_Valid": False,
        "VV_Warnings": "",
    }

    url = (
        f"{base_url}/VariantValidator/variantvalidator/"
        f"{GENOME_BUILD}/{requests.utils.quote(hgvsc, safe='')}/all"
        f"?content-type=application/json"
    )

    max_retries = 5
    data = None
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=timeout)
            if resp.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            resp.raise_for_status()
            data = resp.json()
            break
        except requests.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
                continue
            result["VV_Warnings"] = f"API_ERROR: {e}"
            return result
        except ValueError:
            result["VV_Warnings"] = "API_ERROR: invalid JSON response"
            return result
    else:
        result["VV_Warnings"] = "API_ERROR: max retries exceeded"
        return result

    if data is None:
        result["VV_Warnings"] = "API_ERROR: no response"
        return result

    flag = data.get("flag", "")
    variant_keys = [k for k in data.keys() if k not in ("flag", "metadata")]

    if not variant_keys:
        result["VV_Warnings"] = f"NO_RESULT: flag={flag}"
        return result

    if flag == "intergenic":
        result["VV_Warnings"] = "INTERGENIC"
        return result

    warnings = []

    # Find matching transcript result
    input_transcript = hgvsc.split(":")[0] if ":" in hgvsc else ""
    best_match = None

    for vkey in variant_keys:
        vdata = data[vkey]
        if not isinstance(vdata, dict):
            continue
        vv_hgvsc = vdata.get("hgvs_transcript_variant", "")
        vv_transcript = vv_hgvsc.split(":")[0] if ":" in vv_hgvsc else ""

        input_tx_base = input_transcript.split(".")[0]
        vv_tx_base = vv_transcript.split(".")[0]

        if input_tx_base and vv_tx_base and input_tx_base == vv_tx_base:
            best_match = vdata
            break

    if best_match is None:
        for vkey in variant_keys:
            vdata = data[vkey]
            if isinstance(vdata, dict) and vdata.get("hgvs_transcript_variant"):
                best_match = vdata
                break

    if best_match is None:
        for vkey in variant_keys:
            vdata = data[vkey]
            if isinstance(vdata, dict):
                vw = vdata.get("validation_warnings", [])
                if vw:
                    warnings.extend(vw)
        result["VV_Warnings"] = ("; ".join(warnings) if warnings
                                 else f"NO_MATCH: flag={flag}")
        return result

    vdata = best_match

    vv_hgvsc = vdata.get("hgvs_transcript_variant", "")
    result["VV_HGVSc"] = vv_hgvsc
    result["VV_Transcript"] = vv_hgvsc.split(":")[0] if ":" in vv_hgvsc else ""

    protein = vdata.get("hgvs_predicted_protein_consequence", {})
    if isinstance(protein, dict):
        result["VV_HGVSp"] = protein.get("tlr", "") or protein.get("slr", "")

    pal = vdata.get("primary_assembly_loci", {})
    grch38 = pal.get("grch38", {})
    result["VV_HGVSg"] = grch38.get("hgvs_genomic_description", "")

    vw = vdata.get("validation_warnings", [])
    if vw:
        warnings.extend(vw)

    # Check if VV corrected the HGVS
    if vv_hgvsc and vv_hgvsc != hgvsc:
        input_no_ver = ":".join(p.split(".")[0] if i == 0 else p
                                for i, p in enumerate(hgvsc.split(":")))
        vv_no_ver = ":".join(p.split(".")[0] if i == 0 else p
                             for i, p in enumerate(vv_hgvsc.split(":")))
        if input_no_ver != vv_no_ver:
            warnings.insert(0, f"CORRECTED: {hgvsc} -> {vv_hgvsc}")
        elif hgvsc.split(":")[0] != vv_hgvsc.split(":")[0]:
            warnings.insert(
                0,
                f"TRANSCRIPT_VERSION: {hgvsc.split(':')[0]} -> "
                f"{vv_hgvsc.split(':')[0]}"
            )

    result["VV_Valid"] = True
    result["VV_Warnings"] = "; ".join(warnings)
    return result


def validate_variants(df, base_url, threads, timeout):
    """Validate all variants with HGVSc via parallel queries."""
    mask = ~df["HGVSc"].isin(["-1", "", "nan"]) & df["HGVSc"].notna()
    query_indices = df.index[mask].tolist()
    skip_count = len(df) - len(query_indices)
    log.info("Querying %d variants (%d skipped -- no HGVSc)",
             len(query_indices), skip_count)

    # Build query HGVS using MANE_SELECT RefSeq transcripts where possible
    query_to_indices = {}
    idx_to_query = {}
    no_query_count = 0
    for idx in query_indices:
        hgvsc = str(df.at[idx, "HGVSc"])
        mane = str(df.at[idx, "MANE_SELECT"]) if "MANE_SELECT" in df.columns else ""
        hgvsg = str(df.at[idx, "HGVSg"]) if "HGVSg" in df.columns else ""
        query_hgvs = build_query_hgvs(hgvsc, mane, hgvsg)
        if query_hgvs:
            query_to_indices.setdefault(query_hgvs, []).append(idx)
            idx_to_query[idx] = query_hgvs
        else:
            no_query_count += 1

    unique_hgvsc = list(query_to_indices.keys())
    log.info("Unique HGVS queries: %d (%d could not be converted)",
             len(unique_hgvsc), no_query_count)

    for col in ["VV_HGVSc", "VV_HGVSp", "VV_HGVSg", "VV_Transcript", "VV_Warnings"]:
        df[col] = ""
    df["VV_Valid"] = ""

    results = {}
    completed = 0
    failed = 0
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_hgvsc = {
            executor.submit(query_variant, hgvsc, base_url, timeout): hgvsc
            for hgvsc in unique_hgvsc
        }

        for future in as_completed(future_to_hgvsc):
            hgvsc = future_to_hgvsc[future]
            completed += 1

            try:
                result = future.result()
            except Exception as e:
                result = {
                    "VV_HGVSc": "", "VV_HGVSp": "", "VV_HGVSg": "",
                    "VV_Transcript": "", "VV_Valid": False,
                    "VV_Warnings": f"EXCEPTION: {e}",
                }

            if not result["VV_Valid"]:
                failed += 1

            results[hgvsc] = result

            if completed % 50 == 0 or completed == len(unique_hgvsc):
                elapsed = time.time() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                log.info("  Progress: %d/%d (%.1f/sec, %d failed)",
                         completed, len(unique_hgvsc), rate, failed)

    for query_hgvs, indices in query_to_indices.items():
        result = results.get(query_hgvs, {})
        for idx in indices:
            for col in ["VV_HGVSc", "VV_HGVSp", "VV_HGVSg",
                        "VV_Transcript", "VV_Warnings"]:
                df.at[idx, col] = result.get(col, "")
            df.at[idx, "VV_Valid"] = result.get("VV_Valid", False)

    return df, len(query_indices), len(unique_hgvsc), failed


def main():
    args = parse_args()

    if not os.path.exists(args.input):
        log.error("Input file not found: %s", args.input)
        sys.exit(1)

    outdir = args.outdir or os.path.dirname(args.input)
    os.makedirs(outdir, exist_ok=True)

    if not check_vv_connection(args.vv_url):
        sys.exit(1)

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    log.info("Read %d variants from %s", len(df), args.input)

    df, total_queried, unique_queried, total_failed = validate_variants(
        df, args.vv_url, args.threads, args.timeout,
    )

    # Drop original HGVS columns (VV versions replace them)
    for col in ["HGVSc", "HGVSp", "HGVSg"]:
        if col in df.columns:
            df = df.drop(columns=[col])

    # Reorder columns
    desired_order = [
        "Sample", "Consensus", "Chr", "Start", "End", "Ref", "Alt",
        "Gene", "Consequence",
        "VV_HGVSc", "VV_HGVSp", "VV_HGVSg", "VV_Transcript",
        "VV_Valid", "VV_Warnings",
        "IMPACT", "REF_COUNT", "ALT_COUNT", "Total_Depth", "VAF_pct",
        "COSMIC_ID", "ClinVar", "SIFT", "PolyPhen",
        "gnomAD_exome_AF", "gnomAD_genome_AF", "AF_1KG", "Max_AF", "rsID",
        "MANE_SELECT", "Canonical", "Existing_variation",
        "Dedup_Note", "Filter",
    ]
    final_cols = [c for c in desired_order if c in df.columns]
    extra_cols = [c for c in df.columns if c not in desired_order]
    if extra_cols:
        log.info("Extra columns appended: %s", extra_cols)
    df = df[final_cols + extra_cols]

    # Write output
    basename = os.path.splitext(os.path.basename(args.input))[0]
    basename = basename.replace(".clinical", "").replace(".filtered", "")
    out_path = os.path.join(outdir, f"{basename}.validated.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    log.info("Wrote validated output: %s", out_path)

    # Summary
    total_valid = ((df["VV_Valid"] == True) | (df["VV_Valid"] == "True")).sum()
    total_warnings = df["VV_Warnings"].apply(
        lambda x: bool(str(x).strip())).sum()
    total_corrected = df["VV_Warnings"].str.contains(
        "CORRECTED", na=False).sum()

    log.info("=== VariantValidator Summary ===")
    log.info("  Total variants:       %d", len(df))
    log.info("  Queried (with HGVSc): %d (%d unique)", total_queried, unique_queried)
    log.info("  Validated:            %d", total_valid)
    log.info("  With warnings:        %d", total_warnings)
    log.info("  HGVS corrected:       %d", total_corrected)
    log.info("  Failed:               %d", total_failed)


if __name__ == "__main__":
    main()
