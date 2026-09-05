#!/usr/bin/env python3
"""
patch_validate_hgvs.py

VariantValidator moves from a local Docker container to the public REST API.
Compute nodes on this cluster have no route out, so the step is off in
stage 3 by default and runs post hoc on the login node over each batch's
clinical tables, with a cache so every unique HGVS query is made once.

Edits (exact anchors, each must match once; backups <file>.bak.<ts>):

  scripts/validate_hgvs.py
    - DEFAULT_VV_URL -> https://rest.variantvalidator.org
    - throttle: at least --min-interval seconds between requests (default 1.0)
    - --cache <json>: results keyed on the query HGVS, loaded before querying,
      saved after; transient API errors are not cached
    - User-Agent header on every request
  pipeline/stage3_annotate.sh
    - SKIP_VV defaults from VV_IN_STAGE3 (config, default false)
    - passes --vv-url and --cache from config
  pipeline/config.sh
    - VV_URL, VV_CACHE, VV_IN_STAGE3 appended
  pipeline/validate_hgvs_batch.sh (new)
    - post-hoc runner: validate_hgvs_batch.sh <outdir> [sample ...]

Usage:
    python patch_validate_hgvs.py --root ~/pipelines/anneal [--dry-run]
"""

import argparse
import os
import shutil
import sys
import time

# --------------------------------------------------------------------------
# scripts/validate_hgvs.py
# --------------------------------------------------------------------------

VH_URL_OLD = 'DEFAULT_VV_URL = "http://localhost:5001"\n'
VH_URL_NEW = '''DEFAULT_VV_URL = "https://rest.variantvalidator.org"
USER_AGENT = "anneal/0.3.0 (patkarlab, MRD pipeline)"

# Politeness for the shared public endpoint: at least MIN_INTERVAL seconds
# between requests across all threads. Set from --min-interval in main().
MIN_INTERVAL = 1.0
CACHE_PATH = None
_THROTTLE_LOCK = threading.Lock()
_LAST_REQUEST = [0.0]


def _throttle():
    with _THROTTLE_LOCK:
        wait = MIN_INTERVAL - (time.monotonic() - _LAST_REQUEST[0])
        if wait > 0:
            time.sleep(wait)
        _LAST_REQUEST[0] = time.monotonic()


def load_cache(path):
    """Cache of query HGVS -> result dict. Missing or unreadable -> empty."""
    if not path or not os.path.isfile(path):
        return {}
    try:
        with open(path) as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError) as e:
        log.warning("cache unreadable, starting empty: %s", e)
        return {}


def save_cache(path, cache):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(cache, fh, indent=0, sort_keys=True)
    os.replace(tmp, path)


def is_transient(result):
    w = str(result.get("VV_Warnings", ""))
    return w.startswith("API_ERROR") or w.startswith("EXCEPTION")
'''

VH_IMPORT_OLD = "import time\nfrom concurrent.futures import ThreadPoolExecutor, as_completed\n"
VH_IMPORT_NEW = ("import json\nimport threading\nimport time\n"
                 "from concurrent.futures import ThreadPoolExecutor, as_completed\n")

VH_GET_OLD = "            resp = requests.get(url, timeout=timeout)\n            if resp.status_code == 429:\n"
VH_GET_NEW = ('            _throttle()\n'
              '            resp = requests.get(url, timeout=timeout,\n'
              '                                headers={"User-Agent": USER_AGENT,\n'
              '                                         "Accept": "application/json"})\n'
              '            if resp.status_code == 429:\n')

VH_ARGS_OLD = '''    parser.add_argument("--timeout", type=int, default=120,
                        help="Per-query timeout in seconds (default: 120)")
'''
VH_ARGS_NEW = '''    parser.add_argument("--timeout", type=int, default=120,
                        help="Per-query timeout in seconds (default: 120)")
    parser.add_argument("--cache", default=None,
                        help="JSON cache of results keyed on query HGVS; "
                             "read before querying, updated after")
    parser.add_argument("--min-interval", type=float, default=1.0,
                        help="Minimum seconds between requests to the public "
                             "endpoint (default: 1.0)")
'''

VH_MAIN_OLD = "    if not check_vv_connection(args.vv_url):\n"
VH_MAIN_NEW = ('    global CACHE_PATH, MIN_INTERVAL\n'
               '    CACHE_PATH = args.cache\n'
               '    MIN_INTERVAL = args.min_interval\n'
               '    if not check_vv_connection(args.vv_url):\n')

VH_LOG_OLD = '''    unique_hgvsc = list(query_to_indices.keys())
    log.info("Unique HGVS queries: %d (%d could not be converted)",
             len(unique_hgvsc), no_query_count)
'''
VH_LOG_NEW = '''    unique_hgvsc = list(query_to_indices.keys())
    log.info("Unique HGVS queries: %d (%d could not be converted)",
             len(unique_hgvsc), no_query_count)

    cache = load_cache(CACHE_PATH)
    cached = {h: cache[h] for h in unique_hgvsc if h in cache}
    to_query = [h for h in unique_hgvsc if h not in cache]
    log.info("Cache %s: %d hits, %d to query",
             CACHE_PATH or "(none)", len(cached), len(to_query))
'''

VH_RESULTS_OLD = "    results = {}\n    completed = 0\n    failed = 0\n"
VH_RESULTS_NEW = "    results = dict(cached)\n    completed = 0\n    failed = 0\n"

VH_SUBMIT_OLD = '''            executor.submit(query_variant, hgvsc, base_url, timeout): hgvsc
            for hgvsc in unique_hgvsc
'''
VH_SUBMIT_NEW = '''            executor.submit(query_variant, hgvsc, base_url, timeout): hgvsc
            for hgvsc in to_query
'''

VH_VALID_OLD = '            df.at[idx, "VV_Valid"] = result.get("VV_Valid", False)\n'
VH_VALID_NEW = '            df.at[idx, "VV_Valid"] = str(result.get("VV_Valid", False))\n'

VH_PROG_OLD = '''            if completed % 50 == 0 or completed == len(unique_hgvsc):
                elapsed = time.time() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                log.info("  Progress: %d/%d (%.1f/sec, %d failed)",
                         completed, len(unique_hgvsc), rate, failed)

'''
VH_PROG_NEW = '''            if completed % 50 == 0 or completed == len(to_query):
                elapsed = time.time() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                log.info("  Progress: %d/%d (%.1f/sec, %d failed)",
                         completed, len(to_query), rate, failed)

    if CACHE_PATH:
        fresh = {h: r for h, r in results.items()
                 if h not in cached and not is_transient(r)}
        if fresh:
            cache.update(fresh)
            save_cache(CACHE_PATH, cache)
            log.info("Cache updated: +%d, %d total", len(fresh), len(cache))

'''

# --------------------------------------------------------------------------
# pipeline/stage3_annotate.sh, pipeline/config.sh, batch runner
# --------------------------------------------------------------------------

S3_DEFAULT_OLD = "SKIP_VV=false\n"
S3_DEFAULT_NEW = ('# HGVS validation needs the public VariantValidator API; compute nodes have\n'
                  '# no route out, so it is off unless VV_IN_STAGE3=true in config.sh.\n'
                  '# Run post hoc with pipeline/validate_hgvs_batch.sh instead.\n'
                  'SKIP_VV=true\n'
                  '[ "${VV_IN_STAGE3:-false}" = true ] && SKIP_VV=false\n')

S3_CALL_OLD = '''            python3 "${SCRIPTS_DIR}/validate_hgvs.py" \\
                -i "${CLINICAL_TSV}" \\
                -o "${ANNOTATED_DIR}" 2>&1 || \\
'''
S3_CALL_NEW = '''            python3 "${SCRIPTS_DIR}/validate_hgvs.py" \\
                -i "${CLINICAL_TSV}" \\
                -o "${ANNOTATED_DIR}" \\
                --vv-url "${VV_URL:-https://rest.variantvalidator.org}" \\
                ${VV_CACHE:+--cache "${VV_CACHE}"} 2>&1 || \\
'''

CFG_APPEND = '''
# ---- HGVS validation (VariantValidator public REST API) ----
# Compute nodes have no route out, so this does not run in stage 3 unless
# VV_IN_STAGE3=true. Run it on the login node after a batch:
#   bash pipeline/validate_hgvs_batch.sh <outdir> [sample ...]
# The cache holds patient HGVS strings: gitignored, never committed.
VV_URL="${VV_URL:-https://rest.variantvalidator.org}"
VV_CACHE="${VV_CACHE:-${ANNEAL_ROOT}/vv_cache.json}"
VV_IN_STAGE3="${VV_IN_STAGE3:-false}"
'''

BATCH_RUNNER = '''#!/bin/bash
# =============================================================================
# validate_hgvs_batch.sh -- HGVS validation through the public VariantValidator
#                           API, post hoc, on the login node.
#
# Usage:
#   bash validate_hgvs_batch.sh <outdir> [sample ...]
#
# Runs scripts/validate_hgvs.py on every <sample>.<track>.clinical.tsv under
# <outdir>/<sample>/annotated/ (all samples in <outdir> if none are named),
# sharing one cache so each unique HGVS query hits the API once across the
# batch. About one request per second; a cohort's unique clinical variants
# take minutes, not hours.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
if declare -F activate_conda >/dev/null 2>&1; then
    activate_conda anneal
fi
set -euo pipefail

OUTPUT_DIR="${1:?usage: validate_hgvs_batch.sh <outdir> [sample ...]}"
shift
if [ "$#" -gt 0 ]; then
    SAMPLES=("$@")
else
    SAMPLES=()
    for d in "${OUTPUT_DIR}"/*/annotated; do
        [ -d "${d}" ] && SAMPLES+=("$(basename "$(dirname "${d}")")")
    done
fi

VV_URL="${VV_URL:-https://rest.variantvalidator.org}"
VV_CACHE="${VV_CACHE:-${ANNEAL_ROOT}/vv_cache.json}"
echo "VariantValidator: ${VV_URL}   cache: ${VV_CACHE}   samples: ${#SAMPLES[@]}"

for SAMPLE in "${SAMPLES[@]}"; do
    for track in dcs sscs; do
        tsv="${OUTPUT_DIR}/${SAMPLE}/annotated/${SAMPLE}.${track}.clinical.tsv"
        [ -f "${tsv}" ] || continue
        n=$(($(wc -l < "${tsv}") - 1))
        [ "${n}" -gt 0 ] || { echo "--- ${SAMPLE} ${track}: no clinical variants"; continue; }
        echo "--- ${SAMPLE} ${track}: ${n} variants"
        python3 "${ANNEAL_ROOT}/scripts/validate_hgvs.py" \\
            -i "${tsv}" -o "${OUTPUT_DIR}/${SAMPLE}/annotated" \\
            --vv-url "${VV_URL}" --cache "${VV_CACHE}"
    done
done
echo "done: $(date)"
'''


# --------------------------------------------------------------------------
# machinery
# --------------------------------------------------------------------------

def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        sys.exit("ERROR: %s: anchor matched %d times, expected 1" % (label, n))
    return text.replace(old, new)


def patch_validate(text):
    if "def load_cache(" in text:
        print("  already patched, skipping")
        return text
    text = replace_once(text, VH_IMPORT_OLD, VH_IMPORT_NEW, "imports")
    text = replace_once(text, VH_URL_OLD, VH_URL_NEW, "default URL")
    text = replace_once(text, VH_GET_OLD, VH_GET_NEW, "requests.get")
    text = replace_once(text, VH_ARGS_OLD, VH_ARGS_NEW, "argparse")
    text = replace_once(text, VH_MAIN_OLD, VH_MAIN_NEW, "main globals")
    text = replace_once(text, VH_LOG_OLD, VH_LOG_NEW, "cache load")
    text = replace_once(text, VH_RESULTS_OLD, VH_RESULTS_NEW, "results init")
    text = replace_once(text, VH_SUBMIT_OLD, VH_SUBMIT_NEW, "executor submit")
    text = replace_once(text, VH_PROG_OLD, VH_PROG_NEW, "progress and cache save")
    text = replace_once(text, VH_VALID_OLD, VH_VALID_NEW, "VV_Valid as text")
    text = text.replace("Requires:\n  - Local VariantValidator Docker container running on localhost:5001",
                        "Requires:\n  - Network access to the VariantValidator REST API (login node)")
    return text


def patch_stage3(text):
    if "VV_IN_STAGE3" in text:
        print("  already patched, skipping")
        return text
    text = replace_once(text, S3_DEFAULT_OLD, S3_DEFAULT_NEW, "SKIP_VV default")
    text = replace_once(text, S3_CALL_OLD, S3_CALL_NEW, "validate call")
    return text


def patch_config(text):
    if "VV_IN_STAGE3" in text:
        print("  already patched, skipping")
        return text
    return text.rstrip("\n") + "\n" + CFG_APPEND


def apply(path, fn, dry_run):
    print("== %s" % path)
    if not os.path.isfile(path):
        sys.exit("ERROR: not found: %s" % path)
    original = open(path).read()
    patched = fn(original)
    if patched == original:
        print("  no change")
        return
    if path.endswith(".py"):
        compile(patched, path, "exec")
    if dry_run:
        print("  dry run: would write %d -> %d bytes" % (len(original), len(patched)))
        return
    backup = "%s.bak.%s" % (path, time.strftime("%Y%m%d-%H%M%S"))
    shutil.copy2(path, backup)
    with open(path, "w") as fh:
        fh.write(patched)
    print("  written; backup at %s" % backup)


def write_runner(path, dry_run):
    print("== %s" % path)
    if os.path.exists(path):
        print("  exists, left untouched")
        return
    if dry_run:
        print("  dry run: would create")
        return
    with open(path, "w") as fh:
        fh.write(BATCH_RUNNER)
    os.chmod(path, 0o755)
    print("  created")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    r = args.root
    apply(os.path.join(r, "scripts", "validate_hgvs.py"), patch_validate, args.dry_run)
    apply(os.path.join(r, "pipeline", "stage3_annotate.sh"), patch_stage3, args.dry_run)
    apply(os.path.join(r, "pipeline", "config.sh"), patch_config, args.dry_run)
    write_runner(os.path.join(r, "pipeline", "validate_hgvs_batch.sh"), args.dry_run)


if __name__ == "__main__":
    main()
