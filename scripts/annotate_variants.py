#!/usr/bin/env python3
"""
annotate_variants.py -- Annotate Anneal Rust caller VCFs with VEP + ANNOVAR.

Adapted from targeted-seq-pipeline/13_annotate.py for duplex consensus VCFs
produced by the Rust mpileup variant caller (FORMAT=GT:ALT:TOT:FRAC).
Also handles legacy INFO-only VCFs (DP=N;AF=X;AO=N) via fallback parsing.

Workflow:
  1. Parse Rust caller VCF to extract allelic depths and VAF.
  2. VEP: HGVS nomenclature, consequences, gene symbols, population
     frequencies, MANE Select transcripts.
  3. ANNOVAR: COSMIC, ClinVar, dbSNP, gnomAD genome AF.
  4. Merge all sources on chr:pos:ref:alt into a flat TSV.

Usage:
    python3 annotate_variants.py \\
        --vcf results/SAMPLE/variants/SAMPLE.dcs.vcf \\
        -s SAMPLE \\
        --consensus dcs \\
        -o results/SAMPLE/annotated

Requires:
  - VEP in conda env 'vep': conda run -n vep vep ...
  - ANNOVAR: table_annovar.pl
  - ANNOVAR hg38 databases: refGene, cosmic103, gnomad30_genome,
    clinvar_20220320, avsnp150
"""

import argparse
import csv
import logging
import os
import re
import subprocess
import sys
import time

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.dirname(SCRIPT_DIR)

# Default paths -- override via config or CLI
DEFAULT_REF = os.path.join(
    os.path.expanduser("~"), "targeted-seq-pipeline", "references", "hg38_broad",
    "resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta",
)
VEP_CACHE = os.path.join(
    os.path.expanduser("~"), "targeted-seq-pipeline", "references", "vep_cache",
)
ANNOVAR_DIR = os.path.join(
    os.path.expanduser("~"), "targeted-seq-pipeline", "software", "annovar",
)
ANNOVAR_DB = os.path.join(ANNOVAR_DIR, "humandb")

COLUMNS = [
    "Sample", "Consensus", "Chr", "Start", "End", "Ref", "Alt",
    "Gene", "Consequence", "HGVSc", "HGVSp", "IMPACT",
    "REF_COUNT", "ALT_COUNT", "Total_Depth", "VAF_pct",
    "COSMIC_ID", "ClinVar", "SIFT", "PolyPhen",
    "gnomAD_exome_AF", "gnomAD_genome_AF", "AF_1KG", "Max_AF", "rsID",
    "MANE_SELECT", "Canonical", "HGVSg", "Existing_variation",
]


def parse_args():
    ap = argparse.ArgumentParser(
        description="Annotate Rust caller VCF with VEP + ANNOVAR.",
    )
    ap.add_argument("--vcf", required=True, help="Rust caller VCF")
    ap.add_argument("-s", "--sample-name", required=True, help="Sample name")
    ap.add_argument("--consensus", default="dcs",
                    choices=["dcs", "sscs"],
                    help="Consensus type (default: dcs)")
    ap.add_argument("-r", "--reference", default=DEFAULT_REF,
                    help="Reference FASTA")
    ap.add_argument("-o", "--outdir", default="annotated",
                    help="Output directory")
    ap.add_argument("--vep-fork", type=int, default=4,
                    help="VEP parallel forks (default: 4)")
    ap.add_argument("--vep-cache", default=VEP_CACHE,
                    help="VEP cache directory")
    ap.add_argument("--annovar-dir", default=ANNOVAR_DIR,
                    help="ANNOVAR installation directory")
    ap.add_argument("--skip-vep", action="store_true",
                    help="Skip VEP annotation")
    ap.add_argument("--skip-annovar", action="store_true",
                    help="Skip ANNOVAR annotation")
    return ap.parse_args()


def run(cmd, desc=None, shell=False):
    if desc:
        log.info("%s", desc)
    cmd_str = cmd if shell else " ".join(cmd)
    log.info("  cmd: %s", cmd_str)
    result = subprocess.run(cmd, capture_output=True, text=True, shell=shell)
    if result.returncode != 0:
        log.error("  FAILED (exit %d)", result.returncode)
        for line in (result.stderr or "").strip().splitlines()[-10:]:
            log.error("    %s", line.strip())
    return result.returncode


def run_vep(vcf_in, vcf_out, reference, fork, vep_cache):
    cmd = [
        "conda", "run", "-n", "vep", "vep",
        "--input_file", vcf_in,
        "--output_file", vcf_out,
        "--vcf",
        "--offline",
        "--cache",
        "--dir_cache", vep_cache,
        "--assembly", "GRCh38",
        "--fasta", reference,
        "--fork", str(fork),
        "--force_overwrite",
        "--pick",
        "--everything",
        "--hgvs",
        "--hgvsg",
        "--symbol",
        "--canonical",
        "--mane_select",
    ]
    return run(cmd, desc=f"Running VEP on {os.path.basename(vcf_in)}")


def run_annovar(vcf_in, out_prefix, annovar_dir):
    table_annovar = os.path.join(annovar_dir, "table_annovar.pl")
    annovar_db = os.path.join(annovar_dir, "humandb")

    protocols = []
    operations = []
    db_checks = [
        ("refGene", "g"),
        ("cosmic103", "f"),
        ("gnomad30_genome", "f"),
        ("clinvar_20220320", "f"),
        ("avsnp150", "f"),
    ]
    for db, op in db_checks:
        db_file = os.path.join(annovar_db, f"hg38_{db}.txt")
        if os.path.isfile(db_file):
            protocols.append(db)
            operations.append(op)
        else:
            log.warning("ANNOVAR database not found, skipping: %s", db)

    if not protocols:
        log.error("No ANNOVAR databases available")
        return 1

    cmd = [
        "perl", table_annovar,
        vcf_in,
        annovar_db,
        "-buildver", "hg38",
        "-out", out_prefix,
        "-remove",
        "-protocol", ",".join(protocols),
        "-operation", ",".join(operations),
        "-nastring", ".",
        "-vcfinput",
    ]
    return run(cmd, desc=f"Running ANNOVAR on {os.path.basename(vcf_in)}")


def _get_info_value(info_str, key):
    for field in info_str.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return ""


def parse_rust_vcf(vcf_path):
    """Parse Rust mpileup caller VCF.

    FORMAT fields: GT:ALT:TOT:FRAC
      - ALT = number of alt-supporting reads
      - TOT = total depth (molecular tags)
      - FRAC = VAF as fraction (0.0-1.0)

    Returns dict keyed by chr:pos:ref:alt.
    """
    variants = {}

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            filt = cols[6] if len(cols) > 6 else "."

            alt_count = -1
            tot_depth = -1
            ref_count = -1
            vaf_pct = -1

            # Parse FORMAT fields
            if len(cols) >= 10:
                fmt_keys = cols[8].split(":")
                fmt_vals = cols[9].split(":")
                fmt = dict(zip(fmt_keys, fmt_vals))

                # Rust caller: ALT=alt count, TOT=total depth, FRAC=vaf
                try:
                    alt_count = int(fmt.get("ALT", -1))
                except (ValueError, TypeError):
                    pass
                try:
                    tot_depth = int(fmt.get("TOT", -1))
                except (ValueError, TypeError):
                    pass
                try:
                    vaf_pct = round(float(fmt.get("FRAC", -1)) * 100, 4)
                except (ValueError, TypeError):
                    pass

                # Also handle standard AD field (ref,alt) if present
                ad = fmt.get("AD", "")
                if ad and alt_count < 0:
                    parts = ad.split(",")
                    if len(parts) >= 2:
                        try:
                            ref_count = int(parts[0])
                            alt_count = int(parts[1])
                        except ValueError:
                            pass

                # DP field
                dp = fmt.get("DP", "")
                if dp and tot_depth < 0:
                    try:
                        tot_depth = int(dp)
                    except ValueError:
                        pass

            # Fallback: parse INFO field (for INFO-only VCFs, e.g.
            # older call_variants output with DP=N;AF=X;AO=N)
            if alt_count < 0 or tot_depth < 0:
                info_str = cols[7] if len(cols) > 7 else ""
                if info_str and info_str != ".":
                    info_fields = {}
                    for entry in info_str.split(";"):
                        if "=" in entry:
                            k, v = entry.split("=", 1)
                            info_fields[k] = v

                    if tot_depth < 0 and "DP" in info_fields:
                        try:
                            tot_depth = int(info_fields["DP"])
                        except (ValueError, TypeError):
                            pass
                    if alt_count < 0 and "AO" in info_fields:
                        try:
                            alt_count = int(info_fields["AO"])
                        except (ValueError, TypeError):
                            pass
                    if vaf_pct < 0 and "AF" in info_fields:
                        try:
                            vaf_pct = round(float(info_fields["AF"]) * 100, 4)
                        except (ValueError, TypeError):
                            pass

            # Compute ref_count if we have total and alt
            if ref_count < 0 and tot_depth > 0 and alt_count >= 0:
                ref_count = tot_depth - alt_count

            # Compute VAF if not already set
            if vaf_pct < 0 and alt_count >= 0 and tot_depth > 0:
                vaf_pct = round(alt_count / tot_depth * 100, 4)

            key = f"{chrom}:{pos}:{ref}:{alt}"
            variants[key] = {
                "filter": filt,
                "ref_count": ref_count,
                "alt_count": alt_count,
                "total_depth": tot_depth,
                "vaf_pct": vaf_pct,
            }

    log.info("Parsed %d variants from Rust caller VCF", len(variants))
    return variants


def parse_vep_csq(vep_vcf):
    """Parse VEP VCF output and extract CSQ fields."""
    variants = {}
    csq_fields = []

    with open(vep_vcf) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r'Format: ([^"]+)', line)
                if match:
                    csq_fields = match.group(1).strip().split("|")
                continue
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            info = cols[7]

            csq_data = {}
            for field in info.split(";"):
                if field.startswith("CSQ="):
                    csq_str = field[4:]
                    first_csq = csq_str.split(",")[0]
                    values = first_csq.split("|")
                    for i, val in enumerate(values):
                        if i < len(csq_fields):
                            csq_data[csq_fields[i]] = val
                    break

            key = f"{chrom}:{pos}:{ref}:{alt}"
            variants[key] = csq_data

    log.info("Parsed %d variants from VEP VCF (%d CSQ fields)",
             len(variants), len(csq_fields))
    return variants, csq_fields


def parse_annovar_txt(annovar_txt):
    """Parse ANNOVAR multianno.txt into a dict keyed by chr:pos:ref:alt."""
    variants = {}
    if not os.path.isfile(annovar_txt):
        log.warning("ANNOVAR output not found: %s", annovar_txt)
        return variants

    with open(annovar_txt) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row.get("Chr", "")
            pos = row.get("Start", "")
            ref = row.get("Ref", "")
            alt = row.get("Alt", "")
            key = f"{chrom}:{pos}:{ref}:{alt}"
            variants[key] = row

    log.info("Parsed %d variants from ANNOVAR txt", len(variants))
    return variants


def _clean(val):
    if val is None or val == "" or val == ".":
        return "-1"
    return str(val)


def merge_annotations(vcf_fields, vep_variants, annovar_variants,
                      output_tsv, sample, consensus):
    all_keys = set(vcf_fields.keys()) | set(vep_variants.keys()) | set(annovar_variants.keys())
    log.info("Merging: %d VCF, %d VEP, %d ANNOVAR, %d total unique",
             len(vcf_fields), len(vep_variants), len(annovar_variants), len(all_keys))

    rows = []
    for key in sorted(all_keys):
        parts = key.split(":", 3)
        if len(parts) != 4:
            continue
        chrom, pos, ref, alt = parts

        vcf = vcf_fields.get(key, {})
        vep = vep_variants.get(key, {})
        ann = annovar_variants.get(key, {})

        end = int(pos) + max(len(ref), 1) - 1

        row = {
            "Sample": sample,
            "Consensus": consensus.upper(),
            "Chr": chrom,
            "Start": pos,
            "End": str(end),
            "Ref": ref,
            "Alt": alt,
            "Gene": _clean(vep.get("SYMBOL", ann.get("Gene.refGene", ""))),
            "Consequence": _clean(vep.get("Consequence",
                                          ann.get("ExonicFunc.refGene", ""))),
            "HGVSc": _clean(vep.get("HGVSc", "")),
            "HGVSp": _clean(vep.get("HGVSp", "")),
            "IMPACT": _clean(vep.get("IMPACT", "")),
            "REF_COUNT": _clean(vcf.get("ref_count", -1)),
            "ALT_COUNT": _clean(vcf.get("alt_count", -1)),
            "Total_Depth": _clean(vcf.get("total_depth", -1)),
            "VAF_pct": _clean(vcf.get("vaf_pct", -1)),
            "COSMIC_ID": _clean(ann.get("cosmic103", "")),
            "ClinVar": _clean(ann.get("CLNSIG", ann.get("clinvar_20220320", ""))),
            "SIFT": _clean(vep.get("SIFT", "")),
            "PolyPhen": _clean(vep.get("PolyPhen", "")),
            "gnomAD_exome_AF": _clean(vep.get("gnomADe_AF", "")),
            "gnomAD_genome_AF": _clean(vep.get("gnomADg_AF",
                                                 ann.get("gnomad30_genome", ""))),
            "AF_1KG": _clean(vep.get("AF", "")),
            "Max_AF": _clean(vep.get("MAX_AF", "")),
            "rsID": _clean(ann.get("avsnp150", vep.get("Existing_variation", ""))),
            "MANE_SELECT": _clean(vep.get("MANE_SELECT", "")),
            "Canonical": _clean(vep.get("CANONICAL", "")),
            "HGVSg": _clean(vep.get("HGVSg", "")),
            "Existing_variation": _clean(vep.get("Existing_variation", "")),
        }
        rows.append(row)

    with open(output_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    log.info("Wrote %d annotated variants to %s", len(rows), output_tsv)
    return len(rows)


def main():
    t0 = time.time()
    args = parse_args()
    sample = args.sample_name
    consensus = args.consensus

    log.info("=== Variant Annotation (Anneal) ===")
    log.info("Sample:    %s", sample)
    log.info("Consensus: %s", consensus)
    log.info("VCF:       %s", args.vcf)
    log.info("Reference: %s", args.reference)
    log.info("Output:    %s", args.outdir)

    if not os.path.isfile(args.vcf):
        log.error("VCF not found: %s", args.vcf)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    basename = os.path.splitext(os.path.basename(args.vcf))[0]

    # Step 1: Parse Rust caller VCF
    vcf_fields = parse_rust_vcf(args.vcf)

    # Step 2: VEP
    vep_variants = {}
    if not args.skip_vep:
        if not os.path.isdir(args.vep_cache):
            log.warning("VEP cache not found: %s -- skipping VEP", args.vep_cache)
        else:
            vep_vcf = os.path.join(args.outdir, f"{basename}.vep.vcf")
            rc = run_vep(args.vcf, vep_vcf, args.reference,
                         args.vep_fork, args.vep_cache)
            if rc == 0:
                vep_variants, _ = parse_vep_csq(vep_vcf)
            else:
                log.error("VEP failed -- continuing without VEP annotations")
    else:
        log.info("Skipping VEP (--skip-vep)")

    # Step 3: ANNOVAR
    annovar_variants = {}
    if not args.skip_annovar:
        table_annovar = os.path.join(args.annovar_dir, "table_annovar.pl")
        if not os.path.isfile(table_annovar):
            log.warning("ANNOVAR not found: %s -- skipping", table_annovar)
        else:
            annovar_prefix = os.path.join(args.outdir, basename)
            rc = run_annovar(args.vcf, annovar_prefix, args.annovar_dir)
            annovar_txt = annovar_prefix + ".hg38_multianno.txt"
            if rc == 0:
                annovar_variants = parse_annovar_txt(annovar_txt)
            else:
                log.warning("ANNOVAR failed -- continuing without ANNOVAR")
    else:
        log.info("Skipping ANNOVAR (--skip-annovar)")

    # Step 4: Merge
    output_tsv = os.path.join(args.outdir, f"{basename}.annotated.tsv")
    n = merge_annotations(vcf_fields, vep_variants, annovar_variants,
                          output_tsv, sample, consensus)

    elapsed = time.time() - t0
    log.info("")
    log.info("Output: %s (%d variants)", output_tsv, n)
    log.info("Time:   %.0fs", elapsed)


if __name__ == "__main__":
    main()
