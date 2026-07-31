#!/usr/bin/env python3
"""
Annotate the indel table produced by scan_indels.py.

The indel table is kept separate from the substitution annotation because it
carries provenance the Pisces path does not: supporting read count, forward and
reverse counts, strand fraction, blocklist recurrence across the BNC panel, and
artifact mask status. Merging into the substitution table would either drop
those or force a schema change on both.

Indels are already VCF-anchored in the table (POS is the base before the event,
REF and ALT both carry it), so they are written to a minimal VCF and passed
through the same VEP and ANNOVAR invocations the substitution path uses. The
annotation is then joined back onto the original rows by chrom:pos:ref:alt.

Output keeps every input column and appends Gene, Consequence, HGVSc, HGVSp,
IMPACT, COSMIC_ID, ClinVar, gnomAD_AF and rsID.

    annotate_indels.py --indels S.dcs.indels.tsv --out S.dcs.indels.annotated.tsv \\
        --ref genome.fa --vep-cache /path/vep_cache \\
        --annovar-dir /path/annovar --annovar-db /path/humandb
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

# Conda env bin directories whose perl would shadow VEP's own. VEP loads its
# modules relative to whichever perl is found first on PATH.
VEP_PATH_STRIP = ("envs/targeted-seq/bin", "envs/anneal/bin")
VEP_PERL_VARS = ("PERL5LIB", "PERL_LOCAL_LIB_ROOT", "PERL_MM_OPT", "PERL_MB_OPT")

ANNOT_COLS = ["Gene", "Consequence", "HGVSc", "HGVSp", "IMPACT",
              "COSMIC_ID", "ClinVar", "gnomAD_AF", "rsID"]


def vep_environment():
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join(
        p for p in env.get("PATH", "").split(os.pathsep)
        if not any(s in p for s in VEP_PATH_STRIP))
    for v in VEP_PERL_VARS:
        env.pop(v, None)
    return env


def run(cmd, desc=None, env=None):
    if desc:
        log.info("%s", desc)
    log.info("  cmd: %s", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        log.error("  FAILED (exit %d)", r.returncode)
        for line in (r.stderr or "").strip().splitlines()[-10:]:
            log.error("    %s", line.strip())
    return r.returncode


def read_indels(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    return header, rows


def chrom_sort_key(c):
    c = c.replace("chr", "")
    return (0, int(c)) if c.isdigit() else (1, c)


def write_vcf(rows, path, ref_fai):
    contigs = []
    if os.path.exists(ref_fai):
        for line in open(ref_fai):
            f = line.split("\t")
            if len(f) >= 2:
                contigs.append((f[0], f[1]))
    seen = set()
    uniq = []
    for r in rows:
        k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
        if k not in seen:
            seen.add(k)
            uniq.append(k)
    uniq.sort(key=lambda k: (chrom_sort_key(k[0]), k[1]))

    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for name, length in contigs:
            out.write(f"##contig=<ID={name},length={length}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in uniq:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return len(uniq)


def run_vep(vcf_in, vcf_out, reference, vep_cache, fork):
    cmd = [
        "conda", "run", "-p",
        os.environ.get("VEP_PREFIX", os.path.expanduser("~/miniconda3/envs/vep")),
        "vep",
        "--input_file", vcf_in, "--output_file", vcf_out,
        "--vcf", "--offline", "--cache", "--dir_cache", vep_cache,
        "--assembly", "GRCh38", "--fasta", reference,
        "--fork", str(fork), "--force_overwrite", "--flag_pick",
        "--symbol", "--canonical", "--mane_select", "--hgvs",
    ]
    return run(cmd, desc=f"VEP on {os.path.basename(vcf_in)}", env=vep_environment())


def parse_vep(path):
    """{(chrom,pos,ref,alt): {field: value}} from the picked CSQ entry."""
    out = {}
    if not os.path.exists(path):
        return out
    fields = []
    for line in open(path):
        if line.startswith("##"):
            if "ID=CSQ" in line and "Format:" in line:
                fields = line.split("Format:")[1].strip().rstrip('">').split("|")
            continue
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 8 or "CSQ=" not in c[7]:
            continue
        key = (c[0], int(c[1]), c[3], c[4])
        csq = c[7].split("CSQ=")[1].split(";")[0]
        best = None
        for entry in csq.split(","):
            d = dict(zip(fields, entry.split("|")))
            if d.get("PICK") == "1":
                best = d
                break
            if best is None:
                best = d
        if best:
            out[key] = best
    return out


def run_annovar(vcf_in, out_prefix, annovar_dir, annovar_db):
    table = os.path.join(annovar_dir, "table_annovar.pl")
    protocols, operations = [], []
    for db, op in (("refGene", "g"), ("cosmic103", "f"), ("gnomad211_exome", "f"),
                   ("clinvar_20250721", "f"), ("avsnp151", "f")):
        if os.path.isfile(os.path.join(annovar_db, f"hg38_{db}.txt")):
            protocols.append(db)
            operations.append(op)
        else:
            log.warning("ANNOVAR db missing, skipping: %s", db)
    cmd = ["perl", table, vcf_in, annovar_db, "-buildver", "hg38",
           "-out", out_prefix, "-remove",
           "-protocol", ",".join(protocols), "-operation", ",".join(operations),
           "-nastring", ".", "-vcfinput"]
    return run(cmd, desc=f"ANNOVAR on {os.path.basename(vcf_in)}")


def parse_annovar(path):
    """{(chrom,pos,ref,alt): row} keyed on the original VCF locus, which ANNOVAR
    preserves in the Otherinfo block. Its own Ref/Alt columns are rewritten into
    '-' form for indels and do not match the input."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            # find the original VCF columns in Otherinfo: CHROM POS ID REF ALT
            ovcf = None
            for i in range(len(c) - 5, 0, -1):
                if (c[i].startswith("chr") and i + 4 < len(c)
                        and c[i + 1].isdigit()):
                    ovcf = i
                    break
            if ovcf is None:
                continue
            key = (c[ovcf], int(c[ovcf + 1]), c[ovcf + 3], c[ovcf + 4])
            out[key] = {h: (c[i] if i < len(c) else ".") for h, i in idx.items()}
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indels", required=True, help="scan_indels.py output")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref", required=True, help="reference FASTA (unmasked)")
    ap.add_argument("--vep-cache", required=True)
    ap.add_argument("--annovar-dir", required=True)
    ap.add_argument("--annovar-db", required=True)
    ap.add_argument("--fork", type=int, default=4)
    ap.add_argument("--keep-temp", action="store_true")
    args = ap.parse_args()

    header, rows = read_indels(args.indels)
    if not rows:
        log.warning("no indels in %s, writing header only", args.indels)
        with open(args.out, "w") as out:
            out.write("\t".join(header + ANNOT_COLS) + "\n")
        return

    work = tempfile.mkdtemp(prefix="annot_indels_")
    try:
        vcf = os.path.join(work, "indels.vcf")
        n = write_vcf(rows, vcf, args.ref + ".fai")
        log.info("%d unique indel loci from %d rows", n, len(rows))

        vep_vcf = os.path.join(work, "indels.vep.vcf")
        run_vep(vcf, vep_vcf, args.ref, args.vep_cache, args.fork)
        vep = parse_vep(vep_vcf)
        log.info("VEP annotated %d loci", len(vep))

        av_prefix = os.path.join(work, "indels")
        run_annovar(vcf, av_prefix, args.annovar_dir, args.annovar_db)
        av = parse_annovar(av_prefix + ".hg38_multianno.txt")
        log.info("ANNOVAR annotated %d loci", len(av))

        with open(args.out, "w") as out:
            out.write("\t".join(header + ANNOT_COLS) + "\n")
            for r in rows:
                key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                v = vep.get(key, {})
                a = av.get(key, {})
                extra = [
                    v.get("SYMBOL") or a.get("Gene.refGene", ".") or ".",
                    v.get("Consequence", ".") or ".",
                    v.get("HGVSc", ".") or ".",
                    v.get("HGVSp", ".") or ".",
                    v.get("IMPACT", ".") or ".",
                    a.get("cosmic103", ".") or ".",
                    a.get("clinvar_20250721", ".") or ".",
                    a.get("AF", ".") or ".",
                    a.get("avsnp151", ".") or ".",
                ]
                out.write("\t".join([r.get(h, ".") for h in header] + extra) + "\n")
    finally:
        if not args.keep_temp:
            import shutil
            shutil.rmtree(work, ignore_errors=True)

    log.info("wrote %s (%d rows)", args.out, len(rows))


if __name__ == "__main__":
    main()
