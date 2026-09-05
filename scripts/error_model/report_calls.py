#!/usr/bin/env python3
"""
report_calls.py

Turn a Stage 5 calls table into a readable report: the label column is split
into gene, consequence, protein change, protein id and source; rows are
restricted to protein-altering consequences unless --all is given; the
annotated tables, when supplied, contribute HGVSc, COSMIC, ClinVar, gnomAD and
rsID by chrom/pos/ref/alt; DETECTED rows come first, by VAF.

The calls table itself is untouched; this is a view for reading, and the
calls table remains the record.

Kept by default (VEP terms, matched anywhere in an &-joined list):
    missense_variant, frameshift_variant, stop_gained, stop_lost, start_lost,
    inframe_insertion, inframe_deletion, protein_altering_variant
Dropped by default: intronic, synonymous, splice-site and splice-region,
UTR, upstream/downstream, non-coding.

Usage:
    report_calls.py --calls S.dcs.calls.tsv --out S.dcs.report.tsv \
        [--annotated S.dcs.filtered.tsv] [--indels-annotated S.dcs.indels.annotated.tsv] \
        [--all]
"""

import argparse
import csv
import os
import re
import sys
from urllib.parse import unquote

KEEP = re.compile(r"missense_variant|frameshift_variant|stop_gained|stop_lost|"
                  r"start_lost|inframe_insertion|inframe_deletion|"
                  r"protein_altering_variant")

# annotation columns pulled from the stage 3 tables when present, by any of
# these names (first match wins)
ANNOT_COLUMNS = [
    ("HGVSc",     ["HGVSc"]),
    ("HGVSg",     ["HGVSg"]),
    ("COSMIC",    ["COSMIC_ID", "COSMIC", "cosmic_id", "cosmic"]),
    ("ClinVar",   ["ClinVar", "CLNSIG", "clinvar"]),
    ("gnomAD_AF", ["gnomAD_AF", "gnomAD_exomes_AF", "gnomad_af"]),
    ("rsID",      ["rsID", "Existing_variation", "avsnp"]),
    ("IMPACT",    ["IMPACT"]),
]

CALL_ORDER = {"DETECTED": 0, "not_detected": 1, "not_evaluable": 2}


def split_label(label):
    """'Gene|Consequence|ENSP..:p.X;src=PS' -> gene, consequence, protein, protein_id, source"""
    src = ""
    if ";src=" in label:
        label, src = label.rsplit(";src=", 1)
    parts = label.split("|")
    gene = parts[0] if parts else ""
    consequence = parts[1] if len(parts) > 1 else ""
    protein_id, protein = "", ""
    if len(parts) > 2 and parts[2] not in ("", "-1", "NA"):
        p = unquote(parts[2])
        if ":" in p:
            protein_id, protein = p.split(":", 1)
        else:
            protein = p
    return gene, consequence.replace("&", "; "), protein, protein_id, src


def read_annotation(path, chrom_col, pos_cols, ref_col, alt_col):
    """(chrom, pos, ref, alt) -> {report_col: value} for the columns present."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        if not all(c in cols for c in (chrom_col, ref_col, alt_col)):
            return out
        picks = []
        for report_name, candidates in ANNOT_COLUMNS:
            for c in candidates:
                if c in cols:
                    picks.append((report_name, c))
                    break
        for row in reader:
            vals = {rn: (row.get(c) or "").strip() for rn, c in picks}
            vals = {k: ("" if v in ("-1", ".", "NA") else v) for k, v in vals.items()}
            for pc in pos_cols:
                if pc not in cols:
                    continue
                try:
                    pos = int(row[pc])
                except (TypeError, ValueError):
                    continue
                out.setdefault((row[chrom_col], pos, row[ref_col], row[alt_col]), vals)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--annotated", help="stage 3 substitution table (Chr/Start/End/Ref/Alt ...)")
    ap.add_argument("--indels-annotated", help="stage 3 indel table (chrom/pos/ref/alt ...)")
    ap.add_argument("--all", action="store_true", help="keep every row, not only protein-altering")
    args = ap.parse_args()

    annot = read_annotation(args.annotated, "Chr", ("End", "Start"), "Ref", "Alt")
    for k, v in read_annotation(args.indels_annotated, "chrom", ("pos",), "ref", "alt").items():
        annot.setdefault(k, v)
    annot_names = [rn for rn, _ in ANNOT_COLUMNS]

    with open(args.calls) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols = reader.fieldnames or []
        rows = list(reader)
    for needed in ("label", "chrom", "pos", "ref", "alt", "vaf_pct", "call"):
        if needed not in in_cols:
            sys.exit("ERROR: calls table lacks column %s" % needed)

    passthrough = [c for c in in_cols if c not in ("sample", "label")]
    header = (["sample", "gene", "consequence", "protein_change", "protein_id", "source"]
              + passthrough + annot_names)

    kept, dropped = [], 0
    for r in rows:
        gene, cons, prot, pid, src = split_label(r["label"])
        if not args.all and not KEEP.search(cons):
            dropped += 1
            continue
        try:
            pos = int(r["pos"])
        except ValueError:
            pos = -1
        a = annot.get((r["chrom"], pos, r["ref"], r["alt"]), {})
        rec = [r.get("sample", ""), gene, cons, prot, pid, src]
        rec += [r.get(c, "") for c in passthrough]
        rec += [a.get(n, "") for n in annot_names]
        try:
            vaf = float(r["vaf_pct"])
        except ValueError:
            vaf = -1.0
        kept.append((CALL_ORDER.get(r["call"], 3), -vaf, gene, rec))

    kept.sort(key=lambda t: t[:3])
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for _, _, _, rec in kept:
            w.writerow(rec)

    n_det = sum(1 for t in kept if t[0] == 0)
    sys.stderr.write("%s: %d rows in, %d kept (%d DETECTED), %d dropped as non-protein-altering, "
                     "annotation matched for %d -> %s\n"
                     % (os.path.basename(args.calls), len(rows), len(kept), n_det, dropped,
                        sum(1 for t in kept if any(t[3][len(header) - len(annot_names):])),
                        args.out))


if __name__ == "__main__":
    main()
