#!/usr/bin/env python3
"""
Panel coverage QC for the anneal duplex pipeline.

Reports SSCS and DCS depth gene by gene and exon by exon, with the limit of
detection each region supports. Sample-level consensus metrics are parsed from
the stage 1 stats file rather than recomputed.

Probe coordinates come from BED columns 2-3. Column 4 is a legacy hg19 label;
only the gene and exon name after the semicolon is used.

Outputs:
    {sample}.qc_genes.tsv    one row per gene
    {sample}.qc_exons.tsv    one row per gene/exon
    {sample}.qc_probes.tsv   one row per probe
    {sample}.qc_summary.tsv  one row per sample
    {sample}.qc_genes.png    gene-level coverage, SSCS vs DCS
    {sample}.qc_exons.png    exon-level coverage, SSCS vs DCS

    qc_panel.py --sample S --bed panel.bed \\
        --sscs S.sscs.sc.sorted.bam --dcs S.dcs.sc.sorted.bam \\
        --stats S.stats.txt --outdir qc/ --target-lod 0.1
"""

import argparse
import os
import re
import sys
from collections import OrderedDict

import numpy as np
import pysam

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402

SSCS_C = "#3b6ea5"
DCS_C = "#d98032"
FAIL_C = "#c0392b"


def parse_label(col4):
    """'13:28592540-28592660;FLT3_Ex_20' -> ('FLT3', 'Ex_20')."""
    if not col4 or ";" not in col4:
        return "unknown", ""
    tag = col4.split(";")[-1].strip()
    m = re.match(r"^(.*?)_(Ex_\w+)$", tag)
    if m:
        return m.group(1), m.group(2)
    return tag, ""


def load_probes(bed):
    probes = []
    with open(bed) as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            gene, exon = parse_label(f[3] if len(f) > 3 else "")
            probes.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                           "gene": gene, "exon": exon})
    return probes


def probe_depths(bam, probes, min_bq):
    out = []
    for p in probes:
        try:
            cov = bam.count_coverage(p["chrom"], p["start"], p["end"],
                                     quality_threshold=min_bq)
            per_base = np.array(cov, dtype=np.int64).sum(axis=0)
        except ValueError:
            per_base = np.zeros(p["end"] - p["start"], dtype=np.int64)
        out.append(per_base)
    return out


def parse_stats(path):
    stats = {}
    if not path or not os.path.exists(path):
        return stats
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        k, v = line.split(":", 1)
        v = v.strip()
        if v:
            m = re.match(r"^([\d.]+)", v.replace(",", ""))
            stats[k.strip()] = m.group(1) if m else v
    return stats


def lod_of(depth, min_alt):
    return (100.0 * min_alt / depth) if depth > 0 else float("inf")


def write_tsv(path, rows):
    cols = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def coverage_figure(path, names, sscs, dcs, req_depth, title, target_lod,
                    min_alt, row_h):
    """Horizontal grouped bars, one row per region, worst coverage at top."""
    n = len(names)
    order = np.argsort(dcs)                       # worst first
    names = [names[i] for i in order]
    sscs, dcs = sscs[order], dcs[order]
    fail = dcs < req_depth

    y = np.arange(n)
    h = 0.38
    fig, ax = plt.subplots(figsize=(11, max(4.0, n * row_h)))

    ax.barh(y + h / 2, np.maximum(sscs, 1), height=h, color=SSCS_C,
            label="SSCS", zorder=3)
    ax.barh(y - h / 2, np.maximum(dcs, 1), height=h,
            color=[FAIL_C if f else DCS_C for f in fail],
            label="DCS", zorder=3)

    ax.axvline(req_depth, ls="--", lw=1.2, color="0.25", zorder=4,
               label=f"{target_lod}% LoD needs {req_depth:.0f}x DCS")

    for i, (d, f) in enumerate(zip(dcs, fail)):
        ax.text(max(d, 1) * 1.12, i - h / 2,
                f"{d:.0f}x  ({lod_of(d, min_alt):.2f}%)" if d > 0 else "0x",
                va="center", fontsize=7,
                color=FAIL_C if f else "0.3", zorder=5)

    ax.set_xscale("log")
    ax.set_xlim(left=10, right=max(sscs.max(), req_depth) * 4)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    for lab, f in zip(ax.get_yticklabels(), fail):
        if f:
            lab.set_color(FAIL_C)
            lab.set_fontweight("bold")
    ax.set_xlabel("median depth (log scale)")
    ax.set_title(title, fontsize=12)
    ax.grid(axis="x", ls=":", lw=0.6, alpha=0.6, zorder=0)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--sscs", required=True)
    ap.add_argument("--dcs", required=True)
    ap.add_argument("--stats")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-alt", type=int, default=2,
                    help="alt reads assumed for the LoD calculation")
    ap.add_argument("--target-lod", type=float, default=0.1,
                    help="LoD in percent a region must reach to pass")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    probes = load_probes(args.bed)
    if not probes:
        sys.exit("FATAL: no probes parsed from BED")

    sscs_bam = pysam.AlignmentFile(args.sscs, "rb")
    dcs_bam = pysam.AlignmentFile(args.dcs, "rb")
    s_cov = probe_depths(sscs_bam, probes, args.min_bq)
    d_cov = probe_depths(dcs_bam, probes, args.min_bq)
    sscs_bam.close()
    dcs_bam.close()

    req = 100.0 * args.min_alt / args.target_lod

    # ---- per probe ----
    probe_rows = []
    for p, s, d in zip(probes, s_cov, d_cov):
        sm, dm = float(np.median(s)), float(np.median(d))
        probe_rows.append({
            "chrom": p["chrom"], "start": p["start"], "end": p["end"],
            "gene": p["gene"], "exon": p["exon"],
            "sscs_median": round(sm, 1), "dcs_median": round(dm, 1),
            "dcs_sscs_ratio": round(dm / sm, 4) if sm > 0 else 0.0,
            "lod_pct": round(lod_of(dm, args.min_alt), 4) if dm > 0 else "NA",
            "pass": "PASS" if dm >= req else "FAIL"})

    # ---- per exon and per gene: pool bases, then take the median ----
    def aggregate(keyfn):
        buckets = OrderedDict()
        for p, s, d in zip(probes, s_cov, d_cov):
            k = keyfn(p)
            buckets.setdefault(k, [[], []])
            buckets[k][0].append(s)
            buckets[k][1].append(d)
        rows = []
        for k, (ss, dd) in buckets.items():
            sm = float(np.median(np.concatenate(ss)))
            dm = float(np.median(np.concatenate(dd)))
            rows.append({"region": k, "n_probes": len(ss),
                         "sscs_median": round(sm, 1), "dcs_median": round(dm, 1),
                         "dcs_sscs_ratio": round(dm / sm, 4) if sm > 0 else 0.0,
                         "lod_pct": round(lod_of(dm, args.min_alt), 4) if dm > 0 else "NA",
                         "pass": "PASS" if dm >= req else "FAIL"})
        return rows

    gene_rows = aggregate(lambda p: p["gene"])
    exon_rows = aggregate(lambda p: f"{p['gene']}_{p['exon']}" if p["exon"]
                          else p["gene"])

    write_tsv(os.path.join(args.outdir, f"{args.sample}.qc_probes.tsv"), probe_rows)
    write_tsv(os.path.join(args.outdir, f"{args.sample}.qc_genes.tsv"), gene_rows)
    write_tsv(os.path.join(args.outdir, f"{args.sample}.qc_exons.tsv"), exon_rows)

    # ---- figures ----
    coverage_figure(
        os.path.join(args.outdir, f"{args.sample}.qc_genes.png"),
        [r["region"] for r in gene_rows],
        np.array([r["sscs_median"] for r in gene_rows]),
        np.array([r["dcs_median"] for r in gene_rows]),
        req, f"{args.sample}  coverage by gene", args.target_lod,
        args.min_alt, row_h=0.30)

    coverage_figure(
        os.path.join(args.outdir, f"{args.sample}.qc_exons.png"),
        [r["region"] for r in exon_rows],
        np.array([r["sscs_median"] for r in exon_rows]),
        np.array([r["dcs_median"] for r in exon_rows]),
        req, f"{args.sample}  coverage by exon", args.target_lod,
        args.min_alt, row_h=0.20)

    # ---- summary ----
    stats = parse_stats(args.stats)
    g_fail = [r["region"] for r in gene_rows if r["pass"] == "FAIL"]
    e_fail = [r["region"] for r in exon_rows if r["pass"] == "FAIL"]
    d_all = np.concatenate(d_cov)
    s_all = np.concatenate(s_cov)

    summary = {
        "sample": args.sample,
        "genes": len(gene_rows), "genes_failed": len(g_fail),
        "exons": len(exon_rows), "exons_failed": len(e_fail),
        "probes": len(probe_rows),
        "sscs_median_depth": round(float(np.median(s_all)), 1),
        "dcs_median_depth": round(float(np.median(d_all)), 1),
        "dcs_sscs_ratio": round(float(np.median(d_all) / np.median(s_all)), 4)
                          if np.median(s_all) > 0 else 0.0,
        "panel_median_lod_pct": round(lod_of(float(np.median(d_all)), args.min_alt), 4),
        "target_lod_pct": args.target_lod,
        "singleton_rate_pct": stats.get("Singleton rate", "NA"),
        "sscs_efficiency_pct": stats.get("SSCS efficiency", "NA"),
        "dcs_efficiency_pct": stats.get("DCS efficiency", "NA"),
        "dcs_recovery_pct": stats.get("DCS recovery", "NA"),
        "verdict": "PASS" if not g_fail else f"REVIEW ({len(g_fail)} genes)",
    }
    write_tsv(os.path.join(args.outdir, f"{args.sample}.qc_summary.tsv"), [summary])

    sys.stderr.write(f"\n{args.sample}\n")
    for k in ("genes", "genes_failed", "exons_failed", "sscs_median_depth",
              "dcs_median_depth", "panel_median_lod_pct", "dcs_recovery_pct",
              "verdict"):
        sys.stderr.write(f"  {k:22} {summary[k]}\n")
    if g_fail:
        sys.stderr.write(f"\n  genes below {args.target_lod}% LoD:\n")
        for r in sorted((r for r in gene_rows if r["pass"] == "FAIL"),
                        key=lambda x: x["dcs_median"]):
            sys.stderr.write(f"    {r['region']:14} DCS {r['dcs_median']:>8.0f}x"
                             f"   LoD {r['lod_pct']}%\n")
    sys.stderr.write(f"\nwrote {args.outdir}/{args.sample}.qc_*.tsv and .png\n")


if __name__ == "__main__":
    main()
