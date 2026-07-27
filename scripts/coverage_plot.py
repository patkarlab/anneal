#!/usr/bin/env python3
"""
Coverage plot: SSCS vs DCS median depth per region, split across two stacked
panels so the labels stay readable.

    coverage_plot.py --sample S --bed panel.bed \\
        --sscs S.sscs.sc.sorted.bam --dcs S.dcs.sc.sorted.bam --outdir qc/

--by gene or exon (default exon).
"""

import argparse
import os
import re
from collections import OrderedDict

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402

SSCS_C = "#3b6ea5"
DCS_C = "#d98032"


def region_name(col4, level):
    if not col4 or ";" not in col4:
        return "unknown"
    tag = col4.split(";")[-1].strip()
    tag = tag.split(",")[0]              # probes spanning two annotations
    if level == "exon":
        return tag
    m = re.match(r"^(.*?)_Ex\w*_?\w*$", tag)
    return m.group(1) if m else tag


ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--bed", required=True)
ap.add_argument("--sscs", required=True)
ap.add_argument("--dcs", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--by", choices=["gene", "exon"], default="exon")
ap.add_argument("--min-bq", type=int, default=20)
ap.add_argument("--panels", type=int, default=2)
a = ap.parse_args()

os.makedirs(a.outdir, exist_ok=True)

probes = []
for line in open(a.bed):
    if line.startswith(("#", "track")):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) >= 3:
        probes.append((f[0], int(f[1]), int(f[2]),
                       region_name(f[3] if len(f) > 3 else "", a.by)))

sscs = pysam.AlignmentFile(a.sscs, "rb")
dcs = pysam.AlignmentFile(a.dcs, "rb")

regions = OrderedDict()
for chrom, start, end, name in probes:
    regions.setdefault(name, [[], []])
    for bam, slot in ((sscs, 0), (dcs, 1)):
        try:
            cov = np.array(bam.count_coverage(chrom, start, end,
                                              quality_threshold=a.min_bq),
                           dtype=np.int64).sum(axis=0)
        except ValueError:
            cov = np.zeros(end - start, dtype=np.int64)
        regions[name][slot].append(cov)

sscs.close()
dcs.close()

names = list(regions)
s = np.array([np.median(np.concatenate(regions[n][0])) for n in names])
d = np.array([np.median(np.concatenate(regions[n][1])) for n in names])
s_med, d_med = float(np.median(s)), float(np.median(d))

n_per = int(np.ceil(len(names) / a.panels))
w = 0.4
ymin = max(1, min(d[d > 0].min() if (d > 0).any() else 1, 1) * 0.5)
ymax = s.max() * 2.5

fig, axes = plt.subplots(a.panels, 1,
                         figsize=(max(12, n_per * 0.30), 4.6 * a.panels))
if a.panels == 1:
    axes = [axes]

for pi, ax in enumerate(axes):
    lo, hi = pi * n_per, min((pi + 1) * n_per, len(names))
    if lo >= hi:
        ax.axis("off")
        continue
    sub = names[lo:hi]
    x = np.arange(len(sub))
    ax.bar(x - w / 2, s[lo:hi], w, color=SSCS_C, label="SSCS")
    ax.bar(x + w / 2, d[lo:hi], w, color=DCS_C, label="DCS")
    ax.axhline(s_med, ls="--", lw=1.1, color=SSCS_C)
    ax.axhline(d_med, ls="--", lw=1.1, color=DCS_C)
    ax.set_yscale("log")
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-0.8, n_per - 0.2)
    ax.set_xticks(x)
    ax.set_xticklabels(sub, rotation=90, fontsize=7)
    ax.set_ylabel("median depth")
    ax.grid(axis="y", ls=":", lw=0.6, alpha=0.6)
    if pi == 0:
        ax.legend(fontsize=8, loc="upper right", ncol=2)

fig.suptitle(f"{a.sample}   coverage by {a.by}      "
             f"median SSCS {s_med:,.0f}x      median DCS {d_med:,.0f}x",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])

out = os.path.join(a.outdir, f"{a.sample}.coverage_{a.by}.png")
fig.savefig(out, dpi=140)
print(out)
