#!/usr/bin/env python3
"""
plot_family_sizes.py -- Generate family size distribution plot from Anneal output.

Produces a 4-panel figure:
  A. Full distribution (log-scale y-axis)
  B. Small families 1-20 (linear, highlights singletons)
  C. Cumulative read contribution curve
  D. Summary statistics table

Usage:
    python3 plot_family_sizes.py --input family_sizes.tsv --output family_sizes.png --sample SAMPLE_NAME

For batch overlay of multiple samples:
    python3 plot_family_sizes.py --input-dir /path/to/results/ --output batch_overlay.png --batch
"""

import argparse
import os
import sys

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import numpy as np
except ImportError as e:
    print(f"WARNING: Required library not available ({e}). Skipping plot.", file=sys.stderr)
    sys.exit(0)


def plot_single_sample(input_path, output_path, sample_name):
    """Generate a 4-panel family size distribution plot for one sample."""

    df = pd.read_csv(input_path, sep="\t")
    if "family_size" not in df.columns or "count" not in df.columns:
        print(f"ERROR: Expected columns 'family_size' and 'count' in {input_path}", file=sys.stderr)
        sys.exit(1)

    total_families = df["count"].sum()
    singleton_count = df.loc[df["family_size"] == 1, "count"].values
    singleton_count = singleton_count[0] if len(singleton_count) > 0 else 0
    total_reads = (df["family_size"] * df["count"]).sum()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Anneal 0.1.0 -- Family Size Distribution\n"
        f"({sample_name}, {total_reads/1e6:.0f}M reads, {total_families/1e6:.1f}M families)",
        fontsize=14, fontweight="bold", y=0.98
    )

    # --- Panel A: Full distribution (log scale) ---
    ax = axes[0, 0]
    ax.bar(df["family_size"], df["count"], width=1, color="#2c7bb6", alpha=0.85, edgecolor="none")
    ax.set_yscale("log")
    ax.set_xlabel("Family Size", fontsize=11)
    ax.set_ylabel("Count (log scale)", fontsize=11)
    ax.set_title("A. Full Distribution", fontsize=12, fontweight="bold")
    ax.set_xlim(0, min(200, df["family_size"].max()))
    ax.axvline(x=1, color="red", linestyle="--", alpha=0.5, label="Singletons")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    # --- Panel B: Small families 1-20 (linear) ---
    ax = axes[0, 1]
    small = df[df["family_size"] <= 20].copy()
    colors = ["#d73027" if s == 1 else "#2c7bb6" for s in small["family_size"]]
    ax.bar(small["family_size"], small["count"], color=colors, edgecolor="none", width=0.8)
    ax.set_xlabel("Family Size", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("B. Small Families (1-20)", fontsize=12, fontweight="bold")
    ax.set_xticks(range(1, 21))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.1f}M"))
    ax.grid(axis="y", alpha=0.3)

    singleton_pct = singleton_count / total_families * 100
    ax.annotate(
        f"{singleton_count/1e6:.1f}M\n({singleton_pct:.1f}%)",
        xy=(1, singleton_count), xytext=(4, singleton_count * 0.9),
        fontsize=9, fontweight="bold", color="#d73027",
        arrowprops=dict(arrowstyle="->", color="#d73027", lw=1.2)
    )

    # --- Panel C: Cumulative read contribution ---
    ax = axes[1, 0]
    df["total_reads_contributed"] = df["family_size"] * df["count"]
    df["cumulative_reads"] = df["total_reads_contributed"].cumsum()
    df["cumulative_pct"] = df["cumulative_reads"] / total_reads * 100

    ax.plot(df["family_size"], df["cumulative_pct"], color="#2c7bb6", linewidth=2)
    ax.fill_between(df["family_size"], df["cumulative_pct"], alpha=0.15, color="#2c7bb6")
    ax.set_xlabel("Family Size", fontsize=11)
    ax.set_ylabel("Cumulative % of Total Reads", fontsize=11)
    ax.set_title("C. Cumulative Read Contribution", fontsize=12, fontweight="bold")
    ax.set_xlim(0, min(200, df["family_size"].max()))
    ax.set_ylim(0, 100)
    ax.axhline(y=50, color="gray", linestyle=":", alpha=0.5)
    ax.axhline(y=90, color="gray", linestyle=":", alpha=0.5)
    ax.grid(axis="y", alpha=0.3)

    # --- Panel D: Summary statistics ---
    ax = axes[1, 1]
    ax.axis("off")

    multi = df[df["family_size"] >= 2]
    mean_multi = (
        (multi["family_size"] * multi["count"]).sum() / multi["count"].sum()
        if multi["count"].sum() > 0 else 0
    )
    peak_multi = (
        multi.loc[multi["count"].idxmax(), "family_size"] if len(multi) > 0 else 0
    )
    max_fs = df["family_size"].max()

    table_data = [
        ["Total Reads", f"{total_reads/1e6:.0f}M"],
        ["Total Families", f"{total_families/1e6:.1f}M"],
        ["Singletons", f"{singleton_count/1e6:.1f}M ({singleton_pct:.1f}%)"],
        ["Mean Family (>1)", f"{mean_multi:.0f}"],
        ["Peak Family (>1)", f"{int(peak_multi)}"],
        ["Max Family Size", f"{max_fs}"],
    ]

    table = ax.table(
        cellText=table_data,
        colLabels=["Metric", "Value"],
        loc="center", cellLoc="center", colWidths=[0.45, 0.35]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.6)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#2c7bb6")
            cell.set_text_props(color="white", fontweight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#f0f4f8")

    ax.set_title("D. Summary Statistics", fontsize=12, fontweight="bold", pad=20)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close()
    print(f"Plot saved: {output_path}")


def plot_batch_overlay(input_dir, output_path):
    """Overlay family size distributions from multiple samples."""

    samples = {}
    for sample_name in sorted(os.listdir(input_dir)):
        consensus_dir = os.path.join(input_dir, sample_name, "consensus")
        if not os.path.isdir(consensus_dir):
            continue
        # Look for {sample}.family_sizes.tsv (preferred) or family_sizes.tsv (legacy)
        tsv = os.path.join(consensus_dir, f"{sample_name}.family_sizes.tsv")
        if not os.path.isfile(tsv):
            tsv = os.path.join(consensus_dir, "family_sizes.tsv")
        if os.path.isfile(tsv):
            samples[sample_name] = pd.read_csv(tsv, sep="\t")

    if not samples:
        print(f"ERROR: No family_sizes.tsv files found under {input_dir}", file=sys.stderr)
        sys.exit(1)

    colors = plt.cm.tab20(np.linspace(0, 1, max(len(samples), 2)))

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        f"Anneal 0.1.0 -- Family Size Distribution Overlay ({len(samples)} samples)",
        fontsize=14, fontweight="bold", y=0.98
    )

    # Panel A: overlay distributions (log)
    ax = axes[0, 0]
    for i, (name, df) in enumerate(samples.items()):
        short = name.replace("25NGS", "").replace("-Duplex", "")
        ax.plot(df["family_size"], df["count"], color=colors[i], alpha=0.7, linewidth=1, label=short)
    ax.set_yscale("log")
    ax.set_xlim(0, 200)
    ax.set_xlabel("Family Size")
    ax.set_ylabel("Count (log)")
    ax.set_title("A. All Samples (log scale)", fontweight="bold")
    ax.legend(fontsize=7, ncol=3, loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    # Panel B: overlay distributions (linear, 1-100)
    ax = axes[0, 1]
    for i, (name, df) in enumerate(samples.items()):
        short = name.replace("25NGS", "").replace("-Duplex", "")
        sub = df[(df["family_size"] >= 2) & (df["family_size"] <= 100)]
        ax.plot(sub["family_size"], sub["count"], color=colors[i], alpha=0.7, linewidth=1, label=short)
    ax.set_xlabel("Family Size")
    ax.set_ylabel("Count")
    ax.set_title("B. Multi-read Families (2-100)", fontweight="bold")
    ax.legend(fontsize=7, ncol=3, loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    # Panel C: cumulative read contribution
    ax = axes[1, 0]
    for i, (name, df) in enumerate(samples.items()):
        short = name.replace("25NGS", "").replace("-Duplex", "")
        reads = df["family_size"] * df["count"]
        cum = reads.cumsum() / reads.sum() * 100
        ax.plot(df["family_size"], cum, color=colors[i], alpha=0.7, linewidth=1, label=short)
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 100)
    ax.axhline(y=50, color="gray", linestyle=":", alpha=0.5)
    ax.axhline(y=90, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Family Size")
    ax.set_ylabel("Cumulative % Reads")
    ax.set_title("C. Cumulative Read Contribution", fontweight="bold")
    ax.legend(fontsize=7, ncol=3, loc="lower right")
    ax.grid(axis="y", alpha=0.3)

    # Panel D: summary table
    ax = axes[1, 1]
    ax.axis("off")

    table_data = []
    for i, (name, df) in enumerate(samples.items()):
        short = name.replace("25NGS", "").replace("-Duplex", "")
        total_fam = df["count"].sum()
        sing = df.loc[df["family_size"] == 1, "count"].values
        sing = sing[0] if len(sing) > 0 else 0
        sing_pct = sing / total_fam * 100 if total_fam > 0 else 0
        total_reads = (df["family_size"] * df["count"]).sum()
        multi = df[df["family_size"] >= 2]
        mean_m = (multi["family_size"] * multi["count"]).sum() / multi["count"].sum() if multi["count"].sum() > 0 else 0
        table_data.append([short, f"{total_reads/1e6:.0f}M", f"{total_fam/1e6:.1f}M", f"{sing_pct:.0f}%", f"{mean_m:.0f}"])

    table = ax.table(
        cellText=table_data,
        colLabels=["Sample", "Reads", "Families", "Sing%", "Mean(>1)"],
        loc="center", cellLoc="center", colWidths=[0.2, 0.15, 0.15, 0.12, 0.15]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.35)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#2c7bb6")
            cell.set_text_props(color="white", fontweight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#f0f4f8")
        if col == 0 and row > 0 and row - 1 < len(colors):
            cell.set_text_props(color=colors[row - 1], fontweight="bold")

    ax.set_title("D. Summary Statistics", fontsize=12, fontweight="bold", pad=20)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close()
    print(f"Batch overlay saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot Anneal family size distributions")
    parser.add_argument("--input", help="Single family_sizes.tsv file")
    parser.add_argument("--input-dir", help="Batch results directory (contains sample subdirs)")
    parser.add_argument("--output", required=True, help="Output PNG path")
    parser.add_argument("--sample", default="Sample", help="Sample name for title")
    parser.add_argument("--batch", action="store_true", help="Batch overlay mode (use --input-dir)")

    args = parser.parse_args()

    if args.batch:
        if not args.input_dir:
            parser.error("--batch requires --input-dir")
        plot_batch_overlay(args.input_dir, args.output)
    else:
        if not args.input:
            parser.error("Single-sample mode requires --input")
        plot_single_sample(args.input, args.output, args.sample)


if __name__ == "__main__":
    main()
