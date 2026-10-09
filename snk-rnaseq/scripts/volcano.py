#!/usr/bin/env python3
"""Volcano and MA plots of a PyDESeq2 results table (port extension).

Usage: volcano.py --results <c>.deseq2.results.tsv --title <c> --volcano <pdf>
                  --volcano-png <png> --ma <pdf> [--padj 0.05] [--lfc 1] [--label-top 15]
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

UP, DOWN, NS = "#c0392b", "#2c6fbb", "#b0b0b0"


def classify(res, padj, lfc):
    sig = (res["padj"] < padj) & (res["log2FoldChange"].abs() >= lfc)
    return np.where(sig & (res["log2FoldChange"] > 0), "up", np.where(sig, "down", "ns"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--volcano", required=True)
    ap.add_argument("--volcano-png", required=True)
    ap.add_argument("--ma", required=True)
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--label-top", type=int, default=15)
    args = ap.parse_args()

    res = pd.read_csv(args.results, sep="\t")
    res = res[res["log2FoldChange"].notna()].copy()
    res["label"] = res["gene_name"].fillna(res["gene_id"]).astype(str)
    res["cls"] = classify(res, args.padj, args.lfc)
    colors = {"up": UP, "down": DOWN, "ns": NS}
    n_up, n_down = int((res["cls"] == "up").sum()), int((res["cls"] == "down").sum())

    # Volcano
    v = res[res["padj"].notna()].copy()
    floor = v.loc[v["padj"] > 0, "padj"].min() if (v["padj"] > 0).any() else 1e-300
    v["mlog10"] = -np.log10(v["padj"].clip(lower=floor))
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for cls in ("ns", "down", "up"):
        sub = v[v["cls"] == cls]
        ax.scatter(sub["log2FoldChange"], sub["mlog10"], s=6, c=colors[cls], alpha=0.7, linewidths=0,
                   label={"up": f"up ({n_up})", "down": f"down ({n_down})", "ns": "n.s."}[cls], rasterized=True)
    ax.axhline(-np.log10(args.padj), color="#555555", lw=0.6, ls="--")
    for x in (-args.lfc, args.lfc):
        ax.axvline(x, color="#555555", lw=0.6, ls="--")
    top = v[v["cls"] != "ns"].sort_values("padj").head(args.label_top)
    for _, r in top.iterrows():
        ax.annotate(r["label"], (r["log2FoldChange"], r["mlog10"]), fontsize=6,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(args.title.replace("_vs_", " vs "))
    ax.legend(fontsize=8, frameon=False, markerscale=2)
    fig.tight_layout()
    fig.savefig(args.volcano)
    fig.savefig(args.volcano_png, dpi=200)
    plt.close(fig)

    # MA
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for cls in ("ns", "down", "up"):
        sub = res[res["cls"] == cls]
        ax.scatter(np.log10(sub["baseMean"].clip(lower=1e-1)), sub["log2FoldChange"], s=6, c=colors[cls],
                   alpha=0.7, linewidths=0, rasterized=True)
    ax.axhline(0, color="#555555", lw=0.6)
    ax.set_xlabel("log10 mean of normalized counts")
    ax.set_ylabel("log2 fold change")
    ax.set_title(args.title.replace("_vs_", " vs ") + " (MA)")
    fig.tight_layout()
    fig.savefig(args.ma)
    plt.close(fig)


if __name__ == "__main__":
    main()
