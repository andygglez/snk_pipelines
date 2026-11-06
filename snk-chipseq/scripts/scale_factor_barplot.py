#!/usr/bin/env python3
"""Bar plot of the samplesheet scale factors of one protein (port extension).

One bar per condition: the mean scale factor of its replicates, with an error
bar (SD or SEM) and every replicate drawn as a dot. Also writes the values as a
TSV (condition, rep, scale_factor) next to the plot.

Usage:
  scale_factor_barplot.py --protein H3K9me3 --value <condition>:<rep>:<factor> [...]
      --out plot.pdf --tsv values.tsv [--order c1,c2,...] [--color c1=#5677cd ...]
      [--palette #5677cd,#000000,...] [--title "{protein}"] [--xlabel ""]
      [--ylabel "Scale factor"] [--error sd|sem]
"""

import argparse
import csv
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--value", action="append", required=True, metavar="CONDITION:REP:FACTOR")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--order", default="", help="comma-separated conditions (x axis order)")
    ap.add_argument("--color", action="append", default=[], metavar="CONDITION=COLOR")
    ap.add_argument("--palette", default="", help="comma-separated colours for conditions without --color")
    ap.add_argument("--title", default="{protein}")
    ap.add_argument("--xlabel", default="")
    ap.add_argument("--ylabel", default="Scale factor")
    ap.add_argument("--error", choices=["sd", "sem"], default="sd")
    args = ap.parse_args()

    values = {}  # condition -> [(rep, factor)]
    for v in args.value:
        cond, rep, factor = v.rsplit(":", 2)
        values.setdefault(cond, []).append((rep, float(factor)))

    order = [c for c in args.order.split(",") if c] if args.order else []
    unknown = [c for c in order if c not in values]
    if unknown:
        print(f"warning: conditions {unknown} in the barplot order have no {args.protein} samples; skipped")
    conditions = [c for c in order if c in values] + [c for c in values if c not in order]

    colors = dict(c.split("=", 1) for c in args.color)
    palette = [c for c in args.palette.split(",") if c] or plt.rcParams["axes.prop_cycle"].by_key()["color"]
    free = iter(palette * (len(conditions) // max(len(palette), 1) + 1))
    bar_colors = [colors[c] if c in colors else next(free) for c in conditions]

    means, errors = [], []
    for c in conditions:
        f = [x for _, x in values[c]]
        m = sum(f) / len(f)
        sd = math.sqrt(sum((x - m) ** 2 for x in f) / (len(f) - 1)) if len(f) > 1 else 0.0
        means.append(m)
        errors.append(sd / math.sqrt(len(f)) if args.error == "sem" else sd)

    fig, ax = plt.subplots()
    x = range(len(conditions))
    ax.bar(x, means, color=bar_colors, edgecolor="black", linewidth=1, zorder=1)
    ax.errorbar(x, means, yerr=errors, fmt="none", ecolor="black", elinewidth=2.5, capsize=8,
                capthick=2.5, zorder=2)
    for i, c in enumerate(conditions):
        f = [v for _, v in values[c]]
        n = len(f)
        offsets = [(-0.12 + 0.24 * k / (n - 1)) if n > 1 else 0.0 for k in range(n)]
        for k, (dx, y) in enumerate(zip(offsets, f)):
            # alternate filled and open dots, as in the MS bar plots
            face, edge = ("black", "white") if k % 2 == 0 else ("white", "black")
            ax.scatter(i + dx, y, s=40, facecolor=face, edgecolor=edge, linewidths=1, zorder=3)
    ax.set_xticks(list(x))
    ax.set_xticklabels(conditions)
    ax.set_title(args.title.replace("{protein}", args.protein))
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    ax.set_ylim(bottom=0)
    fig.savefig(args.out)
    plt.close(fig)

    with open(args.tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["condition", "rep", "scale_factor"])
        for c in conditions:
            for rep, f in values[c]:
                w.writerow([c, rep, f])


if __name__ == "__main__":
    main()
