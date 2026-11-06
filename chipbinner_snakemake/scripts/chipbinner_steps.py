#!/usr/bin/env python3
"""One subcommand per Snakefile rule, built on the Python ChIPbinner port (scripts/chipbinner/).

The R scripts of andygglez/snakemake_chipbinner map to these subcommands:
  filter.low.counts.R            -> (not ported: no low-count filtering option)
  chipbinner.norm.bw.R           -> norm-bw
  merge.replicates.R / cp        -> final-bw
  genic.intergenic.plot.R +
  compare.genic.vs.intergenic.py -> compare
  clus.py                        -> clustering
  annotate.clusters.R            -> annotate   (BED instead of .rda)
  density.scatter.plot.R         -> density
  PCA.R                          -> pca / correlation
  (ChIPbinner extras)            -> differential, enrichment, legend, binned-input
"""
import argparse
import json
import logging
import os
import shutil
import string
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from chipbinner import cluster, differential, enrichment, explore, extra_plots, normalize, plots  # noqa: E402
from chipbinner.genomics import read_bins, read_chrom_sizes, write_bed, write_bigwig  # noqa: E402

log = logging.getLogger("chipbinner")
PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf",
           "#7f7f7f", "#bcbd22"]


def save(fig, prefix, formats):
    for fmt in formats:
        path = f"{prefix}.{fmt}"
        if fmt == "html":
            fig.write_html(path, include_plotlyjs=True, config={"displaylogo": False})
        else:
            fig.write_image(path, width=fig.layout.width, height=fig.layout.height,
                            scale=2 if fmt in ("png", "jpeg", "jpg") else 1)


def touch(*paths):
    for p in paths:
        open(p, "w").close()


def _bool(s):
    return str(s).lower() in ("1", "true", "t", "yes", "y")


# --------------------------------------------------------------------------- per replicate
def cmd_binned_input(a):
    """A samplesheet cell that is already a binned-count BED: copy it with integer counts."""
    df = read_bins(a.bed)
    df["score"] = df["score"].round().astype(np.int64)
    write_bed(df, a.out)


def cmd_norm_bw(a):
    df = normalize.norm_bw(a.binned, a.chrom_sizes, a.blacklist, input_binned_file=a.input or None,
                           use_input=bool(a.input), depth_norm=_bool(a.depth_norm), pseudocount=float(a.pseudocount),
                           raw_count_cutoff=float(a.raw_count_cutoff), scaling_factor=float(a.scaling_factor))
    write_bigwig(df, a.out, read_chrom_sizes(a.chrom_sizes))


def cmd_final_bw(a):
    if len(a.reps) == 2:
        write_bigwig(normalize.merge_norm_bw(*a.reps), a.out, read_chrom_sizes(a.chrom_sizes))
        log.info("merging for %s and %s", *a.reps)
    else:
        shutil.copyfile(a.reps[0], a.out)
        log.info("not merging for %s", a.reps[0])


# --------------------------------------------------------------------------- per comparison
def _range(pair, spec, hexbins):
    """Resolve the plot spec (auto / [min, max] / {x, y}) to [min_x, max_x, min_y, max_y].

    auto follows ChIPbinner::genic_intergenic_scatterplot: per axis, the extent of the hex
    centres holding > 10 bins on the scatterplot's own hex grid (plots.auto_range)."""
    auto = plots.auto_range(pair, xbins=int(hexbins))
    if isinstance(spec["range"], str) and spec["range"].lower() == "auto":
        return auto, "auto"
    rng = list(plots.resolve_range(spec["range"]))
    src = "config" if all(v is not None for v in rng) else "config+auto"  # {x: [..], y: auto}
    return [a if v is None else v for v, a in zip(rng, auto)], src


def cmd_compare(a):
    """genic.intergenic.plot.R + compare.genic.vs.intergenic.py: scatterplot, pre_clust, extra panels."""
    spec, g = json.loads(a.plot_spec), json.loads(a.genic_plot)
    pair = cluster.load_pair(a.treated, a.control)
    if len(pair) < 10:
        raise SystemExit(f"{a.prefix}: only {len(pair)} bins scored in both samples")
    rng, src = _range(pair, spec, g["hexbins"])
    with open(f"{a.prefix}.plot_range.yaml", "w") as fh:
        yaml.safe_dump({"comparison": f"{a.t_label}_{a.c_label}", "mark": a.mark, "source": src,
                        "pow": float(spec["pow"]), "x": [round(rng[0], 4), round(rng[1], 4)],
                        "y": [round(rng[2], 4), round(rng[3], 4)]}, fh, sort_keys=False)
    log.info("Pos range x=%s y=%s (%s), pow=%s", rng[:2], rng[2:], src, spec["pow"])
    fig = plots.genic_intergenic_scatterplot(
        a.treated, a.control, gene_bed=a.gene_bed, igr_bed=a.igr_bed, title_of_plot=f"Change in {a.mark}",
        xaxis_label=a.c_label, yaxis_label=a.t_label, plot_range={"x": rng[:2], "y": rng[2:]},
        pow=float(spec["pow"]), hexbins=int(g["hexbins"]), show_scales=_bool(g["show_scales"]),
        show_legend=_bool(g["show_legend"]), legend_pos=g["legend_pos"])
    save(fig, f"{a.prefix}.scatterplot", a.formats)
    cluster.pre_clust(a.treated, a.control, f"{a.prefix}_preclust")
    save(extra_plots.smooth_scatter(pair["x"], pair["y"], a.c_label, a.t_label), f"{a.prefix}_preclust_smoothScatter",
         ["jpeg"])
    save(extra_plots.genic_vs_intergenic(pair, a.gene_bed, a.igr_bed, a.c_label, a.t_label),
         f"{a.prefix}.genic.vs.intergenic", ["png"])


def cmd_legend(a):
    save(extra_plots.genic_legend(), a.out[: -len(".png")], ["png"])


def cmd_clustering(a):
    """clus.py: HDBSCAN(min_cluster_size=minpts, min_samples=minsamps) -> labels .txt + .png + log."""
    mat = np.loadtxt(a.matrix, delimiter=",", ndmin=2)
    t0 = time.time()
    try:
        try:
            import hdbscan

            m = hdbscan.HDBSCAN(min_cluster_size=a.minpts, min_samples=a.minsamps, core_dist_n_jobs=a.cores,
                                cluster_selection_method=a.method).fit(mat)
            impl = "hdbscan"
        except ImportError:
            from sklearn.cluster import HDBSCAN

            m = HDBSCAN(min_cluster_size=a.minpts, min_samples=a.minsamps, n_jobs=a.cores,
                        cluster_selection_method=a.method).fit(mat)
            impl = "sklearn.cluster.HDBSCAN"
        labels, probs = m.labels_.astype(int), m.probabilities_
    except Exception as e:  # e.g. min_samples > number of bins
        labels, probs, impl = np.full(len(mat), -1), None, f"failed: {e}"
    np.savetxt(a.labels_out, labels, fmt="%i")
    with open(a.log, "w") as fh:
        fh.write(f"INFO:root:Starting dbscan ({impl}); minpts={a.minpts} minsamps={a.minsamps} n={len(mat)}\n"
                 f"INFO:root:clusters={len(set(labels) - {-1})} noise={(labels == -1).sum()}\n"
                 f"INFO:root:Hdbscan done in {time.time() - t0:.1f}s\n")
    save(extra_plots.cluster_scatter(mat, labels, probs, a.c_label, a.t_label), a.png_out[: -len(".png")], ["png"])


def cmd_annotate(a):
    """annotate.clusters.R: name clusters A, B, C by mean signal; empty file when too few clusters."""
    labels = np.loadtxt(a.labels, dtype=int, ndmin=1)
    found = len(set(labels) - {-1})
    if found < a.n:
        log.info("%d clusters < %d: empty output", found, a.n)
        touch(a.out)
        return
    if found > a.n:
        log.warning("HDBSCAN found %d clusters, keeping the %d with the lowest mean signal (as ChIPbinner); "
                    "raise minpts/minsamps for fewer clusters", found, a.n)
    write_bed(cluster.annotate_clust(a.matrix, a.pooled, a.labels, a.n), a.out)


def cmd_density(a):
    """density.scatter.plot.R: density-based scatterplot; empty placeholder when there are no clusters."""
    prefix = a.prefix
    if not os.path.exists(a.clusters) or os.path.getsize(a.clusters) == 0:
        touch(*(f"{prefix}.{fmt}" for fmt in a.formats))
        return
    with open(a.range_yaml) as fh:
        rng = yaml.safe_load(fh)
    d = json.loads(a.density_plot)
    fig = plots.density_based_scatterplot(
        a.treated, a.control, a.clusters, gene_bed=a.gene_bed, igr_bed=a.igr_bed, number_of_clusters=a.n,
        title_of_plot=a.mark, xaxis_label=a.c_label, yaxis_label=a.t_label,
        plot_range={"x": rng["x"], "y": rng["y"]}, pow=float(rng["pow"]), hexbins=int(d["hexbins"]),
        show_scales=_bool(d["show_scales"]), include_additional_density_plot=_bool(d["include_additional_density_plot"]),
        filter_extreme_bins=_bool(d["filter_extreme_bins"]), cluster_colors=d["cluster_colors"])
    save(fig, prefix, a.formats)


# --------------------------------------------------------------------------- exploratory
def _colors(conditions):
    cmap = {c: PALETTE[k % len(PALETTE)] for k, c in enumerate(dict.fromkeys(conditions))}
    return [cmap[c] for c in conditions]


def cmd_pca(a):
    fig, _ = explore.plot_pca(a.bws, _colors(a.conditions), a.labels, plot_title=f"{a.mark} · PCA")
    save(fig, a.prefix, a.formats)


def cmd_correlation(a):
    fig, _ = explore.plot_correlation(a.bws, _colors(a.conditions), a.labels, a.method,
                                      plot_title=f"{a.mark} · {a.method} correlation")
    save(fig, a.prefix, a.formats)


# --------------------------------------------------------------------------- ChIPbinner extras
def cmd_differential(a):
    clus = a.clusters if a.clusters and os.path.getsize(a.clusters) > 0 else None
    allb, per, summ, _ = differential.differential_bin_analysis(
        a.treated, a.control, a.blacklist, a.t_label, a.c_label, annotated_clusters=clus,
        bootstrap_value=a.bootstrap, K_value=a.K)
    allb.to_csv(f"{a.prefix}.stats_for_all_bins_no_clustering.csv", index=False)
    if per is not None:
        per.to_csv(f"{a.prefix}.stats_per_cluster.csv", index=False)
    else:
        touch(f"{a.prefix}.stats_per_cluster.csv")
    with open(f"{a.prefix}.ROTS_output_summary.txt", "w") as fh:
        fh.write(summ)


def cmd_enrichment(a):
    letters = list(string.ascii_uppercase[: a.n])
    empty = not os.path.exists(a.clusters) or os.path.getsize(a.clusters) == 0
    found = set() if empty else set(pd.read_csv(a.clusters, sep="\t", header=None)[3].astype(str))
    for letter in letters:
        for region in a.regions:
            stem = f"{a.prefix}.cluster{letter}.{region}"
            agg = None
            if letter in found:
                agg = enrichment.enrich_clust(a.clusters, letter, a.pooled, a.db_dir, a.gene_bed, a.igr_bed,
                                              region=region)
            if agg is None:
                log.info("cluster %s %s: no significant enrichment/depletion", letter, region)
                touch(f"{stem}.csv")
                continue
            agg.to_csv(f"{stem}.csv", index=False)
            save(enrichment.plot_enrichment(agg, a.n_elements, a.cutoff_for_overlap, f"Cluster {letter} · {region}"),
                 stem, a.formats)


# --------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def sp(name, fn):
        p = sub.add_parser(name)
        p.set_defaults(fn=fn)
        return p

    p = sp("binned-input", cmd_binned_input)
    p.add_argument("--bed", required=True)
    p.add_argument("--out", required=True)

    p = sp("norm-bw", cmd_norm_bw)
    p.add_argument("--binned", required=True)
    p.add_argument("--input", default="")
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--blacklist", required=True)
    p.add_argument("--depth-norm", default="true")
    p.add_argument("--pseudocount", default="1e-15")
    p.add_argument("--raw-count-cutoff", default="100")
    p.add_argument("--scaling-factor", default="1")
    p.add_argument("--out", required=True)

    p = sp("final-bw", cmd_final_bw)
    p.add_argument("--reps", nargs="+", required=True)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--out", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--treated", nargs="+")
    common.add_argument("--control", nargs="+")
    common.add_argument("--t-label")
    common.add_argument("--c-label")
    common.add_argument("--mark")
    common.add_argument("--gene-bed")
    common.add_argument("--igr-bed")
    common.add_argument("--prefix")
    common.add_argument("--formats", nargs="+", default=["pdf", "html"])

    def sp2(name, fn):
        p = sub.add_parser(name, parents=[common])
        p.set_defaults(fn=fn)
        return p

    p = sp2("compare", cmd_compare)
    p.add_argument("--plot-spec", required=True, help="JSON {range, pow}")
    p.add_argument("--genic-plot", required=True, help="JSON genic_plot config")

    p = sp("legend", cmd_legend)
    p.add_argument("--out", required=True)

    p = sp2("clustering", cmd_clustering)
    p.add_argument("--matrix", required=True)
    p.add_argument("--minpts", type=int, required=True)
    p.add_argument("--minsamps", type=int, required=True)
    p.add_argument("--method", default="eom")
    p.add_argument("--cores", type=int, default=1)
    p.add_argument("--labels-out", required=True)
    p.add_argument("--png-out", required=True)
    p.add_argument("--log", required=True)

    p = sp("annotate", cmd_annotate)
    p.add_argument("--matrix", required=True)
    p.add_argument("--pooled", required=True)
    p.add_argument("--labels", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--out", required=True)

    p = sp2("density", cmd_density)
    p.add_argument("--clusters", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--range-yaml", required=True)
    p.add_argument("--density-plot", required=True, help="JSON density_plot config")

    for name, fn in (("pca", cmd_pca), ("correlation", cmd_correlation)):
        p = sp2(name, fn)
        p.add_argument("--bws", nargs="+", required=True)
        p.add_argument("--labels", nargs="+", required=True)
        p.add_argument("--conditions", nargs="+", required=True)
        p.add_argument("--method", default="pearson")

    p = sp2("differential", cmd_differential)
    p.add_argument("--blacklist", required=True)
    p.add_argument("--clusters", default="")
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--K", type=int, default=5000)

    p = sp2("enrichment", cmd_enrichment)
    p.add_argument("--clusters", required=True)
    p.add_argument("--pooled", required=True)
    p.add_argument("--db-dir", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--regions", nargs="+", default=["genome_wide", "genic", "intergenic"])
    p.add_argument("--n-elements", type=int, default=10)
    p.add_argument("--cutoff-for-overlap", type=int, default=0)

    a = ap.parse_args()
    for k in ("treated", "control"):  # single-file commands take one bigWig
        v = getattr(a, k, None)
        if isinstance(v, list) and a.cmd not in ("differential",):
            setattr(a, k, v[0])
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    for noisy in ("kaleido", "choreographer", "choreo", "asyncio"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    a.fn(a)


if __name__ == "__main__":
    main()
