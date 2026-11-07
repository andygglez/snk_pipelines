"""One function per workflow step, paths in and paths out.

Each function is one snk-chipseq rule (``chipbinner_cli.py step <name> --json``)
and is also what the standalone :class:`chipbinner.pipeline.Workflow` calls, so
both produce the same files. A step always writes all of its declared outputs:
when there is nothing to compute (too few shared bins, too few HDBSCAN clusters,
no clusters to test) it writes empty placeholder files, as snakemake_chipbinner
does.
"""
from __future__ import annotations

import logging
import os
import shutil
import time

import numpy as np
import pandas as pd
import yaml

from . import binning, cluster, differential as _differential, enrichment as _enrichment, explore, extra_plots
from . import normalize, plots
from .genomics import read_chrom_sizes, write_bed, write_bigwig

log = logging.getLogger("chipbinner")

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf",
           "#7f7f7f", "#bcbd22"]
MIN_SHARED_BINS = 10
GENIC_PLOT = {"hexbins": 150, "show_legend": True, "legend_pos": "right", "show_scales": True}
DENSITY_PLOT = {"hexbins": 150, "show_scales": True, "include_additional_density_plot": True,
                "filter_extreme_bins": True, "cluster_colors": "bright"}


# ----------------------------------------------------------------------------- helpers
def _mkdir_for(*paths):
    for p in paths:
        d = os.path.dirname(p)
        if d:
            os.makedirs(d, exist_ok=True)


def _touch(*paths):
    _mkdir_for(*paths)
    for p in paths:
        open(p, "w").close()


def save_fig(fig, prefix, formats):
    """Write a Plotly figure as <prefix>.<fmt> (html: self-contained; others: kaleido)."""
    _mkdir_for(prefix)
    for fmt in formats:
        path = f"{prefix}.{fmt}"
        if fmt == "html":
            fig.write_html(path, include_plotlyjs=True, config={"displaylogo": False})
        else:
            fig.write_image(path, width=fig.layout.width, height=fig.layout.height,
                            scale=2 if fmt in ("png", "jpeg", "jpg") else 1)


# ----------------------------------------------------------------------------- per replicate
def make_windows(out, chrom_sizes="", bin_size=10000, windows=""):
    """Bin template: the windows BED as given, else bedtools makewindows -w bin_size."""
    win = binning.read_windows(windows) if windows else binning.make_windows(chrom_sizes, int(bin_size))
    _mkdir_for(out)
    write_bed(win, out)


def bin_bam(bam, out, windows="", chrom_sizes="", bin_size=10000, read_bed=""):
    """bam2bed: alignment counts per window (windows BED, else windows from chrom_sizes).

    read_bed: also write the bamtobed read BED there (standalone keep_read_bed)."""
    win = binning.read_windows(windows) if windows else binning.make_windows(chrom_sizes, int(bin_size))
    log.info("bam2bed %s", bam)
    _mkdir_for(out)
    write_bed(binning.bam_to_bins(bam, win, read_bed=read_bed or None), out)


def norm_replicate(binned, out, chrom_sizes, blacklist="", input_binned="", use_input=False, depth_norm=True,
                   pseudocount=1e-15, raw_count_cutoff=100, scaling_factor=1.0):
    """ChIPbinner::norm_bw on one replicate -> bigWig."""
    log.info("norm_bw %s (scaling factor %s)", binned, scaling_factor)
    df = normalize.norm_bw(binned, chrom_sizes, blacklist, input_binned_file=input_binned or None,
                           use_input=bool(use_input) and bool(input_binned), depth_norm=bool(depth_norm),
                           pseudocount=float(pseudocount), raw_count_cutoff=float(raw_count_cutoff),
                           scaling_factor=float(scaling_factor))
    _mkdir_for(out)
    write_bigwig(df, out, read_chrom_sizes(chrom_sizes))


def merge_replicates(bigwigs, out, chrom_sizes):
    """ChIPbinner::merge_norm_bw over any number of replicates (one replicate: copied)."""
    _mkdir_for(out)
    if len(bigwigs) == 1:
        shutil.copyfile(bigwigs[0], out)
        return
    log.info("merge_norm_bw %s", out)
    write_bigwig(normalize.merge_norm_bw(*bigwigs), out, read_chrom_sizes(chrom_sizes))


# ----------------------------------------------------------------------------- per mark
def exploratory(bigwigs, labels, conditions, mark, outdir, formats=("pdf", "html"), pca=True,
                correlation=True, method="pearson"):
    """PCA and correlation plots of every replicate bigWig of one mark (coloured by condition)."""
    os.makedirs(outdir, exist_ok=True)
    cmap = {c: PALETTE[k % len(PALETTE)] for k, c in enumerate(dict.fromkeys(conditions))}
    colors = [cmap[c] for c in conditions]
    if pca:
        fig, _ = explore.plot_pca(list(bigwigs), colors, list(labels), plot_title=f"{mark} · PCA")
        save_fig(fig, os.path.join(outdir, f"{mark}_PCA_plot"), formats)
    if correlation:
        fig, _ = explore.plot_correlation(list(bigwigs), colors, list(labels), method,
                                          plot_title=f"{mark} · {method} correlation")
        save_fig(fig, os.path.join(outdir, f"{mark}_correlation_plot"), formats)


# ----------------------------------------------------------------------------- per comparison
def legend(out):
    """legend.png: key of the genic / intergenic colour scale."""
    save_fig(extra_plots.genic_legend(), out[: -len(".png")], ["png"])


def plot_spec(pair, range="auto", square_auto_range=True):
    """Plot range -> ({"x": [lo, hi], "y": [lo, hi]}, source); auto parts are found from the data."""
    if isinstance(range, str) and range.lower() == "auto":
        return plots.auto_range(pair, square=bool(square_auto_range)), "auto"
    rng = list(plots.resolve_range(range))
    if any(v is None for v in rng):  # partially automatic ({x: auto, y: [..]})
        auto = plots.auto_range(pair, square=False)
        rng = [a if v is None else v for v, a in zip(rng, auto)]
    return rng, "config"


def write_plot_range(path, name, mark, rng, source, pow):
    _mkdir_for(path)
    with open(path, "w") as fh:
        yaml.safe_dump({"comparison": name, "mark": mark, "source": source, "pow": float(pow),
                        "x": [round(rng[0], 4), round(rng[1], 4)], "y": [round(rng[2], 4), round(rng[3], 4)]},
                       fh, sort_keys=False)


def scatter(treated_bw, baseline_bw, gene_bed, igr_bed, prefix, mark, treatment, baseline, range="auto",
            pow=0.75, square_auto_range=True, genic_plot=None, formats=("pdf", "html")):
    """genic_intergenic_scatterplot + genic vs intergenic panels; records the range used.

    Writes <prefix>.plot_range.yaml, <prefix>.scatterplot.<fmt> and
    <prefix>.genic.vs.intergenic.png. Fewer than MIN_SHARED_BINS shared bins:
    empty plot files and a range file with ``source: skipped``.
    """
    name = f"{treatment}_{baseline}"
    plot_files = [f"{prefix}.scatterplot.{f}" for f in formats] + [f"{prefix}.genic.vs.intergenic.png"]
    pair = cluster.load_pair(treated_bw, baseline_bw)
    if len(pair) < MIN_SHARED_BINS:
        log.warning("%s %s: only %d shared bins, skipping", name, mark, len(pair))
        _touch(*plot_files)
        _mkdir_for(prefix)
        with open(f"{prefix}.plot_range.yaml", "w") as fh:
            yaml.safe_dump({"comparison": name, "mark": mark, "source": "skipped", "shared_bins": len(pair)},
                           fh, sort_keys=False)
        return
    rng, src = plot_spec(pair, range, square_auto_range)
    write_plot_range(f"{prefix}.plot_range.yaml", name, mark, rng, src, pow)
    spec = {"x": rng[:2], "y": rng[2:]}
    log.info("%s %s: range x=%s y=%s (%s), pow=%s", name, mark, spec["x"], spec["y"], src, pow)
    scatter_plot(treated_bw, baseline_bw, gene_bed, igr_bed, prefix, mark, treatment, baseline, spec, pow,
                 genic_plot, formats)
    genic_vs_intergenic(pair, gene_bed, igr_bed, prefix, treatment, baseline)


def scatter_plot(treated_bw, baseline_bw, gene_bed, igr_bed, prefix, mark, treatment, baseline, spec, pow,
                 genic_plot=None, formats=("pdf", "html")):
    g = dict(GENIC_PLOT, **(genic_plot or {}))
    fig = plots.genic_intergenic_scatterplot(
        treated_bw, baseline_bw, gene_bed=gene_bed, igr_bed=igr_bed, title_of_plot=f"Change in {mark}",
        xaxis_label=baseline, yaxis_label=treatment, plot_range=spec, pow=float(pow), hexbins=int(g["hexbins"]),
        show_scales=bool(g["show_scales"]), show_legend=bool(g["show_legend"]), legend_pos=g["legend_pos"])
    save_fig(fig, f"{prefix}.scatterplot", formats)


def genic_vs_intergenic(pair, gene_bed, igr_bed, prefix, treatment, baseline):
    save_fig(extra_plots.genic_vs_intergenic(pair, gene_bed, igr_bed, baseline, treatment),
             f"{prefix}.genic.vs.intergenic", ["png"])


def preclust(treated_bw, baseline_bw, prefix, treatment, baseline):
    """ChIPbinner::pre_clust -> <prefix>_preclust_{mat.csv,pooled.bed} + smoothScatter.jpeg.

    Fewer than MIN_SHARED_BINS shared bins: empty files.
    """
    pre = f"{prefix}_preclust"
    _mkdir_for(pre)
    pair = cluster.load_pair(treated_bw, baseline_bw)
    if len(pair) < MIN_SHARED_BINS:
        _touch(f"{pre}_mat.csv", f"{pre}_pooled.bed", f"{pre}_smoothScatter.jpeg")
        return
    cluster.pre_clust(treated_bw, baseline_bw, pre)
    save_fig(extra_plots.smooth_scatter(pair["x"], pair["y"], baseline, treatment), f"{pre}_smoothScatter",
             ["jpeg"])


def _load_matrix(path):
    if os.path.getsize(path) == 0:
        return np.zeros((0, 2))
    return np.loadtxt(path, delimiter=",", ndmin=2)


def run_hdbscan(mat, minpts, minsamps, cluster_selection_method="eom", cores=1):
    """HDBSCAN labels and membership probabilities (hdbscan package, else scikit-learn)."""
    if len(mat) < 2:
        return np.full(len(mat), -1), None, "skipped: fewer than 2 bins"
    try:
        try:
            import hdbscan as _hdb

            m = _hdb.HDBSCAN(min_cluster_size=minpts, min_samples=minsamps, core_dist_n_jobs=int(cores),
                             cluster_selection_method=cluster_selection_method).fit(mat)
            msg = "hdbscan"
        except ImportError:
            from sklearn.cluster import HDBSCAN

            m = HDBSCAN(min_cluster_size=minpts, min_samples=minsamps,
                        cluster_selection_method=cluster_selection_method).fit(mat)
            msg = "sklearn.cluster.HDBSCAN"
        return m.labels_.astype(int), m.probabilities_, msg
    except Exception as e:  # e.g. min_samples > number of bins
        return np.full(len(mat), -1), None, f"failed: {e}"


def hdbscan(matrix, minpts, minsamps, labels_out, png_out, log_out, treatment, baseline,
            cluster_selection_method="eom", cores=1):
    """generate_clust at one (minpts, minsamps) grid point: labels, log and a cluster scatter PNG."""
    t0 = time.time()
    mat = _load_matrix(matrix)
    minpts, minsamps = int(minpts), int(minsamps)
    labels, probs, msg = run_hdbscan(mat, minpts, minsamps, cluster_selection_method, cores)
    _mkdir_for(labels_out, png_out, log_out)
    with open(log_out, "w") as fh:
        fh.write(f"INFO:root:Starting dbscan ({msg}); minpts={minpts} minsamps={minsamps} n={len(mat)}\n"
                 f"INFO:root:clusters={len(set(labels) - {-1})} noise={(labels == -1).sum()}\n"
                 f"INFO:root:Hdbscan done in {time.time() - t0:.1f}s\n")
    log.info("  HDBSCAN %s/%s: %d clusters", minpts, minsamps, len(set(labels) - {-1}))
    np.savetxt(labels_out, labels, fmt="%i")
    if len(mat):
        save_fig(extra_plots.cluster_scatter(mat, labels, probs, baseline, treatment), png_out[: -len(".png")],
                 ["png"])
    else:
        _touch(png_out)


def read_plot_range(path):
    """plot_range.yaml -> (spec {"x": [..], "y": [..]}, pow), or (None, None) when the scatter was skipped."""
    with open(path) as fh:
        r = yaml.safe_load(fh)
    if r.get("source") == "skipped":
        return None, None
    return {"x": list(r["x"]), "y": list(r["y"])}, float(r["pow"])


def density(matrix, pooled, labels, plot_range, n, treated_bw, baseline_bw, gene_bed, igr_bed, annotated_out,
            plot_prefix, mark, treatment, baseline, density_plot=None, formats=("pdf", "html")):
    """annotate_clust (n clusters) + density_based_scatterplot, with the range of the scatter step.

    Fewer than n clusters found (or the scatter was skipped): empty BED and plot files.
    """
    spec, pow = read_plot_range(plot_range)
    density_with_spec(matrix, pooled, labels, spec, pow, n, treated_bw, baseline_bw, gene_bed, igr_bed,
                      annotated_out, plot_prefix, mark, treatment, baseline, density_plot, formats)


def density_with_spec(matrix, pooled, labels, spec, pow, n, treated_bw, baseline_bw, gene_bed, igr_bed,
                      annotated_out, plot_prefix, mark, treatment, baseline, density_plot=None,
                      formats=("pdf", "html")):
    n = int(n)
    lab = np.loadtxt(labels, dtype=int, ndmin=1) if os.path.getsize(labels) else np.zeros(0, int)
    n_found = len(set(lab) - {-1})
    plot_files = [f"{plot_prefix}.{f}" for f in formats]
    if spec is None or n_found < n or len(set(lab)) <= 1:
        log.info("  %s: %d clusters < %s, empty outputs", os.path.basename(annotated_out), n_found, n)
        _touch(annotated_out, *plot_files)
        return
    _mkdir_for(annotated_out)
    write_bed(cluster.annotate_clust(matrix, pooled, labels, n), annotated_out)
    dp = dict(DENSITY_PLOT, **(density_plot or {}))
    fig = plots.density_based_scatterplot(
        treated_bw, baseline_bw, annotated_out, gene_bed=gene_bed, igr_bed=igr_bed, number_of_clusters=n,
        title_of_plot=mark, xaxis_label=baseline, yaxis_label=treatment, plot_range=spec, pow=float(pow),
        hexbins=int(dp["hexbins"]), show_scales=bool(dp["show_scales"]),
        include_additional_density_plot=bool(dp["include_additional_density_plot"]),
        filter_extreme_bins=bool(dp["filter_extreme_bins"]), cluster_colors=dp["cluster_colors"])
    save_fig(fig, plot_prefix, formats)


def _nonempty(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def differential(treated_bws, baseline_bws, blacklist, treatment, baseline, outdir, clusters="", bootstrap=1000,
                 K=5000):
    """differentialBinAnalysis (ROTS). Writes stats_for_all_bins_no_clustering.csv,
    stats_per_cluster.csv (empty without clusters) and ROTS_output_summary.txt."""
    os.makedirs(outdir, exist_ok=True)
    clus = clusters if _nonempty(clusters) else None
    log.info("ROTS %s vs %s (B=%s, K=%s)", treatment, baseline, bootstrap, K)
    allb, per, summ, _ = _differential.differential_bin_analysis(
        list(treated_bws), list(baseline_bws), blacklist, treatment, baseline, annotated_clusters=clus,
        bootstrap_value=int(bootstrap), K_value=int(K))
    allb.to_csv(os.path.join(outdir, "stats_for_all_bins_no_clustering.csv"), index=False)
    per_f = os.path.join(outdir, "stats_per_cluster.csv")
    if per is not None:
        per.to_csv(per_f, index=False)
    else:
        _touch(per_f)
    with open(os.path.join(outdir, "ROTS_output_summary.txt"), "w") as fh:
        fh.write(summ)


def enrichment(clusters, pooled, db_dir, gene_bed, igr_bed, outdir, stem, regions=("genome_wide", "genic",
               "intergenic"), n_elements=10, cutoff_for_overlap=0, formats=("pdf", "html"), skip_done=False):
    """enrich_clust (LOLA) for every cluster letter x region -> <outdir>/<stem>.cluster<L>.<region>.{csv,<fmt>}.

    No clusters: an empty outdir. No significant result: an empty .csv and no plot.
    """
    os.makedirs(outdir, exist_ok=True)
    if not _nonempty(clusters):
        log.warning("enrichment %s: no clusters, skipped", stem)
        return
    letters = sorted(pd.read_csv(clusters, sep="\t", header=None)[3].astype(str).unique())
    for letter in letters:
        for region in regions:
            base = os.path.join(outdir, f"{stem}.cluster{letter}.{region}")
            if skip_done and os.path.exists(f"{base}.csv"):
                continue
            agg = _enrichment.enrich_clust(clusters, letter, pooled, db_dir, gene_bed, igr_bed, region=region)
            if agg is None:
                log.info("  cluster %s %s: no significant enrichment/depletion", letter, region)
                open(f"{base}.csv", "w").close()
                continue
            agg.to_csv(f"{base}.csv", index=False)
            save_fig(_enrichment.plot_enrichment(agg, int(n_elements), int(cutoff_for_overlap),
                                                 f"Cluster {letter} · {region}"), base, formats)


# name used by `chipbinner_cli.py step <name>` -> function
STEPS = {
    "windows": make_windows,
    "bin": bin_bam,
    "norm": norm_replicate,
    "merge": merge_replicates,
    "explore": exploratory,
    "legend": legend,
    "scatter": scatter,
    "preclust": preclust,
    "hdbscan": hdbscan,
    "density": density,
    "differential": differential,
    "enrichment": enrichment,
}

