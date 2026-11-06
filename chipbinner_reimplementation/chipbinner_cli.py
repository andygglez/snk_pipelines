#!/usr/bin/env python3
"""ChIPbinner (padilr1/ChIPbinner) two-condition scatterplots, reimplemented in Python + Plotly.

Subcommands mirror the R functions:
  bin-counts      bedtools bamtobed | intersect -c   BAM -> binned counts BED (snakemake_chipbinner bam2bed)
  filter-low-counts ChIPbinner::filter_low_counts
  apply-norm-factors ChIPbinner::apply_normFactors (DESeq2 | edgeR | depth-norm-only)
  merge-norm-bw   ChIPbinner::merge_norm_bw
  plot-pca        ChIPbinner::plot_PCA
  plot-correlation ChIPbinner::plot_correlation
  differential    ChIPbinner::differentialBinAnalysis (ROTS)
  enrich-clust    ChIPbinner::enrich_clust (LOLA)
  norm-bw         ChIPbinner::norm_bw          binned counts (+ input) -> normalised bigWig
  pre-clust       ChIPbinner::pre_clust        two normalised bigWigs -> <prefix>_mat.csv + _pooled.bed
  generate-clust  ChIPbinner::generate_clust   HDBSCAN on the matrix -> one label per bin
  annotate-clust  ChIPbinner::annotate_clust   name clusters A, B, C ... -> BED (chr start end letter)
  genic-scatter   ChIPbinner::genic_intergenic_scatterplot
  density-scatter ChIPbinner::density_based_scatterplot

Plots: x = wild-type, y = treated. Output: <out_prefix>.{html,png,pdf,svg} (static formats need kaleido).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import json  # noqa: E402

from chipbinner import binning, cluster, differential, enrichment, explore, normalize, plots  # noqa: E402
from chipbinner.genomics import annotation_paths, read_chrom_sizes, write_bed, write_bigwig  # noqa: E402


def _bool(s):
    return str(s).lower() in ("1", "true", "t", "yes", "y")


def _num(s):
    return None if s in (None, "", "None", "null") else float(s)


def _range(s):
    if s in (None, ""):
        return None
    if s.strip().lower() == "auto":
        return "auto"
    if s.strip().startswith("{"):
        return json.loads(s)
    lo, hi = s.split(",")
    return [float(lo), float(hi)]


def add_plot_args(p):
    p.add_argument("--treated", required=True, help="normalised bigWig (or bedGraph) of the treated sample (y)")
    p.add_argument("--wildtype", required=True, help="normalised bigWig (or bedGraph) of the wild-type sample (x)")
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--gene-bed", help="genic regions BED (default: annotations/<assembly>_genic.bed.gz)")
    p.add_argument("--igr-bed", help="intergenic regions BED (default: annotations/<assembly>_intergenic.bed.gz)")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--formats", default="html", help="comma list of html,png,pdf,svg")
    p.add_argument("--title", default="")
    p.add_argument("--xaxis-label")
    p.add_argument("--yaxis-label")
    for a in ("min-x", "max-x", "min-y", "max-y"):
        p.add_argument(f"--{a}", type=_num)
    p.add_argument("--range", dest="plot_range", help='auto | "min,max" | JSON {"x": [..], "y": [..]} '
                   '(overrides --min/--max-x/y)')
    p.add_argument("--axis-title-size", type=int, default=22)
    p.add_argument("--pow", type=float, default=1.25)
    p.add_argument("--show-scales", type=_bool, default=True)
    p.add_argument("--write-table", action="store_true", help="also write <prefix>.bins.tsv.gz")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("norm-bw", help="normalise binned counts")
    p.add_argument("--chip", required=True, help="binned BED, count in column 4")
    p.add_argument("--input", help="binned input/control BED")
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--chrom-sizes")
    p.add_argument("--blacklist")
    p.add_argument("--use-input", type=_bool, default=True)
    p.add_argument("--depth-norm", type=_bool, default=True)
    p.add_argument("--pseudocount", type=float, default=1e-3)
    p.add_argument("--raw-count-cutoff", type=float, default=0)
    p.add_argument("--scaling-factor", type=float, default=1.0)
    p.add_argument("--out", required=True, help=".bw, or .bedGraph[.gz]")

    p = sub.add_parser("pre-clust")
    p.add_argument("--treated", required=True)
    p.add_argument("--wildtype", required=True)
    p.add_argument("--out-prefix", required=True)

    p = sub.add_parser("generate-clust")
    p.add_argument("--matrix", required=True)
    p.add_argument("--minpts", type=int, required=True)
    p.add_argument("--minsamps", type=int, required=True)
    p.add_argument("--cores", type=int, default=1)
    p.add_argument("--out", required=True)

    p = sub.add_parser("annotate-clust")
    p.add_argument("--matrix", required=True)
    p.add_argument("--pooled-bed", required=True)
    p.add_argument("--hdbscan-output", required=True)
    p.add_argument("--number-of-clusters", type=int)
    p.add_argument("--out", required=True, help="BED(.gz): chr start end cluster")

    p = sub.add_parser("genic-scatter")
    add_plot_args(p)
    p.add_argument("--hexbins", type=int, default=150)
    p.add_argument("--show-legend", type=_bool, default=True)
    p.add_argument("--legend-pos", default="left", choices=["left", "right"])

    p = sub.add_parser("density-scatter")
    add_plot_args(p)
    p.add_argument("--annotated-clusters", required=True, help="BED from annotate-clust")
    p.add_argument("--number-of-clusters", type=int)
    p.add_argument("--hexbins", type=int, default=100)
    p.add_argument("--include-additional-density-plot", type=_bool, default=False)
    p.add_argument("--filter-extreme-bins", type=_bool, default=True)
    p.add_argument("--cluster-colors", default="bright", help="bright | alphabet | comma list of colours")
    p.add_argument("--max-points", type=int, default=120000)
    p.add_argument("--show-hulls", action="store_true")

    p = sub.add_parser("bin-counts")
    p.add_argument("--bam", required=True)
    p.add_argument("--windows", help="windows BED; default: --chrom-sizes / --bin-size")
    p.add_argument("--chrom-sizes")
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--bin-size", type=int, default=10000)
    p.add_argument("--read-bed", help="also write the bamtobed read BED")
    p.add_argument("--out", required=True)

    p = sub.add_parser("filter-low-counts")
    p.add_argument("--beds", nargs="+", required=True)
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--blacklist")
    p.add_argument("--cutoff", type=float, default=100)
    p.add_argument("--out-dir", required=True, help="writes <label>_filt.bed")

    p = sub.add_parser("apply-norm-factors")
    p.add_argument("--norm-method", required=True, choices=["DESeq2", "edgeR", "depth-norm-only"])
    p.add_argument("--treated", nargs="+", required=True)
    p.add_argument("--wildtype", nargs="+", required=True)
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--blacklist")
    p.add_argument("--out-dir", required=True, help="writes <label>_sizeFactorNorm.bed")

    p = sub.add_parser("merge-norm-bw")
    p.add_argument("--rep1", required=True)
    p.add_argument("--rep2", required=True)
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--out", required=True)

    for nm in ("plot-pca", "plot-correlation"):
        p = sub.add_parser(nm)
        p.add_argument("--bigwigs", nargs="+", required=True)
        p.add_argument("--colors", nargs="+", required=True)
        p.add_argument("--sample-labels", nargs="+")
        p.add_argument("--title", default="")
        p.add_argument("--out-prefix", required=True)
        p.add_argument("--formats", default="html")
        if nm == "plot-correlation":
            p.add_argument("--method", default="pearson", choices=["pearson", "spearman", "kendall"])

    p = sub.add_parser("differential")
    p.add_argument("--treated", nargs="+", required=True)
    p.add_argument("--wildtype", nargs="+", required=True)
    p.add_argument("--treated-label", required=True)
    p.add_argument("--wildtype-label", required=True)
    p.add_argument("--annotated-clusters")
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--blacklist")
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--K", type=int, required=True)
    p.add_argument("--out-dir", required=True)

    p = sub.add_parser("enrich-clust")
    p.add_argument("--annotated-clusters", required=True)
    p.add_argument("--query-cluster", required=True)
    p.add_argument("--pooled-bed", required=True)
    p.add_argument("--db-dir", required=True, help="LOLA regionDB/<genome> folder")
    p.add_argument("--region", default="genome_wide", choices=["genome_wide", "genic", "intergenic"])
    p.add_argument("--genome-assembly", default="hg38", choices=["hg38", "mm10"])
    p.add_argument("--n-elements", type=int, default=10)
    p.add_argument("--cutoff-for-overlap", type=int, default=0)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--formats", default="html")

    a = ap.parse_args(argv)

    def save(fig, prefix, formats):
        for f in formats.split(","):
            if f == "html":
                fig.write_html(f"{prefix}.html", include_plotlyjs=True)
            else:
                fig.write_image(f"{prefix}.{f}", width=fig.layout.width, height=fig.layout.height)
        print("written:", ", ".join(f"{prefix}.{f}" for f in formats.split(",")))

    ann = annotation_paths(getattr(a, "genome_assembly", "hg38"))
    if a.cmd == "bin-counts":
        win = (binning.read_windows(a.windows) if a.windows
               else binning.make_windows(a.chrom_sizes or ann["chrom_sizes"], a.bin_size))
        df = binning.bam_to_bins(a.bam, win, read_bed=a.read_bed)
        write_bed(df, a.out)
        print(f"{len(df)} windows, {df['count'].sum()} read-window hits -> {a.out}")
        return
    if a.cmd == "filter-low-counts":
        os.makedirs(a.out_dir, exist_ok=True)
        for lab, df in normalize.filter_low_counts(a.beds, a.blacklist or ann["blacklist"], a.cutoff).items():
            write_bed(df, os.path.join(a.out_dir, f"{lab}_filt.bed"))
            print(f"{lab}: {len(df)} bins kept")
        return
    if a.cmd == "apply-norm-factors":
        os.makedirs(a.out_dir, exist_ok=True)
        res, table = normalize.apply_norm_factors(a.norm_method, a.treated, a.wildtype, a.blacklist or ann["blacklist"])
        for lab, df in res.items():
            write_bed(df, os.path.join(a.out_dir, f"{lab}_sizeFactorNorm.bed"))
        print(table.to_string(index=False))
        return
    if a.cmd == "merge-norm-bw":
        write_bigwig(normalize.merge_norm_bw(a.rep1, a.rep2), a.out, read_chrom_sizes(ann["chrom_sizes"]))
        print("->", a.out)
        return
    if a.cmd in ("plot-pca", "plot-correlation"):
        if a.cmd == "plot-pca":
            fig, tab = explore.plot_pca(a.bigwigs, a.colors, a.sample_labels, plot_title=a.title)
            save(fig, f"{a.out_prefix}_PCA_plot", a.formats)
        else:
            fig, tab = explore.plot_correlation(a.bigwigs, a.colors, a.sample_labels, a.method, plot_title=a.title)
            save(fig, f"{a.out_prefix}_correlation_plot", a.formats)
        return
    if a.cmd == "differential":
        os.makedirs(a.out_dir, exist_ok=True)
        allb, per, summ, _ = differential.differential_bin_analysis(
            a.treated, a.wildtype, a.blacklist or ann["blacklist"], a.treated_label, a.wildtype_label,
            annotated_clusters=a.annotated_clusters, bootstrap_value=a.bootstrap, K_value=a.K)
        allb.to_csv(os.path.join(a.out_dir, "stats_for_all_bins_no_clustering.csv"), index=False)
        if per is not None:
            per.to_csv(os.path.join(a.out_dir, "stats_per_cluster.csv"), index=False)
        open(os.path.join(a.out_dir, "ROTS_output_summary.txt"), "w").write(summ)
        print(summ)
        return
    if a.cmd == "enrich-clust":
        agg = enrichment.enrich_clust(a.annotated_clusters, a.query_cluster, a.pooled_bed, a.db_dir,
                                      ann["gene"], ann["igr"], region=a.region)
        if agg is None:
            print("No significant enrichment or depletion.")
            return
        agg.to_csv(f"{a.out_prefix}.csv", index=False)
        save(enrichment.plot_enrichment(agg, a.n_elements, a.cutoff_for_overlap,
                                        f"Cluster {a.query_cluster} · {a.region}"), a.out_prefix, a.formats)
        return
    if a.cmd == "norm-bw":
        ann = annotation_paths(a.genome_assembly)
        cs = a.chrom_sizes or ann["chrom_sizes"]
        df = normalize.norm_bw(a.chip, cs, a.blacklist or ann["blacklist"], input_binned_file=a.input,
                               use_input=a.use_input, depth_norm=a.depth_norm, pseudocount=a.pseudocount,
                               raw_count_cutoff=a.raw_count_cutoff, scaling_factor=a.scaling_factor)
        if a.out.endswith((".bw", ".bigwig", ".bigWig")):
            write_bigwig(df, a.out, read_chrom_sizes(cs))
        else:
            write_bed(df, a.out)
        print(f"{len(df)} bins -> {a.out}")
    elif a.cmd == "pre-clust":
        pair = cluster.pre_clust(a.treated, a.wildtype, a.out_prefix)
        print(f"{len(pair)} bins -> {a.out_prefix}_mat.csv, {a.out_prefix}_pooled.bed")
    elif a.cmd == "generate-clust":
        lab = cluster.generate_clust(a.matrix, a.minpts, a.minsamps, a.cores, a.out)
        print(f"{len(set(lab) - {-1})} clusters, {(lab == -1).sum()} noise bins -> {a.out}")
    elif a.cmd == "annotate-clust":
        df = cluster.annotate_clust(a.matrix, a.pooled_bed, a.hdbscan_output, a.number_of_clusters)
        write_bed(df, a.out)
        print(df["name"].value_counts().sort_index().to_dict(), "->", a.out)
    else:
        common = dict(genome_assembly=a.genome_assembly, gene_bed=a.gene_bed, igr_bed=a.igr_bed,
                      out_prefix=a.out_prefix, formats=a.formats.split(","), title_of_plot=a.title,
                      xaxis_label=a.xaxis_label, yaxis_label=a.yaxis_label, min_x=a.min_x, max_x=a.max_x,
                      min_y=a.min_y, max_y=a.max_y, pow=a.pow, hexbins=a.hexbins,
                      show_scales=a.show_scales, write_table=a.write_table, axis_title_size=a.axis_title_size,
                      plot_range=_range(a.plot_range))
        if a.cmd == "genic-scatter":
            plots.genic_intergenic_scatterplot(a.treated, a.wildtype, show_legend=a.show_legend,
                                               legend_pos=a.legend_pos, **common)
        else:
            cc = a.cluster_colors if a.cluster_colors in ("bright", "alphabet") else a.cluster_colors.split(",")
            plots.density_based_scatterplot(
                a.treated, a.wildtype, a.annotated_clusters, number_of_clusters=a.number_of_clusters,
                include_additional_density_plot=a.include_additional_density_plot,
                filter_extreme_bins=a.filter_extreme_bins, cluster_colors=cc, max_points=a.max_points,
                show_hulls=a.show_hulls, **common)
        print("written:", ", ".join(f"{a.out_prefix}.{f}" for f in a.formats.split(",")))


if __name__ == "__main__":
    main()
