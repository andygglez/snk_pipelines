"""Config-driven ChIPbinner workflow producing the snakemake_chipbinner output tree.

Layout (``<outdir>`` = chipbinner.outdir; names as in andygglez/snakemake_chipbinner):

    <outdir>/bam2bed/<sample>.binned.bed                 bin counts (BAM inputs)
    <outdir>/bam2bed/<sample>.bed                         read BED (only with binning.keep_read_bed)
    <outdir>/bam2bed/<sample>.binned_filt.bed            filter_low_counts (when enabled)
    <outdir>/{MS,notMS}bigwigs/<sample>.binned[_filt].bw  norm_bw per replicate
    <outdir>/{MS,notMS}bigwigs_final/<cond>.<mark>.bw    merge_norm_bw (or the single replicate)
    <outdir>/comparisons/<T>_<C>/
        <T>_<C>_<mark>.scatterplot.{pdf,html}             genic_intergenic_scatterplot
        legend.png
        <T>_<C>_<mark>.plot_range.yaml                    range used (auto-detected or from config)
        <T>_<C>_<mark>_preclust_mat.csv / _preclust_pooled.bed / _preclust_smoothScatter.jpeg
        <T>_<C>_<mark>.genic.vs.intergenic.png
        clustering/clus.<T>_<C>_<mark>_clustered.<i>.<j>.{txt,png}, log.<i>.<j>.txt
        clustering/<T>_<C>_<mark>.<n>clusters.<i>.<j>.annotated_clusters.bed      (R: .rda)
        clustering/<T>_<C>_<mark>_<n>clusters.<i>.<j>.plots_density_scatterplot.{pdf,html}
        differential/  (ROTS; needs 2 replicates per condition)   [port extension of the tree]
        enrichment/    (LOLA-style region DB)                      [port extension of the tree]
    <outdir>/exploratory/<mark>_{PCA,correlation}_plot.{pdf,html}  [port extension of the tree]

Empty placeholder files are written where snakemake_chipbinner writes them
(annotation / density plot when HDBSCAN finds too few clusters).
"""
from __future__ import annotations

import logging
import os
import shutil
import time

import numpy as np
import pandas as pd
import yaml

from . import binning, cluster, differential, enrichment, explore, extra_plots, normalize, plots
from .genomics import annotation_paths, read_chrom_sizes, write_bed, write_bigwig

log = logging.getLogger("chipbinner")

DEFAULTS = {
    "outdir": "outdir",
    "genome_assembly": "hg38",
    "annotations": {"gene_bed": "", "igr_bed": "", "blacklist": "", "chrom_sizes": ""},
    "samplesheet": "",
    "sample_folder": "",
    "binning": {"bin_size": 10000, "windows": "", "keep_read_bed": False},
    "filter_low_counts": False,
    "MS_normalize": False,
    "MS_coefficients": {},
    "merge_replicates": True,
    "norm": {"use_input": False, "depth_norm": True, "raw_count_cutoff": 100, "pseudocount": 1e-15},
    "comparisons": [],
    "plot_defaults": {"range": "auto", "pow": 0.75, "square_auto_range": True},
    "genic_plot": {"hexbins": 150, "show_legend": True, "legend_pos": "right", "show_scales": True},
    "density_plot": {"hexbins": 150, "show_scales": True, "include_additional_density_plot": True,
                     "filter_extreme_bins": True, "cluster_colors": "bright"},
    "clustering": {"minpts": [100, 400, 700, 1000], "minsamps": [100, 400, 700, 1000],
                   "number_of_clusters": [2, 3], "cluster_selection_method": "eom", "cores": 1},
    "formats": ["pdf", "html"],
    "pca": {"enabled": True},
    "correlation": {"enabled": True, "method": "pearson"},
    "differential": {"enabled": False, "bootstrap": 1000, "K": 5000,
                     "clusters": {"minpts": 100, "minsamps": 100, "number_of_clusters": 3}},
    "enrichment": {"enabled": False, "db_dir": "", "regions": ["genome_wide", "genic", "intergenic"],
                   "n_elements": 10, "cutoff_for_overlap": 0,
                   "clusters": {"minpts": 100, "minsamps": 100, "number_of_clusters": 3}},
}
PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf",
           "#7f7f7f", "#bcbd22"]


def _merge(base, over):
    out = dict(base)
    for k, v in (over or {}).items():
        out[k] = _merge(base[k], v) if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return out


def load_config(path: str) -> dict:
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    cfg = _merge(DEFAULTS, raw.get("chipbinner", raw))
    cfg["_config_dir"] = os.path.dirname(os.path.abspath(path))
    return cfg


def _stem(name: str) -> str:
    base = os.path.basename(str(name))
    for ext in (".bam", ".binned.bed.gz", ".binned.bed", ".bed.gz", ".bed", ".bedGraph.gz", ".bedGraph",
                ".bw", ".bigwig", ".bigWig"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def _kind(path: str) -> str:
    p = str(path)
    if p.endswith(".bam"):
        return "bam"
    if p.endswith((".bw", ".bigwig", ".bigWig")):
        return "bigwig"
    return "bed"


def parse_comparison(c) -> dict:
    if isinstance(c, str):  # snakemake_chipbinner style "<treatment>_<control>"
        t, ctl = c.split("_", 1)
        c = {"treatment": t, "control": ctl}
    return dict(c)


class Workflow:
    def __init__(self, cfg: dict, force: bool = False):
        self.cfg = cfg
        self.force = force
        self.out = cfg["outdir"]
        ann = annotation_paths(cfg["genome_assembly"])
        a = cfg["annotations"]
        self.gene_bed = a.get("gene_bed") or ann["gene"]
        self.igr_bed = a.get("igr_bed") or ann["igr"]
        self.blacklist = a.get("blacklist") or ann["blacklist"]
        self.chrom_sizes = a.get("chrom_sizes") or ann["chrom_sizes"]
        self.bw_folder = "MSbigwigs" if cfg["MS_normalize"] else "notMSbigwigs"
        self.formats = list(cfg["formats"])
        self.sheet = self._read_sheet()
        ms = cfg["MS_coefficients"]
        if isinstance(ms, str) and ms:
            with open(self._path(ms)) as fh:
                ms = yaml.safe_load(fh)
        self.ms = ms or {}

    # ------------------------------------------------------------------ helpers
    def _path(self, p):
        if not p or os.path.isabs(p):
            return p
        for base in (self.cfg.get("sample_folder") or "", os.getcwd(), self.cfg["_config_dir"]):
            q = os.path.join(base, p) if base else p
            if os.path.exists(q):
                return q
        return os.path.join(self.cfg.get("sample_folder") or "", p)

    def _read_sheet(self):
        sh = pd.read_csv(self._path(self.cfg["samplesheet"]), dtype=str).fillna("")
        need = {"Rep1", "Condition", "Mark"}
        if not need <= set(sh.columns):
            raise ValueError(f"samplesheet needs columns {sorted(need)} (+ optional Rep2, Input1, Input2)")
        for c in ("Rep2", "Input1", "Input2"):
            if c not in sh:
                sh[c] = ""
        return sh

    def _done(self, *paths):
        return (not self.force) and all(os.path.exists(p) for p in paths)

    def _mk(self, *parts):
        d = os.path.join(self.out, *parts)
        os.makedirs(d, exist_ok=True)
        return d

    def _save(self, fig, prefix, formats=None):
        for fmt in formats or self.formats:
            path = f"{prefix}.{fmt}"
            if fmt == "html":
                fig.write_html(path, include_plotlyjs=True, config={"displaylogo": False})
            else:
                fig.write_image(path, width=fig.layout.width, height=fig.layout.height,
                                scale=2 if fmt in ("png", "jpeg", "jpg") else 1)

    # ------------------------------------------------------------------ per replicate
    def binned(self, rep: str) -> str:
        """Raw binned counts BED for a BAM / binned-BED replicate (bam2bed)."""
        src = self._path(rep)
        key = _stem(rep)
        d = self._mk("bam2bed")
        out = os.path.join(d, f"{key}.binned.bed")
        if self._done(out):
            return out
        if _kind(src) == "bam":
            b = self.cfg["binning"]
            if b.get("windows"):
                win = binning.read_windows(self._path(b["windows"]))
            else:
                win = binning.make_windows(self.chrom_sizes, int(b["bin_size"]))
            log.info("bam2bed %s", key)
            counts = binning.bam_to_bins(src, win, read_bed=os.path.join(d, f"{key}.bed") if b.get("keep_read_bed") else None)
            write_bed(counts, out)
        else:
            from .genomics import read_bins

            df = read_bins(src)
            df["score"] = df["score"].round().astype(np.int64)
            write_bed(df, out)
        return out

    def replicate_bw(self, rep: str, inp: str = "") -> str:
        src = self._path(rep)
        key = _stem(rep)
        d = self._mk(self.bw_folder)
        if _kind(src) == "bigwig":  # already normalised: used as is
            out = os.path.join(d, f"{key}.binned.bw")
            if not self._done(out):
                shutil.copyfile(src, out)
            return out
        bed = self.binned(rep)
        ibed = self.binned(inp) if inp else None
        cutoff = self.cfg["filter_low_counts"]
        if cutoff not in (False, None, "", 0, "no"):
            fout = bed.replace(".binned.bed", ".binned_filt.bed")
            if not self._done(fout):
                log.info("filter_low_counts %s (cutoff %s)", key, cutoff)
                res = normalize.filter_low_counts([bed], self.blacklist, float(cutoff))
                write_bed(list(res.values())[0], fout)
            bed = fout
        label = os.path.basename(bed)[: -len(".bed")]
        out = os.path.join(d, f"{label}.bw")
        if self._done(out):
            return out
        sf = 1.0
        if self.cfg["MS_normalize"]:
            if key not in self.ms:
                raise KeyError(f"no MS coefficient for {key} in MS_coefficients")
            sf = float(self.ms[key])
        n = self.cfg["norm"]
        log.info("norm_bw %s (scaling factor %s)", key, sf)
        df = normalize.norm_bw(bed, self.chrom_sizes, self.blacklist, input_binned_file=ibed,
                               use_input=bool(n["use_input"]) and bool(ibed), depth_norm=bool(n["depth_norm"]),
                               pseudocount=float(n["pseudocount"]), raw_count_cutoff=float(n["raw_count_cutoff"]),
                               scaling_factor=sf)
        write_bigwig(df, out, read_chrom_sizes(self.chrom_sizes))
        return out

    def replicate_bws(self, row) -> list[str]:
        reps = [(row["Rep1"], row["Input1"])]
        if row["Rep2"] and self.cfg["merge_replicates"]:
            reps.append((row["Rep2"], row["Input2"] or row["Input1"]))
        return [self.replicate_bw(r, i) for r, i in reps]

    def final_bw(self, cond: str, mark: str) -> str:
        rows = self.sheet[(self.sheet["Condition"] == cond) & (self.sheet["Mark"] == mark)]
        if rows.empty:
            raise KeyError(f"no samplesheet row for condition {cond}, mark {mark}")
        row = rows.iloc[0]
        d = self._mk(f"{self.bw_folder}_final")
        out = os.path.join(d, f"{cond}.{mark}.bw")
        reps = self.replicate_bws(row)
        if self._done(out):
            return out
        if len(reps) == 2:
            log.info("merge_norm_bw %s.%s", cond, mark)
            write_bigwig(normalize.merge_norm_bw(*reps), out, read_chrom_sizes(self.chrom_sizes))
        else:
            shutil.copyfile(reps[0], out)
        return out

    # ------------------------------------------------------------------ per comparison
    def _plot_params(self, comp: dict, mark: str) -> dict:
        p = dict(self.cfg["plot_defaults"])
        p.update({k: comp[k] for k in ("range", "pow") if k in comp})
        p.update((comp.get("overrides") or {}).get(mark, {}))
        return p

    def compare(self, comp: dict, mark: str):
        t, c = comp["treatment"], comp["control"]
        name = f"{t}_{c}"
        cdir = self._mk("comparisons", name)
        pre = os.path.join(cdir, f"{name}_{mark}")
        bw_t, bw_c = self.final_bw(t, mark), self.final_bw(c, mark)
        pp = self._plot_params(comp, mark)
        pair = cluster.load_pair(bw_t, bw_c)
        if len(pair) < 10:
            log.warning("%s %s: only %d shared bins, skipping", name, mark, len(pair))
            return
        # ---- plot range: auto-detected or from config; recorded next to the plots
        if isinstance(pp["range"], str) and pp["range"].lower() == "auto":
            rng = plots.auto_range(pair, square=bool(pp.get("square_auto_range", True)))
            src = "auto"
        else:
            r = plots.resolve_range(pp["range"])
            rng, src = list(r), "config"
            auto = None
            if any(v is None for v in rng):  # partially automatic ({x: auto, y: [..]})
                auto = plots.auto_range(pair, square=False)
                rng = [a if v is None else v for v, a in zip(rng, auto)]
        with open(f"{pre}.plot_range.yaml", "w") as fh:
            yaml.safe_dump({"comparison": name, "mark": mark, "source": src, "pow": float(pp["pow"]),
                            "x": [round(rng[0], 4), round(rng[1], 4)], "y": [round(rng[2], 4), round(rng[3], 4)]},
                           fh, sort_keys=False)
        spec = {"x": rng[:2], "y": rng[2:]}
        log.info("%s %s: range x=%s y=%s (%s), pow=%s", name, mark, spec["x"], spec["y"], src, pp["pow"])

        # ---- genic/intergenic scatterplot (+ legend.png)
        g = self.cfg["genic_plot"]
        if not self._done(f"{pre}.scatterplot.{self.formats[0]}"):
            fig = plots.genic_intergenic_scatterplot(
                bw_t, bw_c, gene_bed=self.gene_bed, igr_bed=self.igr_bed, title_of_plot=f"Change in {mark}",
                xaxis_label=c, yaxis_label=t, plot_range=spec, pow=float(pp["pow"]), hexbins=int(g["hexbins"]),
                show_scales=bool(g["show_scales"]), show_legend=bool(g["show_legend"]),
                legend_pos=g["legend_pos"])
            self._save(fig, f"{pre}.scatterplot")
            self._save(extra_plots.genic_legend(), os.path.join(cdir, "legend"), ["png"])

        # ---- pre_clust (+ smoothScatter) and the genic vs intergenic panels
        mat_f, bed_f = f"{pre}_preclust_mat.csv", f"{pre}_preclust_pooled.bed"
        if not self._done(mat_f, bed_f):
            cluster.pre_clust(bw_t, bw_c, f"{pre}_preclust")
            self._save(extra_plots.smooth_scatter(pair["x"], pair["y"], c, t), f"{pre}_preclust_smoothScatter", ["jpeg"])
        if not self._done(f"{pre}.genic.vs.intergenic.png"):
            self._save(extra_plots.genic_vs_intergenic(pair, self.gene_bed, self.igr_bed, c, t),
                       f"{pre}.genic.vs.intergenic", ["png"])

        # ---- HDBSCAN grid, annotation, density-based scatterplots
        cl = self.cfg["clustering"]
        kdir = self._mk("comparisons", name, "clustering")
        dp = self.cfg["density_plot"]
        mat = np.loadtxt(mat_f, delimiter=",", ndmin=2)
        for i in cl["minpts"]:
            for j in cl["minsamps"]:
                stem = f"{name}_{mark}_clustered"
                lab_f = os.path.join(kdir, f"clus.{stem}.{i}.{j}.txt")
                if not self._done(lab_f):
                    labels, probs = self._hdbscan(mat, int(i), int(j), cl, os.path.join(kdir, f"log.{i}.{j}.txt"))
                    np.savetxt(lab_f, labels, fmt="%i")
                    self._save(extra_plots.cluster_scatter(mat, labels, probs, c, t),
                               os.path.join(kdir, f"clus.{stem}.{i}.{j}"), ["png"])
                labels = np.loadtxt(lab_f, dtype=int, ndmin=1)
                n_found = len(set(labels) - {-1})
                for n in cl["number_of_clusters"]:
                    ann = os.path.join(kdir, f"{name}_{mark}.{n}clusters.{i}.{j}.annotated_clusters.bed")
                    dens = os.path.join(kdir, f"{name}_{mark}_{n}clusters.{i}.{j}.plots_density_scatterplot")
                    if self._done(ann, f"{dens}.{self.formats[0]}"):
                        continue
                    if n_found < int(n) or len(set(labels)) <= 1:
                        log.info("  minpts %s minsamps %s: %d clusters < %s, empty outputs", i, j, n_found, n)
                        open(ann, "w").close()
                        for fmt in self.formats:
                            open(f"{dens}.{fmt}", "w").close()
                        continue
                    if n_found > int(n):
                        # R's annotate_clust keeps the n clusters with the lowest mean signal and silently
                        # drops the rest; say so, because the kept clusters can then be tiny.
                        log.warning("  minpts %s minsamps %s: HDBSCAN found %d clusters, keeping the %s with the "
                                    "lowest mean signal (as ChIPbinner); raise minpts/minsamps for fewer clusters",
                                    i, j, n_found, n)
                    write_bed(cluster.annotate_clust(mat_f, bed_f, lab_f, int(n)), ann)
                    fig = plots.density_based_scatterplot(
                        bw_t, bw_c, ann, gene_bed=self.gene_bed, igr_bed=self.igr_bed, number_of_clusters=int(n),
                        title_of_plot=mark, xaxis_label=c, yaxis_label=t, plot_range=spec, pow=float(pp["pow"]),
                        hexbins=int(dp["hexbins"]), show_scales=bool(dp["show_scales"]),
                        include_additional_density_plot=bool(dp["include_additional_density_plot"]),
                        filter_extreme_bins=bool(dp["filter_extreme_bins"]), cluster_colors=dp["cluster_colors"])
                    self._save(fig, dens)

        if self.cfg["differential"]["enabled"]:
            self.differential(comp, mark, name, kdir)
        if self.cfg["enrichment"]["enabled"]:
            self.enrichment(comp, mark, name, kdir, bed_f)

    def _hdbscan(self, mat, minpts, minsamps, cl, logfile):
        t0 = time.time()
        try:
            import hdbscan

            m = hdbscan.HDBSCAN(min_cluster_size=minpts, min_samples=minsamps,
                                core_dist_n_jobs=int(cl.get("cores", 1)),
                                cluster_selection_method=cl["cluster_selection_method"]).fit(mat)
            labels, probs = m.labels_.astype(int), m.probabilities_
            msg = "hdbscan"
        except ImportError:
            from sklearn.cluster import HDBSCAN

            m = HDBSCAN(min_cluster_size=minpts, min_samples=minsamps,
                        cluster_selection_method=cl["cluster_selection_method"]).fit(mat)
            labels, probs = m.labels_.astype(int), m.probabilities_
            msg = "sklearn.cluster.HDBSCAN"
        except Exception as e:  # e.g. min_samples > number of bins
            labels, probs, msg = np.full(len(mat), -1), None, f"failed: {e}"
        with open(logfile, "w") as fh:
            fh.write(f"INFO:root:Starting dbscan ({msg}); minpts={minpts} minsamps={minsamps} n={len(mat)}\n"
                     f"INFO:root:clusters={len(set(labels) - {-1})} noise={(labels == -1).sum()}\n"
                     f"INFO:root:Hdbscan done in {time.time() - t0:.1f}s\n")
        log.info("  HDBSCAN %s/%s: %d clusters", minpts, minsamps, len(set(labels) - {-1}))
        return labels, probs

    def _selected_clusters(self, name, mark, kdir, sel):
        f = os.path.join(kdir, f"{name}_{mark}.{sel['number_of_clusters']}clusters.{sel['minpts']}."
                               f"{sel['minsamps']}.annotated_clusters.bed")
        return f if os.path.exists(f) and os.path.getsize(f) > 0 else None

    def _rep_bws(self, cond, mark):
        row = self.sheet[(self.sheet["Condition"] == cond) & (self.sheet["Mark"] == mark)].iloc[0]
        return [self.replicate_bw(r, i) for r, i in ((row["Rep1"], row["Input1"]), (row["Rep2"], row["Input2"])) if r]

    def differential(self, comp, mark, name, kdir):
        dcfg = self.cfg["differential"]
        t, c = comp["treatment"], comp["control"]
        tb, cb = self._rep_bws(t, mark), self._rep_bws(c, mark)
        if len(tb) < 2 or len(cb) < 2:
            log.warning("differential %s %s: needs 2 replicates per condition, skipped", name, mark)
            return
        ddir = self._mk("comparisons", name, "differential")
        if self._done(os.path.join(ddir, "stats_for_all_bins_no_clustering.csv")):
            return
        clus = self._selected_clusters(name, mark, kdir, dcfg["clusters"])
        log.info("ROTS %s %s (B=%s, K=%s)", name, mark, dcfg["bootstrap"], dcfg["K"])
        allb, per, summ, _ = differential.differential_bin_analysis(
            tb, cb, self.blacklist, t, c, annotated_clusters=clus, bootstrap_value=int(dcfg["bootstrap"]),
            K_value=int(dcfg["K"]))
        allb.to_csv(os.path.join(ddir, "stats_for_all_bins_no_clustering.csv"), index=False)
        if per is not None:
            per.to_csv(os.path.join(ddir, "stats_per_cluster.csv"), index=False)
        with open(os.path.join(ddir, "ROTS_output_summary.txt"), "w") as fh:
            fh.write(summ)

    def enrichment(self, comp, mark, name, kdir, pooled):
        ecfg = self.cfg["enrichment"]
        clus = self._selected_clusters(name, mark, kdir, ecfg["clusters"])
        if not clus or not ecfg["db_dir"]:
            log.warning("enrichment %s %s: no clusters or no db_dir, skipped", name, mark)
            return
        edir = self._mk("comparisons", name, "enrichment")
        letters = sorted(pd.read_csv(clus, sep="\t", header=None)[3].astype(str).unique())
        for letter in letters:
            for region in ecfg["regions"]:
                stem = os.path.join(edir, f"{name}_{mark}.cluster{letter}.{region}")
                if self._done(f"{stem}.csv"):
                    continue
                agg = enrichment.enrich_clust(clus, letter, pooled, self._path(ecfg["db_dir"]), self.gene_bed,
                                              self.igr_bed, region=region)
                if agg is None:
                    log.info("  cluster %s %s: no significant enrichment/depletion", letter, region)
                    open(f"{stem}.csv", "w").close()
                    continue
                agg.to_csv(f"{stem}.csv", index=False)
                self._save(enrichment.plot_enrichment(agg, int(ecfg["n_elements"]), int(ecfg["cutoff_for_overlap"]),
                                                      f"Cluster {letter} · {region}"), stem)

    # ------------------------------------------------------------------ exploratory
    def exploratory(self):
        for mark in sorted(self.sheet["Mark"].unique()):
            rows = self.sheet[self.sheet["Mark"] == mark]
            bws, labs = [], []
            for _, r in rows.iterrows():
                for rep, inp in ((r["Rep1"], r["Input1"]), (r["Rep2"], r["Input2"])):
                    if rep:
                        bws.append(self.replicate_bw(rep, inp))
                        labs.append(f"{r['Condition']}_{_stem(rep)}")
            if len(bws) < 2:
                continue
            edir = self._mk("exploratory")
            conds = list(dict.fromkeys(rows["Condition"]))
            cmap = {cnd: PALETTE[k % len(PALETTE)] for k, cnd in enumerate(conds)}
            colors = [cmap[next(cnd for cnd in conds if l.startswith(cnd + "_"))] for l in labs]
            if self.cfg["pca"]["enabled"] and not self._done(os.path.join(edir, f"{mark}_PCA_plot.{self.formats[0]}")):
                fig, _ = explore.plot_pca(bws, colors, labs, plot_title=f"{mark} · PCA")
                self._save(fig, os.path.join(edir, f"{mark}_PCA_plot"))
            if self.cfg["correlation"]["enabled"] and not self._done(os.path.join(edir, f"{mark}_correlation_plot.{self.formats[0]}")):
                fig, _ = explore.plot_correlation(bws, colors, labs, self.cfg["correlation"]["method"],
                                                  plot_title=f"{mark} · {self.cfg['correlation']['method']} correlation")
                self._save(fig, os.path.join(edir, f"{mark}_correlation_plot"))

    # ------------------------------------------------------------------ run
    def run(self):
        os.makedirs(self.out, exist_ok=True)
        for _, row in self.sheet.iterrows():
            self.final_bw(row["Condition"], row["Mark"])
        self.exploratory()
        marks = {cond: set(self.sheet.loc[self.sheet["Condition"] == cond, "Mark"]) for cond in self.sheet["Condition"]}
        for comp in map(parse_comparison, self.cfg["comparisons"]):
            common = sorted(marks.get(comp["treatment"], set()) & marks.get(comp["control"], set()))
            wanted = comp.get("marks", "all")
            if wanted not in (None, "all"):
                common = [m for m in common if m in wanted]
            for mark in common:
                log.info("comparison %s_%s, mark %s", comp["treatment"], comp["control"], mark)
                self.compare(comp, mark)
