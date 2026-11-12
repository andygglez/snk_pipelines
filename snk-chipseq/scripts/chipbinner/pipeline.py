"""Config-driven ChIPbinner workflow producing the snakemake_chipbinner output tree.

Layout (``<outdir>`` = chipbinner.outdir; names as in andygglez/snakemake_chipbinner):

    <outdir>/bam2bed/<sample>.binned.bed                 bin counts (BAM inputs)
    <outdir>/bam2bed/<sample>.bed                         read BED (only with binning.keep_read_bed)
    <outdir>/{MS,notMS}bigwigs/<sample>.binned.bw         norm_bw per replicate
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

Each step is a function of :mod:`chipbinner.steps`; snk-chipseq runs the same
functions as Snakemake rules. This class adds the samplesheet, the skip-if-done
logic and the snakemake_chipbinner file layout.
"""
from __future__ import annotations

import logging
import os
import shutil

import numpy as np
import pandas as pd
import yaml

from . import cluster, steps
from .genomics import annotation_paths, write_bed

log = logging.getLogger("chipbinner")

DEFAULTS = {
    "outdir": "outdir",
    "genome_assembly": "hg38",
    "annotations": {"gene_bed": "", "igr_bed": "", "blacklist": "", "chrom_sizes": ""},
    "samplesheet": "",
    "sample_folder": "",
    "binning": {"bin_size": 10000, "windows": "", "keep_read_bed": False},
    "MS_normalize": False,
    "MS_coefficients": {},
    "merge_replicates": True,
    "norm": {"use_input": False, "depth_norm": True, "raw_count_cutoff": 100, "pseudocount": 1e-15},
    "comparisons": [],
    "plot_defaults": {"range": "auto", "pow": 0.75, "square_auto_range": True, "auto_quantiles": [0, 1]},
    "genic_plot": {"hexbins": 150, "show_legend": False, "legend_pos": "right", "show_scales": True},
    "density_plot": {"hexbins": 150, "show_scales": True, "include_additional_density_plot": True,
                     "filter_extreme_bins": True, "cluster_colors": "bright"},
    "clustering": {"enabled": True, "minpts": [100, 400, 700, 1000], "minsamps": [100, 400, 700, 1000],
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
PALETTE = steps.PALETTE


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
            steps.bin_bam(src, out, windows=self._path(b["windows"]) if b.get("windows") else "",
                          chrom_sizes=self.chrom_sizes, bin_size=int(b["bin_size"]),
                          read_bed=os.path.join(d, f"{key}.bed") if b.get("keep_read_bed") else "")
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
        ibed = self.binned(inp) if inp else ""
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
        steps.norm_replicate(bed, out, self.chrom_sizes, self.blacklist, input_binned=ibed,
                             use_input=n["use_input"], depth_norm=n["depth_norm"], pseudocount=n["pseudocount"],
                             raw_count_cutoff=n["raw_count_cutoff"], scaling_factor=sf)
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
        if not self._done(out):
            steps.merge_replicates(reps, out, self.chrom_sizes)
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
        if len(pair) < steps.MIN_SHARED_BINS:
            log.warning("%s %s: only %d shared bins, skipping", name, mark, len(pair))
            return
        # ---- plot range: auto-detected or from config; recorded next to the plots
        rng, src = steps.plot_spec(pair, pp["range"], pp.get("square_auto_range", True),
                                   pp.get("auto_quantiles", (0.0, 1.0)))
        steps.write_plot_range(f"{pre}.plot_range.yaml", name, mark, rng, src, pp["pow"])
        spec = {"x": rng[:2], "y": rng[2:]}
        log.info("%s %s: range x=%s y=%s (%s), pow=%s", name, mark, spec["x"], spec["y"], src, pp["pow"])

        # ---- genic/intergenic scatterplot (+ legend.png)
        if not self._done(f"{pre}.scatterplot.{self.formats[0]}"):
            steps.scatter_plot(bw_t, bw_c, self.gene_bed, self.igr_bed, pre, mark, t, c, spec, pp["pow"],
                               self.cfg["genic_plot"], self.formats)
            steps.legend(os.path.join(cdir, "legend.png"))

        # ---- genic vs intergenic panels
        if not self._done(f"{pre}.genic.vs.intergenic.png"):
            steps.genic_vs_intergenic(pair, self.gene_bed, self.igr_bed, pre, t, c)

        cl = self.cfg["clustering"]
        kdir = os.path.join(self.out, "comparisons", name, "clustering")
        if not cl.get("enabled", True):  # stop after the scatterplots
            if self.cfg["differential"]["enabled"]:
                self.differential(comp, mark, name, kdir)
            return

        # ---- pre_clust (+ smoothScatter)
        mat_f, bed_f = f"{pre}_preclust_mat.csv", f"{pre}_preclust_pooled.bed"
        if not self._done(mat_f, bed_f):
            steps.preclust(bw_t, bw_c, pre, t, c)

        # ---- HDBSCAN grid, annotation, density-based scatterplots
        kdir = self._mk("comparisons", name, "clustering")
        for i in cl["minpts"]:
            for j in cl["minsamps"]:
                stem = f"{name}_{mark}_clustered"
                lab_f = os.path.join(kdir, f"clus.{stem}.{i}.{j}.txt")
                if not self._done(lab_f):
                    steps.hdbscan(mat_f, i, j, lab_f, os.path.join(kdir, f"clus.{stem}.{i}.{j}.png"),
                                  os.path.join(kdir, f"log.{i}.{j}.txt"), t, c,
                                  cl["cluster_selection_method"], cl.get("cores", 1))
                for n in cl["number_of_clusters"]:
                    ann = os.path.join(kdir, f"{name}_{mark}.{n}clusters.{i}.{j}.annotated_clusters.bed")
                    dens = os.path.join(kdir, f"{name}_{mark}_{n}clusters.{i}.{j}.plots_density_scatterplot")
                    if self._done(ann, f"{dens}.{self.formats[0]}"):
                        continue
                    steps.density_with_spec(mat_f, bed_f, lab_f, spec, pp["pow"], n, bw_t, bw_c, self.gene_bed,
                                            self.igr_bed, ann, dens, mark, t, c, self.cfg["density_plot"],
                                            self.formats)

        if self.cfg["differential"]["enabled"]:
            self.differential(comp, mark, name, kdir)
        if self.cfg["enrichment"]["enabled"]:
            self.enrichment(comp, mark, name, kdir, bed_f)

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
        steps.differential(tb, cb, self.blacklist, t, c, ddir, clusters=clus or "",
                           bootstrap=dcfg["bootstrap"], K=dcfg["K"])

    def enrichment(self, comp, mark, name, kdir, pooled):
        ecfg = self.cfg["enrichment"]
        clus = self._selected_clusters(name, mark, kdir, ecfg["clusters"])
        if not clus or not ecfg["db_dir"]:
            log.warning("enrichment %s %s: no clusters or no db_dir, skipped", name, mark)
            return
        steps.enrichment(clus, pooled, self._path(ecfg["db_dir"]), self.gene_bed, self.igr_bed,
                         self._mk("comparisons", name, "enrichment"), f"{name}_{mark}", regions=ecfg["regions"],
                         n_elements=ecfg["n_elements"], cutoff_for_overlap=ecfg["cutoff_for_overlap"],
                         formats=self.formats, skip_done=not self.force)

    # ------------------------------------------------------------------ exploratory
    def exploratory(self):
        for mark in sorted(self.sheet["Mark"].unique()):
            rows = self.sheet[self.sheet["Mark"] == mark]
            bws, labs, conds = [], [], []
            for _, r in rows.iterrows():
                for rep, inp in ((r["Rep1"], r["Input1"]), (r["Rep2"], r["Input2"])):
                    if rep:
                        bws.append(self.replicate_bw(rep, inp))
                        labs.append(f"{r['Condition']}_{_stem(rep)}")
                        conds.append(r["Condition"])
            if len(bws) < 2:
                continue
            edir = self._mk("exploratory")
            do_pca = self.cfg["pca"]["enabled"] and not self._done(
                os.path.join(edir, f"{mark}_PCA_plot.{self.formats[0]}"))
            do_cor = self.cfg["correlation"]["enabled"] and not self._done(
                os.path.join(edir, f"{mark}_correlation_plot.{self.formats[0]}"))
            if do_pca or do_cor:
                steps.exploratory(bws, labs, conds, mark, edir, self.formats, pca=do_pca, correlation=do_cor,
                                  method=self.cfg["correlation"]["method"])

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
