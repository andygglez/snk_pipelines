#!/usr/bin/env python3
"""Differential expression with PyDESeq2 (port extension, not in nf-core/rnaseq).

Reads a tximport gene table (gene_id, gene_name, <sample>...), rounds it to
integer counts, fits DESeq2 with design ~condition and writes one results table
per comparison (treatment vs baseline). A model is fitted once per distinct
baseline, with that baseline as the reference level, so LFC shrinkage can use
the treatment coefficient directly.

Outputs (in --outdir):
  deseq2.normalized_counts.tsv, deseq2.vst.tsv, deseq2.size_factors.tsv,
  deseq2.pca.pdf, deseq2.pca.tsv
  <treatment>_vs_<baseline>/<treatment>_vs_<baseline>.deseq2.results.tsv
  <treatment>_vs_<baseline>/<treatment>_vs_<baseline>.deseq2.significant.tsv
"""

import argparse
import csv
import os

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def read_conditions(samplesheet):
    conditions = {}
    with open(samplesheet, newline="") as fh:
        for row in csv.DictReader(fh):
            conditions[row["sample"].strip()] = row["condition"].strip()
    return conditions


def size_factors(dds):
    if "size_factors" in dds.obs:
        return dds.obs["size_factors"].to_numpy()
    return np.asarray(dds.obsm["size_factors"])


def fit(counts, metadata, baseline, args, inference):
    meta = metadata.copy()
    levels = [baseline] + sorted(c for c in meta["condition"].unique() if c != baseline)
    meta["condition"] = pd.Categorical(meta["condition"], categories=levels)
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition",
                       refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    return dds


def shrink_coeff(dds, treatment):
    cols = list(dds.varm["LFC"].columns)
    for name in (f"condition[T.{treatment}]", f"condition_{treatment}_vs_"):
        hits = [c for c in cols if c.startswith(name)]
        if hits:
            return hits[0]
    raise ValueError(f"No LFC coefficient for {treatment!r} in {cols}")


def pca_plot(vst, metadata, pdf, tsv):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top = vst.loc[:, vst.var(axis=0).sort_values(ascending=False).index[:500]]
    x = top.to_numpy() - top.to_numpy().mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    pcs = u * s
    var = s ** 2 / np.sum(s ** 2) * 100
    pd.DataFrame({"sample": top.index, "condition": metadata.loc[top.index, "condition"],
                  "PC1": pcs[:, 0], "PC2": pcs[:, 1] if pcs.shape[1] > 1 else 0.0}
                 ).to_csv(tsv, sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 5))
    for cond in sorted(metadata["condition"].unique()):
        idx = [i for i, smp in enumerate(top.index) if metadata.loc[smp, "condition"] == cond]
        ax.scatter(pcs[idx, 0], pcs[idx, 1] if pcs.shape[1] > 1 else np.zeros(len(idx)), label=cond, s=40)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)" if len(var) > 1 else "PC2")
    ax.set_title("PCA of VST counts (top 500 variable genes)")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(pdf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--samplesheet", required=True, help="CSV with sample (sample ID) and condition columns")
    ap.add_argument("--comparison", action="append", required=True, metavar="BASELINE:TREATMENT")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-counts", type=int, default=10)
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--shrink", action="store_true")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    table = pd.read_csv(args.counts, sep="\t")
    gene_names = table.set_index("gene_id")["gene_name"]
    raw = table.drop(columns=["gene_name"]).set_index("gene_id")
    counts = raw.round().astype(int).T  # samples x genes
    counts = counts.loc[:, counts.sum(axis=0) >= args.min_counts]

    conditions = read_conditions(args.samplesheet)
    missing = [s for s in counts.index if s not in conditions]
    if missing:
        raise ValueError(f"Samples {missing} of {args.counts} are not in {args.samplesheet}")
    metadata = pd.DataFrame({"condition": [conditions[s] for s in counts.index]}, index=counts.index)

    comparisons = [tuple(c.split(":", 1)) for c in args.comparison]
    n = metadata["condition"].value_counts()
    for b, t in comparisons:
        for cond in (b, t):
            if n.get(cond, 0) < 2:
                raise ValueError(f"Condition {cond!r} has {n.get(cond, 0)} quantified samples in "
                                 f"{args.counts}; DESeq2 needs at least 2.")

    inference = DefaultInference(n_cpus=args.threads)
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda name: os.path.join(args.outdir, name)  # noqa: E731

    models = {}
    for b in dict.fromkeys(b for b, _ in comparisons):
        models[b] = fit(counts, metadata, b, args, inference)

    # Run-level tables come from the first model (they do not depend on the reference level).
    dds = models[comparisons[0][0]]
    norm = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names).T
    norm.insert(0, "gene_name", gene_names.reindex(norm.index).values)
    norm.to_csv(out("deseq2.normalized_counts.tsv"), sep="\t", index_label="gene_id")
    pd.DataFrame({"sample": dds.obs_names, "size_factor": size_factors(dds)}).to_csv(
        out("deseq2.size_factors.tsv"), sep="\t", index=False)
    dds.vst(use_design=False)
    vst = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)
    vst_out = vst.T.copy()
    vst_out.insert(0, "gene_name", gene_names.reindex(vst_out.index).values)
    vst_out.to_csv(out("deseq2.vst.tsv"), sep="\t", index_label="gene_id")
    pca_plot(vst, metadata, out("deseq2.pca.pdf"), out("deseq2.pca.tsv"))

    for b, t in comparisons:
        name = f"{t}_vs_{b}"
        stats = DeseqStats(models[b], contrast=["condition", t, b], alpha=args.padj,
                           inference=inference, quiet=True)
        stats.summary()
        res = stats.results_df.copy()
        res["log2FoldChange_unshrunk"] = res["log2FoldChange"]
        if args.shrink:
            stats.lfc_shrink(coeff=shrink_coeff(models[b], t))
            res["log2FoldChange"] = stats.results_df["log2FoldChange"]
            res["lfcSE"] = stats.results_df["lfcSE"]
        res.insert(0, "gene_name", gene_names.reindex(res.index).values)
        res = res.sort_values(["padj", "pvalue"], na_position="last")
        d = os.path.join(args.outdir, name)
        os.makedirs(d, exist_ok=True)
        res.to_csv(os.path.join(d, f"{name}.deseq2.results.tsv"), sep="\t", index_label="gene_id")
        sig = res[(res["padj"] < args.padj) & (res["log2FoldChange"].abs() >= args.lfc)]
        sig.to_csv(os.path.join(d, f"{name}.deseq2.significant.tsv"), sep="\t", index_label="gene_id")
        print(f"{name}: {len(sig)} significant genes (padj < {args.padj}, |log2FC| >= {args.lfc}; "
              f"{(sig['log2FoldChange'] > 0).sum()} up, {(sig['log2FoldChange'] < 0).sum()} down)")


if __name__ == "__main__":
    main()
