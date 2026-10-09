#!/usr/bin/env python3
"""Preranked GSEA of a PyDESeq2 results table with GSEApy (port extension).

Genes are ranked by the Wald statistic. Gene IDs are mapped to --gene-id
(gene_name or gene_id), genes without a statistic are dropped, and duplicate
names keep the entry with the largest |stat|.

Outputs (in --outdir): gsea.rnk, gsea.results.tsv, gsea.top_terms.pdf and
enrichment_plots/<n>_<term>.pdf for the --top-terms best terms by FDR.
"""

import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")
import gseapy as gp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

COLUMNS = ["Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes"]


def placeholder_pdf(path, text):
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.axis("off")
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=10)
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--gmt", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--gene-id", choices=["gene_name", "gene_id"], default="gene_name")
    ap.add_argument("--min-size", type=int, default=15)
    ap.add_argument("--max-size", type=int, default=500)
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--top-terms", type=int, default=10)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "enrichment_plots"), exist_ok=True)
    out = lambda name: os.path.join(args.outdir, name)  # noqa: E731

    res = pd.read_csv(args.results, sep="\t")
    res = res[res["stat"].notna()].copy()
    res["id"] = res[args.gene_id].fillna(res["gene_id"]).astype(str)
    res = res.reindex(res["stat"].abs().sort_values(ascending=False).index).drop_duplicates("id")
    rnk = res.set_index("id")["stat"].sort_values(ascending=False)
    rnk.to_csv(out("gsea.rnk"), sep="\t", header=False)

    try:
        pre = gp.prerank(rnk=rnk, gene_sets=args.gmt, min_size=args.min_size, max_size=args.max_size,
                         permutation_num=args.permutations, seed=args.seed, threads=args.threads,
                         outdir=None, no_plot=True, verbose=False)
        table = pre.res2d
    except (LookupError, ValueError) as err:  # no gene set passes the size filters
        print(f"WARNING: GSEA found no testable gene sets: {err}")
        pre, table = None, pd.DataFrame(columns=COLUMNS)

    table = table.drop(columns=["Name"], errors="ignore")
    if len(table):
        for col in ("NES", "NOM p-val", "FDR q-val", "FWER p-val"):
            table[col] = pd.to_numeric(table[col])
        table = table.sort_values(["FDR q-val", "NOM p-val"])
    table.to_csv(out("gsea.results.tsv"), sep="\t", index=False)

    if pre is None or not len(table):
        placeholder_pdf(out("gsea.top_terms.pdf"), "No gene set passed the size filters.")
        return

    top = table.head(args.top_terms)
    try:
        gp.dotplot(top, column="FDR q-val", title="Preranked GSEA", cutoff=1.01, top_term=args.top_terms,
                   size=6, figsize=(6, max(3, 0.4 * len(top) + 1.5)), ofname=out("gsea.top_terms.pdf"))
    except Exception as err:  # noqa: BLE001  (plotting must not fail the analysis)
        print(f"WARNING: dotplot failed: {err}")
        placeholder_pdf(out("gsea.top_terms.pdf"), "Dot plot failed; see gsea.results.tsv.")

    for i, term in enumerate(top["Term"], start=1):
        safe = re.sub(r"[^\w.-]+", "_", str(term))[:80]
        pre.plot(terms=[term], ofname=os.path.join(args.outdir, "enrichment_plots", f"{i:02d}_{safe}.pdf"))
        plt.close("all")


if __name__ == "__main__":
    main()
