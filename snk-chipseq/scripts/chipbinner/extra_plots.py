"""Companion figures written by ChIPbinner / snakemake_chipbinner, in Plotly.

* smooth_scatter       pre_clust's smoothScatter (<prefix>_smoothScatter.jpeg)
* cluster_scatter      clus.py's HDBSCAN scatter (clus.<name>.<minpts>.<minsamps>.png)
* genic_vs_intergenic  compare.genic.vs.intergenic.py (<comp>_<mark>.genic.vs.intergenic.png)
* genic_legend         genic_intergenic_scatterplot's legend.png
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .genomics import read_regions
from .plots import ARROW, FONT, _genic_legend_shapes

CLUSTER_DEEP = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3", "#8C8C8C"]


def _square_layout(fig, w, h, title, size=14):
    fig.update_layout(template="simple_white", width=w, height=h, font=dict(family=FONT, size=size),
                      title=dict(text=f"<b>{title}</b>", x=0.5) if title else None)


def smooth_scatter(x, y, xlab, ylab, nbins=128):
    """graphics::smoothScatter equivalent: binned 2-D density, white -> blue, darkest = densest."""
    h, xe, ye = np.histogram2d(x, y, bins=nbins)
    from scipy.ndimage import gaussian_filter

    z = gaussian_filter(h.T, 1.2)
    z = np.where(z > 0, z, np.nan)
    blues = [[0, "#ffffff"], [0.15, "#c6dbef"], [0.4, "#6baed6"], [0.7, "#2171b5"], [1, "#08306b"]]
    fig = go.Figure(go.Heatmap(x=(xe[:-1] + xe[1:]) / 2, y=(ye[:-1] + ye[1:]) / 2, z=np.sqrt(z),
                               colorscale=blues, showscale=False, hoverinfo="skip"))
    _square_layout(fig, 560, 560, "")
    fig.update_xaxes(title_text=xlab, mirror=True, ticks="outside", title_font_size=18)
    fig.update_yaxes(title_text=ylab, mirror=True, ticks="outside", title_font_size=18)
    return fig


def cluster_scatter(mat: np.ndarray, labels: np.ndarray, probabilities: np.ndarray | None, xlab, ylab,
                    max_points=150000, seed=0):
    """clus.py figure: bins coloured by HDBSCAN cluster (seaborn 'deep'), desaturated by
    membership probability, noise grey. x = wild type, y = treated."""
    n = len(labels)
    idx = np.arange(n)
    if n > max_points:
        idx = np.sort(np.random.default_rng(seed).choice(n, max_points, replace=False))
    fig = go.Figure()
    lab = labels[idx]
    for c in sorted(set(lab)):
        m = idx[lab == c]
        col = "#808080" if c < 0 else CLUSTER_DEEP[c % len(CLUSTER_DEEP)]
        op = 0.25 if probabilities is None or c < 0 else None
        fig.add_trace(go.Scattergl(
            x=mat[m, 1], y=mat[m, 0], mode="markers", name="noise" if c < 0 else f"cluster {c}",
            marker=dict(size=4, color=col, opacity=op if op is not None else 0.35,
                        line=dict(width=0))))
    _square_layout(fig, 640, 560, "")
    fig.update_xaxes(title_text=xlab, mirror=True, ticks="outside", title_font_size=18)
    fig.update_yaxes(title_text=ylab, mirror=True, ticks="outside", title_font_size=18)
    fig.update_layout(legend=dict(itemsizing="constant"))
    return fig


def _contained(bins: pd.DataFrame, regions: pd.DataFrame) -> np.ndarray:
    """bedtools intersect -a bins -b regions -wa -f 1: bin fully inside one region."""
    out = np.zeros(len(bins), bool)
    for c, sub in regions.groupby("chr"):
        m = (bins["chr"] == c).to_numpy()
        if not m.any():
            continue
        o = np.argsort(sub["start"].to_numpy())
        rs = sub["start"].to_numpy()[o]
        re_ = np.maximum.accumulate(sub["end"].to_numpy()[o])
        j = np.searchsorted(rs, bins.loc[m, "start"].to_numpy(), side="right") - 1
        ok = j >= 0
        hit = np.zeros(m.sum(), bool)
        hit[ok] = re_[j[ok]] >= bins.loc[m, "end"].to_numpy()[ok]
        out[m] = hit
    return out


def genic_vs_intergenic(pair: pd.DataFrame, gene_bed: str, igr_bed: str, xlab, ylab,
                        max_points=200000, seed=0):
    """compare.genic.vs.intergenic.py: all bins | + intergenic (red) | + genic (blue), y = x line.

    Genic / intergenic = bins fully contained in a gene-body / intergenic region.
    """
    gen = _contained(pair, read_regions(gene_bed))
    igr = _contained(pair, read_regions(igr_bed))
    n = len(pair)
    idx = np.arange(n)
    if n > max_points:
        idx = np.sort(np.random.default_rng(seed).choice(n, max_points, replace=False))
    x, y = pair["x"].to_numpy()[idx], pair["y"].to_numpy()[idx]
    g, i = gen[idx], igr[idx]
    fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.07,
                        subplot_titles=["All bins", f"Intergenic bins ({igr.sum():,})", f"Genic bins ({gen.sum():,})"])
    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    for col in (1, 2, 3):
        fig.add_trace(go.Scattergl(x=x, y=y, mode="markers", name="All bins", showlegend=col == 1,
                                   marker=dict(size=2, color="silver", opacity=0.15)), row=1, col=col)
        if col == 2:
            fig.add_trace(go.Scattergl(x=x[i], y=y[i], mode="markers", name="Intergenic bins",
                                       marker=dict(size=2, color="#d62728", opacity=0.12)), row=1, col=col)
        if col == 3:
            fig.add_trace(go.Scattergl(x=x[g], y=y[g], mode="markers", name="Genic bins",
                                       marker=dict(size=2, color="#1f4fd6", opacity=0.12)), row=1, col=col)
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color="black", width=1),
                                 showlegend=False, hoverinfo="skip"), row=1, col=col)
        fig.update_xaxes(title_text=f"<b>{xlab}</b>", title_font_size=18, mirror=True, ticks="outside", row=1, col=col)
        fig.update_yaxes(title_text=f"<b>{ylab}</b>", title_font_size=18, mirror=True, ticks="outside", row=1, col=col)
    _square_layout(fig, 1500, 500, "")
    fig.update_layout(legend=dict(itemsizing="constant", orientation="h", y=-0.2, x=0.5, xanchor="center"))
    return fig


def genic_legend(size=260):
    """The standalone 2-D colour x opacity key (legend.png)."""
    fig = go.Figure()
    shapes, ann = _genic_legend_shapes((0.18, 0.95), (0.18, 0.95), np.array([255.0, 255, 255]), "#222", 14)
    fig.update_layout(width=size, height=size, margin=dict(l=10, r=10, t=10, b=10), shapes=shapes,
                      annotations=ann, xaxis=dict(visible=False), yaxis=dict(visible=False),
                      plot_bgcolor="white", paper_bgcolor="white")
    return fig
