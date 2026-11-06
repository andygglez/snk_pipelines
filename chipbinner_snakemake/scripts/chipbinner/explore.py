"""Ports of ChIPbinner::plot_PCA and plot_correlation (Plotly)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .genomics import read_bins, sample_label
from .normalize import _common_bins

FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"


def signal_matrix(bigwigs, sample_labels=None):
    frames = _common_bins([read_bins(f) for f in bigwigs])
    mat = pd.DataFrame({lab: f["score"].to_numpy(float) for lab, f in
                        zip(sample_labels or [sample_label(f) for f in bigwigs], frames)})
    return mat.dropna()


def pca(mat: pd.DataFrame):
    """stats::prcomp(t(signal), center = TRUE, scale. = TRUE).

    Bins with zero variance across samples cannot be scaled (prcomp stops);
    they are dropped here.
    """
    x = mat.to_numpy().T
    sd = x.std(axis=0, ddof=1)
    keep = sd > 0
    z = (x[:, keep] - x[:, keep].mean(axis=0)) / sd[keep]
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    scores = u * s
    var = s ** 2 / (x.shape[0] - 1)
    return scores, var / var.sum()


def plot_pca(bigwigs, colors, sample_labels=None, plot_title="", font_size=14, point_size=12):
    if len(bigwigs) < 2:
        raise ValueError("At least two bigWig files are needed for the PCA plot.")
    if len(colors) != len(bigwigs):
        raise ValueError("The number of colors must match the number of files.")
    mat = signal_matrix(bigwigs, sample_labels)
    scores, prop = pca(mat)
    labs = list(mat.columns)
    fig = go.Figure()
    fig.add_hline(y=0, line=dict(color="rgba(0,0,0,0.12)"))
    fig.add_vline(x=0, line=dict(color="rgba(0,0,0,0.12)"))
    for k, lab in enumerate(labs):
        fig.add_trace(go.Scatter(x=[scores[k, 0]], y=[scores[k, 1] if scores.shape[1] > 1 else 0],
                                 mode="markers", name=lab,
                                 marker=dict(size=point_size, color=colors[k], line=dict(color="white", width=1))))
    pc = [f"PC{i + 1} ({p * 100:.2g}%)" for i, p in enumerate(prop)]
    fig.update_layout(template="simple_white", width=760, height=560, font=dict(family=FONT, size=font_size),
                      title=dict(text=f"<b>{plot_title}</b>", x=0.5), legend_title_text="Samples",
                      xaxis=dict(title=pc[0], mirror=True, ticks="outside"),
                      yaxis=dict(title=pc[1] if len(pc) > 1 else "PC2", mirror=True, ticks="outside"))
    return fig, pd.DataFrame(scores, index=labs, columns=pc)


def plot_correlation(bigwigs, colors, sample_labels=None, correlation_method="pearson", plot_title="",
                     key_title="", font_size=13):
    """Correlation heatmap with hierarchical clustering (gplots::heatmap.2 defaults:
    Euclidean distance between rows of the correlation matrix, complete linkage)."""
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist

    if len(bigwigs) < 2:
        raise ValueError("At least two bigWig files are needed for the correlation plot.")
    if len(colors) != len(bigwigs):
        raise ValueError("The number of colors must match the number of files.")
    mat = signal_matrix(bigwigs, sample_labels)
    corr = mat.corr(method=correlation_method)
    labs = list(corr.columns)
    Z = linkage(pdist(corr.to_numpy()), method="complete")
    dn = dendrogram(Z, no_plot=True, labels=labs)
    order = [labs.index(x) for x in dn["ivl"]]
    c = corr.to_numpy()[np.ix_(order, order)]
    ol = [labs[i] for i in order]
    n = len(labs)
    fig = go.Figure()
    # heatmap on x/y; dendrograms on x2 (top) and y2 (left); colour strips as shapes
    fig.add_trace(go.Heatmap(z=c, x=list(range(n)), y=list(range(n)), colorscale="RdBu_r", zmid=None,
                             colorbar=dict(title=key_title or correlation_method, len=0.5, y=0.3),
                             text=np.round(c, 3), texttemplate="%{text}", hovertemplate=
                             "%{customdata[0]} vs %{customdata[1]}: %{z:.3f}<extra></extra>",
                             customdata=np.dstack(np.meshgrid(ol, ol)), xaxis="x", yaxis="y"))
    for xs, ys in zip(dn["icoord"], dn["dcoord"]):
        xs = (np.array(xs) - 5) / 10  # leaf positions 5, 15, ... -> 0, 1, ...
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color="#444", width=1.2),
                                 xaxis="x2", yaxis="y2", hoverinfo="skip", showlegend=False))
        fig.add_trace(go.Scatter(x=ys, y=xs, mode="lines", line=dict(color="#444", width=1.2),
                                 xaxis="x3", yaxis="y3", hoverinfo="skip", showlegend=False))
    col = {lab: colors[i] for i, lab in enumerate(labs)}
    for k, lab in enumerate(ol):
        fig.add_shape(type="rect", xref="x", yref="paper", x0=k - 0.5, x1=k + 0.5, y0=0.765, y1=0.79,
                      fillcolor=col[lab], line=dict(width=0))
        fig.add_shape(type="rect", xref="paper", yref="y", x0=0.185, x1=0.21, y0=k - 0.5, y1=k + 0.5,
                      fillcolor=col[lab], line=dict(width=0))
    rng = [-0.5, n - 0.5]
    fig.update_layout(
        template="simple_white", width=300 + 70 * n, height=300 + 70 * n, font=dict(family=FONT, size=font_size),
        title=dict(text=f"<b>{plot_title}</b>", x=0.5), showlegend=False,
        xaxis=dict(domain=[0.22, 0.86], range=rng, tickvals=list(range(n)), ticktext=ol, tickangle=-45,
                   showline=False, ticks=""),
        yaxis=dict(domain=[0, 0.755], range=rng, tickvals=list(range(n)), ticktext=ol, side="right",
                   showline=False, ticks="", anchor="x"),
        xaxis2=dict(domain=[0.22, 0.86], range=rng, visible=False, anchor="y2"),
        yaxis2=dict(domain=[0.80, 0.97], visible=False, anchor="x2"),
        xaxis3=dict(domain=[0.0, 0.18], visible=False, autorange="reversed", anchor="y3"),
        yaxis3=dict(domain=[0, 0.755], range=rng, visible=False, anchor="x3"))
    return fig, corr
