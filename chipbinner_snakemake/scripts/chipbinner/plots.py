"""Plotly ports of ChIPbinner::genic_intergenic_scatterplot and density_based_scatterplot.

Both compare two conditions bin by bin: x = wild-type score, y = treated score.

* Genic/intergenic panel - bins are hex-binned (hexbin::hexbin port); each hex is
  coloured by the mean genic class of its bins (+1 genic only, -1 intergenic only,
  0 both; RdYlBu, blue = genic) and its opacity grows with count^pow, with a
  2-D colour x opacity legend inset.
* Density panel (optional) - MASS::kde2d density as filled iso-bands (magma).
* Density-based clusters panel - the same density as a white haze, bins coloured
  by their HDBSCAN cluster (A, B, C ... from annotate_clust) and labelled at the
  cluster centroid.
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import diverging, hex_to_rgb, sequential, unlabel_rgb

from .cluster import assign_clusters, genic_class, load_pair
from .genomics import annotation_paths, read_regions, sample_label
from .stats import expand_range, hex_geometry, hexbin, kde2d, pretty, quantile_filter

FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"
ARROW = "►"
# pals::alphabet (Green-Armytage), ChIPbinner's cluster palette
ALPHABET = ["#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919", "#005C31", "#2BCE48", "#FFCC99",
            "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380", "#FFA405", "#FFA8BB",
            "#426600", "#FF0010", "#5EF1F2", "#00998F", "#E0FF66", "#740AFF", "#990000", "#FFFF80",
            "#FFE100", "#FF5005"]
# Default: bright, colour-blind-friendly hues that stay legible on a black panel
BRIGHT = ["#FFB000", "#33BBEE", "#EE3377", "#8EE06B", "#B58CFF", "#FF7F3F", "#2EE6C5", "#FFE45C",
          "#FF9DD4", "#7FA6FF"] + ALPHABET
NA_COLOR = "#7f7f7f"


# --------------------------------------------------------------------------- colours
def _to_rgb(c):
    return np.array(unlabel_rgb(c) if c.startswith("rgb") else hex_to_rgb(c), float)


def color_ramp(colors, n):
    """grDevices::colorRampPalette(colors)(n) -> (n, 3) array of RGB 0-255."""
    rgb = np.array([_to_rgb(c) for c in colors])
    pos = np.linspace(0, 1, len(rgb))
    t = np.linspace(0, 1, n)
    return np.stack([np.interp(t, pos, rgb[:, k]) for k in range(3)], axis=1)


RDYLBU50 = color_ramp(diverging.RdYlBu, 50)  # pals::brewer.rdylbu(50): red (intergenic) -> blue (genic)
MAGMA = sequential.Magma


def _rgb_str(a):
    return "rgb(%d,%d,%d)" % tuple(np.clip(np.rint(a), 0, 255))


# --------------------------------------------------------------------------- data
def _labels(treated, wildtype, xaxis_label, yaxis_label):
    return (xaxis_label or sample_label(wildtype)), (yaxis_label or sample_label(treated))


def _annotations(genome_assembly, gene_bed, igr_bed):
    p = annotation_paths(genome_assembly)
    return gene_bed or p["gene"], igr_bed or p["igr"]


def _limits(hexcells, min_x, max_x, min_y, max_y):
    """Axis limits: user values, else the range of hex centres holding > 10 bins
    (ChIPbinner's default for min_x/min_y); then expanded by 5 % on each side."""
    dense = hexcells[hexcells["count"] > 10]
    if dense.empty:
        dense = hexcells
    lo_x = dense["x"].min() if min_x is None else float(min_x)
    hi_x = dense["x"].max() if max_x is None else float(max_x)
    lo_y = dense["y"].min() if min_y is None else float(min_y)
    hi_y = dense["y"].max() if max_y is None else float(max_y)
    return expand_range(lo_x, hi_x), expand_range(lo_y, hi_y)


def resolve_range(spec):
    """Plot range spec -> (min_x, max_x, min_y, max_y); None entries mean automatic.

    spec: None / "auto" | [min, max] (both axes) | {"x": [min, max] | "auto", "y": ...}.
    """
    if spec is None or (isinstance(spec, str) and spec.lower() == "auto"):
        return (None, None, None, None)
    if isinstance(spec, dict):
        x = resolve_range(spec.get("x", "auto"))
        y = resolve_range(spec.get("y", "auto"))
        return (x[0], x[1], y[2], y[3])
    lo, hi = (float(v) for v in spec)
    if not lo < hi:
        raise ValueError(f"range min must be < max: {spec}")
    return (lo, hi, lo, hi)


def auto_range(pair, xbins=150, quantiles=(0.001, 0.999), min_count=10, square=False):
    """Automatic plot range, as ChIPbinner::genic_intergenic_scatterplot decides it.

    R: bins trimmed to the 0.1-99.9 % quantiles of x and y, hexbin(x, y, xbins = hexbins),
    then min_x = min(hex x centre with count > 10) and min_y likewise, each axis on its
    own; max_x / max_y are user inputs in R, so the automatic max applies the same rule
    (max hex centre with count > 10). Pass the plot's own ``hexbins`` as ``xbins``.
    square=True (not ChIPbinner) would force one shared range on both axes.
    Returns [min_x, max_x, min_y, max_y] before ChIPbinner's 5 % expand_range."""
    pdat = quantile_filter(pair, *quantiles)
    cells, _ = hexbin(pdat["x"].to_numpy(), pdat["y"].to_numpy(), xbins=xbins)
    dense = cells[cells["count"] > min_count]
    if dense.empty:
        dense = cells
    r = []
    for ax in ("x", "y"):
        lo, hi = float(dense[ax].min()), float(dense[ax].max())
        if not lo < hi:  # too few bins for > min_count hexes to span this axis (tiny data / fine grid)
            lo, hi = float(cells[ax].min()), float(cells[ax].max())
            warnings.warn(f"auto range: hexes with > {min_count} bins do not span {ax}; "
                          f"using every occupied hex ({lo:.3g}..{hi:.3g})")
        r += [lo, hi]
    if square:
        lo, hi = min(r[0], r[2]), max(r[1], r[3])
        r = [lo, hi, lo, hi]
    return r


def hex_table(pdat, xbins, pow, alpha_range):
    """Hex cells with count, mean genic class (prop) and opacity (ggplot scale_alpha)."""
    cells, cid = hexbin(pdat["x"].to_numpy(), pdat["y"].to_numpy(), xbins=xbins)
    prop = pd.Series(pdat["r"].to_numpy(float)).groupby(cid).mean()
    cells["prop"] = prop.reindex(cells["cell"]).to_numpy()
    t = cells["count"].to_numpy(float) ** pow
    span = t.max() - t.min()
    a = alpha_range[0] + ((t - t.min()) / span if span > 0 else 1.0) * (alpha_range[1] - alpha_range[0])
    cells["alpha"] = np.clip(a, 0, 1)
    return cells


# --------------------------------------------------------------------------- traces
def _hex_traces(cells, bg_rgb, xa, ya, n_alpha=32):
    """Filled hexagon polygons, grouped by (colour, opacity) so each group is one trace.

    Colours are pre-composited over the panel background, so polygons can get a
    hairline outline of their own colour (no seams between neighbouring hexes).
    """
    geom = hex_geometry(cells.attrs["dx"], cells.attrs["dy"])
    ci = np.rint((cells["prop"].to_numpy() + 1) / 2 * 49).astype(int)
    ai = np.rint(cells["alpha"].to_numpy() * (n_alpha - 1)).astype(int)
    keep = ai > 0
    traces = []
    groups = pd.DataFrame({"c": ci[keep], "a": ai[keep]}, index=cells.index[keep]).groupby(["c", "a"]).groups
    for (c, a), idx in groups.items():
        sub = cells.loc[idx]
        alpha = a / (n_alpha - 1)
        col = _rgb_str(RDYLBU50[c] * alpha + bg_rgb * (1 - alpha))
        vx = sub["x"].to_numpy()[:, None] + geom[:, 0][None, :]
        vy = sub["y"].to_numpy()[:, None] + geom[:, 1][None, :]
        nan = np.full((len(sub), 1), np.nan)
        xs = np.hstack([vx, vx[:, :1], nan]).ravel()
        ys = np.hstack([vy, vy[:, :1], nan]).ravel()
        traces.append(go.Scatter(x=xs, y=ys, mode="lines", fill="toself", fillcolor=col,
                                 line=dict(color=col, width=0.4), hoverinfo="skip",
                                 showlegend=False, xaxis=xa, yaxis=ya))
    pct = (cells["prop"].to_numpy() + 1) / 2 * 100
    traces.append(go.Scatter(
        x=cells["x"], y=cells["y"], mode="markers", xaxis=xa, yaxis=ya, showlegend=False,
        marker=dict(size=7, opacity=0), customdata=np.c_[cells["count"], cells["prop"], pct],
        hovertemplate="x %{x:.2f}, y %{y:.2f}<br>%{customdata[0]:,} bins<br>"
                      "genic score %{customdata[1]:.2f} (%{customdata[2]:.0f}% genic)<extra></extra>"))
    return traces


def _diagonal(lim, xa, ya, color="rgba(160,160,160,0.6)"):
    lo = min(lim[0][0], lim[1][0])
    hi = max(lim[0][1], lim[1][1])
    return go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color=color, width=1),
                      hoverinfo="skip", showlegend=False, xaxis=xa, yaxis=ya)


def _genic_legend_shapes(dom_x, dom_y, bg_rgb, text_color, font_size=11):
    """2-D legend (9 colours x 9 opacities) drawn as paper-coordinate rectangles."""
    len_ = 9
    pal = color_ramp(diverging.RdYlBu, len_)
    alphas = np.round(np.linspace(0, 255, len_ + 1)[1:]) / 255
    (x0, x1), (y0, y1) = dom_x, dom_y
    w = (x1 - x0) / len_
    h = (y1 - y0) / len_
    shapes = []
    for r in range(len_):          # rows: % genic (bottom red -> top blue)
        for c in range(len_):      # cols: number of bins (left faint -> right opaque)
            col = _rgb_str(pal[r] * alphas[c] + bg_rgb * (1 - alphas[c]))
            shapes.append(dict(type="rect", xref="paper", yref="paper", layer="above",
                               x0=x0 + c * w, x1=x0 + (c + 1) * w, y0=y0 + r * h, y1=y0 + (r + 1) * h,
                               fillcolor=col, line=dict(width=0)))
    shapes.append(dict(type="rect", xref="paper", yref="paper", x0=x0, x1=x1, y0=y0, y1=y1,
                       line=dict(color=text_color, width=0.6), layer="above"))
    ann = [dict(text=f"# of bins {ARROW}", x=(x0 + x1) / 2, y=y0, xref="paper", yref="paper",
                yanchor="top", showarrow=False, font=dict(size=font_size, color=text_color, family=FONT)),
           dict(text=f"% genic {ARROW}", x=x0, y=(y0 + y1) / 2, xref="paper", yref="paper",
                xanchor="right", textangle=-90, showarrow=False,
                font=dict(size=font_size, color=text_color, family=FONT))]
    return shapes, ann


def _density_bands(pdat):
    gx, gy, z = kde2d(pdat["x"].to_numpy(), pdat["y"].to_numpy(), n=100)
    brks = pretty(z, 30)
    return gx, gy, z, brks


def _contour(gx, gy, z, brks, colorscale, xa, ya, colorbar=None):
    return go.Contour(
        x=gx, y=gy, z=z.T, xaxis=xa, yaxis=ya, zmin=brks[0], zmax=brks[-1],
        contours=dict(start=brks[0], end=brks[-1], size=brks[1] - brks[0], coloring="fill", showlines=False),
        line=dict(width=0), colorscale=colorscale, showscale=colorbar is not None, colorbar=colorbar,
        hovertemplate="x %{x:.2f}, y %{y:.2f}<br>density %{z:.3g}<extra></extra>")


# --------------------------------------------------------------------------- layout
class _Panels:
    """Square panels laid out in a row, with facet-style title strips."""

    def __init__(self, n, panel=440, gap=96, ml=104, mr=36, mt=118, mb=96, strip=26):
        self.n, self.panel, self.gap, self.strip = n, panel, gap, strip
        self.width = ml + n * panel + (n - 1) * gap + mr
        self.height = mt + panel + mb
        self.margin = dict(l=ml, r=mr, t=mt, b=mb)
        self.area_w = n * panel + (n - 1) * gap

    def dom_x(self, k):
        a = k * (self.panel + self.gap) / self.area_w
        return [a, a + self.panel / self.area_w]

    def frac_x(self, px):
        return px / self.area_w

    def frac_y(self, px):
        return px / self.panel

    def axis_names(self, k):
        return ("x", "y") if k == 0 else (f"x{k + 1}", f"y{k + 1}")


def _axis(dom, anchor, lim, title, show_scales, axis_color, size, tick_size=16):
    return dict(domain=dom, anchor=anchor, range=list(lim),
                title=dict(text=f"<b>{title}</b>" if title else "", font=dict(size=size), standoff=12),
                tickfont=dict(size=tick_size),
                showgrid=False, zeroline=False, showline=True, mirror=True, linecolor=axis_color,
                linewidth=1, ticks="outside" if show_scales else "", showticklabels=show_scales,
                tickcolor=axis_color, ticklen=4, constrain="domain")


def _strip(lay, dom, text, size, frac):
    lay.shapes += (dict(type="rect", xref="paper", yref="paper", x0=dom[0], x1=dom[1], y0=1.0,
                        y1=1.0 + frac, fillcolor="black", line=dict(width=1, color="black"),
                        layer="above"),)
    lay.annotations += (dict(text=f"<b>{text}</b>", x=(dom[0] + dom[1]) / 2, y=1.0 + frac / 2,
                             xref="paper", yref="paper", showarrow=False,
                             font=dict(color="white", size=size, family=FONT)),)


def _base_layout(p, title, subtitle, font_size):
    lay = go.Layout(
        width=p.width, height=p.height, margin=p.margin, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family=FONT, size=font_size, color="#222"), hovermode="closest",
        title=dict(text=f"<b>{title}</b>" + (f"<br><span style='font-size:{font_size}px;color:#666'>"
                                              f"{subtitle}</span>" if subtitle else ""),
                   x=0.5, xanchor="center", y=1 - 34 / p.height, yanchor="middle", yref="container",
                   font=dict(size=font_size + 7)),
        showlegend=False, shapes=(), annotations=())
    return lay


def _save(fig, out_prefix, formats, scale=2):
    if not out_prefix:
        return []
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    written = []
    for fmt in formats:
        path = f"{out_prefix}.{fmt}"
        if fmt == "html":
            fig.write_html(path, include_plotlyjs=True, full_html=True,
                           config={"displaylogo": False, "toImageButtonOptions": {"format": "svg"}})
        else:
            fig.write_image(path, width=fig.layout.width, height=fig.layout.height,
                            scale=scale if fmt == "png" else 1)
        written.append(path)
    return written


# --------------------------------------------------------------------------- public API
def genic_intergenic_scatterplot(treated_samp_norm_bw, wildtype_samp_norm_bw, genome_assembly="hg38",
                                 gene_bed=None, igr_bed=None, out_prefix=None, formats=("html",),
                                 title_of_plot="", xaxis_label=None, yaxis_label=None,
                                 min_x=None, max_x=None, min_y=None, max_y=None, pow=1.25, hexbins=150,
                                 show_scales=True, show_legend=True, legend_pos="left",
                                 alpha_range=(0.01, 5), quantiles=(0.001, 0.999), panel_size=560,
                                 font_size=14, axis_title_size=22, tick_font_size=16, plot_range=None,
                                 write_table=False):
    """Port of genic_intergenic_scatterplot(): one white panel of genic/intergenic hexbins.

    Defaults follow ChIPbinner (hexbins 150, pow 1.25, alpha range 0.01-5 clipped
    at 1, 0.1 %/99.9 % quantile trimming). ``max_x``/``max_y`` are optional here
    (default: as for min, the extent of hexes holding > 10 bins).
    """
    gene_bed, igr_bed = _annotations(genome_assembly, gene_bed, igr_bed)
    xl, yl = _labels(treated_samp_norm_bw, wildtype_samp_norm_bw, xaxis_label, yaxis_label)
    pair = load_pair(treated_samp_norm_bw, wildtype_samp_norm_bw)
    pair["r"] = genic_class(pair, gene_bed, igr_bed)
    pdat = quantile_filter(pair, *quantiles)
    cells = hex_table(pdat, hexbins, pow, alpha_range)
    if plot_range is not None:
        min_x, max_x, min_y, max_y = resolve_range(plot_range)
    lim = _limits(cells, min_x, max_x, min_y, max_y)

    bg = np.array([255.0, 255, 255])
    p = _Panels(1, panel=panel_size, mt=92)
    lay = _base_layout(p, title_of_plot or f"{yl} vs {xl}", None, font_size)
    lay.xaxis = _axis(p.dom_x(0), "y", lim[0], f"{xl} {ARROW}", show_scales, "black", axis_title_size, tick_font_size)
    lay.yaxis = _axis([0, 1], "x", lim[1], f"{yl} {ARROW}", show_scales, "black", axis_title_size, tick_font_size)
    fig = go.Figure(layout=lay)
    fig.add_traces(_hex_traces(cells, bg, "x", "y"))
    fig.add_trace(_diagonal(lim, "x", "y", "rgba(128,128,128,0.8)"))
    if show_legend:
        x0, x1 = (0.04, 0.30) if legend_pos == "left" else (0.70, 0.96)
        shapes, ann = _genic_legend_shapes((x0 + 0.03, x1), (0.70, 0.96), bg, "#222", font_size - 2)
        fig.layout.shapes += tuple(shapes)
        fig.layout.annotations += tuple(ann)
    _save(fig, out_prefix, formats)
    if write_table and out_prefix:
        pdat.to_csv(f"{out_prefix}.bins.tsv.gz", sep="\t", index=False)
    return fig


def density_based_scatterplot(treated_samp_norm_bw, wildtype_samp_norm_bw, annotated_clusters,
                              genome_assembly="hg38", gene_bed=None, igr_bed=None, number_of_clusters=None,
                              out_prefix=None, formats=("html",), title_of_plot="", xaxis_label=None,
                              yaxis_label=None, min_x=None, max_x=None, min_y=None, max_y=None, pow=1.25,
                              hexbins=100, show_scales=True, include_additional_density_plot=False,
                              filter_extreme_bins=True, cluster_colors="bright", point_size=2,
                              point_opacity=0.12, alpha_range=(0.01, 3), max_points=120000, show_hulls=False, panel_size=460,
                              font_size=14, axis_title_size=22, tick_font_size=16, plot_range=None, seed=0,
                              write_table=False):
    """Port of density_based_scatterplot(): [density |] genic/intergenic | density-based clusters.

    ``annotated_clusters`` is the BED written by ``annotate-clust`` (chr, start,
    end, cluster letter). Defaults follow ChIPbinner (hexbins 100, pow 1.25,
    alpha range 0.01-1, 1 %/99 % trimming when ``filter_extreme_bins``).
    ``alpha_range`` defaults to (0.01, 3), clipped at 1 (ChIPbinner: (0.01, 1),
    which leaves most hexes near-black on the black panel).
    ``max_points`` subsamples the cluster scatter for display only (centroids,
    density and the bin table use every bin). ``cluster_colors``: "bright"
    (default, legible on black), "alphabet" (ChIPbinner's pals::alphabet) or a list.
    """
    gene_bed, igr_bed = _annotations(genome_assembly, gene_bed, igr_bed)
    xl, yl = _labels(treated_samp_norm_bw, wildtype_samp_norm_bw, xaxis_label, yaxis_label)
    pair = load_pair(treated_samp_norm_bw, wildtype_samp_norm_bw)
    pair["r"] = genic_class(pair, gene_bed, igr_bed)
    clusters = read_regions(annotated_clusters, name_col=3)
    clusters["name"] = clusters["name"].astype(str)
    letters = sorted(clusters["name"].unique())
    if number_of_clusters is not None:
        letters = letters[:int(number_of_clusters)]
        clusters = clusters[clusters["name"].isin(letters)]
    pair["cluster"] = assign_clusters(pair, clusters)
    pdat = quantile_filter(pair, 0.01, 0.99) if filter_extreme_bins else pair

    cells = hex_table(pdat, hexbins, pow, alpha_range)
    if plot_range is not None:
        min_x, max_x, min_y, max_y = resolve_range(plot_range)
    lim = _limits(cells, min_x, max_x, min_y, max_y)
    gx, gy, z, brks = _density_bands(pdat)

    palette = ALPHABET if cluster_colors == "alphabet" else BRIGHT if cluster_colors == "bright" else list(cluster_colors)
    cclr = dict(zip(letters, palette))

    kinds = (["density"] if include_additional_density_plot else []) + ["genic", "clusters"]
    p = _Panels(len(kinds), panel=panel_size)
    n_clu = pdat["cluster"].notna().sum()
    lay = _base_layout(p, title_of_plot or f"{yl} vs {xl}", None, font_size)
    fig = go.Figure(layout=lay)
    black = np.array([0.0, 0, 0])
    strip_titles = {"density": title_of_plot or "Bin density", "genic": "Intergenic vs genic ratio",
                    "clusters": "Density-based clusters"}

    for k, kind in enumerate(kinds):
        xa, ya = p.axis_names(k)
        dom = p.dom_x(k)
        ax_x = _axis(dom, ya, lim[0], f"{xl} {ARROW}", show_scales, "black", axis_title_size, tick_font_size)
        ax_y = _axis([0, 1], xa, lim[1], f"{yl} {ARROW}", show_scales, "black", axis_title_size, tick_font_size)
        fig.layout["xaxis" if k == 0 else f"xaxis{k + 1}"] = ax_x
        fig.layout["yaxis" if k == 0 else f"yaxis{k + 1}"] = ax_y
        # black panel background
        fig.add_shape(type="rect", xref=f"{xa} domain", yref=f"{ya} domain", x0=0, x1=1, y0=0, y1=1,
                      fillcolor="black", line=dict(width=0), layer="below")
        _strip(fig.layout, dom, strip_titles[kind], font_size, p.frac_y(p.strip))

        if kind == "density":
            cb = dict(orientation="h", x=dom[0] + p.frac_x(16), y=1 - p.frac_y(14), xanchor="left",
                      yanchor="top", len=p.frac_x(110), lenmode="fraction", thickness=9,
                      showticklabels=False, ticks="", outlinecolor="white", outlinewidth=1,
                      title=dict(text=f"# of bins {ARROW}", side="top",
                                 font=dict(color="white", size=font_size - 2)))
            fig.add_trace(_contour(gx, gy, z, brks, MAGMA, xa, ya, colorbar=cb))
            fig.add_trace(_diagonal(lim, xa, ya))
        elif kind == "genic":
            fig.add_trace(_diagonal(lim, xa, ya))
            fig.add_traces(_hex_traces(cells, black, xa, ya))
            w = dom[1] - dom[0]
            shapes, ann = _genic_legend_shapes((dom[0] + 0.09 * w, dom[0] + 0.33 * w), (0.70, 0.94),
                                               black, "white", font_size - 2)
            fig.layout.shapes += tuple(shapes)
            fig.layout.annotations += tuple(ann)
        else:
            haze = [[0, "rgba(255,255,255,0)"], [1, "rgba(255,255,255,0.92)"]]
            fig.add_trace(_contour(gx, gy, z, brks, haze, xa, ya))
            fig.add_trace(_diagonal(lim, xa, ya))
            rng = np.random.default_rng(seed)
            show = pdat
            if max_points and len(pdat) > max_points:
                show = pdat.iloc[np.sort(rng.choice(len(pdat), max_points, replace=False))]
            groups = [(None, "unclustered", NA_COLOR, point_opacity)] + [(c, c, cclr[c], point_opacity) for c in letters]
            for key, name, col, op in groups:
                sub = show[show["cluster"].isna()] if key is None else show[show["cluster"] == key]
                if sub.empty:
                    continue
                fig.add_trace(go.Scattergl(
                    x=sub["x"], y=sub["y"], mode="markers", name=name, xaxis=xa, yaxis=ya,
                    marker=dict(size=point_size, color=col, opacity=op, line=dict(width=0)),
                    customdata=np.c_[sub["chr"], sub["start"], sub["end"]],
                    hovertemplate=f"<b>{name}</b> %{{customdata[0]}}:%{{customdata[1]:,}}-%{{customdata[2]:,}}"
                                  f"<br>{xl} %{{x:.2f}}, {yl} %{{y:.2f}}<extra></extra>"))
            cen = pdat.dropna(subset=["cluster"]).groupby("cluster")[["x", "y"]].mean()
            for c in letters:
                if c not in cen.index:
                    continue
                if show_hulls:
                    from scipy.spatial import ConvexHull
                    pts = pdat.loc[pdat["cluster"] == c, ["x", "y"]].to_numpy()
                    hv = pts[ConvexHull(pts).vertices]
                    fig.add_trace(go.Scatter(x=np.r_[hv[:, 0], hv[:1, 0]], y=np.r_[hv[:, 1], hv[:1, 1]],
                                             mode="lines", line=dict(color=cclr[c], width=1.2),
                                             hoverinfo="skip", showlegend=False, xaxis=xa, yaxis=ya))
                fig.add_annotation(x=cen.loc[c, "x"], y=cen.loc[c, "y"], xref=xa, yref=ya, text=f"<b>{c}</b>",
                                   showarrow=False, bgcolor="white", bordercolor=cclr[c], borderwidth=2,
                                   borderpad=3, font=dict(color="black", size=font_size + 1),
                                   hovertext=f"cluster {c}: {(pdat['cluster'] == c).sum():,} bins")
            # small colour key for the clusters (top-left of the panel)
            w = dom[1] - dom[0]
            for i, c in enumerate(letters):
                fig.add_annotation(x=dom[0] + 0.04 * w, y=0.95 - i * 0.065, xref="paper", yref="paper",
                                   xanchor="left", showarrow=False,
                                   text=f"<span style='color:{cclr[c]}'>●</span> {c}  "
                                        f"<span style='color:#bbb'>{(pdat['cluster'] == c).sum():,}</span>",
                                   font=dict(color="white", size=font_size - 1))

    _save(fig, out_prefix, formats)
    if write_table and out_prefix:
        pdat.to_csv(f"{out_prefix}.bins.tsv.gz", sep="\t", index=False)
    return fig
