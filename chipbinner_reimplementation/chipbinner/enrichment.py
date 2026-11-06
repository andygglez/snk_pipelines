"""Port of ChIPbinner::enrich_clust (LOLA::runLOLA + qvalue::qvalue) with a Plotly lollipop plot.

Region databases use the LOLA directory layout, so LOLA Core / ChIPbinner
databases work unchanged::

    <db_dir>/<collection>/regions/<set>.bed[.gz]
    <db_dir>/<collection>/index.txt     (optional, tab-separated, column "filename")

runLOLA, per region set: support = query bins overlapping the set,
b = universe bins overlapping the set - support, c = |query| - support,
d = |universe| - support - b - c; one-sided Fisher exact test ("greater" for
enrichment, "less" for depletion) with the conditional-MLE odds ratio (as R's
fisher.test); qValue = qvalue::qvalue(p) (Storey, smoother pi0).
"""
from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd

from .genomics import RegionSet, read_regions

SIG_CUTS = [(1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*")]


def load_region_db(db_dir: str, collections=None) -> list[dict]:
    sets = []
    for coll in sorted(os.listdir(db_dir)):
        cdir = os.path.join(db_dir, coll, "regions")
        if not os.path.isdir(cdir) or (collections and coll not in collections):
            continue
        for f in sorted(glob.glob(os.path.join(cdir, "*.bed")) + glob.glob(os.path.join(cdir, "*.bed.gz"))):
            sets.append({"collection": coll, "filename": os.path.basename(f), "regions": read_regions(f)})
    if not sets:
        raise ValueError(f"no region sets found under {db_dir}/<collection>/regions/")
    return sets


# --------------------------------------------------------------------------- qvalue
def _smooth_spline_df(x, y, df=3.0):
    """Cubic smoothing spline with ~df effective degrees of freedom (R smooth.spline(x, y, df))."""
    from scipy.interpolate import make_smoothing_spline

    def eff_df(lam):
        eye = np.eye(len(x))
        return sum(make_smoothing_spline(x, eye[i], lam=lam)(x[i]) for i in range(len(x)))

    lo, hi = -12.0, 6.0  # log10 lambda bisection: df decreases with lambda
    for _ in range(50):
        mid = (lo + hi) / 2
        if eff_df(10 ** mid) > df:
            lo = mid
        else:
            hi = mid
    return make_smoothing_spline(x, y, lam=10 ** ((lo + hi) / 2))


def qvalue(p: np.ndarray, lambdas=np.arange(0.05, 0.951, 0.05)) -> np.ndarray:
    """qvalue::qvalue(p)$qvalues with pi0.method = "smoother" (smooth.log.pi0 = FALSE)."""
    p = np.asarray(p, float)
    m = len(p)
    if m == 0:
        return p
    if p.max() < lambdas.max():
        # qvalue::pi0est stops here; LOLA catches it and reports no qValue
        raise ValueError("maximum p-value is smaller than lambda range")
    pi0s = np.array([np.mean(p >= l) / (1 - l) for l in lambdas])
    pi0 = float(min(_smooth_spline_df(lambdas, pi0s, df=3)(lambdas)[-1], 1.0))
    if pi0 <= 0:
        raise ValueError("pi0 <= 0 in qvalue")
    o = np.argsort(-p, kind="stable")  # decreasing p
    i = np.arange(m, 0, -1)
    qs = pi0 * np.minimum(1.0, np.minimum.accumulate(p[o] * m / i))
    out = np.empty(m)
    out[o] = qs
    return out


# --------------------------------------------------------------------------- runLOLA
def run_lola(query: pd.DataFrame, universe: pd.DataFrame, db: list[dict], direction="enrichment") -> pd.DataFrame:
    from scipy.stats import fisher_exact
    from scipy.stats.contingency import odds_ratio

    alt = "less" if direction == "depletion" else "greater"
    rows = []
    for s in db:
        rs = RegionSet(s["regions"])
        support = int(rs.overlaps_any(query).sum())
        b = int(rs.overlaps_any(universe).sum()) - support
        c = len(query) - support
        d = len(universe) - support - b - c
        if b < 0 or c < 0 or d < 0:
            raise ValueError("query bins must be a subset of the universe")
        tab = np.array([[support, b], [c, d]])
        pv = fisher_exact(tab, alternative=alt).pvalue
        orr = odds_ratio(tab, kind="conditional").statistic
        rows.append(dict(collection=s["collection"], filename=s["filename"], pValue=pv, oddsRatio=orr,
                         support=support, b=b, c=c, d=d, size=len(s["regions"])))
    t = pd.DataFrame(rows)
    try:
        t["qValue"] = qvalue(t["pValue"].to_numpy())
    except ValueError:
        pass  # LOLA warns and leaves qValue out
    t["pValueLog"] = -np.log10(t["pValue"] + 1e-322)
    t["rnkSup"] = t["support"].rank(ascending=False, method="min")
    t["rnkPV"] = t["pValueLog"].rank(ascending=False, method="min")
    t["rnkOR"] = t["oddsRatio"].rank(ascending=False, method="min")
    t["maxRnk"] = t[["rnkSup", "rnkPV", "rnkOR"]].max(axis=1)
    t["meanRnk"] = t[["rnkSup", "rnkPV", "rnkOR"]].mean(axis=1).round(3)
    return t.sort_values(["pValueLog", "meanRnk"], ascending=[False, True]).reset_index(drop=True)


def _signif(q):
    for cut, sym in SIG_CUTS:
        if q <= cut:
            return sym
    return "ns"


def _summarise(t: pd.DataFrame, kind: str, exclude: str | None):
    if exclude:
        t = t[t["filename"] != exclude]
    if t.empty or "qValue" not in t:
        return None
    t = t.copy()
    t["reg"] = t["filename"].str.replace(".bed.gz", "", regex=False).str.replace(".bed", "", regex=False)
    t["sig"] = [_signif(q) for q in t["qValue"]]
    pos = t.loc[t["oddsRatio"] > 0, "oddsRatio"]
    t.loc[~(t["oddsRatio"] > 0), "oddsRatio"] = pos.min() if len(pos) else 0
    out = t.assign(overlap=t["support"])[["oddsRatio", "reg", "qValue", "overlap", "sig"]]
    out["type"] = kind
    return None if (out["sig"] != "ns").sum() == 0 else out


def enrich_clust(annotated_clusters: str, query_cluster: str, pooled_bed_file: str, db_dir: str,
                 gene_bed: str, igr_bed: str, region: str = "genome_wide", collections=None):
    """Port of ChIPbinner::enrich_clust. Returns the aggregated table (or None).

    region: genome_wide (all pooled bins as background), genic or intergenic
    (query and background restricted to genic-only / intergenic-only bins; the
    gene.bed.gz / intergenic.bed.gz set itself is then dropped, as in R).
    """
    clus = read_regions(annotated_clusters, name_col=3)
    q = clus[clus["name"].astype(str) == str(query_cluster)][["chr", "start", "end"]].reset_index(drop=True)
    uni = read_regions(pooled_bed_file)
    db = load_region_db(db_dir, collections)
    g, i = RegionSet.from_bed(gene_bed), RegionSet.from_bed(igr_bed)
    exclude = None
    if region == "genic":
        uni = uni[g.overlaps_any(uni) & ~i.overlaps_any(uni)]
        q = q[g.overlaps_any(q) & ~i.overlaps_any(q)]
        exclude = "gene.bed.gz"
    elif region == "intergenic":
        uni = uni[i.overlaps_any(uni) & ~g.overlaps_any(uni)]
        q = q[i.overlaps_any(q) & ~g.overlaps_any(q)]
        exclude = "intergenic.bed.gz"
    elif region != "genome_wide":
        raise ValueError("region must be genome_wide, genic or intergenic")
    enr = _summarise(run_lola(q, uni, db, "enrichment"), "enrichment", exclude)
    dep = _summarise(run_lola(q, uni, db, "depletion"), "depletion", exclude)
    parts = [x for x in (enr, dep) if x is not None]
    return pd.concat(parts, ignore_index=True) if parts else None


def plot_enrichment(agg: pd.DataFrame, n_elements=10, cutoff_for_overlap=0, title=""):
    """Lollipop plot of the top significant sets per direction (ChIPbinner's ggplot)."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    kinds = [k for k in ("enrichment", "depletion") if k in set(agg["type"])]
    fig = make_subplots(rows=1, cols=len(kinds), subplot_titles=[k.capitalize() for k in kinds],
                        horizontal_spacing=0.28)
    for j, k in enumerate(kinds, start=1):
        sub = (agg[(agg["type"] == k) & (agg["overlap"] > cutoff_for_overlap) & (agg["sig"] != "ns")]
               .sort_values("oddsRatio", ascending=False).head(n_elements))
        sizes = 14 + 18 * (sub["overlap"] - sub["overlap"].min()) / max(1, np.ptp(sub["overlap"]))
        for _, r in sub.iterrows():
            fig.add_trace(go.Scatter(x=[0, r["oddsRatio"]], y=[r["reg"]] * 2, mode="lines",
                                     line=dict(color="#bbb", width=2), hoverinfo="skip", showlegend=False),
                          row=1, col=j)
        fig.add_trace(go.Scatter(
            x=sub["oddsRatio"], y=sub["reg"], mode="markers+text", text=sub["sig"], textposition="top center",
            marker=dict(size=sizes, color="#ff8c00", line=dict(color="white", width=1)), showlegend=False,
            customdata=np.c_[sub["qValue"], sub["overlap"]],
            hovertemplate="%{y}<br>odds ratio %{x:.3g}<br>FDR %{customdata[0]:.3g}<br>"
                          "overlap %{customdata[1]}<extra></extra>"), row=1, col=j)
        fig.update_yaxes(autorange="reversed", showgrid=True, gridcolor="#ddd", row=1, col=j)
        fig.update_xaxes(title_text="Odds ratio", showgrid=True, gridcolor="#eee", griddash="dash",
                         rangemode="tozero", row=1, col=j)
    fig.update_layout(template="simple_white", width=420 * len(kinds) + 200,
                      height=max(360, 60 + 34 * min(n_elements, len(agg))),
                      title=dict(text=f"<b>{title}</b>", x=0.5), font=dict(size=13))
    return fig
