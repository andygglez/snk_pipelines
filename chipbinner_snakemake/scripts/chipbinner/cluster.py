"""Ports of ChIPbinner::pre_clust, generate_clust (inst/python/clus.py) and annotate_clust."""
from __future__ import annotations

import string

import numpy as np
import pandas as pd

from .genomics import RegionSet, match_bins, read_bins


def paired_bins(treated: pd.DataFrame, wildtype: pd.DataFrame) -> pd.DataFrame:
    """Bins scored in both samples, as built at the top of every ChIPbinner scatterplot.

    R: r <- Reduce(overlap, list(treated[score != 0], wildtype[score != 0])), i.e.
    the non-zero treated bins that overlap a non-zero wild-type bin (treated
    coordinates). Returns chr/start/end plus ``x`` (wild-type score) and ``y``
    (treated score), the orientation of every ChIPbinner scatterplot.
    """
    t = treated[treated["score"] != 0].reset_index(drop=True)
    w = wildtype[wildtype["score"] != 0].reset_index(drop=True)
    j = match_bins(t, w)
    ok = j >= 0
    out = t.loc[ok, ["chr", "start", "end"]].reset_index(drop=True)
    out["x"] = w["score"].to_numpy(float)[j[ok]]
    out["y"] = t.loc[ok, "score"].to_numpy(float)
    return out


def load_pair(treated_samp_norm_bw: str, wildtype_samp_norm_bw: str) -> pd.DataFrame:
    return paired_bins(read_bins(treated_samp_norm_bw), read_bins(wildtype_samp_norm_bw))


def genic_class(bins: pd.DataFrame, gene_bed: str, igr_bed: str) -> np.ndarray:
    """+1 genic only, -1 intergenic only, 0 both/neither (ChIPbinner's ``olap``)."""
    g = RegionSet.from_bed(gene_bed).overlaps_any(bins)
    i = RegionSet.from_bed(igr_bed).overlaps_any(bins)
    return np.where(g & ~i, 1, np.where(~g & i, -1, 0)).astype(np.int8)


def pre_clust(treated_samp_norm_bw: str, wildtype_samp_norm_bw: str, out_prefix: str):
    """Write ``<prefix>_mat.csv`` (treated, wildtype; no header) and ``<prefix>_pooled.bed``.

    Same files as ChIPbinner::pre_clust, so they can be fed to R's annotate_clust
    too. (R's base-graphics smoothScatter is drawn by extra_plots.smooth_scatter.)
    """
    pair = load_pair(treated_samp_norm_bw, wildtype_samp_norm_bw)
    pair[["y", "x"]].to_csv(f"{out_prefix}_mat.csv", header=False, index=False, float_format="%.15g")
    pair[["chr", "start", "end"]].to_csv(f"{out_prefix}_pooled.bed", sep="\t", header=False, index=False)
    return pair


def generate_clust(matrix_file: str, minpts: int, minsamps: int, cores: int = 1,
                   out_file: str | None = None) -> np.ndarray:
    """HDBSCAN(min_cluster_size=minpts, min_samples=minsamps) on the pre_clust matrix.

    Uses the ``hdbscan`` package (what ChIPbinner installs) when importable,
    otherwise scikit-learn's HDBSCAN (same algorithm, EOM cluster selection).
    Labels: -1 = noise.
    """
    d = np.loadtxt(matrix_file, delimiter=",", ndmin=2)
    try:
        import hdbscan

        labels = hdbscan.HDBSCAN(min_cluster_size=minpts, min_samples=minsamps,
                                 core_dist_n_jobs=cores).fit(d).labels_
    except ImportError:
        from sklearn.cluster import HDBSCAN

        labels = HDBSCAN(min_cluster_size=minpts, min_samples=minsamps, n_jobs=cores).fit(d).labels_
    labels = labels.astype(int)
    if out_file:
        np.savetxt(out_file, labels, fmt="%i")
    return labels


def annotate_clust(matrix_file: str, pooled_bed_file: str, hdbscan_output_file: str,
                   number_of_clusters: int | None = None) -> pd.DataFrame:
    """Name HDBSCAN clusters A, B, C ... by increasing mean signal 0.5*(treated + wildtype).

    Returns the bins of the first ``number_of_clusters`` clusters (all clusters
    when None) as chr/start/end/name. ChIPbinner accepts 2..25 clusters. As in R,
    when HDBSCAN found more clusters than requested only the lowest-signal ones
    are kept (on a full genome with small minpts these can be a few hundred bins).
    """
    mat = pd.read_csv(matrix_file, header=None, names=["X1", "X2"])
    bed = pd.read_csv(pooled_bed_file, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chr", "start", "end"], dtype={0: str})
    mat["clu"] = np.loadtxt(hdbscan_output_file, dtype=int, ndmin=1)
    ctr = (mat[mat["clu"] != -1].groupby("clu")[["X1", "X2"]].mean()
           .assign(mu=lambda d: 0.5 * (d["X1"] + d["X2"])).sort_values("mu", kind="stable"))
    n = len(ctr) if number_of_clusters is None else int(number_of_clusters)
    if not 1 < n < 26 or n > len(ctr):
        raise ValueError(f"number_of_clusters must be 2..25 and <= clusters found ({len(ctr)}); got {n}")
    parts = []
    for letter, clu in zip(string.ascii_uppercase, ctr.index[:n]):
        sub = bed[(mat["clu"] == clu).to_numpy()].copy()
        sub["name"] = letter
        parts.append(sub)
    return pd.concat(parts, ignore_index=True)


def assign_clusters(bins: pd.DataFrame, clusters: pd.DataFrame) -> np.ndarray:
    """Cluster letter of each bin (first matching cluster in A, B, C order; None otherwise)."""
    lab = np.full(len(bins), None, dtype=object)
    for letter in sorted(clusters["name"].unique()):
        hit = RegionSet(clusters[clusters["name"] == letter]).overlaps_any(bins)
        lab[(lab == None) & hit] = letter  # noqa: E711 (elementwise)
    return lab
