"""Port of ChIPbinner::differentialBinAnalysis and of ROTS (Bioconductor ROTS::ROTS, two groups, unpaired).

ROTS (reproducibility-optimised test statistic, Suomi et al. 2017):
  * 2B stratified bootstrap resamplings (columns resampled within each group)
    and 2B label permutations;
  * for each a1 in ssq = (0..20)/100, (11..50)/50, (6..25)/5 the statistic
    d / (a1 + s) (and the plain d, "slr"), top-list reproducibility between
    bootstrap pairs (b, b + B) for every top-list size N < K; the same for the
    permutation pairs;
  * pick (a1, a2, k) maximising Z = (R_boot - R_perm) / sd(R_boot);
  * p-values and FDR from the permutation null of the selected statistic.
This is the ROTS algorithm, re-implemented; resampling uses numpy's RNG, so
results agree with R statistically, not draw for draw.

Differences in implementation (not in the method): the bootstrap statistics are
streamed (no nbins x 2B matrices), and the top-list overlap is computed for all
N at once from max(rank1, rank2) instead of N separate scans.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .genomics import RegionSet, read_bins, read_regions, sample_label
from .normalize import _common_bins, _drop_blacklist

SSQ = np.r_[np.arange(0, 21) / 100, np.arange(11, 51) / 50, np.arange(6, 26) / 5]
N_GRID = np.r_[np.arange(1, 21) * 5, np.arange(11, 51) * 10, np.arange(21, 41) * 25, np.arange(11, 1001) * 100]


def _test_statistic(X, Y):
    """ROTS testStatistic(paired = FALSE) for two groups: d = mean(Y) - mean(X), pooled SE s."""
    mX, mY = X.mean(axis=1), Y.mean(axis=1)
    sX = ((X - mX[:, None]) ** 2).sum(axis=1)
    sY = ((Y - mY[:, None]) ** 2).sum(axis=1)
    nX, nY = X.shape[1], Y.shape[1]
    d = mY - mX
    s = np.sqrt(((sX + sY) / (nX + nY - 2)) * (1 / nX + 1 / nY))
    return d, s


def _overlap_curve(r1, r2, N):
    """NeedForSpeed calculateOverlap: for each N, share of the top-N by r1 whose r2 is in r2's top N.

    Ties follow the C++ code: r1 ties broken by r2 (both descending); "in the top N"
    means r2 >= the N-th largest r2, i.e. fewer than N values strictly greater.
    """
    order = np.lexsort((r2, r1))[::-1]
    r2o = r2[order]
    asc = np.sort(r2)
    n_greater = len(r2) - np.searchsorted(asc, r2o, side="right")
    m = np.maximum(np.arange(len(r2)), n_greater)
    return np.searchsorted(np.sort(m), N, side="left") / N


def rots(data: np.ndarray, groups, B: int = 1000, K: int | None = None, seed: int = 123,
         progress: bool = False):
    """ROTS(data, groups, B, K, log = TRUE) for two unpaired groups.

    ``groups`` are labels; like ROTS (character groups -> factor levels), group 1
    is the alphabetically first label. Returns a dict with d, logfc, pvalue, FDR,
    a1, a2, k, R, Z, ztable and the ordered labels.
    """
    rng = np.random.default_rng(seed)
    groups = np.asarray(groups).astype(str)
    levels = sorted(set(groups))
    if len(levels) != 2:
        raise ValueError("ROTS port supports exactly two groups")
    order = np.argsort([levels.index(g) for g in groups], kind="stable")
    data = np.asarray(data, float)[:, order]
    cl = np.array([levels.index(g) + 1 for g in groups[order]])
    g1, g2 = np.flatnonzero(cl == 1), np.flatnonzero(cl == 2)
    if len(g1) < 2 or len(g2) < 2:
        raise ValueError("At least two replicates per condition are needed.")
    n = data.shape[0]
    if K is None:
        K = n // 4
    K = min(K, n)
    N = N_GRID[N_GRID < K]
    if len(N) == 0:
        raise ValueError(f"K={K} too small (needs > 5)")
    logfc = data[:, g1].mean(axis=1) - data[:, g2].mean(axis=1)

    def boot_cols():
        cols = np.empty(len(cl), int)
        cols[g1] = rng.choice(g1, len(g1), replace=True)
        cols[g2] = rng.choice(g2, len(g2), replace=True)
        return cols

    boots = [boot_cols() for _ in range(2 * B)]
    perms = [rng.permutation(len(cl)) for _ in range(2 * B)]

    def stat(cols):
        return _test_statistic(data[:, cols[g1]], data[:, cols[g2]])

    nrow = len(SSQ) + 1
    ov = np.zeros((B, nrow, len(N)))
    ovp = np.zeros((B, nrow, len(N)))
    np_err = np.seterr(divide="ignore", invalid="ignore")  # s = 0 -> Inf, as in R
    for b in range(B):
        D1, S1 = stat(boots[b])
        D2, S2 = stat(boots[b + B])
        P1, PS1 = stat(perms[b])
        P2, PS2 = stat(perms[b + B])
        for i, a in enumerate(SSQ):
            ov[b, i] = _overlap_curve(np.abs(D1 / (S1 + a)), np.abs(D2 / (S2 + a)), N)
            ovp[b, i] = _overlap_curve(np.abs(P1 / (PS1 + a)), np.abs(P2 / (PS2 + a)), N)
        ov[b, -1] = _overlap_curve(np.abs(D1), np.abs(D2), N)
        ovp[b, -1] = _overlap_curve(np.abs(P1), np.abs(P2), N)
        if progress and (b + 1) % max(1, B // 20) == 0:
            print(f"  ROTS bootstrap {b + 1}/{B}", flush=True)
    np.seterr(**np_err)
    repro = ov.mean(axis=0)
    repro_p = ovp.mean(axis=0)
    repro_sd = ov.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (repro - repro_p) / repro_sd
    fin = np.isfinite(z)
    zmax = z[fin].max()
    # which(arr.ind = TRUE) scans column-major; take the first hit
    hits = [(r, c) for c in range(z.shape[1]) for r in range(z.shape[0]) if fin[r, c] and z[r, c] == zmax]
    r_sel, c_sel = hits[0]
    a1, a2 = (float(SSQ[r_sel]), 1.0) if r_sel < len(SSQ) else (1.0, 0.0)
    k = int(N[c_sel])

    d0, s0 = _test_statistic(data[:, g1], data[:, g2])
    d = d0 / (a1 + a2 * s0)
    absd = np.abs(d)
    ordd = np.argsort(-absd, kind="stable")
    a_sorted = absd[ordd]
    # p-value: share of all permuted |d| >= observed |d| (pooled over the 2B permutations)
    count_ge = np.zeros(n)
    A = np.empty((n, 2 * B), dtype=np.float32)
    ranks = np.arange(1, n + 1)
    a_asc = np.sort(a_sorted)
    for i, p in enumerate(perms):
        pd_, ps_ = stat(p)
        pr = np.sort(np.abs(pd_ / (a1 + a2 * ps_)))
        count_ge += len(pr) - np.searchsorted(pr, a_sorted, side="left")
        # biggerN(a, a.rand): #(a.rand > a) + (a %in% a.rand)
        greater = len(pr) - np.searchsorted(pr, a_sorted, side="right")
        pos = np.searchsorted(pr, a_sorted, side="left")
        isin = (pos < len(pr)) & (pr[np.minimum(pos, len(pr) - 1)] == a_sorted)
        A[:, i] = (greater + isin) / ranks
    pval = np.empty(n)
    pval[ordd] = count_ge / (n * 2 * B)
    fdr_sorted = np.minimum(np.median(A, axis=1), 1.0)
    fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
    fdr = np.empty(n)
    fdr[ordd] = fdr_sorted
    del a_asc
    return dict(d=d, logfc=logfc, pvalue=pval, FDR=fdr, a1=a1, a2=a2, k=k, R=float(repro[r_sel, c_sel]),
                Z=float(zmax), ztable=pd.DataFrame(z, index=[*map(str, SSQ), "slr"], columns=N),
                levels=levels, B=B)


def rots_summary(res, coords, fdr=0.05) -> str:
    """Text like ROTS summary(rots.out, fdr = 0.05) (first lines, as ChIPbinner saves)."""
    o = np.argsort(-np.abs(res["d"]), kind="stable")
    sel = [i for i in o if res["FDR"][i] <= fdr]
    lines = ["ROTS results: ", "", f"Number of resamplings:  {res['B']} ", "",
             f"a1:                     {res['a1']} ", f"a2:                     {res['a2']} ",
             f"Top list size:          {res['k']} ", f"Reproducibility value:  {res['R']:.6g} ",
             f"Z-score:                {res['Z']:.6g} ", "", f"{len(sel)} rows satisfy the condition."]
    if sel:
        lines.append("Row\tROTS-statistic\tpvalue\tFDR\tbin")
        for i in sel[:10]:
            lines.append(f"{i + 1}\t{res['d'][i]:.6g}\t{res['pvalue'][i]:.6g}\t{res['FDR'][i]:.6g}\t{coords[i]}")
    return "\n".join(lines) + "\n"


def differential_bin_analysis(treated_sample_bigWigFiles, wildtype_sample_bigWigFiles, blacklist,
                              treated_condition_label, wildtype_condition_label, annotated_clusters=None,
                              bootstrap_value=1000, K_value=None, seed=123, progress=True):
    """Port of ChIPbinner::differentialBinAnalysis.

    Returns (all_bins, per_cluster, summary_text). Columns follow ChIPbinner:
    all bins -> genomic_coord, chr, start, end, logFC, FDR; per cluster ->
    genomic_coord, seqnames, start, end, strand, logFC, FDR, cluster. ``start``
    is 1-based (GRanges convention) as in the R output, and genomic_coord is
    "chr.start.end". Extra columns: pvalue, ROTS_statistic and
    log2FC_treated_vs_wildtype (treated - wildtype mean; ROTS's logFC is
    mean(first group alphabetically) - mean(second)).
    """
    bl = RegionSet.from_bed(blacklist)
    files = list(treated_sample_bigWigFiles) + list(wildtype_sample_bigWigFiles)
    frames = _common_bins([_drop_blacklist(read_bins(f), bl) for f in files])
    mat = np.column_stack([f["score"].to_numpy(float) for f in frames])
    labels = [sample_label(f) for f in files]
    # ChIPbinner orders the columns wildtype first, then treated
    nt = len(treated_sample_bigWigFiles)
    cols = list(range(nt, len(files))) + list(range(nt))
    mat = mat[:, cols]
    groups = [wildtype_condition_label] * (len(files) - nt) + [treated_condition_label] * nt
    ok = np.all(np.isfinite(mat), axis=1)
    bins = frames[0].loc[ok, ["chr", "start", "end"]].reset_index(drop=True)
    mat = mat[ok]
    coords = (bins["chr"] + "." + (bins["start"] + 1).astype(str) + "." + bins["end"].astype(str)).to_numpy()
    res = rots(mat, groups, B=bootstrap_value, K=K_value, seed=seed, progress=progress)
    nw = len(files) - nt
    all_bins = pd.DataFrame({
        "genomic_coord": coords, "chr": bins["chr"], "start": bins["start"] + 1, "end": bins["end"],
        "logFC": res["logfc"], "FDR": res["FDR"], "pvalue": res["pvalue"], "ROTS_statistic": res["d"],
        "log2FC_treated_vs_wildtype": mat[:, nw:].mean(axis=1) - mat[:, :nw].mean(axis=1)})
    summary = rots_summary(res, coords)
    per_cluster = None
    if annotated_clusters:
        clus = read_regions(annotated_clusters, name_col=3)
        parts = []
        for letter in sorted(clus["name"].astype(str).unique()):
            hit = RegionSet(clus[clus["name"].astype(str) == letter]).overlaps_any(bins)
            sub = all_bins[hit].rename(columns={"chr": "seqnames"}).copy()
            sub.insert(4, "strand", "*")
            sub["cluster"] = letter
            parts.append(sub)
        per_cluster = pd.concat(parts, ignore_index=True)[
            ["genomic_coord", "seqnames", "start", "end", "strand", "logFC", "FDR", "cluster",
             "pvalue", "ROTS_statistic", "log2FC_treated_vs_wildtype"]]
    return all_bins, per_cluster, summary, labels
