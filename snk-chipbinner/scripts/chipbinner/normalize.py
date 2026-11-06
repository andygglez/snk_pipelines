"""Port of ChIPbinner::norm_bw (R/norm_bw.R)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .genomics import RegionSet, match_bins, read_bins, read_chrom_sizes


def norm_bw(immunoprecipitated_binned_file: str,
            chrom_sizes: str,
            blacklist: str,
            input_binned_file: str | None = None,
            use_input: bool = True,
            depth_norm: bool = True,
            pseudocount: float = 1e-3,
            raw_count_cutoff: float = 0,
            scaling_factor: float = 1.0) -> pd.DataFrame:
    """Normalise a binned ChIP count file, optionally against its input.

    Returns a DataFrame chr/start/end/score with one row per retained bin.

    Score (same formulas as ChIPbinner, including its log bases):
      depth_norm & input : log2((ChIP*s / (libChIP/1e6) + pc) / (input / (libInput/1e6) + pc))
      depth_norm, no input: ln(ChIP*s / libChIP * 1e6 + pc)
      no depth, input     : log2((ChIP*s + pc) / (input + pc))
      no depth, no input  : ln(ChIP*s + pc)
    Library sizes are summed over bins not overlapping the blacklist. Bins are
    kept when they lie on a chromosome of ``chrom_sizes`` (chrM excluded), do not
    overlap the blacklist and max(raw ChIP, raw input) > ``raw_count_cutoff``.

    Differences from the R code: ChIPbinner drops chrM from the bin template but
    not from the score vectors, which misaligns them when chrM bins are present;
    here chrM is dropped consistently.
    """
    bl = RegionSet.from_bed(blacklist)
    sizes = read_chrom_sizes(chrom_sizes)
    chip = read_bins(immunoprecipitated_binned_file)
    use_input = bool(use_input and input_binned_file)

    chip_lib = chip.loc[~bl.overlaps_any(chip), "score"].sum()
    if use_input:
        ctrl = read_bins(input_binned_file)
        ctrl_lib = ctrl.loc[~bl.overlaps_any(ctrl), "score"].sum()
        # subsetByOverlaps in both directions -> keep bins present in both files
        j = match_bins(chip, ctrl)
        chip = chip[j >= 0].reset_index(drop=True)
        ctrl_score = ctrl["score"].to_numpy(float)[j[j >= 0]]
    print(f"ChIP sample={immunoprecipitated_binned_file} depth={chip_lib:.0f}"
          + (f"; input={input_binned_file} depth={ctrl_lib:.0f}" if use_input else "; not using input"))

    keep = chip["chr"].isin(sizes.keys()).to_numpy() & (chip["chr"] != "chrM").to_numpy()
    raw = chip["score"].to_numpy(float)
    mx = np.maximum(raw, ctrl_score) if use_input else raw
    keep &= (mx > raw_count_cutoff) & ~bl.overlaps_any(chip)

    c = raw[keep] * scaling_factor
    if depth_norm:
        if use_input:
            score = np.log2((c / (chip_lib / 1e6) + pseudocount)
                            / (ctrl_score[keep] / (ctrl_lib / 1e6) + pseudocount))
        else:
            score = np.log(c / chip_lib * 1e6 + pseudocount)
    else:
        if use_input:
            score = np.log2((c + pseudocount) / (ctrl_score[keep] + pseudocount))
        else:
            score = np.log(c + pseudocount)

    out = chip.loc[keep, ["chr", "start", "end"]].copy()
    out["score"] = score
    return out[np.isfinite(out["score"])].reset_index(drop=True)


# --------------------------------------------------------------------------- shared helpers
def _common_bins(frames: list[pd.DataFrame]) -> list[pd.DataFrame]:
    """Restrict every frame to the bins present in all of them, rows aligned to frames[0].

    Equivalent to ChIPbinner's Reduce(subsetByOverlaps) over the samples
    (one bin grid shared by all files).
    """
    base = frames[0].reset_index(drop=True)
    idx = [np.arange(len(base))]
    keep = np.ones(len(base), bool)
    for f in frames[1:]:
        j = match_bins(base, f)
        keep &= j >= 0
        idx.append(j)
    out = []
    for f, j in zip(frames, idx):
        sub = f.iloc[j[keep]].reset_index(drop=True)
        sub[["chr", "start", "end"]] = base.loc[keep, ["chr", "start", "end"]].to_numpy()
        out.append(sub)
    return out


def _drop_blacklist(df: pd.DataFrame, bl: RegionSet) -> pd.DataFrame:
    return df[~bl.overlaps_any(df)].reset_index(drop=True)


# --------------------------------------------------------------------------- filter_low_counts
def filter_low_counts(sample_bedfiles: list[str], blacklist: str, cutoff: float = 100) -> dict:
    """Port of ChIPbinner::filter_low_counts.

    Drops blacklisted bins, keeps bins present in every sample, then keeps bins
    whose maximum raw count across the samples is > ``cutoff``. Returns
    {label: DataFrame(chr, start, end, score)}; label = file name without its
    last extension (ChIPbinner writes ``<label>_filt.bed``).
    R writes 1-based starts into the BED (as.data.frame of a GRanges); here the
    BED stays 0-based.
    """
    from .genomics import sample_label

    bl = RegionSet.from_bed(blacklist)
    frames = _common_bins([_drop_blacklist(read_bins(f), bl) for f in sample_bedfiles])
    mx = np.max(np.column_stack([f["score"].to_numpy(float) for f in frames]), axis=1)
    keep = mx > cutoff
    return {sample_label(p): f[keep].reset_index(drop=True) for p, f in zip(sample_bedfiles, frames)}


# --------------------------------------------------------------------------- merge_norm_bw
def merge_norm_bw(rep1: str, rep2: str) -> pd.DataFrame:
    """Port of ChIPbinner::merge_norm_bw: mean score of two replicates over their shared bins."""
    a, b = _common_bins([read_bins(rep1), read_bins(rep2)])
    out = a[["chr", "start", "end"]].copy()
    out["score"] = (a["score"].to_numpy(float) + b["score"].to_numpy(float)) / 2
    return out


# --------------------------------------------------------------------------- apply_normFactors
def deseq2_size_factors(counts: np.ndarray) -> np.ndarray:
    """DESeq2::estimateSizeFactors (median-of-ratios, rows with a zero are ignored)."""
    with np.errstate(divide="ignore"):
        logc = np.log(counts)
    loggeo = logc.mean(axis=1)
    ok = np.isfinite(loggeo)
    sf = np.empty(counts.shape[1])
    for j in range(counts.shape[1]):
        m = ok & (counts[:, j] > 0)
        sf[j] = np.exp(np.median(logc[m, j] - loggeo[m]))
    return sf


def _rank_avg(x):
    from scipy.stats import rankdata
    return rankdata(x, method="average")


def tmm_factors(counts: np.ndarray, logratio_trim=0.3, sum_trim=0.05, do_weighting=True,
                a_cutoff=-1e10) -> np.ndarray:
    """edgeR::calcNormFactors(method = "TMM") (edgeR 3.x reference-column rule)."""
    x = counts.astype(float)
    lib = x.sum(axis=0)
    f75 = np.quantile(x, 0.75, axis=0) / lib
    ref = int(np.argmax(np.sqrt(x).sum(axis=0))) if np.median(f75) < 1e-20 else int(np.argmin(np.abs(f75 - f75.mean())))
    out = np.ones(x.shape[1])
    r, nr = x[:, ref], lib[ref]
    for j in range(x.shape[1]):
        o, no = x[:, j], lib[j]
        with np.errstate(divide="ignore", invalid="ignore"):
            logr = np.log2((o / no) / (r / nr))
            abse = (np.log2(o / no) + np.log2(r / nr)) / 2
            v = (no - o) / no / o + (nr - r) / nr / r
        fin = np.isfinite(logr) & np.isfinite(abse) & (abse > a_cutoff)
        logr, abse, v = logr[fin], abse[fin], v[fin]
        if len(logr) == 0 or np.max(np.abs(logr)) < 1e-6:
            continue
        n = len(logr)
        lo_l = np.floor(n * logratio_trim) + 1
        hi_l = n + 1 - lo_l
        lo_s = np.floor(n * sum_trim) + 1
        hi_s = n + 1 - lo_s
        rl, ra = _rank_avg(logr), _rank_avg(abse)
        keep = (rl >= lo_l) & (rl <= hi_l) & (ra >= lo_s) & (ra <= hi_s)
        if do_weighting:
            f = np.nansum(logr[keep] / v[keep]) / np.nansum(1 / v[keep])
        else:
            f = np.nanmean(logr[keep])
        out[j] = 2 ** (0 if np.isnan(f) else f)
    return out / np.exp(np.mean(np.log(out)))


def apply_norm_factors(norm_method: str, treated_sample_bedfiles: list[str],
                       wildtype_sample_bedfiles: list[str], blacklist: str) -> tuple[dict, pd.DataFrame]:
    """Port of ChIPbinner::apply_normFactors.

    norm_method: "DESeq2" (counts / median-of-ratios size factor), "edgeR"
    (TMM-normalised CPM) or "depth-norm-only" (CPM). Returns
    ({label: DataFrame(chr, start, end, score)}, factors table); ChIPbinner
    writes ``<label>_sizeFactorNorm.bed``.
    """
    from .genomics import sample_label

    bl = RegionSet.from_bed(blacklist)
    files = list(wildtype_sample_bedfiles) + list(treated_sample_bedfiles)
    frames = _common_bins([_drop_blacklist(read_bins(f), bl) for f in files])
    counts = np.column_stack([f["score"].to_numpy(float) for f in frames])
    labels = [sample_label(f) for f in files]
    if norm_method == "DESeq2":
        sf = deseq2_size_factors(counts)
        norm = counts / sf
        table = pd.DataFrame({"sample": labels, "sizeFactor": sf})
    elif norm_method in ("edgeR", "depth-norm-only"):
        nf = tmm_factors(counts)
        lib = counts.sum(axis=0)
        eff = lib * (nf if norm_method == "edgeR" else 1.0)
        norm = counts / eff * 1e6
        table = pd.DataFrame({"sample": labels, "lib.size": lib, "norm.factors": nf})
    else:
        raise ValueError("norm_method must be DESeq2, edgeR or depth-norm-only")
    out = {}
    for k, lab in enumerate(labels):
        df = frames[0][["chr", "start", "end"]].copy()
        df["score"] = norm[:, k]
        out[lab] = df
    return out, table
