"""Genomic I/O and interval helpers (replaces rtracklayer / GenomicRanges).

All coordinates are 0-based half-open (BED convention). Two intervals
[s1, e1) and [s2, e2) overlap when s1 < e2 and s2 < e1, which is what
GenomicRanges::overlapsAny() computes on the equivalent 1-based ranges.
"""
from __future__ import annotations

import gzip
import os

import numpy as np
import pandas as pd

BIGWIG_EXT = (".bw", ".bigwig", ".bigWig")


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_bins(path: str, score_col: int = 3) -> pd.DataFrame:
    """Read binned scores as a DataFrame with columns chr, start, end, score.

    Accepts bigWig (one interval per bin, as written by ChIPbinner's norm_bw or
    by :func:`chipbinner.normalize.norm_bw`) or BED/bedGraph (optionally gzipped),
    taking the score from column ``score_col`` (0-based; ChIPbinner binned BEDs
    keep the raw count in column 4, i.e. index 3, the BED "name" field).
    """
    if str(path).endswith(BIGWIG_EXT):
        import pyBigWig

        bw = pyBigWig.open(str(path))
        parts = []
        for chrom in bw.chroms():
            iv = bw.intervals(chrom)
            if not iv:
                continue
            a = np.asarray(iv, dtype=float)
            parts.append(pd.DataFrame({"chr": chrom, "start": a[:, 0].astype(np.int64),
                                       "end": a[:, 1].astype(np.int64), "score": a[:, 2]}))
        bw.close()
        df = pd.concat(parts, ignore_index=True)
    else:
        df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, score_col],
                         dtype={0: str}, low_memory=False)
        df.columns = ["chr", "start", "end", "score"]
        df = df[~df["chr"].str.startswith(("track", "browser"))]
        df["start"] = df["start"].astype(np.int64)
        df["end"] = df["end"].astype(np.int64)
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
    return df.reset_index(drop=True)


def read_regions(path: str, name_col: int | None = None) -> pd.DataFrame:
    """Read a BED file of regions (chr, start, end[, name])."""
    names = ["chr", "start", "end"] + (["name"] if name_col is not None else [])
    if not path or not os.path.exists(path) or _is_empty(path):  # e.g. no blacklist
        return pd.DataFrame({c: pd.Series(dtype=str if c in ("chr", "name") else np.int64) for c in names})
    cols = [0, 1, 2] + ([name_col] if name_col is not None else [])
    df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=cols, dtype={0: str})
    df.columns = names
    return df


def _is_empty(path: str) -> bool:
    """True for a file with no data lines (0 bytes, or only blank / comment lines)."""
    with _open_text(path) as fh:
        return not any(line.strip() and not line.startswith("#") for line in fh)


def read_chrom_sizes(path: str) -> dict:
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["chr", "len"], dtype={0: str})
    return dict(zip(df["chr"], df["len"].astype(np.int64)))


class RegionSet:
    """Merged, per-chromosome sorted interval set with a vectorised overlapsAny()."""

    def __init__(self, regions: pd.DataFrame):
        self.idx = {}
        for chrom, sub in regions.groupby("chr", sort=False):
            s = sub["start"].to_numpy(np.int64)
            e = sub["end"].to_numpy(np.int64)
            o = np.argsort(s, kind="stable")
            s, e = s[o], e[o]
            # merge overlapping / touching intervals so ends become monotonic
            e_run = np.maximum.accumulate(e)
            new = np.ones(len(s), bool)
            new[1:] = s[1:] > e_run[:-1]
            gid = np.cumsum(new) - 1
            ms = s[new]
            me = np.zeros(len(ms), np.int64)
            np.maximum.at(me, gid, e)
            self.idx[chrom] = (ms, me)

    @classmethod
    def from_bed(cls, path: str) -> "RegionSet":
        return cls(read_regions(path))

    def overlaps_any(self, bins: pd.DataFrame) -> np.ndarray:
        out = np.zeros(len(bins), bool)
        chrs = bins["chr"].to_numpy()
        starts = bins["start"].to_numpy(np.int64)
        ends = bins["end"].to_numpy(np.int64)
        for chrom in pd.unique(chrs):
            if chrom not in self.idx:
                continue
            m = chrs == chrom
            ms, me = self.idx[chrom]
            j = np.searchsorted(ms, ends[m], side="left") - 1  # last region starting before bin end
            ok = j >= 0
            hit = np.zeros(m.sum(), bool)
            hit[ok] = me[j[ok]] > starts[m][ok]
            out[m] = hit
        return out


def match_bins(a: pd.DataFrame, b: pd.DataFrame) -> np.ndarray:
    """For every bin of ``a`` return the index (into ``b``) of the overlapping bin of ``b``.

    Equivalent to ``to(findOverlaps(a, b))`` for two tilings of the genome on the
    same grid (one hit per bin; a 1-bp coordinate shift between the two files, as
    produced by ChIPbinner's R objects, is tolerated). -1 means no overlap.
    """
    out = np.full(len(a), -1, np.int64)
    a_chr = a["chr"].to_numpy()
    b_chr = b["chr"].to_numpy()
    a_s = a["start"].to_numpy(np.int64)
    a_e = a["end"].to_numpy(np.int64)
    b_s_all = b["start"].to_numpy(np.int64)
    b_e_all = b["end"].to_numpy(np.int64)
    for chrom in pd.unique(a_chr):
        ma = np.flatnonzero(a_chr == chrom)
        mb = np.flatnonzero(b_chr == chrom)
        if len(mb) == 0:
            continue
        o = mb[np.argsort(b_s_all[mb], kind="stable")]
        bs, be = b_s_all[o], b_e_all[o]
        j = np.searchsorted(bs, a_e[ma], side="left") - 1
        ok = j >= 0
        ok[ok] = be[j[ok]] > a_s[ma][ok]
        out[ma[ok]] = o[j[ok]]
    return out


def write_bed(df: pd.DataFrame, path: str, cols=None):
    cols = cols or list(df.columns)
    df[cols].to_csv(path, sep="\t", header=False, index=False,
                    compression="gzip" if str(path).endswith(".gz") else None)


def write_bigwig(df: pd.DataFrame, path: str, chrom_sizes: dict | None = None):
    """Write chr/start/end/score bins to a bigWig (chromosomes in chrom_sizes order)."""
    import pyBigWig

    if chrom_sizes is None:
        chrom_sizes = df.groupby("chr", sort=False)["end"].max().to_dict()
    present = set(df["chr"])
    chroms = [c for c in chrom_sizes if c in present]
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(c, int(chrom_sizes[c])) for c in chroms])
    for c in chroms:
        sub = df[df["chr"] == c].sort_values("start")
        bw.addEntries([c] * len(sub), sub["start"].tolist(), ends=sub["end"].tolist(),
                      values=sub["score"].astype(float).tolist())
    bw.close()


def sample_label(path: str) -> str:
    """basename(tools::file_path_sans_ext(path)) as used by ChIPbinner for labels."""
    base = os.path.basename(str(path))
    return base.rsplit(".", 1)[0] if "." in base else base


def annotation_paths(genome_assembly: str, annot_dir: str | None = None) -> dict:
    """Paths of the genic / intergenic / blacklist / chrom sizes files shipped with the port."""
    d = (annot_dir or os.environ.get("CHIPBINNER_ANNOTATIONS")
         or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "annotations"))
    return {
        "gene": os.path.join(d, f"{genome_assembly}_genic.bed.gz"),
        "igr": os.path.join(d, f"{genome_assembly}_intergenic.bed.gz"),
        "blacklist": os.path.join(d, f"{genome_assembly}_blacklist.bed.gz"),
        "chrom_sizes": os.path.join(d, f"{genome_assembly}.chrom.sizes"),
    }
