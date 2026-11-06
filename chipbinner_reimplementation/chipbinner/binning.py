"""BAM -> binned counts, as in the snakemake_chipbinner ``bam2bed`` rule.

That rule runs ``bedtools bamtobed -i x.bam | sort`` and then
``bedtools intersect -a windows.bed -b reads.bed -sorted -c``: every mapped
alignment record (each mate, secondary/supplementary and duplicates included -
bamtobed does not filter them) counts once in every window it overlaps by >= 1 bp.
The alignment span is the reference span of the CIGAR (M/D/N/=/X), like bamtobed.
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

from .genomics import read_chrom_sizes


def make_windows(chrom_sizes: str, bin_size: int, chroms=None) -> pd.DataFrame:
    """bedtools makewindows -g chrom.sizes -w bin_size (last window truncated)."""
    sizes = read_chrom_sizes(chrom_sizes)
    parts = []
    for c, n in sizes.items():
        if chroms is not None and c not in chroms:
            continue
        s = np.arange(0, n, bin_size, dtype=np.int64)
        parts.append(pd.DataFrame({"chr": c, "start": s, "end": np.minimum(s + bin_size, n)}))
    return pd.concat(parts, ignore_index=True)


def read_windows(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"],
                     dtype={0: str}, comment="#")
    return df


def bam_to_bins(bam: str, windows: pd.DataFrame, read_bed: str | None = None,
                skip_flags: int = 0x4, min_mapq: int = 0) -> pd.DataFrame:
    """Count alignments overlapping each window; returns chr/start/end/count (window order).

    ``skip_flags`` defaults to unmapped only (bamtobed behaviour). ``read_bed``
    optionally writes the bamtobed-style read BED (chr, start, end, name, mapq, strand).
    """
    import pysam

    out = np.zeros(len(windows), dtype=np.int64)
    wchr = windows["chr"].to_numpy()
    fh = None
    if read_bed:
        fh = gzip.open(read_bed, "wt") if read_bed.endswith(".gz") else open(read_bed, "w")
    with pysam.AlignmentFile(bam, "rb") as bf:
        per_chrom = {}
        for rec in bf.fetch(until_eof=True):
            if rec.flag & skip_flags or rec.is_unmapped or rec.mapping_quality < min_mapq:
                continue
            c = rec.reference_name
            st, en = rec.reference_start, rec.reference_end
            lst = per_chrom.setdefault(c, ([], []))
            lst[0].append(st)
            lst[1].append(en)
            if fh is not None:
                name = rec.query_name + ("/1" if rec.is_read1 else "/2" if rec.is_read2 else "")
                fh.write(f"{c}\t{st}\t{en}\t{name}\t{rec.mapping_quality}\t{'-' if rec.is_reverse else '+'}\n")
    if fh is not None:
        fh.close()
    for c, (st, en) in per_chrom.items():
        m = np.flatnonzero(wchr == c)
        if len(m) == 0:
            continue
        st = np.sort(np.asarray(st, np.int64))
        en = np.sort(np.asarray(en, np.int64))
        ws = windows["start"].to_numpy(np.int64)[m]
        we = windows["end"].to_numpy(np.int64)[m]
        # reads with start < window end, minus reads that end at/before window start
        out[m] = np.searchsorted(st, we, side="left") - np.searchsorted(en, ws, side="right")
    res = windows[["chr", "start", "end"]].copy()
    res["count"] = out
    return res
