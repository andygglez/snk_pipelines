#!/usr/bin/env python3
"""One-off converter: ChIPbinner's R data objects (.rda) -> plain BED / bigWig files.

ChIPbinner ships its genome annotations (genic, intergenic, blacklist, chrom sizes)
and its example normalised bigWigs / annotated clusters as R objects. This script
parses them with the pure-Python `rdata` package (no R code is executed) so the
Python port never needs R.

GRanges coordinates are 1-based closed; they are written as 0-based half-open BED
(start - 1, end), exactly what rtracklayer::export.bed would write.

Usage:
    rda_to_bed.py sysdata  <sysdata.rda> <outdir>
    rda_to_bed.py granges  <obj.rda>     <out.bed.gz|out.bw> [--chrom-sizes f]
    rda_to_bed.py clusters <obj.rda>     <out.bed.gz>
    rda_to_bed.py regiondb <db.rda>      <out_dir>   (LOLA regionDB -> <out_dir>/<collection>/regions/*.bed.gz)
"""
import argparse
import gzip
import sys
import warnings

import numpy as np
import pandas as pd
import rdata

warnings.filterwarnings("ignore", message="Missing constructor")


def load(path):
    return rdata.conversion.convert(rdata.parser.parse_file(path))


def granges_to_df(g):
    """Convert a GRanges (as SimpleNamespace) into a BED-like DataFrame."""
    vals = np.asarray(g.seqnames.values.astype(str))
    lens = np.asarray(g.seqnames.lengths, dtype=np.int64)
    chrom = np.repeat(vals, lens)
    start1 = np.asarray(g.ranges.start, dtype=np.int64)
    width = np.asarray(g.ranges.width, dtype=np.int64)
    df = pd.DataFrame({"chr": chrom, "start": start1 - 1, "end": start1 - 1 + width})
    for k, v in getattr(g.elementMetadata, "listData", {}).items():
        df[str(k)] = np.asarray(v)
    return df


def write_bed(df, out):
    df.to_csv(out, sep="\t", header=False, index=False,
              compression="gzip" if out.endswith(".gz") else None)


def write_bw(df, out, chrom_sizes=None):
    import pyBigWig
    df = df.sort_values(["chr", "start"])
    if chrom_sizes is not None:
        cs = pd.read_csv(chrom_sizes, sep="\t", header=None, names=["chr", "len"])
        sizes = dict(zip(cs.chr, cs.len))
    else:
        sizes = df.groupby("chr").end.max().to_dict()
    chroms = [c for c in sizes if c in set(df.chr)]
    bw = pyBigWig.open(out, "w")
    bw.addHeader([(c, int(sizes[c])) for c in chroms])
    for c in chroms:
        sub = df[df.chr == c]
        bw.addEntries([c] * len(sub), sub.start.tolist(), ends=sub.end.tolist(),
                      values=sub.score.astype(float).tolist())
    bw.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["sysdata", "granges", "clusters", "regiondb"])
    ap.add_argument("rda")
    ap.add_argument("out")
    ap.add_argument("--chrom-sizes")
    a = ap.parse_args()
    objs = load(a.rda)
    if a.mode == "sysdata":
        names = {"gene": "genic", "igr": "intergenic", "bl": "blacklist"}
        for asm in ("hg38", "mm10"):
            for key, label in names.items():
                df = granges_to_df(objs[f"{asm}_{key}"])[["chr", "start", "end"]]
                write_bed(df, f"{a.out}/{asm}_{label}.bed.gz")
                print(f"{asm}_{label}: {len(df)} regions", file=sys.stderr)
            cs = objs[f"{asm}_chrom_sizes"]
            cs = pd.DataFrame(cs) if not isinstance(cs, pd.DataFrame) else cs
            cs.to_csv(f"{a.out}/{asm}.chrom.sizes", sep="\t", header=False, index=False)
    elif a.mode == "granges":
        df = granges_to_df(list(objs.values())[0])
        if a.out.endswith(".bw"):
            write_bw(df, a.out, a.chrom_sizes)
        else:
            write_bed(df, a.out)
        print(f"{len(df)} bins", file=sys.stderr)
    elif a.mode == "regiondb":
        import os
        db = list(objs.values())[0]
        anno = db["regionAnno"]
        grl = db["regionGRL"]
        flat = granges_to_df(grl.unlistData)
        ends = np.asarray(grl.partitioning.end, dtype=np.int64)
        starts = np.r_[0, ends[:-1]]
        for k, (_, row) in enumerate(anno.iterrows()):
            d = os.path.join(a.out, str(row["collection"]), "regions")
            os.makedirs(d, exist_ok=True)
            sub = flat.iloc[starts[k]:ends[k]][["chr", "start", "end"]]
            write_bed(sub, os.path.join(d, str(row["filename"])))
            print(row["collection"], row["filename"], len(sub), file=sys.stderr)
        for coll, sub in anno.groupby("collection"):
            sub.drop(columns=[c for c in sub.columns if c == "size"]).to_csv(
                os.path.join(a.out, str(coll), "index.txt"), sep="\t", index=False)
    else:
        clus = list(objs.values())[0]
        parts = []
        for lab, g in clus.items():
            df = granges_to_df(g)[["chr", "start", "end"]]
            df["name"] = str(lab)
            parts.append(df)
        df = pd.concat(parts)
        write_bed(df, a.out)
        print(df.name.value_counts().sort_index().to_dict(), file=sys.stderr)


if __name__ == "__main__":
    main()
