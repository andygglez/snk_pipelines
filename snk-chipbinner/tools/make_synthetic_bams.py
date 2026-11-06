#!/usr/bin/env python3
"""Make small SYNTHETIC ChIP BAMs for testing the workflow end to end (not real data).

hg38 chr21 + chr22, 10 kb bins. Per-bin read rates mimic an NSD1-KO-like
H3K36me2 change: WT has signal on genic and intergenic bins; KO keeps genic
signal and loses most intergenic signal. Two replicates per condition with
different depths, 100 bp single-end reads, coordinate-sorted and indexed.

    make_synthetic_bams.py <out_dir> [--seed 1]
"""
import argparse
import os
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chipbinner.binning import make_windows  # noqa: E402
from chipbinner.genomics import RegionSet, annotation_paths, read_chrom_sizes  # noqa: E402

CHROMS = ["chr21", "chr22"]
DESIGN = {  # condition: (genic rate, intergenic rate, other rate)
    "WT": (60, 45, 30),
    "KO": (55, 9, 20),
}
DEPTH = {"rep1": 1.0, "rep2": 0.7}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    ann = annotation_paths("hg38")
    sizes = read_chrom_sizes(ann["chrom_sizes"])
    win = make_windows(ann["chrom_sizes"], 10000, chroms=set(CHROMS))
    g = RegionSet.from_bed(ann["gene"]).overlaps_any(win)
    i = RegionSet.from_bed(ann["igr"]).overlaps_any(win)
    cls = np.where(g & ~i, 0, np.where(i & ~g, 1, 2))
    # shared per-bin "biology" (so replicates correlate), per-sample noise on top
    base = rng.gamma(4.0, 0.25, len(win))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": c, "LN": int(sizes[c])} for c in CHROMS]}
    for cond, rates in DESIGN.items():
        for rep, depth in DEPTH.items():
            lam = np.array(rates)[cls] * base * depth * rng.gamma(20.0, 0.05, len(win))
            n = rng.poisson(lam)
            path = os.path.join(a.out_dir, f"{cond}_{rep}_H3K36me2.bam")
            with pysam.AlignmentFile(path, "wb", header=header) as out:
                k = 0
                for c in CHROMS:
                    m = (win["chr"] == c).to_numpy()
                    st = win["start"].to_numpy()[m]
                    en = win["end"].to_numpy()[m]
                    pos = np.concatenate([rng.integers(s, max(s + 1, e - 100), size=x)
                                          for s, e, x in zip(st, en, n[m])])
                    pos.sort()
                    tid = CHROMS.index(c)
                    for p in pos:
                        r = pysam.AlignedSegment()
                        r.query_name = f"r{k}"
                        r.flag = 16 if k % 2 else 0
                        r.reference_id = tid
                        r.reference_start = int(p)
                        r.mapping_quality = 42
                        r.cigartuples = [(0, 100)]
                        r.query_sequence = "A" * 100
                        r.query_qualities = pysam.qualitystring_to_array("I" * 100)
                        out.write(r)
                        k += 1
            pysam.index(path)
            print(f"{path}: {n.sum()} reads")


if __name__ == "__main__":
    main()
