#!/usr/bin/env python
"""GC phasing tracks for cooltools eigs-cis, at several resolutions in one pass.

  gc_tracks.py --fasta F --fai F.fai --canonical RE --res 1000 5000 ... \
               --out 1000.gc 5000.gc ...

Port extension used by `pipeline references` (no nf-core counterpart).

Bins follow `cooltools genome binnify`: [i*res, min((i+1)*res, length)) over the
canonical chromosomes (names matching --canonical), in natural order. GC is
(G+C)/(bin length), both cases counted, Ns counted in the length: this is
bioframe.frac_gc(..., mapped_only=False), so all-N bins give 0. (The cooltools
0.5.1 CLI `genome gc` cannot do this: its --mapped-only flag defaults to True,
so it always divides by the ACGT count and gives NaN for all-N bins.)

Output: TSV with header "chrom start end GC", GC as %.6f; use it as
`cooltools eigs-cis --phasing-track <file>::GC`.

Each chromosome is read once through the .fai offsets (memory ~3x its length).
"""

import argparse
import re
import sys

import numpy as np

from make_arms import natural_key


def read_fai(path):
    fai = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fai[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return fai


def gc_mask(fh, length, offset, linebases, linewidth):
    """Boolean array: True where the base is G/C (any case)."""
    if length == 0:
        return np.zeros(0, dtype=bool)
    full, rest = divmod(length, linebases)
    fh.seek(offset)
    raw = fh.read(full * linewidth + rest)
    seq = raw.replace(b"\n", b"").replace(b"\r", b"")
    if len(seq) != length:
        sys.exit("gc_tracks.py: sequence length mismatch (%d vs %d); "
                 "is the .fai from this FASTA?" % (len(seq), length))
    up = np.frombuffer(seq, dtype=np.uint8) & 0xDF  # a-z -> A-Z
    return (up == ord("G")) | (up == ord("C"))


def bin_counts(mask, res):
    """G+C count and length of each bin of size res."""
    n = len(mask)
    full = n // res
    counts = mask[:full * res].reshape(full, res).sum(axis=1, dtype=np.int64)
    lengths = np.full(full, res, dtype=np.int64)
    if n % res:
        counts = np.append(counts, mask[full * res:].sum(dtype=np.int64))
        lengths = np.append(lengths, n % res)
    return counts, lengths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--canonical", required=True, help="regex of the chromosomes kept")
    ap.add_argument("--res", required=True, nargs="+", type=int)
    ap.add_argument("--out", required=True, nargs="+")
    args = ap.parse_args(argv)
    if len(args.res) != len(args.out):
        ap.error("--res and --out need the same number of values")
    if any(r <= 0 for r in args.res):
        ap.error("resolutions must be positive")

    fai = read_fai(args.fai)
    canonical = re.compile(args.canonical)
    chroms = sorted((c for c in fai if canonical.search(c)), key=natural_key)
    if not chroms:
        sys.exit("gc_tracks.py: no chromosome matches --canonical %r" % args.canonical)

    outs = [open(p, "w") for p in args.out]
    for out in outs:
        out.write("chrom\tstart\tend\tGC\n")
    with open(args.fasta, "rb") as fh:
        for chrom in chroms:
            length, offset, linebases, linewidth = fai[chrom]
            mask = gc_mask(fh, length, offset, linebases, linewidth)
            for res, out in zip(args.res, outs):
                counts, lengths = bin_counts(mask, res)
                starts = np.arange(len(counts), dtype=np.int64) * res
                ends = starts + lengths
                frac = counts / lengths
                out.write("".join("%s\t%d\t%d\t%.6f\n" % (chrom, s, e, g)
                                  for s, e, g in zip(starts, ends, frac)))
            del mask
    for out in outs:
        out.close()


if __name__ == "__main__":
    main()
