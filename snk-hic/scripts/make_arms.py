#!/usr/bin/env python
"""Chromosome arms and the cooltools view frame, for any UCSC assembly.

  make_arms.py arms --sizes S --cytoband C --canonical RE --out arms.bed
  make_arms.py view --arms arms.bed --exclude RE --out view_df.bed

Port extension used by `pipeline references` (no nf-core counterpart).

arms: the canonical chromosomes (names matching --canonical), in natural order
(chr1, chr2, ..., chr10, ..., chrX, chrY). A chromosome with `acen` bands in the
UCSC cytoBand table is split at the end of its p-side acen band (the p11.1/q11.1
boundary; the midpoint of the acen span when no band name starts with "p") into
<chrom>_p and <chrom>_q. A chromosome without acen bands stays whole, named
after itself (mouse and fly cytoBand tables have none). bioframe's
fetch_centromeres only covers human assemblies, hence the cytoBand parsing.

view: the arms minus the chromosomes matching --exclude.

Output: 4-column BED (chrom, start, end, name), no header.
"""

import argparse
import gzip
import re
import sys


def natural_key(chrom):
    base = chrom[3:] if chrom.startswith("chr") else chrom
    m = re.match(r"(\d+)(.*)$", base)
    if m:
        return (0, int(m.group(1)), m.group(2))
    return (1, 0, base)


def read_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                sizes[f[0]] = int(f[1])
    return sizes


def read_acen(path):
    """chrom -> list of (start, end, band name) of its acen bands."""
    opener = gzip.open if path.endswith(".gz") else open
    acen = {}
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[4] == "acen":
                acen.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3]))
    return acen


def centromere(bands):
    p_bands = [b for b in bands if b[2].startswith("p")]
    if p_bands:
        return max(b[1] for b in p_bands)
    return (min(b[0] for b in bands) + max(b[1] for b in bands)) // 2


def arms(sizes, acen, canonical):
    rows = []
    for chrom in sorted((c for c in sizes if canonical.search(c)), key=natural_key):
        length = sizes[chrom]
        bands = acen.get(chrom)
        mid = centromere(bands) if bands else None
        if mid is not None and 0 < mid < length:
            rows.append((chrom, 0, mid, chrom + "_p"))
            rows.append((chrom, mid, length, chrom + "_q"))
        else:
            rows.append((chrom, 0, length, chrom))
    return rows


def write_bed(rows, path):
    with open(path, "w") as out:
        for r in rows:
            out.write("%s\t%d\t%d\t%s\n" % r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    a = sub.add_parser("arms")
    a.add_argument("--sizes", required=True, help="chrom sizes (chrom<TAB>length)")
    a.add_argument("--cytoband", required=True, help="UCSC cytoBand.txt(.gz)")
    a.add_argument("--canonical", required=True, help="regex of the chromosomes kept")
    a.add_argument("--out", required=True)
    v = sub.add_parser("view")
    v.add_argument("--arms", required=True)
    v.add_argument("--exclude", required=True, help="regex of the chromosomes dropped")
    v.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "arms":
        rows = arms(read_sizes(args.sizes), read_acen(args.cytoband), re.compile(args.canonical))
        if not rows:
            sys.exit("make_arms.py: no chromosome matches --canonical %r" % args.canonical)
        write_bed(rows, args.out)
    elif args.cmd == "view":
        drop = re.compile(args.exclude)
        rows = []
        with open(args.arms) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4 and not drop.search(f[0]):
                    rows.append((f[0], int(f[1]), int(f[2]), f[3]))
        write_bed(rows, args.out)
    else:
        ap.error("choose a command: arms or view")


if __name__ == "__main__":
    main()
