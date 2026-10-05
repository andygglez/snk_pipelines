#!/usr/bin/env python3
"""RustQC output handling from nf-core/rnaseq 3.27.0 (Groovy in the pipeline).

  published_name  conf/modules/rustqc.config publishDir saveAs: maps
                  '<sample>/<tool>/<file>' to its path under <aligner>/rustqc/
                  (None = not published)
  mqc_files       workflows/rnaseq: the mqcKeep filter on the RUSTQC outputs
                  that go to MultiQC, in emit order

Usage: rustqc_publish.py mqc <rustqc sample dir>   (prints MultiQC files)
       rustqc_publish.py ie <rustqc sample dir>    (prints infer_experiment.txt)
"""

import os
import re
import sys

RSEQC_SUBDIR_TOOLS = ("junction_annotation", "junction_saturation", "inner_distance", "read_duplication")


EMITTED = ("dupradar/", "featurecounts/", "preseq/", "samtools/", "rseqc/", "qualimap/")


def published_name(filename):
    stripped = re.sub(r"^[^/]+/", "", filename, count=1)
    if not stripped.startswith(EMITTED):  # only RUSTQC's emitted outputs are published
        return None
    base = [t for t in stripped.split("/") if t][-1]

    if stripped.startswith("qualimap/"):
        sample = filename.split("/")[0]
        return re.sub(r"^qualimap/", f"qualimap/{sample}/", stripped, count=1)

    if stripped.startswith("featurecounts/"):
        if base.endswith(".featureCounts.biotype.tsv.summary"):
            return "featurecounts/" + base.replace(".biotype.tsv.summary", ".tsv.summary")
        if base.endswith(".featureCounts.tsv.summary"):
            return None
        return "featurecounts/" + base

    if stripped.startswith("samtools/"):
        return "samtools_stats/" + base

    if stripped.startswith("dupradar/"):
        if "Boxplot" in base:
            return "dupradar/box_plot/" + base
        if "ExpDens" in base and "Curve_mqc" not in base:
            return "dupradar/scatter_plot/" + base
        if "expressionHist" in base:
            return "dupradar/histogram/" + base
        if "dupMatrix" in base:
            return "dupradar/gene_data/" + base
        if "intercept_slope" in base:
            return "dupradar/intercepts_slope/" + base
        return "dupradar/" + base

    if stripped.startswith("rseqc/"):
        rseqc_path = re.sub(r"^rseqc/", "", stripped, count=1)
        tool = [t for t in rseqc_path.split("/") if t][0]
        if tool in RSEQC_SUBDIR_TOOLS:
            if base.endswith(".r"):
                return f"rseqc/{tool}/rscript/{base}"
            if base.endswith(".png") or base.endswith(".svg"):
                return f"rseqc/{tool}/plot/{base}"
            if base.endswith(".xls"):
                return f"rseqc/{tool}/xls/{base}"
            if base.endswith(".bed"):
                return f"rseqc/{tool}/bed/{base}"
            if base.endswith(".junction_annotation.log"):
                return f"rseqc/{tool}/log/{base}"
            if base.endswith(".txt"):
                return f"rseqc/{tool}/txt/{base}"
        return stripped

    return stripped


_KEEP = re.compile(r"(?i)\.(txt|tsv|xls|log|stats|flagstat|idxstats|html)$")


def _keep(name):
    if name.endswith(".featureCounts.tsv.summary"):
        return False
    return bool(_KEEP.search(name)) or "_mqc." in name


def _files(top, recursive):
    if not os.path.isdir(top):
        return []
    if not recursive:
        return sorted(os.path.join(top, f) for f in os.listdir(top) if os.path.isfile(os.path.join(top, f)))
    out = []
    for root, _dirs, files in os.walk(top):
        out += [os.path.join(root, f) for f in files]
    return sorted(out)


def mqc_files(sample_dir):
    """RUSTQC dupradar, featurecounts, preseq, samtools, rseqc, qualimap outputs kept for MultiQC."""
    out = []
    for sub, recursive in (("dupradar", False), ("featurecounts", False), ("preseq", False),
                           ("samtools", False), ("rseqc", True), ("qualimap", True)):
        out += [f for f in _files(os.path.join(sample_dir, sub), recursive) if _keep(os.path.basename(f))]
    return out


def infer_experiment(sample_dir):
    return [f for f in _files(os.path.join(sample_dir, "rseqc"), True) if f.endswith(".infer_experiment.txt")][:1]


if __name__ == "__main__":
    mode, sample_dir = sys.argv[1], sys.argv[2]
    for f in (mqc_files if mode == "mqc" else infer_experiment)(sample_dir):
        print(f)
