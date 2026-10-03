#!/usr/bin/env python3
"""Build the MultiQC side files that nf-core/rnaseq 3.27.0 creates in Groovy.

Mirrors subworkflows/local/multiqc_rnaseq (merged-report branch):
  fail_trimmed_samples_mqc.tsv     samples with <= min_trimmed_reads after trimming
  fail_mapped_samples_mqc.tsv      samples below min_mapped_reads % uniquely mapped
  strand_check_summary_mqc.json    strandCheckSummaryYaml
  strand_check_composition_mqc.json strandCheckCompositionYaml
  multiqc_sample_merge.yml         multiqcSampleMergeYaml
  name_replacement.txt             multiqcNameReplacements
Files that nf-core would not create (no rows) are not written.
"""

import argparse
import csv
import json
import os
import re
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rnaseq_strand import (  # noqa: E402
    classify_strand,
    java_double_str,
    java_float_str,
    round_one_decimal,
)


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def load_asset(path):
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return {k: v for k, v in data.items() if not str(k).startswith("_")}


def simple_name(path):
    return os.path.basename(path).split(".")[0]


def inference_certainty(analysis):
    if not analysis:
        return None
    fwd, rev = analysis["forwardFragments"], analysis["reverseFragments"]
    stranded = fwd + rev
    if stranded == 0:
        return None
    if analysis["inferred_strandedness"] == "forward":
        return fwd / stranded * 100
    if analysis["inferred_strandedness"] == "reverse":
        return rev / stranded * 100
    return None


def summary_cells(provided, status, salmon, rseqc):
    salmon = salmon or {}
    rseqc = rseqc or {}
    return {
        "provided": provided,
        "salmon_inferred": salmon.get("inferred_strandedness") or "-",
        "salmon_pct": round_one_decimal(inference_certainty(salmon)),
        "salmon_s": round_one_decimal(salmon.get("forwardFragments")),
        "salmon_a": round_one_decimal(salmon.get("reverseFragments")),
        "salmon_u": round_one_decimal(salmon.get("unstrandedFragments")),
        "rseqc_inferred": rseqc.get("inferred_strandedness") or "-",
        "rseqc_pct": round_one_decimal(inference_certainty(rseqc)),
        "rseqc_s": round_one_decimal(rseqc.get("forwardFragments")),
        "rseqc_a": round_one_decimal(rseqc.get("reverseFragments")),
        "rseqc_u": round_one_decimal(rseqc.get("unstrandedFragments")),
        "status": status,
    }


def composition(analysis):
    return {
        "Sense": round_one_decimal(analysis["forwardFragments"]),
        "Antisense": round_one_decimal(analysis["reverseFragments"]),
        "Unstranded": round_one_decimal(analysis["unstrandedFragments"]),
    }


def sample_merge_yaml(rows):
    pe_ids = sorted({r["sample"] for r in rows if r.get("fastq_2")})
    if not pe_ids:
        return "table_sample_merge: {}\n"

    def pattern(sid, read):
        esc = re.sub(r"[\\^$.|?*+()\[\]{}/]", lambda m: "\\" + m.group(0), sid).replace("'", "''")
        return f"    - type: regex\n      pattern: '(?<=^{esc})_{read}$'"

    r1 = "\n".join(pattern(s, 1) for s in pe_ids)
    r2 = "\n".join(pattern(s, 2) for s in pe_ids)
    return f'table_sample_merge:\n  "Read 1":\n{r1}\n  "Read 2":\n{r2}\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--trim-status", nargs="*", default=[])
    ap.add_argument("--map-status", nargs="*", default=[])
    ap.add_argument("--strand-status", nargs="*", default=[])
    ap.add_argument("--infer-experiment", nargs="*", default=[],
                    help="<sample>.infer_experiment.txt files; empty when infer_experiment is not run")
    ap.add_argument("--aligner-display-name", default="Aligned reads")
    ap.add_argument("--min-trimmed-reads", type=float, required=True)
    ap.add_argument("--stranded-threshold", type=float, required=True)
    ap.add_argument("--unstranded-threshold", type=float, required=True)
    ap.add_argument("--summary-asset", required=True)
    ap.add_argument("--composition-asset", required=True)
    ap.add_argument("--status-header", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    def out(name):
        return os.path.join(args.outdir, name)

    # fail_trimmed_samples_mqc.tsv (header kept once, as collectFile keepHeader)
    failed = []
    for f in args.trim_status:
        st = load_json(f)
        if st["num_reads"] is not None and st["num_reads"] <= args.min_trimmed_reads:
            failed.append(f"{st['sample']}\t{java_double_str(st['num_reads'])}")
    if failed:
        with open(out("fail_trimmed_samples_mqc.tsv"), "w") as fh:
            fh.write("Sample\tReads after trimming\n" + "".join(r + "\n" for r in sorted(failed)))

    # fail_mapped_samples_mqc.tsv
    with open(args.status_header) as fh:
        status_header = fh.read()
    failed = []
    for f in args.map_status:
        st = load_json(f)
        if st["pass"] is False:
            failed.append(f"{st['sample']}\t{java_float_str(st['percent_mapped'])}")
    if failed:
        with open(out("fail_mapped_samples_mqc.tsv"), "w") as fh:
            fh.write(status_header + f"Sample\t{args.aligner_display_name} (%)\n")
            fh.write("".join(r + "\n" for r in sorted(failed)))

    # Strandedness check rows: [sample, provided, status, salmon, rseqc]
    strand = {st["sample"]: st for st in map(load_json, args.strand_status)}
    rows = []
    if args.infer_experiment:
        for f in args.infer_experiment:
            sample = os.path.basename(f)[: -len(".infer_experiment.txt")]
            provided, status, salmon, rseqc = classify_strand(
                strand[sample], f, args.stranded_threshold, args.unstranded_threshold
            )
            rows.append((sample, provided, status, salmon, rseqc))
    else:
        for sample, st in strand.items():
            if st.get("salmon_strand_analysis"):
                rows.append((sample, "auto", "-", st["salmon_strand_analysis"], None))
    rows.sort(key=lambda r: r[0])

    if rows:
        summary_static = load_asset(args.summary_asset)
        header_keys = list(summary_static["headers"].keys())
        data = {}
        for sample, provided, status, salmon, rseqc in rows:
            raw = summary_cells(provided, status, salmon, rseqc)
            unknown = set(raw) - set(header_keys)
            if unknown:
                sys.exit(f"strand_check_summary.yaml headers do not declare columns: {unknown}")
            data[sample] = {k: raw[k] for k in header_keys if raw.get(k) is not None}
        with open(out("strand_check_summary_mqc.json"), "w") as fh:
            json.dump({**summary_static, "data": data}, fh, indent=4)

        comp_static = load_asset(args.composition_asset)
        for k in ("parent_id", "parent_name", "parent_description"):
            if k in summary_static:
                comp_static[k] = summary_static[k]
        rseqc_data = {s: composition(r) for s, _, _, _, r in rows if r}
        salmon_data = {s: composition(sa) for s, _, _, sa, _ in rows if sa}
        datasets, labels = [], []
        if rseqc_data:
            datasets.append(rseqc_data)
            labels.append("RSeQC")
        if salmon_data:
            datasets.append(salmon_data)
            labels.append("Salmon")
        pconfig = dict(comp_static["pconfig"])
        if len(datasets) > 1:
            pconfig["data_labels"] = [{"name": lab, "ylab": pconfig.get("ylab")} for lab in labels]
        comp_static["pconfig"] = pconfig
        comp_static["data"] = datasets[0] if len(datasets) == 1 else datasets
        with open(out("strand_check_composition_mqc.json"), "w") as fh:
            json.dump(comp_static, fh, indent=4)

    # Samplesheet-derived files: first run of each sample
    first_run = {}
    with open(args.samplesheet, newline="") as fh:
        for row in csv.DictReader(fh):
            first_run.setdefault(str(row["sample"]).strip(), row)
    sheet_rows = [{"sample": s, **r} for s, r in first_run.items()]

    with open(out("multiqc_sample_merge.yml"), "w") as fh:
        fh.write(sample_merge_yaml(sheet_rows))

    mappings = []
    for r in sheet_rows:
        paired = bool(r.get("fastq_2"))
        suffixes = ["_1", "_2"] if paired else [""]
        name1 = simple_name(r["fastq_1"])
        if name1 != r["sample"]:
            mappings.append(f"{name1}\t{r['sample']}{suffixes[0]}")
            if paired:
                mappings.append(f"{simple_name(r['fastq_2'])}\t{r['sample']}{suffixes[1]}")
    if mappings:
        with open(out("name_replacement.txt"), "w") as fh:
            fh.write("".join(m + "\n" for m in sorted(mappings)))


if __name__ == "__main__":
    main()
