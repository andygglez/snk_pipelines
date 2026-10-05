"""Sample-status logic that nf-core/rnaseq 3.27.0 runs in Groovy (no container).

Ported functions, with the nf-core source they mirror:
  calculate_strandedness      fastq_qc_trim_filter_setstrandedness: calculateStrandedness
  salmon_strandedness         fastq_qc_trim_filter_setstrandedness: getSalmonInferredStrandedness
  inferexperiment_strandedness utils_nfcore_rnaseq_pipeline: getInferexperimentStrandedness
  trimgalore_reads_after_filtering fastq_fastqc_umitools_trimgalore: getTrimGaloreReadsAfterFiltering
  star_percent_mapped         align_star: getStarPercentMapped
  bowtie2_percent_mapped      align_bowtie2: getBowtie2PercentMapped
  fastp_reads_after_filtering fastq_fastqc_umitools_fastp: getFastpReadsAfterFiltering
  seqkit_read_length          fastq_remove_rrna: getReadLengthFromSeqkitStats

Groovy's `toFloat()` truncates to 32-bit precision before the value is promoted
to double, so f32() is applied at the same points to reproduce nf-core's numbers.
"""

import json
import math
import re
import struct


def f32(value):
    return struct.unpack("f", struct.pack("f", float(value)))[0]


def _java_str(value, precision):
    """Java Float/Double.toString: shortest repr, scientific outside [1e-3, 1e7)."""
    if value == 0:
        return "0.0"
    for p in range(1, precision + 1):
        s = f"{value:.{p}g}"
        if (f32(float(s)) if precision == 9 else float(s)) == value:
            break
    mant = float(s)
    if 1e-3 <= abs(mant) < 1e7:
        out = f"{mant:f}".rstrip("0")
        return out + "0" if out.endswith(".") else out
    m, e = f"{mant:.{p - 1}e}".split("e")
    m = m.rstrip("0") if "." in m else m + ".0"
    m = m + "0" if m.endswith(".") else m
    return f"{m}E{int(e)}"


def java_float_str(value):
    return _java_str(f32(value), 9)


def java_double_str(value):
    return _java_str(float(value), 17)


def calculate_strandedness(forward, reverse, unstranded, stranded_threshold=0.8, unstranded_threshold=0.1):
    total = forward + reverse + unstranded
    total_stranded = forward + reverse
    strandedness = "undetermined"
    if total_stranded > 0:
        forward_prop = forward / float(total_stranded)
        reverse_prop = reverse / float(total_stranded)
        if forward_prop >= stranded_threshold:
            strandedness = "forward"
        elif reverse_prop >= stranded_threshold:
            strandedness = "reverse"
        elif abs(forward_prop - reverse_prop) <= unstranded_threshold:
            strandedness = "unstranded"
    return {
        "inferred_strandedness": strandedness,
        "forwardFragments": forward / float(total) * 100,
        "reverseFragments": reverse / float(total) * 100,
        "unstrandedFragments": unstranded / float(total) * 100,
    }


def salmon_strandedness(json_file, stranded_threshold=0.8, unstranded_threshold=0.1):
    with open(json_file) as fh:
        lib = json.load(fh)
    unstranded = sum(lib.get(k) or 0 for k in ("IU", "U", "MU"))
    n_assigned = float(lib.get("num_assigned_fragments") or 0)
    bias = float(lib.get("strand_mapping_bias") or 0.0)
    forward = bias * n_assigned
    reverse = (1 - bias) * n_assigned
    return calculate_strandedness(forward, reverse, unstranded, stranded_threshold, unstranded_threshold)


_IE_PATTERNS = {
    "unstranded": re.compile(r"Fraction of reads failed to determine:\s([\d\.]+)"),
    "forward": re.compile(r'Fraction of reads explained by "\++,--":\s([\d\.]+)'),
    "reverse": re.compile(r'Fraction of reads explained by "\+-,-\+":\s([\d\.]+)'),
    "pe_forward": re.compile(r'Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)'),
    "pe_reverse": re.compile(r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)'),
}


def inferexperiment_strandedness(txt_file, stranded_threshold=0.8, unstranded_threshold=0.1):
    vals = {"forward": 0, "reverse": 0, "unstranded": 0}
    with open(txt_file) as fh:
        for line in fh:
            for key, pat in _IE_PATTERNS.items():
                m = pat.search(line)
                if m:
                    vals[key.replace("pe_", "")] = f32(m.group(1)) * 100
    return calculate_strandedness(
        vals["forward"], vals["reverse"], vals["unstranded"], stranded_threshold, unstranded_threshold
    )


def trimgalore_reads_after_filtering(report):
    total = filtered = 0.0
    with open(report) as fh:
        for line in fh:
            m = re.search(r"([\d\.]+)\ssequences processed in total", line)
            if m:
                total = f32(m.group(1))
            m = re.search(r"shorter than the length cutoff[^:]+:\s*([\d\.]+)", line)
            if m:
                filtered = f32(m.group(1))
    return total - filtered


def star_percent_mapped(log_final):
    percent = 0.0
    with open(log_final) as fh:
        for line in fh:
            m = re.search(r"Uniquely mapped reads %\s*\|\s*([\d\.]+)%", line)
            if m:
                percent = f32(m.group(1))
    return percent


def bowtie2_percent_mapped(log):
    percent = 0.0
    with open(log) as fh:
        for line in fh:
            m = re.search(r"(\d+\.\d+)% overall alignment rate", line)
            if m:
                percent = f32(m.group(1))
    return percent


def fastp_reads_after_filtering(json_file):
    with open(json_file) as fh:
        return int(json.load(fh)["summary"]["after_filtering"]["total_reads"])


def seqkit_read_length(stats_file):
    """Mean avg_len over the files in a `seqkit stats --tabular` table (Math.round)."""
    with open(stats_file) as fh:
        lines = fh.read().splitlines()
    if len(lines) < 2:
        return 100
    header = lines[0].split("\t")
    if "avg_len" not in header:
        return 100
    idx = header.index("avg_len")
    lens = [f32(line.split("\t")[idx]) for line in lines[1:]]
    return int(math.floor(f32(sum(lens) / len(lens)) + 0.5))


def classify_strand(strand_status, infer_experiment_txt, stranded_threshold, unstranded_threshold):
    """utils_nfcore_rnaseq_pipeline: classifyStrand -> (provided, status, salmon, rseqc)."""
    rseqc = inferexperiment_strandedness(infer_experiment_txt, stranded_threshold, unstranded_threshold)
    salmon = strand_status.get("salmon_strand_analysis")
    status = "fail"
    if salmon:
        provided = "auto"
        if salmon["inferred_strandedness"] == rseqc["inferred_strandedness"] != "undetermined":
            status = "pass"
    else:
        provided = strand_status["strandedness"]
        if provided == rseqc["inferred_strandedness"]:
            status = "pass"
    return provided, status, salmon, rseqc


def round_one_decimal(value):
    """Java Math.round(v * 10) / 10.0d (round half up)."""
    return None if value is None else math.floor(value * 10 + 0.5) / 10.0
