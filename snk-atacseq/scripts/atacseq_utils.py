"""Logic that nf-core/atacseq 2.1.2 runs in Groovy (no container).

Ported functions, with the nf-core source they mirror:
  trimgalore_reads_after_filtering  fastq_fastqc_umitools_trimgalore: getTrimGaloreReadsAfterFiltering
  read_validated_samplesheet        input_check: create_fastq_channel (on the
                                    output of bin/check_samplesheet.py)

Groovy's `toFloat()` truncates to 32-bit precision, so f32() is applied at the
same points to reproduce nf-core's comparison against min_trimmed_reads.
"""

import csv
import re
import struct


def f32(value):
    return struct.unpack("f", struct.pack("f", float(value)))[0]


def trimgalore_reads_after_filtering(report):
    total = filtered = 0.0
    with open(report) as fh:
        for line in fh:
            m = re.search(r"([\d\.]+)\ssequences processed in total", line)
            if m:
                total = f32(m.group(1))
            m = re.search(r"shorter than the length cutoff[^:]+:\s([\d\.]+)", line)
            if m:
                filtered = f32(m.group(1))
    return total - filtered


def strip_run(sample_id):
    """meta.id - ~/_T\\d+$/ (library id -> merged-library id)."""
    return re.sub(r"_T\d+$", "", sample_id, count=1)


def strip_rep(sample_id):
    """meta.id - ~/_REP\\d+$/ (merged-library id -> merged-replicate id)."""
    return re.sub(r"_REP\d+$", "", sample_id, count=1)


def read_validated_samplesheet(path, seq_center=""):
    """Rows of samplesheet.valid.csv as create_fastq_channel metadata."""
    libs = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            lib = row["sample"]
            rg = f"@RG\\tID:{lib}\\tSM:{strip_run(lib)}\\tPL:ILLUMINA\\tLB:{lib}\\tPU:1"
            if seq_center:
                rg += f"\\tCN:{seq_center}"
            libs.append({
                "id": lib,
                "single_end": row["single_end"] in ("1", "true", "True"),
                "control": row["control"],
                "read_group": f"'{rg}'",
                "fastq_1": row["fastq_1"],
                "fastq_2": row["fastq_2"],
            })
    return libs
