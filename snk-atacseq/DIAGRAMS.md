# snk-atacseq — pipeline diagrams

Snakemake port of **nf-core/atacseq 2.1.2**. These diagrams show the flow of commands and data
through `Snakefile`, and how each key in `config/config.yaml` changes those commands. They are
written in Mermaid. You can view them on GitHub, in VS Code (with a Mermaid preview extension),
or by opening `DIAGRAMS.html` in a browser.

**How to read them**

- `«key»` is a value from `config/config.yaml`. `«args.X»` is the fixed tool argument string under `args:` (nf-core `ext.args`).
- Dashed orange boxes run only under some settings. Pink hexagons are **gates** (checkpoints that drop libraries or samples).
- Levels: a **library** (`_T<m>`) is one samplesheet row. A **merged library** (`.mLb.`) is one replicate. A **merged replicate** (`.mRp.`) is all replicates of a sample. Peaks are called at both the mLb and the mRp level.
- Reads → filtered BAM is the same as in snk-chipseq (see `../snk-chipseq/DIAGRAMS.md`, figure 2b), with the ATAC-specific differences shown here.

---

## 1. Analysis overview

```mermaid
flowchart TD
    classDef input fill:#dbeafe,stroke:#2563eb,color:#111
    classDef step fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef opt fill:#fff7ed,stroke:#ea580c,color:#111,stroke-dasharray:5 4
    classDef gate fill:#fce7f3,stroke:#db2777,color:#111
    classDef out fill:#dcfce7,stroke:#16a34a,color:#111

    SS[/"samplesheet: sample, fastq_1, fastq_2, replicate<br/>(+ control, control_replicate if with_control)"/]:::input
    REF[/"FASTA + GTF/GFF, blacklist, mito_name"/]:::input
    GR["reference prep: gene BED, TSS BED, sizes, autosomes,<br/>include regions (minus blacklist and chrM), index, khmer gsize"]:::step
    TG["FastQC · Trim Galore"]:::opt
    G1{{"gate: min_trimmed_reads"}}:::gate
    AL["align (bwa / bowtie2 / chromap / STAR) · sort"]:::step
    MLB["merged library: merge runs · MarkDuplicates · filter → mLb.clN BAM"]:::out
    QC["Preseq · Picard metrics · deepTools profiles · fingerprint"]:::opt
    MRP["merged replicate: merge clN BAMs · MarkDuplicates REMOVE → mRp.clN BAM<br/>unless skip_merge_replicates"]:::opt
    BW["CPM bigWigs (mLb and mRp)"]:::out
    PK["MACS2 callpeak (mLb and mRp)"]:::out
    G2{{"gate: non-empty peak file"}}:::gate
    AN["FRiP · counts · HOMER · QC plots"]:::opt
    CS["consensus peaks · featureCounts · DESeq2 QC (per level)"]:::opt
    AQ["ataqv + mkarv report (mLb peaks)"]:::opt
    OUT["IGV session · MultiQC"]:::out

    REF --> GR
    SS --> TG --> G1 --> AL --> MLB --> QC & BW & PK
    MLB --> MRP --> BW & PK
    PK --> G2 --> AN & CS & AQ
    BW & AN & CS & AQ & QC --> OUT
```

---

## 2. Command flow — ATAC-specific steps

```mermaid
flowchart TD
    classDef step fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef opt fill:#fff7ed,stroke:#ea580c,color:#111,stroke-dasharray:5 4
    classDef gate fill:#fce7f3,stroke:#db2777,color:#111
    classDef out fill:#dcfce7,stroke:#16a34a,color:#111

    subgraph REF["Reference extras"]
        TSS["tss_extract<br/>awk gene BED → 1-bp TSS BED<br/>unless «tss_bed»"]:::opt
        AUT["get_autosomes<br/>get_autosomes.py fai → autosomes.txt"]:::step
        BLR["genome_blacklist_regions<br/>complementBed of «blacklist» (or whole genome)<br/>| awk drop «mito_name» unless «keep_mito»"]:::step
    end

    MKD[/"mLb.mkD BAM"/]
    FLT["bamtools_filter (same flags as chipseq)<br/>-F 0x0400 unless «keep_dups», -q 1 unless «keep_multi_map»<br/>-L include_regions (no chrM)"]:::step
    CLN[/"mLb.clN BAM"/]

    subgraph DT["deepTools — unless «skip_plot_profile»"]
        CMS["computematrix_scale_regions<br/>computeMatrix «args.computematrix_scale_regions» -R gene BED"]:::opt
        CMR["computematrix_reference_point<br/>computeMatrix «args.computematrix_reference_point» -R TSS BED"]:::opt
        PP["plotProfile (scale-regions matrix)"]:::opt
        PH["plotHeatmap (reference-point matrix)"]:::opt
    end
    FPR["deeptools_plotfingerprint<br/>«args.plotfingerprint» --numberOfSamples «fingerprint_bins»<br/>[+ control if «with_control»]"]:::opt

    subgraph REP["Merged replicate — unless «skip_merge_replicates»; groups with more than one replicate"]
        MRS["picard_mergesamfiles_replicate<br/>MergeSamFiles «args.picard_mergesamfiles» on mLb.clN BAMs"]:::opt
        MRD["picard_markduplicates_replicate<br/>MarkDuplicates «args.picard_markduplicates_replicate»<br/>(REMOVE_DUPLICATES true) → mRp.clN BAM"]:::opt
    end

    GC["bedtools_genomecov → ucsc_bedgraphtobigwig<br/>scale 1e6/mapped [SE: -fs «fragment_size»] — mLb and mRp"]:::out
    MC["macs2_callpeak — mLb and mRp<br/>macs2 callpeak «args.macs2»<br/>[--broad --broad-cutoff «broad_cutoff» unless «narrow_peak»]<br/>[--bdg --SPMR if «save_macs_pileup»]<br/>mLb only: [--pvalue «macs_pvalue»] [--qvalue «macs_fdr»]<br/>--gsize «macs_gsize» or khmer -k «read_length»<br/>[--control if «with_control»]"]:::out
    G2{{"peak_summary checkpoint per level<br/>drop empty peak files"}}:::gate
    FR["frip_score · multiqc_custom_peaks<br/>intersectBed «args.frip»"]:::step
    HO["homer_annotatepeaks «args.homer_annotatepeaks»<br/>unless «skip_peak_annotation»"]:::opt
    PQ["plot_macs2_qc · plot_homer_annotatepeaks<br/>unless «skip_peak_qc»"]:::opt
    CON["macs2_consensus per level (more than one peak file)<br/>mergeBed · macs2_merged_expand.py<br/>mLb only: --min_replicates «min_reps_consensus»<br/>unless «skip_consensus_peaks»"]:::opt
    FC["featurecounts «args.featurecounts» on mLb.clN BAMs"]:::opt
    DQ["deseq2_qc_library / deseq2_qc_replicate<br/>«args.deseq2_qc» [--vst TRUE if «deseq2_vst»]<br/>unless «skip_deseq2_qc»"]:::opt
    AQ["ataqv<br/>ataqv «args.ataqv» [--mitochondrial-reference-name «ataqv_mito_reference» or «mito_name»]<br/>--peak-file mLb peaks --tss-file TSS --autosomal-reference-file NA mkD.bam"]:::opt
    MK["ataqv_mkarv<br/>mkarv → interactive HTML<br/>both unless «skip_ataqv»"]:::opt

    BLR --> FLT
    MKD --> FLT --> CLN
    CLN --> CMS & CMR & FPR & GC & MC
    TSS --> CMR & AQ
    CMS --> PP
    CMR --> PH
    CLN --> MRS --> MRD --> GC & MC
    MC --> G2 --> FR & HO & CON
    HO --> PQ
    CON --> FC --> DQ
    G2 --> AQ
    MKD --> AQ
    AUT --> AQ --> MK
```

---

## 3. How config switches shape the run

```mermaid
flowchart TD
    classDef key fill:#fef9c3,stroke:#ca8a04,color:#111
    classDef rule fill:#f3f4f6,stroke:#6b7280,color:#111

    WC{"«with_control»?"}:::key
    WC -->|"no (default)"| W1["every merged library calls peaks alone<br/>samplesheet without control columns"]:::rule
    WC -->|yes| W2["samples with a control call peaks vs it<br/>(mLb and mRp), fingerprint includes it"]:::rule
    MR{"«skip_merge_replicates»?"}:::key
    MR -->|no| M1["mRp BAMs, bigWigs, peaks, consensus, DESeq2"]:::rule
    MR -->|yes| M2["mLb level only"]:::rule
    MT{"«keep_mito»?"}:::key
    MT -->|no| T1["«mito_name» removed from include regions"]:::rule
    NP{"«narrow_peak»?"}:::key
    NP -->|no| N1["broad peaks, outputs under macs2/broad_peak/"]:::rule
    NP -->|yes| N2["narrow peaks + summits, macs2/narrow_peak/"]:::rule
    AL{"«aligner»"}:::key --> A1["bwa / bowtie2 / chromap (--Tn5-shift) / star<br/>index built unless «*_index»"]:::rule
    GS{"«macs_gsize» set?"}:::key -->|no| G1["khmer -k «read_length»"]:::rule
    SK["«skip_qc» «skip_fastqc» «skip_preseq» (default true)<br/>«skip_picard_metrics» «skip_plot_profile» «skip_plot_fingerprint»<br/>«skip_ataqv» «skip_peak_annotation» «skip_peak_qc»<br/>«skip_consensus_peaks» «skip_deseq2_qc» «skip_igv» «skip_multiqc»"]:::key --> S1["each turns its step off"]:::rule
```

---

## 4. Config key reference (ATAC-specific; shared keys behave as in snk-chipseq)

| Key | Rule(s) | Effect on the command |
|---|---|---|
| `with_control` | samplesheet check, macs2_callpeak, plotFingerprint | Read control columns and call peaks vs controls. |
| `tss_bed` | tss_extract, computeMatrix reference-point, ataqv | Given TSS BED; otherwise derived from the gene BED. |
| `mito_name`, `keep_mito` | genome_blacklist_regions | Drop the mitochondrial contig from the include regions unless kept. |
| `ataqv_mito_reference` | ataqv | `--mitochondrial-reference-name` (falls back to `mito_name`). |
| `skip_merge_replicates` | merged-replicate rules | Turn off the mRp level. |
| `macs_fdr`, `macs_pvalue` | macs2_callpeak | `--qvalue` / `--pvalue`, merged-library level only. |
| `min_reps_consensus` | macs2_consensus | `--min_replicates`, merged-library consensus only. |
| `skip_ataqv` | ataqv, ataqv_mkarv | Turn ataqv off. |
| `skip_preseq` | preseq_lcextrap | Default `true` here. |
| `args.picard_markduplicates_library` / `_replicate` | picard_markduplicates_* | Duplicates are marked (mLb) or removed (mRp). |
| `args.computematrix_scale_regions` / `_reference_point` | deeptools_computematrix_* | Gene-body and TSS profiles. |
| `args.chromap` | chromap | Includes `--Tn5-shift`. |
| `args.macs2`, `args.ataqv` | macs2_callpeak, ataqv | `--keep-dup all --nomodel`; `--ignore-read-groups`. |
| other keys (`fasta`, `gtf`, `aligner`, trimming, `keep_dups`, `keep_multi_map`, `fragment_size`, `narrow_peak`, `save_*`, `resources`, ...) | — | Same as `snk-chipseq/DIAGRAMS.md`, section 4. |

To regenerate the real rule DAG: `snakemake --configfile config/test.yaml --rulegraph mermaid-js > rulegraph.mmd`
(rules after the checkpoints appear only once those have run).
