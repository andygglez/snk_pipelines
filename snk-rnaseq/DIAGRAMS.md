# snk-rnaseq — pipeline diagrams

Snakemake port of **nf-core/rnaseq 3.27.0**, using the default route: `aligner: star_salmon` and
`trimmer: trimgalore`, plus BBSplit, `additional_fasta` and the optional Salmon pseudo-alignment.
These diagrams show the flow of commands and data through `Snakefile`, and how each key in
`config/config.yaml` changes those commands. They are written in Mermaid. You can view them on
GitHub, in VS Code (with a Mermaid preview extension), or by opening `DIAGRAMS.html` in a browser.

**How to read them**

- `«key»` inside a command is a value that comes straight from `config/config.yaml`.
- Dashed orange boxes are rules that only run under some settings (the condition is on the box).
- Pink hexagons are **gates**: aggregate checkpoints that drop samples which fail a threshold.
- Each rule runs in the image given by `config/containers.yaml`.
- `outdir` holds the files nf-core publishes; `workdir` holds intermediates that nf-core keeps in `work/`.

---

## 1. Analysis overview

```mermaid
flowchart TD
    classDef input fill:#dbeafe,stroke:#2563eb,color:#111
    classDef step fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef opt fill:#fff7ed,stroke:#ea580c,color:#111,stroke-dasharray:5 4
    classDef gate fill:#fce7f3,stroke:#db2777,color:#111
    classDef out fill:#dcfce7,stroke:#16a34a,color:#111

    SS[/"samplesheet<br/>sample, fastq_1, fastq_2, strandedness"/]:::input
    REF[/"FASTA + GTF or GFF<br/>optional: transcripts, gene BED, indices"/]:::input

    subgraph G["1 · Reference preparation"]
        GREF["GTF from GFF · GTF filter · spike-ins<br/>gene BED · transcript FASTA · fai/sizes"]:::step
        GIDX["STAR index · Salmon index · BBSplit index"]:::step
    end

    subgraph R["2 · Read processing"]
        CAT["merge runs (cat)"]:::step
        LINT["fq lint"]:::opt
        FQC["FastQC raw"]:::opt
        TG["Trim Galore"]:::opt
        G1{{"gate: reads after trimming<br/>≥ min_trimmed_reads"}}:::gate
        BBS["BBSplit: remove other genomes"]:::opt
        STR["strandedness: given, or inferred<br/>by Salmon on 1M reads"]:::step
    end

    subgraph A["3 · Alignment and quantification"]
        STAR["STAR 2-pass → genome BAM + transcriptome BAM"]:::step
        G2{{"gate: % mapped<br/>≥ min_mapped_reads"}}:::gate
        SORT["samtools sort / index / stats"]:::step
        MD["Picard MarkDuplicates"]:::opt
        SQ["Salmon quant on transcriptome BAM"]:::step
        TXI["tx2gene → tximport → SummarizedExperiment"]:::out
        DSQ["DESeq2 QC: PCA, sample distances"]:::opt
        PS["Salmon pseudo-alignment<br/>(pseudo_aligner: salmon)"]:::opt
    end

    subgraph Q["4 · Post-alignment QC and tracks"]
        ST["StringTie"]:::opt
        BQ["RSeQC · Qualimap · dupRadar · biotype counts · Preseq"]:::opt
        BW["bigWig coverage (+ strand-specific)"]:::opt
    end

    MQC["MultiQC report"]:::out

    REF --> GREF --> GIDX
    SS --> CAT --> LINT --> TG --> G1 --> BBS --> STR
    CAT --> FQC
    GIDX --> STAR & SQ & PS & BBS & STR
    STR --> STAR --> G2 --> SORT --> MD
    STAR --> SQ --> TXI --> DSQ
    STR --> PS --> TXI
    MD --> ST & BQ & BW
    FQC & TG & STAR & MD & BQ & DSQ --> MQC
```

The `min_trimmed_reads` gate drops a sample from every later step. The `min_mapped_reads`
gate only drops it from MarkDuplicates onward: Salmon quantification still includes every
sample that passed trimming (as in nf-core).

---

## 2. Command flow (what each rule runs)

### 2a · Reference preparation (PREPARE_GENOME)

```mermaid
flowchart TD
    classDef input fill:#dbeafe,stroke:#2563eb,color:#111
    classDef step fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef opt fill:#fff7ed,stroke:#ea580c,color:#111,stroke-dasharray:5 4
    classDef out fill:#dcfce7,stroke:#16a34a,color:#111

    FA[/"«fasta»"/]:::input
    GTF[/"«gtf»"/]:::input
    GFF[/"«gff»"/]:::input
    ADD[/"«additional_fasta»"/]:::input
    TX[/"«transcript_fasta»"/]:::input
    GB[/"«gene_bed»"/]:::input

    GZ["gunzip — any .gz reference<br/>gzip -cd"]:::opt
    LF["link_fasta<br/>put FASTA in genome dir next to its .fai"]:::step
    GFR["gffread_gtf<br/>gffread --keep-exon-attrs -F -T<br/>only if no «gtf»"]:::opt
    GF["custom_gtffilter<br/>gtffilter.py: keep lines on FASTA contigs with transcript_id<br/>[--skip_transcript_id_check if «skip_gtf_transcript_filter»]<br/>unless «skip_gtf_filter»"]:::opt
    CAF["custom_catadditionalfasta<br/>fasta2gtf.py — append spike-ins to FASTA + GTF<br/>biotype = gene_type if «gencode» else «featurecounts_group_type»"]:::opt
    G2B["eautils_gtf2bed<br/>gtf2bed.pl → gene BED<br/>unless «gene_bed»"]:::opt
    PTF["preprocess_transcripts_fasta_gencode<br/>cut -d '|' -f1<br/>«transcript_fasta» + «gencode»"]:::opt
    GFT["gffread_transcripts<br/>gffread -g FASTA -w transcripts.fasta<br/>if «gffread_transcript_fasta»"]:::opt
    MTF["make_transcripts_fasta<br/>rsem-prepare-reference --gtf<br/>default when no «transcript_fasta»"]:::opt
    FAI["samtools_faidx<br/>samtools faidx · cut -f1,2 → .sizes"]:::out

    UT["untar — index given as .tar.gz<br/>tar -xavf --strip-components"]:::opt
    BBI["bbsplit_index<br/>bbsplit.sh ref_primary=FASTA ref_NAME=... build=1<br/>refs from «bbsplit_fasta_list»<br/>unless «skip_bbsplit» or «bbsplit_index»"]:::opt
    SGG["star_genomegenerate<br/>STAR --runMode genomeGenerate --sjdbGTFfile<br/>--genomeSAindexNbases from genome size<br/>--limitGenomeGenerateRAM from resources<br/>unless «star_index»"]:::opt
    SIX["salmon_index<br/>decoys = genome contigs, gentrome = tx + genome<br/>salmon index [--gencode] -k «pseudo_aligner_kmer_size» «extra_salmon_index_args»<br/>unless «salmon_index»"]:::opt

    FA --> GZ --> LF
    FA --> LF
    GFF --> GFR
    GTF --> GF
    GFR --> GF
    LF --> GF --> CAF
    ADD --> CAF
    CAF --> G2B & MTF & GFT & FAI
    GF --> G2B
    GB -.-> G2B
    TX --> PTF
    FAI --> GFT
    LF --> FAI
    CAF --> BBI & SGG & SIX
    MTF & GFT & PTF --> SIX
    TX --> SIX
    UT -.-> SGG & SIX & BBI
```

Without `additional_fasta`, the filtered GTF and the linked FASTA feed every downstream rule
directly (the `custom_catadditionalfasta` box is skipped). `save_reference: true` writes all of
these to `outdir/genome/` instead of `workdir/genome/`.

### 2b · Read processing (FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS)

```mermaid
flowchart TD
    classDef input fill:#dbeafe,stroke:#2563eb,color:#111
    classDef step fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef opt fill:#fff7ed,stroke:#ea580c,color:#111,stroke-dasharray:5 4
    classDef gate fill:#fce7f3,stroke:#db2777,color:#111

    SS[/"«input» samplesheet"/]:::input
    CAT["cat_fastq_se / cat_fastq_pe<br/>cat runs [| gzip]<br/>if several runs, plain FASTQ, or «save_merged_fastq»"]:::opt
    L1["fq_lint raw<br/>fq lint «extra_fqlint_args»<br/>unless «skip_linting»"]:::opt
    F1["fastqc raw<br/>fastqc --threads --memory<br/>unless «skip_fastqc» or «skip_qc»"]:::opt
    TG["trimgalore_pe / trimgalore_se<br/>trim_galore --fastqc_args '-t N' «extra_trimgalore_args»<br/>--cores N [--paired] --gzip<br/>unless «skip_trimming»"]:::opt
    TS["trim_status<br/>reads left after trimming"]:::step
    G1{{"trim_summary checkpoint<br/>pass = reads ≥ «min_trimmed_reads»"}}:::gate
    L2["fq_lint trimmed"]:::opt
    BB["bbsplit_pe / bbsplit_se<br/>bbsplit.sh ambiguous2=all maxindel=150000<br/>keep primary reads; other genomes kept if «save_bbsplit_reads»<br/>unless «skip_bbsplit»"]:::opt
    L3["fq_lint bbsplit"]:::opt
    F2["fastqc filtered<br/>only when BBSplit ran"]:::opt
    SUB["fq_subsample<br/>fq subsample --record-count 1000000 --seed 1<br/>only for strandedness: auto"]:::opt
    SST["salmon_strand<br/>salmon quant --skipQuant → lib_format_counts.json"]:::opt
    STR["strandedness<br/>auto → forward / reverse / unstranded using<br/>«stranded_threshold» and «unstranded_threshold»<br/>(undetermined → unstranded)"]:::step

    SS --> CAT --> L1 --> TG --> TS --> G1
    CAT --> F1
    G1 --> L2 --> BB --> L3
    BB --> F2
    L3 --> SUB --> SST --> STR
    SS -->|"given strandedness"| STR
```

Each `fq_lint` report is an input of the next stage, so a lint failure stops that sample, as
in nf-core. With `skip_trimming` or `skip_bbsplit`, the reads skip that box and the gate.

### 2c · Alignment and quantification

```mermaid
flowchart TD
    classDef step fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef opt fill:#fff7ed,stroke:#ea580c,color:#111,stroke-dasharray:5 4
    classDef gate fill:#fce7f3,stroke:#db2777,color:#111
    classDef out fill:#dcfce7,stroke:#16a34a,color:#111

    RD[/"filtered reads + strandedness"/]
    SA["star_align<br/>STAR --quantMode TranscriptomeSAM --twopassMode Basic<br/>--outSAMtype BAM Unsorted --outFilterMultimapNmax 20 ...<br/>[--outReadsUnmapped Fastx if «save_unaligned»]<br/>«extra_star_align_args» (same flag overrides the default)<br/>[--sjdbGTFfile GTF unless «star_ignore_sjdbgtf»]<br/>--outSAMattrRGline ID SM [PL «seq_platform»] [CN «seq_center»]"]:::step
    MS["map_status<br/>% uniquely mapped from Log.final.out"]:::step
    G2{{"sample_summary checkpoint<br/>pass = % mapped ≥ «min_mapped_reads»"}}:::gate
    SO["samtools_sort → samtools_index_sorted<br/>samtools sort --reference FASTA"]:::step
    SST["samtools_stats / flagstat / idxstats<br/>on .sorted.bam and .markdup.sorted.bam"]:::out
    MD["picard_markduplicates → samtools_index_markdup<br/>picard MarkDuplicates --REMOVE_DUPLICATES false<br/>unless «skip_markduplicates»"]:::opt
    SQB["salmon_quant_bam<br/>salmon quant -t transcripts -a toTranscriptome.bam<br/>--libType «salmon_quant_libtype» or from strandedness<br/>«extra_salmon_quant_args»"]:::step
    SQR["salmon_quant_reads<br/>salmon quant --index -1 -2 --libType ... «extra_salmon_quant_args»<br/>if «pseudo_aligner» = salmon, unless «skip_pseudo_alignment»"]:::opt
    T2G["custom_tx2gene<br/>tx2gene.py id=«gtf_group_features» extra=«gtf_extra_attributes»"]:::step
    TXI["tximport<br/>tximport.r → salmon.merged.gene/transcript TPM, counts, lengths"]:::out
    SE["summarizedexperiment<br/>gene and transcript .rds"]:::out
    DQ["deseq2_qc<br/>deseq2_qc.r --count_col 3 [--vst TRUE if «deseq2_vst»]<br/>unless «skip_deseq2_qc» or «skip_qc»"]:::opt

    RD --> SA --> MS --> G2 --> SO --> MD
    SO & MD --> SST
    SA -->|"toTranscriptome.bam"| SQB
    RD --> SQR
    SQB -->|"star_salmon/"| T2G
    SQR -->|"salmon/"| T2G
    T2G --> TXI --> SE
    TXI --> DQ
```

tx2gene, tximport, SummarizedExperiment and DESeq2 QC run once per quantifier directory:
`star_salmon/` always, and also `salmon/` with pseudo-alignment.

### 2d · BAM QC, transcript assembly and coverage

```mermaid
flowchart TD
    classDef step fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef opt fill:#fff7ed,stroke:#ea580c,color:#111,stroke-dasharray:5 4
    classDef out fill:#dcfce7,stroke:#16a34a,color:#111

    FB[/"final BAM<br/>markdup.sorted.bam (or sorted.bam if «skip_markduplicates»)"/]

    STG["stringtie<br/>stringtie -e -G GTF -A -C -b ballgown [--fr | --rf] -v<br/>unless «skip_stringtie»"]:::opt
    PRE["preseq_lcextrap<br/>preseq lc_extrap -bam -seed 1 [-pe]<br/>only if «skip_preseq» is false"]:::opt
    FC["featurecounts<br/>featureCounts -B -C -g BIOTYPE -t «featurecounts_feature_type» -s 0/1/2<br/>BIOTYPE = gene_type if «gencode» else «featurecounts_group_type»<br/>unless «skip_biotype_qc» or biotype missing from GTF"]:::opt
    MCB["multiqc_custom_biotype<br/>mqc_features_stat.py"]:::opt
    QS["samtools_sort_qualimap<br/>samtools sort -n"]:::opt
    QM["qualimap_rnaseq<br/>qualimap rnaseq --sorted -p strand-specific-...<br/>unless «skip_qualimap»"]:::opt
    DR["dupradar<br/>dupradar.r<br/>unless «skip_dupradar»"]:::opt
    RQ["rseqc_* — one rule per module in «rseqc_modules»<br/>bam_stat · infer_experiment · inner_distance · junction_annotation<br/>junction_saturation · read_distribution · read_duplication · tin<br/>unless «skip_rseqc»"]:::opt
    GC["bedtools_genomecov<br/>bedtools genomecov -split -bg<br/>+ -du -strand ± for forward/reverse libraries"]:::opt
    BC["ucsc_bedclip<br/>bedClip"]:::opt
    BW["ucsc_bedgraphtobigwig<br/>bedGraphToBigWig → .bigWig (.forward/.reverse)<br/>unless «skip_bigwig»"]:::out
    MI["multiqc_rnaseq_inputs<br/>fail tables, strand check, sample merge config<br/>uses «min_trimmed_reads», thresholds"]:::step
    MQ["multiqc<br/>multiqc --config assets + «multiqc_config»<br/>[--title «multiqc_title»] [logo «multiqc_logo»]<br/>unless «skip_multiqc»"]:::out

    FB --> STG & PRE & FC & QS & DR & RQ & GC
    FC --> MCB
    QS --> QM
    GC --> BC --> BW
    MCB & QM & DR & RQ & PRE & MI --> MQ
```

`skip_qc: true` turns off every QC box here (and FastQC and DESeq2 QC), but not StringTie or bigWigs.

---

## 3. How config switches shape the run

Colours: yellow = config key, blue = value derived in the Snakefile, grey = rules affected.

### 3a · Which reference files are built

```mermaid
flowchart LR
    classDef key fill:#fef9c3,stroke:#ca8a04,color:#111
    classDef rule fill:#f3f4f6,stroke:#6b7280,color:#111

    A{"«gtf» set?"}:::key -->|yes| A1["use it (gunzip if .gz)"]:::rule
    A -->|no| A2["gffread_gtf from «gff»"]:::rule
    B{"«skip_gtf_filter»?"}:::key -->|no| B1["custom_gtffilter"]:::rule
    C{"«additional_fasta»?"}:::key -->|yes| C1["custom_catadditionalfasta<br/>new FASTA + GTF used everywhere"]:::rule
    D{"«transcript_fasta»?"}:::key -->|yes| D1["use it; + gencode header fix if «gencode»"]:::rule
    D -->|no| E{"«gffread_transcript_fasta»?"}:::key
    E -->|yes| E1["gffread_transcripts"]:::rule
    E -->|no| E2["make_transcripts_fasta (RSEM)"]:::rule
    F{"«gene_bed»?"}:::key -->|no| F1["eautils_gtf2bed"]:::rule
    I{"«star_index» / «salmon_index» / «bbsplit_index»"}:::key -->|"empty"| I1["build the index"]:::rule
    I -->|".tar.gz"| I2["untar"]:::rule
    I -->|"directory"| I3["use as is"]:::rule
```

### 3b · Which samples and read stages go forward

```mermaid
flowchart TD
    classDef key fill:#fef9c3,stroke:#ca8a04,color:#111
    classDef rule fill:#f3f4f6,stroke:#6b7280,color:#111
    classDef gate fill:#fce7f3,stroke:#db2777,color:#111

    R0["raw / merged reads"]:::rule
    K1{"«skip_trimming»?"}:::key
    TG["Trim Galore + gate «min_trimmed_reads»"]:::gate
    K2{"«skip_bbsplit»?<br/>(default true)"}:::key
    BB["BBSplit with «bbsplit_fasta_list» or «bbsplit_index»"]:::rule
    K3{"sample strandedness"}:::key
    SA["Salmon inference: «stranded_threshold» / «unstranded_threshold»"]:::rule
    ST["strandedness used by Salmon --libType, StringTie,<br/>featureCounts -s, Qualimap -p, dupRadar, bigWig strands"]:::rule
    K4{"«pseudo_aligner» = salmon<br/>and not «skip_pseudo_alignment»?"}:::key
    PS["also run salmon_quant_reads → outdir/salmon/"]:::rule
    STAR["STAR + gate «min_mapped_reads»"]:::gate
    K5{"«skip_markduplicates»?"}:::key
    FB1["final BAM = markdup.sorted.bam"]:::rule
    FB2["final BAM = sorted.bam"]:::rule

    R0 --> K1
    K1 -->|no| TG --> K2
    K1 -->|yes| K2
    K2 -->|no| BB --> K3
    K2 -->|yes| K3
    K3 -->|auto| SA --> ST
    K3 -->|"forward / reverse / unstranded"| ST
    ST --> STAR --> K5
    ST --> K4 -->|yes| PS
    K5 -->|no| FB1
    K5 -->|yes| FB2
```

### 3c · Which QC tools run

```mermaid
flowchart LR
    classDef key fill:#fef9c3,stroke:#ca8a04,color:#111
    classDef rule fill:#f3f4f6,stroke:#6b7280,color:#111

    Q{"«skip_qc»?"}:::key -->|yes| NONE["no FastQC, BAM QC or DESeq2 QC"]:::rule
    Q -->|no| T["QC tool set"]:::rule
    T --> P["«skip_preseq» false → preseq<br/>(off by default)"]:::rule
    T --> BQ["«skip_biotype_qc» false and biotype present → featureCounts"]:::rule
    T --> QM["«skip_qualimap» false → qualimap"]:::rule
    T --> DR["«skip_dupradar» false → dupRadar"]:::rule
    T --> RS["«skip_rseqc» false → each module in «rseqc_modules»"]:::rule
    T --> DE["«skip_deseq2_qc» false → deseq2_qc"]:::rule
    T --> FQ["«skip_fastqc» false → fastqc raw + filtered"]:::rule
    O["Independent of skip_qc"]:::key --> S1["«skip_stringtie» → stringtie"]:::rule
    O --> S2["«skip_bigwig» → coverage tracks"]:::rule
    O --> S3["«skip_linting» → fq_lint"]:::rule
    O --> S4["«skip_multiqc» → multiqc"]:::rule
```

---

## 4. Config key reference

| Key | Rule(s) | Effect on the command |
|---|---|---|
| `input` | all per-sample rules, summarizedexperiment | Samplesheet. Rows with the same `sample` are runs, concatenated by `cat_fastq_*`. `strandedness: auto` triggers Salmon inference. |
| `outdir` / `workdir` | all | Published files / intermediates (nf-core `work/`). |
| `fasta` | link_fasta, samtools_faidx, indices, sort/stats | Genome; `.gz` is unzipped first. Required. |
| `gtf` / `gff` | gffread_gtf, custom_gtffilter, everything using the GTF | One is required. GFF is converted with `gffread --keep-exon-attrs -F -T`. |
| `skip_gtf_filter` | custom_gtffilter | Turns the filter off. |
| `skip_gtf_transcript_filter` | custom_gtffilter | Adds `--skip_transcript_id_check`. |
| `additional_fasta` | custom_catadditionalfasta | Appends spike-ins (e.g. ERCC) to the genome and GTF. |
| `transcript_fasta` | salmon_index, salmon_quant_bam | Given transcriptome; otherwise built by RSEM or gffread. |
| `gffread_transcript_fasta` | gffread_transcripts | Build transcripts with gffread instead of `rsem-prepare-reference`. |
| `gencode` | salmon_index, preprocess_transcripts_fasta_gencode, featurecounts | `--gencode` for the index, keeps only the first pipe-separated field of transcript headers, biotype attribute becomes `gene_type`. |
| `gene_bed` | RSeQC rules | Given BED12; otherwise `eautils_gtf2bed`. |
| `star_index`, `salmon_index`, `bbsplit_index` | star_align, salmon_*, bbsplit_* | Directory → used; `.tar.gz` → `untar`; empty → built. |
| `bbsplit_fasta_list` | bbsplit_index | CSV `name,fasta`; each becomes `ref_<name>=` for `bbsplit.sh build=1`. |
| `save_reference` | reference rules | Write `genome/` and `genome/index/` under `outdir`. |
| `gtf_group_features`, `gtf_extra_attributes` | custom_tx2gene | `id` and `extra` columns of the tx2gene table. |
| `featurecounts_feature_type`, `featurecounts_group_type` | featurecounts, custom_catadditionalfasta | `featureCounts -t` and `-g` (biotype QC). |
| `skip_linting`, `extra_fqlint_args` | fq_lint | Turn lint off / arguments to `fq lint`. |
| `trimmer`, `skip_trimming` | trimgalore_* | Only `trimgalore` is ported; skip passes reads straight through. |
| `extra_trimgalore_args` | trimgalore_* | Merged into the Trim Galore arguments; R2 options are dropped for single-end. |
| `min_trimmed_reads` | trim_summary, multiqc_rnaseq_inputs | Samples with fewer reads after trimming are dropped. |
| `save_trimmed`, `save_merged_fastq` | trimgalore_*, cat_fastq_* | Publish trimmed / merged FASTQ. `save_merged_fastq` also forces `cat_fastq` for single-run samples. |
| `skip_bbsplit`, `save_bbsplit_reads` | bbsplit_* | Run BBSplit (default off) / keep the reads assigned to the other genomes. |
| `stranded_threshold`, `unstranded_threshold` | strandedness | Salmon-based calls for `auto` samples. |
| `aligner` | — | Only `star_salmon` is ported (others are rejected at start-up). |
| `extra_star_align_args` | star_align | Merged with the default STAR flags; the same flag replaces the default. |
| `star_ignore_sjdbgtf` | star_align | Leaves out `--sjdbGTFfile`. |
| `seq_platform`, `seq_center` | star_align | `PL:` / `CN:` in `--outSAMattrRGline` (samplesheet columns win). |
| `save_unaligned` | star_align | `--outReadsUnmapped Fastx`; unmapped reads go to `star_salmon/unmapped/`. |
| `save_align_intermeds` | star_align, samtools_sort | Publish the unsorted STAR BAMs and the sorted BAM. |
| `min_mapped_reads` | sample_summary | Samples below this % uniquely mapped skip MarkDuplicates, StringTie, BAM QC and bigWigs. |
| `skip_markduplicates` | picard_markduplicates | Final BAM becomes `sorted.bam`. |
| `salmon_quant_libtype` | salmon_quant_* | Fixed `--libType`; empty → from strandedness (ISF/ISR/IU or SF/SR/U). |
| `extra_salmon_quant_args` | salmon_quant_* | Appended to `salmon quant`. |
| `pseudo_aligner`, `skip_pseudo_alignment`, `pseudo_aligner_kmer_size`, `extra_salmon_index_args` | salmon_quant_reads, salmon_index | Run Salmon on reads too; `-k` and extra flags for `salmon index`. |
| `skip_stringtie` | stringtie | Turns StringTie off. |
| `skip_qc` and `skip_fastqc`, `skip_preseq`, `skip_biotype_qc`, `skip_qualimap`, `skip_dupradar`, `skip_rseqc`, `skip_deseq2_qc` | QC rules | Turn off each QC tool (or all of them). Preseq is off by default. |
| `rseqc_modules` | rseqc_* | Comma list; one rule per module. |
| `deseq2_vst` | deseq2_qc | `--vst TRUE` (otherwise rlog). |
| `skip_bigwig` | bedtools_genomecov, ucsc_* | Turns the coverage tracks off. |
| `skip_multiqc`, `multiqc_config`, `multiqc_title`, `multiqc_logo` | multiqc | Report on/off, extra `--config`, `--title`, custom logo. |
| `resources.<label>` | every rule | `threads`/`mem_mb` per nf-core label. Memory also sets STAR `--limitGenomeGenerateRAM`, BBSplit/Picard `-Xmx`, Qualimap and FastQC memory. |
| `with_umi`, `remove_ribo_rna`, `skip_alignment`, `skip_quantification_merge`, `stringtie_ignore_gtf`, `contaminant_screening`, `use_rustqc`, `bam_csi_index` | — | **Not ported.** Must stay at their defaults; the Snakefile refuses to start otherwise. |

---

## 5. Exact rule graph

To regenerate the real rule DAG for a config (no jobs are run), from inside `snk-rnaseq/`:

```bash
snakemake --configfile config/test.yaml --rulegraph mermaid-js > rulegraph.mmd
```

Rules after the `trim_summary` and `sample_summary` checkpoints only appear once those
checkpoints have run (for example after a real test run).
