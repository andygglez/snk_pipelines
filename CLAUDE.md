# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

This repository ports nf-core Nextflow pipelines to Snakemake, **one pipeline at a time**. For each pipeline:

1. Download the nf-core pipeline source.
2. Download every Singularity image it uses into the port's own `containers/` folder.
3. Reimplement the pipeline's functionality in Snakemake. Every step must run from the downloaded Singularity images.

Finish one pipeline before you start the next. The nf-core source is the reference implementation. Read its `main.nf`, `workflows/`, `subworkflows/`, `modules/`, `conf/modules.config`, `nextflow.config` and `bin/` to find the exact commands, arguments and container used at each step.

## Directory layout

Use `hic` as the example pipeline:

```
pipelines/
├── nfcore-hic/               # Downloaded nf-core/hic source (reference only; do not edit)
└── snk-hic/                  # Snakemake port
    ├── Snakefile             # All rules, declared explicitly
    ├── config/
    │   ├── config.yaml       # Every parameter that controls the commands in Snakefile
    │   └── containers.yaml   # Path to the Singularity image used by each rule
    ├── containers/           # This pipeline's Singularity images (.img / .sif)
    └── scripts/              # Custom scripts called by rules
```

Naming: `nfcore-<pipeline>` holds the source and `snk-<pipeline>` holds the port. Both sit at the top level. Each port keeps its own images in `snk-<pipeline>/containers/`. Images are not shared between pipelines.

## Downloading a pipeline

Use the nf-core tools CLI, sending the images straight to the port's `containers/` folder:

```bash
mkdir -p snk-<pipeline>/containers
export NXF_SINGULARITY_CACHEDIR=$PWD/snk-<pipeline>/containers
nf-core pipelines download nf-core/<pipeline> \
    --revision <version> \
    --outdir nfcore-<pipeline> \
    --compress none \
    --container-system singularity \
    --container-cache-utilisation amend
```

Note the downloaded revision in the header of `snk-<pipeline>/config/config.yaml` so the port stays pinned to it.

## Snakemake port conventions

- **One `Snakefile`.** Declare every rule explicitly in `snk-<pipeline>/Snakefile`. Don't split rules into `rules/*.smk` includes.
- **All parameters in `config/config.yaml`.** This covers tool arguments, thresholds, reference paths, sample sheets, resolutions, threads and memory. Don't hardcode values in the Snakefile that nf-core exposes as `params` or sets in `conf/modules.config`.
- **Container paths in `config/containers.yaml`.** Keep container paths in this file, not in `config.yaml`. Map each tool or rule to its image path, which defaults to the file in `containers/`. Users can repoint images by editing this file alone. Every rule sets a `container:` directive read from this mapping, for example:

  ```yaml
  # config/containers.yaml
  bowtie2: containers/depot.galaxyproject.org-singularity-bowtie2-2.4.4--py39hbb4e92a_0.img
  ```
  ```python
  # Snakefile
  configfile: "config/config.yaml"
  configfile: "config/containers.yaml"
  ...
  rule bowtie2_align:
      container: config["bowtie2"]
  ```

  Use the same image that the nf-core module uses for that step.
- **Binaries vs. scripts:**
  - Call standard tools (STAR, bowtie2, samtools, cooler, MultiQC, etc.) directly in the rule's `shell:` block, using the binary inside the rule's container. Don't wrap them in scripts.
  - Custom logic goes in `snk-<pipeline>/scripts/`, for example nf-core `bin/` scripts or bespoke formatting such as building Hi-C pileups in a specific format. Call it from the rule with `script:` or `shell:`, and run it inside the container that provides its interpreter and dependencies.
- **Functional parity.** Each nf-core process should map to an equivalent rule with the same command semantics and outputs. When a step is skipped or changed, record that in a comment in the Snakefile.

## Running a port

`pipeline.py` (symlinked as `~/.local/bin/pipeline`) wraps the commands below: `pipeline run <name> -c my.yaml -t N [-n] [-B dir] [-- snakemake args]` runs from the current directory, passes `config.yaml` + `containers.yaml` + the user config as `--configfile` (later wins), and refuses to start while an image is missing. Pipelines are named rnaseq, chipseq, atacseq, hic and wgbs (wgbs = `snk-methylseq`; mapping in `NAMES`). `pipeline list` shows each pipeline with an image-download progress bar, `pipeline config <name> -p <dir>` copies the default config, and `pipeline pull <name>` downloads missing images (one at a time) from the `container_urls:` section of `containers.yaml`. Pull on a login node, because compute nodes are offline. When adding a port, fill `container_urls:` with the nf-core module URLs. `pipeline run` also binds automatically the workdir, the port dir, the directories of the path-valued config keys listed per port in `BIND_KEYS` (plus the files listed in `MANIFEST_KEYS`), and the samplesheet's `fastq_1`/`fastq_2` paths, symlink targets included; `-B` adds extras, `--no-auto-bind` disables it. When adding a port or a new path parameter, add its keys to `BIND_KEYS` in `pipeline.py`.

Run from inside `snk-<pipeline>/`, binding any paths the containers need. Both config files are loaded by the Snakefile's `configfile:` directives. To override either one, pass `--configfile`.

```bash
# Snakemake >= 8
snakemake --cores <N> \
    --software-deployment-method apptainer \
    --apptainer-args "--bind <data_dir>,<ref_dir>"

# Snakemake 7
snakemake --cores <N> --use-singularity \
    --singularity-args "--bind <data_dir>,<ref_dir>"

# Dry run / DAG check
snakemake -n
```

## Testing on this machine

This laptop is small (4 cores, ~10 GB RAM). **Don't run heavy computations locally, or the laptop will crash.** When you test locally:

- Prefer `snakemake -n` (dry runs) and DAG/lint checks.
- Real runs must use tiny test data only (for example, a few thousand read pairs and a single small chromosome), with `--cores 1` or `--cores 2` and low memory settings.
- Don't build full-genome indexes, align full datasets, or pull and convert large images in parallel.
- Send anything heavier to the cluster (see the slurmer agent and the remote-slurm skill).
- **Local container runs don't work.** Ubuntu 26.04 restricts unprivileged user namespaces through AppArmor (`kernel.apparmor_restrict_unprivileged_userns = 1`), and the conda Apptainer has no AppArmor profile. Every `apptainer exec` fails with `Could not write info to setgroups: Permission denied`. Locally, stick to dry runs and lint, and run every real test on Rorqual (see below). Don't change AppArmor or sysctl settings yourself; that's the user's call.
- Snakemake and Apptainer live in the `nfcore-snk` conda env (`/home/andy/programs/anaconda3/envs/nfcore-snk/bin`); they aren't on the default PATH.
- For a local test run, cap memory so a runaway job can't take the laptop down: `systemd-run --user --scope -p MemoryMax=5G -p MemorySwapMax=0 nice -n 19 snakemake --cores 2 ...`

## Test profile for each port

- Each port has a `config/test.yaml` that mirrors the nf-core `conf/test.config`, with the same test-dataset URLs, the same params and small `resources`. Run it with `--configfile config/test.yaml`, so outputs go to `results_test/`.
- Remote inputs (http/ftp) are fetched by a `fetch_remote` rule. Compute nodes are offline, so the cluster test uses copies of the config (`config/test_local.yaml`, `config/samplesheet_local.csv`) that point at pre-downloaded files with absolute paths. Leave the originals unchanged. The job count drops by the number of `fetch_remote` jobs (snk-hic: 41 → 38).
- Match nf-core's output layout (directory names and file names under `outdir`), so results can be compared file by file.

## Running tests on Rorqual (slurmer agent)

Delegate every real run to the `slurmer` agent (remote-slurm skill). Some steps need the user first:

1. **SSH socket.** If `remote_ctl.sh status` exits with code 1, the ControlMaster socket is dead. Ask the user to run `ssh andygg98@rorqual.alliancecan.ca`, approve Duo, and keep that session open.
2. **Allocation.** If it exits with code 2, there's no `claude_session` job. Ask the user to start one manually, for example `salloc --job-name=claude_session --account=def-majewski --cpus-per-task=8 --mem=24G --time=2:00:00`, and wait until it's RUNNING. Exit code 3 means it's still PENDING.

Cluster setup that worked for snk-hic:

- Layout: everything sits under `/scratch/andygg98/claude_jobs/`, with the port in `snk-<pipeline>/`, the nf-core source in `nfcore-<pipeline>/`, and logs in `logs/job_<tag>.log`. Copy the port without `.snakemake/` and `results*/`.
- Compute nodes have **no internet**. Download test data and pip wheels locally, then rsync them over.
- Snakemake: the venv is `/scratch/andygg98/claude_jobs/venv_snk` with `snakemake==9.27.0`, installed from locally downloaded wheels (`pip install --no-index --find-links wheels`), because the Alliance wheelhouse only goes up to 9.17.1. Reuse that venv.
- Apptainer: `module load apptainer` (1.3.5, StdEnv/2023). Add `--apptainer-args "--bind /scratch/andygg98/claude_jobs/<dir>"`.
- `remote_ctl.sh` quirks:
  - `deploy` takes one path and rsyncs it into `claude_jobs/`, so stage a directory locally and deploy that.
  - `run-async` runs `bash <script>` from `claude_jobs/`, so launch workflows through a small wrapper `.sh`.
  - `run-async` hangs and never prints its TAG, but the job does launch. Run it in the background and find the tag from the `logs/job_*.log` name.

## Validating a port against nf-core

A port counts as validated once its test-profile outputs match an nf-core run of the pinned revision on the same data. snk-hic (2.1.0) passed this check: every science output was identical.

- **nf-core reference run** (on Rorqual, offline):
  - `module load nextflow` (24.04.4 worked) and `module load apptainer`, with `NXF_OFFLINE=true`.
  - Point `NXF_SINGULARITY_CACHEDIR` and `NXF_APPTAINER_CACHEDIR` at the port's `containers/` to reuse the same images.
  - Use `-profile test,singularity` with `--input`/`--fasta` pointing at the local test data.
  - Put a small `-c` config outside the pipeline dir, setting `singularity.runOptions = '--bind /scratch/andygg98'` and `singularity.envWhitelist = 'MPLCONFIGDIR'`.
  - Never edit the nf-core source.
- **Matplotlib warning:** containers with Matplotlib (for example cooltools) print a warning when `HOME` isn't writable. In nf-core, that text lands in `versions.yml` and breaks `CUSTOM_DUMPSOFTWAREVERSIONS`. Fix it by exporting `MPLCONFIGDIR=/tmp/mplcfg` and whitelisting that variable.
- **What to compare**, giving a verdict of identical or differs for each:
  - every stats/QC text file (`diff`)
  - the main tabular outputs (line counts plus a sorted diff)
  - binary matrices (for example `cooler info` nnz/sum/nbins, pixels, and balanced weights compared by max abs difference)
  - every downstream text output
  - the file inventory, meaning files that are only in one tree or are named differently

## Port status

- **snk-hic** (nf-core/hic 2.1.0): validated on Rorqual on 2026-10-02, with outputs identical to nf-core. Still open: a few output names and paths differ from nf-core:
  - distance decay: `<sample>.<res>_distcount.*` → `<sample>_distcount.*`
  - FastQC: `<sample>_0_{1,2}_fastqc.*` → `<sample>_{1,2}_fastqc.*`
  - mapped-pairs BAM: `hicpro/mapping/<sample>/` → `hicpro/mapping/`
  - validated samplesheet: `pipeline_info/samplesheet.valid.csv` → `samplesheet/samplesheet.valid.csv`
- **snk-rnaseq** (nf-core/rnaseq 3.27.0): validated on Rorqual on 2026-10-02 against an offline nf-core run (nextflow/26.04.4; the pipeline needs ≥25.10.4; the nf-schema 2.7.2 plugin comes from registry.nextflow.io and goes into `NXF_PLUGINS_DIR`). All quantification and count tables, bigWigs and samtools stats are identical. Remaining differences are run-to-run noise: BAM tie order (STAR threads), BBSplit thread nondeterminism, and RSeQC junction_saturation sampling. Still open: no software-versions output (`pipeline_info`, MultiQC versions table), and several options are not ported yet (the Snakefile rejects them at start-up). Intermediates go to `workdir` (`work/` or `work_test/`), which also has to be left out when copying to the cluster.
- **snk-chipseq** (nf-core/chipseq 2.1.0): validated on Rorqual on 2026-10-02 against an offline nf-core run. All four aligners are ported; iGenomes is not. Tested against nf-core: bwa and bowtie2. chromap and STAR are only dry-run checked (not to be tested on ChIP-seq data, per the user). The file inventory is identical (644 files). All BAMs (records), bigWigs, MACS3 peaks/xls, consensus BED/SAF, samtools stats, FRiP, SPP and featureCounts counts are identical. The remaining differences are timestamps (Picard, Trim Galore, FastQC), row or column order, float noise of about 1e-15, and Perl-hash tie-breaking: HOMER picks between equally near genes in a few rows, and Trim Galore's runner-up adapter differs when both counts are 0. Row/column order differs in computeMatrix, IGV, macs3_peak.summary and plotFingerprint (raw.txt rows only). In featureCounts, nf-core orders the columns by work-dir path and the port orders them by sample name.
  - Extra scenarios, 2026-10-02 (`config/test_{A,B,C,D}.yaml` overlay `test_local.yaml`; nf-core uses `nfcore-chipseq/params_X.yaml`; compared with `compare_scenario.py`). Every scenario has the same file inventory as nf-core and no unexplained differences:
    - A: bowtie2, narrow peaks, all `save_*` options. Intermediate BAMs and unmapped FASTQs have the same records in a different order (multithreading); the final BAMs are identical.
    - B: single-end, gzipped FASTA and GTF, a blacklist, a pre-built bwa index directory, preseq.
    - D: a `.tar.gz` bwa index, a replicate with two runs, a library below `min_trimmed_reads` (dropped), a single-replicate antibody (no consensus) and a 0-byte peak file. All behave as in nf-core.
    - C: as D but with a GFF. It runs in the port, but nf-core 2.1.0 crashes with `--gff` (GFFREAD is called with 1 of its 2 inputs), so it can't be compared.
    - preseq: on the small INPUT test libraries it fails identically in both pipelines ("max count before zero ..."), so A runs with preseq off.
    - DESeq2: the PCA signs can flip, because the featureCounts column order differs; the values are identical by sample name.
  - Cluster quirks:
    - Lmod exports a `which` shell function that breaks `which` inside containers (phantompeakqualtools). The Snakefile's `shell.prefix` and the nf-core wrapper both `unset -f which`.
    - The `nextflow` module is a self-contained build that ignores non-core plugins. The nf-core run uses the standard 24.04.4 launcher at `claude_jobs/nfcore-chipseq/nextflow`, with nf-validation 1.1.3 in `nxf_home/plugins`.
  - The images are on Rorqual only; the local `snk-chipseq/containers/` is empty.
- **snk-atacseq** (nf-core/atacseq 2.1.2): validated on Rorqual on 2026-10-02 against offline nf-core runs, in three scenarios: `test`, `test_controls` (as in nf-core, this profile leaves `with_control` false) and `with_control` (`test_controls` + `--with_control`; config `test_with_control.yaml`). Every scenario has the same file inventory as nf-core (727 / 749 / 717 files).
  - Identical: all BAMs, MACS2 peaks/xls/gappedPeak, consensus BED/SAF/boolean, bigWigs and scale factors, FRiP and peak counts, samtools stats, Picard metrics and ataqv JSON (apart from timestamps), and featureCounts counts per sample.
  - Remaining differences:
    - timestamps: Picard, ataqv, Trim Galore, FastQC
    - row order: computeMatrix, plotFingerprint raw.txt, macs2_peak summary, IGV list
    - column order: featureCounts and the HOMER summary plot table. nf-core orders featureCounts by task completion, so the DESeq2 PCA signs flip and the distance matrix is permuted.
    - HOMER: equally near genes are tie-broken differently in 1–7 rows per file
    - merged-replicate BAMs: the same records, but `PG:Z:MarkDuplicates.A/.B` can swap. nf-core sorts MergeSamFiles inputs by staged path, so the order is effectively random.
  - Not ported: iGenomes and a `.tar.gz` chromap index (nf-core 2.1.2 itself crashes on it: `UNTAR.out` is undefined). Software versions are not produced. Unlike nf-core, a preseq failure is not ignored (skip_preseq defaults to true).
  - Snakemake 9: params that read a checkpoint must list the checkpoint's output as an input. They also must not read `input.<name>` from an unpack()/checkpoint input function (they raise AttributeError at run time), so compute those params from the wildcards instead.
  - nf-core reference run: atacseq asks for nf-validation `latest`, so pin `nf-validation@1.1.3` in the `-c` config to run offline. It also needs the multiqc 1.14 image (CUSTOM_DUMPSOFTWAREVERSIONS).
  - Images (34, including multiqc 1.14 for the nf-core run) are on Rorqual only; the local `snk-atacseq/containers/` is empty. Test data is reused from `claude_jobs/chipseq_testdata`.
- **snk-methylseq** (nf-core/methylseq 4.2.0): validated on Rorqual on 2026-10-03 against offline nf-core runs: the `test` profile plus nine scenarios (`config/scenarios/{A..I}.yaml` overlay `test_local.yaml`; nf-core uses `nfcore-methylseq/params_X.yaml`; compared with `compare_methylseq.py`):
  - A: bwameth, Qualimap, Preseq. B: bwamem + TAPS (rastair), save_reference, save_align_intermeds. C: bismark_hisat, cytosine_report, unmapped. D: skip_deduplication, nomeseq, Qualimap, save_reference/intermeds. E: rrbs, em_seq, clip_r1, comprehensive. F: bwameth targeted sequencing, all_contexts, collecthsmetrics. G: replicate runs (CAT_FASTQ), skip_trimming, bismark targeted. H: Bowtie2 index `.tar.gz`, Preseq. I: bwameth index `.tar.gz`, rrbs.
  - Every scenario has the same file inventory as nf-core, apart from `multiqc_software_versions.txt`. All BAMs, BAIs, methylation calls/coverage/bedGraphs, M-bias, splitting reports, coverage2cytosine, MethylDackel, rastair, samtools stats, Qualimap data and indices are identical.
  - Remaining differences: absolute paths (Bismark alignment reports, MultiQC sources, picard `.dict`/interval list `UR:`), timestamps and run times (Bismark HTML reports, Picard metrics, Trim Galore, FastQC zip entries, Qualimap HTML), Perl hash order in Trim Galore's adapter-detection line, sample order in the Bismark summary (nf-core collects in completion order), and the missing MultiQC versions/methods sections.
  - nf-core quirks kept on purpose (documented in the Snakefile header): PICARD_COLLECTHSMETRICS never runs (empty ch_gzi); coverage2cytosine is published under `bismark/` even for bismark_hisat; targeted bedGraphs go to `methyldackel/` for every aligner except bismark; Qualimap and BWA_INDEX publish their `versions.yml`; a failing Preseq is ignored (it fails on the tiny test libraries in both pipelines); `.tar.gz` indices are published to `untar/`.
  - MultiQC plot export (Kaleido) times out in the multiqc 1.32 container on Rorqual in both pipelines, so `multiqc_plots/` stays empty in both.
  - Not ported: iGenomes, Parabricks GPU aligner, software versions / pipeline_info.
  - nf-core reference run: needs Nextflow ≥ 25.04 (`module load nextflow/26.04.4`) with the nf-schema 2.5.1 plugin in `nfcore-methylseq/plugins` (`NXF_PLUGINS_DIR`). Set `NXF_SYNTAX_PARSER=v1`, because the 4.2.0 release includes `conf/aws/batch/nextflow.config`, which it doesn't ship. Pass `--igenomes_ignore --igenomes_base <local dir>`, because nf-schema can't check the S3 path offline. Set `executor.cpus`/`memory` in the `-c` config, because the srun step otherwise exposes 1 CPU. Wrappers on scratch: `run_methylseq_test.sh`, `run_nfcore_methylseq_test.sh`, `run_methylseq_scenarios.sh`, `run_compare_methylseq_all.sh`.
  - Images (14) and test data (`claude_jobs/methylseq_testdata`) are on Rorqual only; the local `snk-methylseq/containers/` is empty.
