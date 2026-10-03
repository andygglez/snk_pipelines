# Snakemake ports of nf-core pipelines

This repository reimplements nf-core Nextflow pipelines in Snakemake. Every step
runs inside the same Singularity image that the nf-core module uses, so the
outputs can be compared with nf-core file by file.

| Pipeline  | Port folder      | nf-core reference              | Analysis                         |
|-----------|------------------|--------------------------------|----------------------------------|
| `rnaseq`  | `snk-rnaseq/`    | nf-core/rnaseq 3.27.0          | RNA-seq (STAR + Salmon)          |
| `chipseq` | `snk-chipseq/`   | nf-core/chipseq 2.1.0          | ChIP-seq peaks                   |
| `atacseq` | `snk-atacseq/`   | nf-core/atacseq 2.1.2          | ATAC-seq peaks                   |
| `hic`     | `snk-hic/`       | nf-core/hic 2.1.0              | Hi-C contact maps                |
| `wgbs`    | `snk-methylseq/` | nf-core/methylseq 4.2.0        | Bisulfite / EM-seq methylation   |

The `nfcore-<pipeline>/` folders hold the downloaded nf-core source. They are
reference only and are never edited.

## Contents

1. [Setup](#setup)
2. [Workflow at a glance](#workflow-at-a-glance)
3. [The `pipeline` command](#the-pipeline-command)
4. [Container specifications](#container-specifications)
5. [Configuring a pipeline](#configuring-a-pipeline)
6. [Running each pipeline](#running-each-pipeline)
7. [Running on an offline cluster (Rorqual)](#running-on-an-offline-cluster-rorqual)
8. [Repository layout](#repository-layout)

---

## Setup

You need:

- **Snakemake ≥ 8** (tested with 9.27.0) and **Apptainer** (or Singularity).
- **Python 3 with PyYAML** to run `pipeline.py`. The Python in Snakemake's
  environment already has it.

Make `pipeline` callable from any directory:

```bash
ln -s /path/to/pipelines/pipeline.py ~/.local/bin/pipeline   # ~/.local/bin must be on PATH
```

`pipeline` finds Snakemake on your `PATH`, so activate its environment first.
If it lives somewhere else, set `PIPELINE_SNAKEMAKE=/path/to/snakemake` or pass
`--snakemake`.

---

## Workflow at a glance

```
1. pipeline list                        which pipelines exist, which images are downloaded
2. pipeline pull <name>                 download missing images     (needs internet)
3. pipeline config <name> -p <dir>      copy the default config.yaml
4. edit <dir>/<name>.config.yaml        set input, outdir, references, options
5. pipeline run <name> -c <cfg> -n      dry run: check the job graph
6. pipeline run <name> -c <cfg> -t N    real run                    (works offline)
```

Steps 1–2 only need doing once per machine or cluster. Steps 3–6 are repeated
for each analysis. Run step 2 on a machine with internet access, such as a
cluster login node. After that the pipeline never needs the network for its
images.

---

## The `pipeline` command

`pipeline.py` is a thin front end over Snakemake. It always takes the
pipeline's official name: `rnaseq`, `chipseq`, `atacseq`, `hic` or `wgbs`.
(`methylseq` also works as an alias for `wgbs`.)

### `pipeline list`

Lists the available pipelines, each with a bar showing how many of its images
have been downloaded:

```
rnaseq
  ████████████████████████████████  30/30 images  100%

chipseq
  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░   0/34 images    0%
...
```

A pipeline can only run when its bar is full (100%).

### `pipeline pull <name>`

Downloads every image of `<name>` that is missing, one at a time, into
`snk-<folder>/containers/`. Images already present are skipped, so it is safe
to re-run after an interrupted download.

```bash
pipeline pull rnaseq           # download what is missing
pipeline pull rnaseq -n        # only list what would be downloaded
pipeline pull rnaseq -c my.yaml  # use image paths overridden in my.yaml
```

Downloads go to a `.part` file first and are only renamed when complete, so a
failed download never leaves a broken image behind.

### `pipeline config <name> -p <path>`

Copies the pipeline's default `config.yaml`, with every option documented, so
you can edit it.

```bash
pipeline config rnaseq -p ~/projects/exp1/          # -> ~/projects/exp1/rnaseq.config.yaml
pipeline config hic -p ~/projects/exp2/local.hic.yaml   # explicit file name
pipeline config wgbs                                # -> ./wgbs.config.yaml
```

It refuses to overwrite an existing file unless you pass `-f/--force`.

### `pipeline run <name> -c <config> [options]`

Runs the pipeline with your config.

```bash
pipeline run rnaseq -c rnaseq.config.yaml -t 12
```

| Option | Meaning |
|---|---|
| `-c, --configfile` | Your config YAML (required). |
| `-t, --threads` | Total cores Snakemake may use (`--cores`). Default 1. |
| `-d, --workdir` | Directory to run in. Default: the current directory. |
| `-B, --bind` | Extra host path to make visible inside the containers. Repeatable. |
| `-n, --dry-run` | Show the jobs without running anything. |
| `--no-check` | Start even if some images are missing. |
| `--snakemake` | Path to the snakemake executable. |
| `-- <args>` | Anything after `--` is passed to Snakemake unchanged. |

What `run` does:

1. **Checks the images.** If any image is missing, it stops before Snakemake
   starts and tells you to run `pipeline pull`. A dry run only warns.
2. **Runs in your directory, not in the repository.** Outputs (`outdir`,
   `workdir`, `.snakemake/`) are written to the run directory, and relative
   paths in your config are resolved from there.
3. **Layers the configs.** It passes, in order, the port's `config.yaml`, its
   `containers.yaml`, then your file. Later files win, so your config only
   needs the keys you want to change.
4. **Uses Apptainer.** It runs Snakemake with `--sdm apptainer` and binds the
   run directory, plus every `-B` path, into the containers. Data or references
   that live outside the run directory must be bound with `-B`.

It prints the full Snakemake command before running it, so the call can be
reproduced by hand.

Examples:

```bash
pipeline run chipseq -c chipseq.config.yaml -t 8 -B /scratch/me/genomes
pipeline run hic -c hic.config.yaml -n                      # dry run
pipeline run wgbs -c wgbs.config.yaml -t 4 -- --rerun-incomplete --keep-going
```

---

## Container specifications

Each pipeline lists its images in **`snk-<folder>/config/containers.yaml`**.
This file has two sections.

**`containers:`** gives the image each rule runs in. Paths are relative to the
port folder (or absolute):

```yaml
containers:
  # FASTQC
  fastqc: containers/depot.galaxyproject.org-singularity-fastqc-0.12.1--hdfd78af_0.img
  # STAR_GENOMEGENERATE, STAR_ALIGN
  star: containers/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-26-...-data.img
```

Each key is used by one or more rules through a `container:` directive in the
Snakefile, and the comment names the nf-core processes it serves. The images
are exactly the ones the nf-core modules pin, under the file names Nextflow
gives them. That means a folder of images downloaded by
`nf-core pipelines download` can be reused as-is.

**`container_urls:`** gives the URL each image comes from. It has the same keys
as `containers:`:

```yaml
container_urls:
  fastqc: https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0
  star: https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/.../data
```

Only `pipeline pull` reads this section, and it doesn't affect runs.

Because `containers:` holds plain file paths, Snakemake never downloads images
on its own. It runs `apptainer exec <file>` directly. That is why a run works
on offline compute nodes once the images have been pulled.

**Using a different image** for one step: either edit `containers.yaml`, or
override the key in your own config (it is merged over `containers.yaml`):

```yaml
# my.yaml
containers:
  fastqc: /shared/images/fastqc-0.12.1.sif
```

Pass the same `-c my.yaml` to `pipeline pull` so it checks the overridden paths.

| Pipeline | Images |
|---|---|
| rnaseq | 30 |
| chipseq | 34 |
| atacseq | 33 |
| hic | 14 |
| wgbs | 14 |

---

## Configuring a pipeline

### How the configuration is layered

Each port ships a default **`config/config.yaml`** that mirrors nf-core's
`nextflow.config` and `conf/modules.config` for the pinned release. Every tool
argument, threshold, reference path, option and resource limit lives there, and
nothing is hardcoded in the Snakefile.

When you run `pipeline run <name> -c my.yaml`, the configs are merged in this
order, and a later file overrides an earlier one:

```
snk-<folder>/config/config.yaml       defaults (nf-core params)
snk-<folder>/config/containers.yaml   image paths
my.yaml                               your settings
```

Nested sections such as `resources:` are merged key by key. So you can either
edit a full copy made with `pipeline config`, or write a short file holding
only what differs from the defaults:

```yaml
# rnaseq.config.yaml (minimal)
input: samplesheet.csv
outdir: results
fasta: /scratch/me/genomes/GRCh38.fa.gz
gtf: /scratch/me/genomes/gencode.v44.gtf.gz
gencode: true
resources:
  process_high: {threads: 12, mem_mb: 64000}
```

### Reading `config.yaml`

Every option is annotated with a `###` marker that gives the accepted values,
sometimes followed by a short explanation:

```yaml
aligner: "bwa"         ### ["bwa", "bowtie2", "chromap", "star"]
macs_gsize: ""         ### <integer, or "" to compute from read_length>  # Effective genome size...
fasta: ""              ### <path or URL: .fa/.fasta(.gz)>  # Genome FASTA (.gz allowed). Required.
```

- `["a", "b"]`: choose one of the listed values.
- `<...>`: a free value of the described type.
- `""` usually means "not set". For references, an empty value means the file
  or index is built from `fasta`. Pass an existing index to skip that step.

### What every config needs

| Key | Meaning |
|---|---|
| `input` | Samplesheet CSV (format per pipeline below). |
| `outdir` | Published results, with the same layout as nf-core's `--outdir`. |
| `workdir` | Intermediate files that nf-core keeps in its `work/` (all pipelines except hic). |
| `fasta` | Genome FASTA (`.gz` allowed). Required. |
| `resources` | `{threads, mem_mb}` per process label (see below). |

Paths may be local or http(s)/ftp URLs. Remote files are downloaded by a
`fetch_remote` rule, which needs internet. Use local paths on offline compute
nodes.

### Resources

`resources:` mirrors nf-core's process labels (`process_low`,
`process_medium`, `process_high`, …; for hic: `low`, `medium`, `high`, …). Each
label sets the threads and memory for the jobs that carry it:

```yaml
resources:
  process_medium: {threads: 6,  mem_mb: 36864}
  process_high:   {threads: 12, mem_mb: 73728}
```

Snakemake never gives one job more threads than `-t`. Memory is only enforced
by a cluster scheduler, or locally if you add `-- --resources mem_mb=<N>`.
Lower these values on small machines.

---

## Running each pipeline

All pipelines follow the same steps. Only the samplesheet and the key options
differ.

```bash
mkdir -p ~/projects/exp1 && cd ~/projects/exp1
pipeline config <name> -p .          # then edit <name>.config.yaml
pipeline run <name> -c <name>.config.yaml -n
pipeline run <name> -c <name>.config.yaml -t 8
```

Each port also has a test profile, `snk-<folder>/config/test.yaml`, that
mirrors nf-core's `conf/test.config` (the same small public dataset and
params). It is a good first run to check an installation.

### rnaseq

Samplesheet:

```csv
sample,fastq_1,fastq_2,strandedness
WT_REP1,WT_REP1_R1.fastq.gz,WT_REP1_R2.fastq.gz,auto
```

- `strandedness`: `forward`, `reverse`, `unstranded`, or `auto` (inferred
  with Salmon). Several rows with the same `sample` are runs of one sample and
  are concatenated.
- Key options: `fasta` plus `gtf` (or `gff`), `gencode`, `star_index` /
  `salmon_index` (built if empty), `pseudo_aligner`, `skip_bbsplit` /
  `bbsplit_fasta_list`, `additional_fasta`.
- Ported route: `aligner: star_salmon` with Trim Galore. Options that are not
  ported yet are rejected at start-up with a clear message.

### chipseq

Samplesheet:

```csv
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
SPT5_T0,SPT5_T0_R1.fastq.gz,SPT5_T0_R2.fastq.gz,1,SPT5,SPT5_INPUT,1
SPT5_INPUT,INPUT_R1.fastq.gz,INPUT_R2.fastq.gz,1,,,
```

- Leave `fastq_2` empty for single-end data. Control (input) rows leave
  `antibody` and `control` empty.
- Key options: `aligner` (`bwa`, `bowtie2`, `chromap`, `star`), `narrow_peak`,
  `macs_gsize` *or* `read_length`, `blacklist`, `gtf`/`gff`, prebuilt
  `<aligner>_index`.

### atacseq

Samplesheet:

```csv
sample,fastq_1,fastq_2,replicate
OSMOTIC_STRESS_T0,T0_R1.fastq.gz,T0_R2.fastq.gz,1
```

- Add `control,control_replicate` columns and set `with_control: true` to use
  controls.
- Key options: `aligner`, `narrow_peak`, `macs_gsize` *or* `read_length`,
  `blacklist`, `gtf`/`gff`, prebuilt `<aligner>_index`.

### hic

Samplesheet:

```csv
sample,fastq_1,fastq_2
SRR4292758,SRR4292758_R1.fastq.gz,SRR4292758_R2.fastq.gz
```

- Several rows per sample are lanes or runs, each processed as one chunk and
  then merged.
- Key options: `digestion` (`hindiii`, `mboi`, `dpnii`, `arima`), or set
  `restriction_site` and `ligation_site` directly. Also `bin_size` (list of
  resolutions in bp), `bwt2_index` and `chromosome_size` (built if empty).
- hic has no `workdir`. Everything goes under `outdir`.

### wgbs

Samplesheet:

```csv
sample,fastq_1,fastq_2,genome
SRR389222_sub1,SRR389222_sub1.fastq.gz,,
```

- Several rows with the same `sample` are concatenated. Leave `fastq_2` empty
  for single-end data.
- Key options: `aligner` (`bismark`, `bismark_hisat`, `bwameth`, `bwamem`),
  `bismark_index` / `bwameth_index` (built if empty), library presets such as
  `rrbs`, `em_seq`, `pbat`, `nomeseq`, and TAPS with `bwamem`.

---

## Running on an offline cluster (Rorqual)

Compute nodes have no internet, so do everything that downloads on a login
node:

```bash
# login node (internet)
module load apptainer
source /scratch/<user>/claude_jobs/venv_snk/bin/activate
pipeline list
pipeline pull rnaseq                     # download images once
# copy input data and references to /scratch beforehand; use local paths in the config

# compute node (salloc / sbatch)
module load apptainer
source /scratch/<user>/claude_jobs/venv_snk/bin/activate
cd /scratch/<user>/exp1
pipeline run rnaseq -c rnaseq.config.yaml -t $SLURM_CPUS_PER_TASK -B /scratch/<user>
```

- Bind `/scratch/<user>` (or whichever tree holds data and references). Only
  the run directory is bound automatically.
- Config paths must be local, because URLs trigger `fetch_remote`, which fails
  offline.
- If a run stops early with "image(s) missing", go back to a login node and run
  `pipeline pull`.

---

## Repository layout

```
pipelines/
├── pipeline.py               # the `pipeline` command
├── README.md
├── nfcore-<pipeline>/        # nf-core source, reference only (never edited)
└── snk-<pipeline>/           # Snakemake port
    ├── Snakefile             # all rules
    ├── config/
    │   ├── config.yaml       # defaults for every parameter
    │   ├── containers.yaml   # image per rule + download URLs
    │   ├── test.yaml         # nf-core test profile
    │   └── samplesheet.csv   # test samplesheet
    ├── containers/           # Singularity images (filled by `pipeline pull`)
    ├── scripts/              # nf-core bin/ scripts and custom helpers
    └── assets/               # nf-core assets (MultiQC configs, blacklists, …); not in hic
```
