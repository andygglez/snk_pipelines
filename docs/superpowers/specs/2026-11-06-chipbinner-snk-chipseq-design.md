# chipbinner downstream of snk-chipseq: design

Date: 2026-11-06. Status: approved in brainstorming, not yet implemented.

## Goal

Run the Python port of padilr1/ChIPbinner (`chipbinner_reimplementation/`) as a port extension of
snk-chipseq (nf-core/chipseq 2.1.0). It has no nf-core counterpart. It consumes snk-chipseq's
final BAMs and is configured by a new `chipbinner:` section of `snk-chipseq/config/config.yaml`.

## Decisions

| topic | decision |
|---|---|
| input | final BAMs `{ML}/<id>.final.bam` (+ `.bai`); chipbinner bins and normalises them itself |
| annotations | always snk-chipseq's: chrom.sizes (`custom_getchromsizes`), gene BED, intergenic BED, blacklist; no shipped hg38/mm10 files |
| container | custom `snk-chipseq/containers/chipbinner.def` → `chipbinner.sif`, built on Rorqual |
| rules | one rule per step (12 rules) |
| code | moved (not copied) to `snk-chipseq/scripts/`; `chipbinner_reimplementation/` keeps validation, tools and examples |
| `filter_low_counts` | removed from the software entirely |
| versions | not in `versions.yaml` (not an nf-core process) |

## 1. Config section

Placed after the heatmap / `scale_factor_barplot` block of `config.yaml`. The `comparisons` entry is
live, not commented out, so it doubles as the example; it matches the default port samplesheet
(T0/T15, SPT5), so a default run includes chipbinner.

```yaml
# =============================================================================
# chipbinner: port extension, no nf-core counterpart. Python port of
# padilr1/ChIPbinner (two-condition 10 kb bin scatterplots + HDBSCAN clusters),
# run on the final BAMs (<id>.final.bam). Port-format samplesheets only.
# Off while `comparisons` is empty. Outputs: <aligner>/merged_library/chipbinner/
# =============================================================================
chipbinner:
  # One entry per comparison. Remove all entries (comparisons: []) to turn chipbinner off.
  comparisons:
    - baseline: T0          # condition on the x axis (samplesheet `condition`)
      treatment: T15        # condition on the y axis
      proteins: all         # all | [SPT5, ...]: ChIP proteins present in both conditions
      range: auto           # auto (extent of hexes with > 10 bins) | [min, max] (both axes)
                            # | {x: [min, max] | auto, y: [min, max] | auto}
      pow: 0.75             # exponent of the hexbin count shading
      overrides: {}         # per protein, e.g. {SPT5: {range: [-5, 1], pow: 1.25}}

  binning:
    bin_size: 10000          # bedtools makewindows -w over the chipseq chrom.sizes
    windows: ""              # or a windows BED (overrides bin_size)
  MS_normalize: false        # true: each replicate's samplesheet scale_factor is its
                             # MS coefficient -> MSbigwigs/ (empty scale_factor = 1)
  norm:                      # ChIPbinner::norm_bw
    use_input: false         # true: log2(ChIP/input), input = the control snk-chipseq assigns
    depth_norm: true
    raw_count_cutoff: 100
    pseudocount: 1.0e-15

  plot_defaults: {range: auto, pow: 0.75, square_auto_range: true}
  genic_plot:   {hexbins: 150, show_legend: true, legend_pos: right, show_scales: true}
  density_plot: {hexbins: 150, show_scales: true, include_additional_density_plot: true,
                 filter_extreme_bins: true, cluster_colors: bright}   # bright | alphabet
  clustering:                # HDBSCAN grid: one job per (minpts, minsamps)
    minpts: [100, 400, 700, 1000]
    minsamps: [100, 400, 700, 1000]
    number_of_clusters: [2, 3]
    cluster_selection_method: eom
  formats: [pdf, html]

  pca: {enabled: true}
  correlation: {enabled: true, method: pearson}
  differential: {enabled: false, bootstrap: 1000, K: 5000,          # ROTS; needs >= 2 reps per condition
                 clusters: {minpts: 100, minsamps: 100, number_of_clusters: 3}}
  enrichment:   {enabled: false, db_dir: "", regions: [genome_wide, genic, intergenic],
                 n_elements: 10, cutoff_for_overlap: 0,
                 clusters: {minpts: 100, minsamps: 100, number_of_clusters: 3}}
```

Compared with `chipbinner_reimplementation/config/chipbinner.yaml`:

- Removed because snk-chipseq supplies them: `outdir`, `genome_assembly`, `annotations`,
  `samplesheet`, `sample_folder`, `MS_coefficients` (replaced by the samplesheet `scale_factor`).
- Removed (YAGNI): `merge_replicates` (always the mean of all N replicates), `binning.keep_read_bed`,
  `clustering.cores`.
- Removed from the software: `filter_low_counts`.
- Renamed: `treatment`/`control` → `treatment`/`baseline` (as in snk-rnaseq `comparisons`);
  `marks` → `proteins`.
- Threads/memory go to the existing `resources:` section as `chipbinner_<rule>` entries.
- `containers.yaml`: `chipbinner: containers/chipbinner.sif`.

`test.yaml` (nf-core-format samplesheet) leaves chipbinner off, because chipbinner runs only with
port samplesheets.

## 2. Connection to snk-chipseq

**Activation.** chipbinner outputs are added to `all_outputs()` only when `NEW_SHEET` is true and
`chipbinner.comparisons` is non-empty.

**Samples.** The chipbinner sample table is built in memory from `SAMPLES`; no chipbinner
samplesheet is written by the user.

| chipbinner concept | snk-chipseq source |
|---|---|
| Condition, Mark | `condition`, `protein` of ChIP rows (input/IgG rows excluded) |
| replicates | every `<condition>_<protein>_<rep>` sample (any number) |
| input of a replicate | `SAMPLES[id]["control"]` (assigned by `read_port_samplesheet`) |
| MS coefficient | `SAMPLES[id]["scale_factor"]` (used only when `MS_normalize`) |
| BAM | `{ML}/<id>.final.bam` + `.bai` |

**Annotations:**
- chrom.sizes comes from `custom_getchromsizes`.
- The gene BED is config `gene_bed`, otherwise built by gtf2bed.
- The intergenic BED is config `intergenic_bed`, otherwise built by the existing `intergenic_bed`
  rule. That rule must now also be defined when chipbinner is active.
- The blacklist is `BLACKLIST` (gunzipped). When it is unset, chipbinner gets an empty file.

All of these are declared as rule inputs.

**Outputs.** `CB = {ML}/chipbinner` (`<outdir>/<aligner>/merged_library/chipbinner/`). Binned BEDs go
to `{W}/chipbinner/` (not published).

## 3. Rules

`BWD` = `MSbigwigs` when `MS_normalize`, otherwise `notMSbigwigs`. `<cmp>` = `<treatment>_<baseline>`.
`wildcard_constraints` limit `cmp`, `protein`, `cond` and `id` to known values with `alternation()`,
because condition names may contain `_`.

| rule | scope | inputs | outputs |
|---|---|---|---|
| `chipbinner_windows` | once | chrom.sizes or `binning.windows` | `{W}/chipbinner/windows.bed` |
| `chipbinner_bin` | ChIP replicate and its control | `<id>.final.bam`, `.bai`, windows | `{W}/chipbinner/bam2bed/<id>.binned.bed` |
| `chipbinner_norm` | ChIP replicate | binned BED (+ control's when `use_input`), chrom.sizes, blacklist | `CB/BWD/<id>.binned.bw` |
| `chipbinner_merge` | condition × protein | its replicate bigWigs | `CB/BWD_final/<cond>.<protein>.bw` (mean; 1 replicate = copy) |
| `chipbinner_explore` | protein with ≥ 2 replicate bigWigs | those bigWigs | `CB/exploratory/<protein>_{PCA,correlation}_plot.<fmt>` (each only if enabled) |
| `chipbinner_legend` | cmp | none | `CB/comparisons/<cmp>/legend.png` |
| `chipbinner_scatter` | cmp × protein | 2 final bigWigs, gene/intergenic BEDs | `<P>.scatterplot.<fmt>`, `<P>.plot_range.yaml`, `<P>.genic.vs.intergenic.png` |
| `chipbinner_preclust` | cmp × protein | 2 final bigWigs | `<P>_preclust_mat.csv`, `<P>_preclust_pooled.bed`, `<P>_preclust_smoothScatter.jpeg` |
| `chipbinner_hdbscan` | × (minpts i, minsamps j) | mat.csv | `<K>/clus.<cmp>_<protein>_clustered.<i>.<j>.{txt,png}`, `<K>/log.<cmp>_<protein>.<i>.<j>.txt` |
| `chipbinner_density` | × n | mat, pooled BED, labels, plot_range.yaml, final bigWigs, BEDs | `<K>/<cmp>_<protein>.<n>clusters.<i>.<j>.annotated_clusters.bed`, `<K>/<cmp>_<protein>_<n>clusters.<i>.<j>.plots_density_scatterplot.<fmt>` |
| `chipbinner_differential` | cmp × protein, if enabled | replicate bigWigs of both conditions, blacklist, selected annotated BED | `CB/comparisons/<cmp>/<protein>/differential/{stats_for_all_bins_no_clustering.csv,stats_per_cluster.csv,ROTS_output_summary.txt}` |
| `chipbinner_enrichment` | cmp × protein, if enabled | selected annotated BED, pooled BED, gene/intergenic BEDs, `db_dir` | `directory(CB/comparisons/<cmp>/<protein>/enrichment)` |

Paths and placeholders:
- `<P>` = `CB/comparisons/<cmp>/<cmp>_<protein>` and `<K>` = `CB/comparisons/<cmp>/clustering`.
- "Selected annotated BED" is the `chipbinner_density` output for `differential.clusters` or
  `enrichment.clusters`. It may be empty, in which case differential runs without clusters and
  writes an empty `stats_per_cluster.csv`, and enrichment writes an empty directory.
- Empty placeholder outputs are kept as in the current Workflow: too few shared bins, or fewer
  HDBSCAN clusters than `n`.

Departures from the snakemake_chipbinner tree, all needed for several proteins per comparison:
- `legend.png` has its own rule.
- HDBSCAN logs are named `log.<cmp>_<protein>.<i>.<j>.txt`.
- `differential/` and `enrichment/` live under `<cmp>/<protein>/`.

**Code layout:**
- `snk-chipseq/scripts/chipbinner/`: the library, moved from `chipbinner_reimplementation/chipbinner/`.
- `snk-chipseq/scripts/chipbinner/steps.py` (new): one function per rule step, extracted from
  `Workflow` (auto-range + scatter, HDBSCAN with log and PNG, annotate + density with placeholders,
  exploratory, differential, enrichment).
- `snk-chipseq/scripts/chipbinner_cli.py`: one subcommand per rule, each a thin wrapper over
  `steps.py` or the library. `merge-norm-bw` accepts N bigWigs; the `filter-low-counts` subcommand is
  removed.
- `chipbinner/pipeline.py` (`Workflow`) calls `steps.py`.
- `chipbinner_reimplementation/chipbinner_pipeline.py` imports the package from `snk-chipseq/scripts/`.
  Remove `filter_low_counts` from its `config/test_synthetic.yaml`.

Each rule runs `python3 {SCRIPTS}/chipbinner_cli.py <subcommand> ...` with `container: config["chipbinner"]`.

## 4. Container, errors, testing

**Container:**
- `snk-chipseq/containers/chipbinner.def` bootstraps from `docker://mambaorg/micromamba`. It installs
  one pinned conda-forge/bioconda env with python, numpy, pandas, scipy, scikit-learn, hdbscan,
  plotly, python-kaleido, pyBigWig, pysam and pyyaml, at the versions the reimplementation was
  validated with. The Chromium that kaleido needs goes into the image, because compute nodes are
  offline.
- Build on a Rorqual login node. If an unprivileged `apptainer build` is not allowed there, ask the
  user before switching to the fallback: a sandbox built from a conda-pack tarball.
- `containers.yaml` gets a `container_defs:` section (`chipbinner: containers/chipbinner.def`).
  `pipeline.py pull` builds those with `apptainer build` instead of pulling a URL.
- Rules export `MPLCONFIGDIR=/tmp/mplcfg` and point `HOME` to a writable temporary directory for
  kaleido/Chromium.

**Start-up checks** (fail before any job, with a message naming the entry):
- a comparison's baseline or treatment is not a samplesheet condition
- a listed protein is missing from either condition, or no protein is shared
- `differential.enabled` while a compared condition × protein has fewer than 2 replicates
- `enrichment.enabled` with an empty `db_dir`
- `norm.use_input` while a compared ChIP replicate has no control
- `MS_normalize` with a non-positive scale_factor (empty means 1)

**Tests** (all on Rorqual; nothing runs locally, not even a dry run):

1. Refactor parity. Run the standalone `chipbinner_pipeline.py` on the synthetic BAMs with the code
   before and after the move (both configs without `filter_low_counts`), then diff every output.
   They must be identical.
2. Port run. `snk-chipseq/config/scenarios/chipbinner_test.yaml` overlays `test_local.yaml` with
   `input: config/samplesheet_port_local.csv`. This is a new file: a copy of the port-format
   `config/samplesheet.csv` (T0/T15 SPT5, two replicates each, pooled inputs, with scale_factor
   values filled in) whose FASTQ URLs are replaced by absolute paths to pre-downloaded files on
   Rorqual. The existing `samplesheet_local.csv` is in nf-core format and stays unchanged. Rorqual's
   `claude_jobs/` was purged, so the test FASTQs must be downloaded locally and rsynced again, and
   the images pulled again. Its settings are:
   - `bin_size: 1000`
   - `minpts`/`minsamps: [20, 50]`, `number_of_clusters: [2, 3]`
   - `use_input: true`, `MS_normalize: true`
   - `differential.enabled: true` with `bootstrap: 100`
   - enrichment off (there is no yeast region DB)
3. Port vs standalone. Run `Workflow` on the same final BAMs with an equivalent standalone config.
   Compare:
   - bigWigs via `bigWigToBedGraph`
   - CSV/BED/label/plot_range files with `diff`
   - plots by existence only

   Verdict: identical, or differs (with an explanation).
4. No regression. `test.yaml` (nf-core sheet) runs with chipbinner off and keeps its file inventory.
