# chipbinner_snakemake

[andygglez/snakemake_chipbinner](https://github.com/andygglez/snakemake_chipbinner), rebuilt on the Python reimplementation of [padilr1/ChIPbinner](https://github.com/padilr1/ChIPbinner) (`scripts/chipbinner/`; development copy and validation in `../chipbinner_reimplementation/`). It keeps the old pipeline's rules, config sections and `outdir/` tree. The R scripts are replaced by `scripts/chipbinner_steps.py`, one subcommand per rule, and every step runs in a container.

## Running

Standalone pipeline: it is deliberately not registered with the repo's `pipeline` command (the folder has no `snk-` prefix). Once it is validated it will be integrated downstream of snk-chipseq.

```bash
cd chipbinner_snakemake
# edit config/config.yaml (organism, comparisons, ...) and config/sample_list.csv
snakemake --cores 8 --software-deployment-method apptainer --apptainer-args "--bind <data_dir>"
```

Build the images once, on a node with internet (the Rorqual login node):

```bash
apptainer build containers/chipbinner.sif containers/chipbinner.def
apptainer pull containers/depot.galaxyproject.org-singularity-bedtools-2.30.0--hc088bd4_0.img \
    https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0
```

`config/sample_list.csv` keeps the old columns, `Rep1,Rep2,Condition,Mark`, and accepts optional `Input1,Input2` (used with `general.use_input`). Each cell is a BAM, as before, a binned-count BED, or an already normalised bigWig. Relative cells are looked up in `inputs.sample_folder`. MS coefficients (`references.MS_yaml`) are keyed by the file name without `.bam`, as before.

## Rules (old → new)

| snakemake_chipbinner | here | output |
|---|---|---|
| `bam2bed` (bedtools) | `binning_windows`, `bam2bed` (bedtools) | `outdir/bam2bed/<sample>.{bed,binned.bed}` |
| `filter_low_counts` | removed (no low-count filtering option) | — |
| `bed2bw` | same | `outdir/{MS,notMS}bigwigs/<sample>.binned.bw` |
| `merged_bw_files` / `notmerged_bw_files` | same | `outdir/{MS,notMS}bigwigs_final/<condition>.<mark>.bw` |
| `compare` | `compare`, `legend` | `comparisons/<T>_<C>/<T>_<C>_<mark>.scatterplot.{pdf,html}`, `_preclust_{mat.csv,pooled.bed,smoothScatter.jpeg}`, `.genic.vs.intergenic.png`, `.plot_range.yaml`, `legend.png` |
| `clustering` (clus.py loop + annotate.clusters.R) | `clustering` (one job per minpts × minsamps), `annotate_clusters` | `clustering/clus.<T>_<C>_<mark>_clustered.<i>.<j>.{txt,png}`, `<T>_<C>_<mark>.<n>clusters.<i>.<j>.annotated_clusters.bed` |
| `density_plot` | same | `clustering/<T>_<C>_<mark>_<n>clusters.<i>.<j>.plots_density_scatterplot.{pdf,html}` |
| (PCA.R, never wired in) | `pca`, `correlation` | `outdir/exploratory/<mark>_{PCA,correlation}_plot.*` |
| — | `differential` (ROTS), `enrichment` (LOLA), off by default | `comparisons/<T>_<C>/{differential,enrichment}/<T>_<C>_<mark>.*` |

## Plot ranges

`comparison_parameters` keeps the old `"<min>_<max>_<pow>"` strings, keyed `<T>_<C>_<mark>`. As in snakemake_chipbinner, min/max apply to both axes. Each entry can also be `auto` or `{range: auto | [min, max] | {x: [..], y: auto}, pow: …}`. A comparison without an entry uses `plot_defaults`, or the `range`/`pow` set on the comparison entry, which is `auto` by default.

`auto` decides the limits the way `ChIPbinner::genic_intergenic_scatterplot` does:
1. **Trim:** keep the bins scored in both samples, trimmed to the 0.1–99.9 % quantiles of x and y.
2. **Hex grid:** run `hexbin(x, y, xbins = genic_plot.hexbins)`.
3. **Per axis:** each axis runs from the lowest to the highest hex centre holding more than 10 bins. ChIPbinner applies this rule to `min_x`/`min_y` only; it has no default for `max_x`/`max_y` (they are required arguments), so the port uses the same rule for the maximum.
4. **Margin:** each axis gets ChIPbinner's 5 % `expand_range` on both sides.

x and y are set independently, as in R. The values used (before the 5 % margin) go to `<T>_<C>_<mark>.plot_range.yaml`. The density plots read that file, just as snakemake_chipbinner passed the same min/max to `density_based_scatterplot`, whose limits have no automatic default.

Plots carry only their title; there is no "N bins …" subtitle.

## Deliberate changes from snakemake_chipbinner

- **Clusters:** annotated clusters are written as BED, not `.rda`. Every grid pair's density plot is a declared output; the old rule declared only `100.100`.
- **Clustering logs:** they are `log.<T>_<C>_<mark>.<i>.<j>.txt`, because several marks share the clustering folder.
- **Samplesheet:** one condition may have two replicates while another has one. Condition names may contain `_`, because comparisons are split on known conditions.
- **Reruns:** nothing needs deleting before a rerun, because Snakemake tracks every output.
- **Formats:** plots are written in every `formats` entry (default `pdf` + interactive `html`).
- **Clustering defaults:** `minpts`/`minsamps` default to the old `seq 100 300 1000` grid, and `number_of_clusters` to `[2, 3]`. When HDBSCAN finds more clusters than requested, ChIPbinner keeps the *n* with the lowest mean signal; on whole genomes with small `minpts` these can be tiny, and the rule logs a warning.

## Test

`config/test.yaml` runs every rule, including MS normalisation, merging, both range styles, differential and enrichment. It uses synthetic BAMs (hg38 chr21 + chr22, made by `tools/make_synthetic_bams.py`) and ChIPbinner's downsampled region DB, both placed under `test_data/` (not in git):

```bash
snakemake --cores 2 --configfile config/test.yaml --software-deployment-method apptainer   # -> results_test/outdir
```
