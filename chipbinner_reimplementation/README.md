# ChIPbinner scatterplots, reimplemented in Python + Plotly

This folder ports the two-condition scatterplots of [padilr1/ChIPbinner](https://github.com/padilr1/ChIPbinner) (R), plus the steps that feed them, to Python. Each subcommand of `chipbinner_cli.py` matches an R function:

| subcommand | R function | output |
|---|---|---|
| `norm-bw` | `norm_bw` | normalised bigWig/bedGraph (log2 ChIP/input, CPM-scaled) |
| `pre-clust` | `pre_clust` | `<prefix>_mat.csv` (treated, wildtype) + `<prefix>_pooled.bed` |
| `generate-clust` | `generate_clust` (HDBSCAN) | one label per bin |
| `annotate-clust` | `annotate_clust` | BED `chr start end A/B/C…` (ordered by mean signal) |
| `genic-scatter` | `genic_intergenic_scatterplot` | hexbins coloured by % genic and shaded by bin count |
| `density-scatter` | `density_based_scatterplot` | [magma density \|] genic/intergenic \| density-based clusters |
| `bin-counts` | (`bedtools bamtobed \| intersect -c`, snakemake_chipbinner `bam2bed`) | BAM → binned-count BED |
| `filter-low-counts` | `filter_low_counts` | binned BED without low-count bins |
| `apply-norm-factors` | `apply_normFactors` | DESeq2 / edgeR TMM / depth-only scaling factors |
| `merge-norm-bw` | `merge_norm_bw` | mean of replicate bigWigs |
| `plot-pca` / `plot-correlation` | `plot_PCA` / `plot_correlation` | PCA and correlation heatmap |
| `differential` | `differentialBinAnalysis` (ROTS) | per-bin and per-cluster statistics |
| `enrich-clust` | `enrich_clust` (LOLA + qvalue) | Fisher enrichment/depletion table + lollipop plot |

Every plot puts the wild type on x and the treated sample on y. Plots are written as interactive HTML (zoom, plus hover showing bin counts, % genic and genomic coordinates) and, through kaleido, as PNG/PDF/SVG.

`./run_example.sh` rebuilds `results_example/` from `example_data/` (the full Cal27 H3K36me2 NSD1-KO vs WT bigWigs and annotated clusters that ship with ChIPbinner).

## Workflow driver (future `chipbinner:` section of snk-chipseq)

`chipbinner_pipeline.py <config.yaml> [--force] [--outdir DIR]` runs the whole chain from a config and writes the output tree of [andygglez/snakemake_chipbinner](https://github.com/andygglez/snakemake_chipbinner). `config/chipbinner.yaml` documents every key; its defaults mirror that pipeline. Finished outputs are skipped unless `--force` is given.

- **Samplesheet:** `Rep1,Rep2,Condition,Mark`, as in snakemake_chipbinner, plus optional `Input1,Input2`. A cell can be a BAM, a binned-count BED or an already normalised bigWig.
- **Comparisons:** one entry per `{treatment, control}`, with x = control and y = treatment. The old `"<T>_<C>"` strings are accepted too.
  - Each entry takes `range: auto | [min, max] | {x: …, y: …}` and `pow`. Per-mark `overrides` replace either one for that mark.
  - `auto` uses the extent of the hexes holding more than 10 bins, with the same limits on both axes by default. The range actually used is written to `<T>_<C>_<mark>.plot_range.yaml`.
- **Axis text:** axis titles and tick labels are enlarged on every plot.

```
<outdir>/bam2bed/<sample>.binned[_filt].bed           BAM inputs only
<outdir>/{MS,notMS}bigwigs/<sample>.binned[_filt].bw   norm_bw per replicate
<outdir>/{MS,notMS}bigwigs_final/<cond>.<mark>.bw     merged replicates
<outdir>/comparisons/<T>_<C>/
    <T>_<C>_<mark>.scatterplot.{pdf,html}, legend.png, <T>_<C>_<mark>.plot_range.yaml
    <T>_<C>_<mark>_preclust_{mat.csv,pooled.bed,smoothScatter.jpeg}, <T>_<C>_<mark>.genic.vs.intergenic.png
    clustering/clus.<T>_<C>_<mark>_clustered.<minpts>.<minsamps>.{txt,png}, log.<minpts>.<minsamps>.txt
    clustering/<T>_<C>_<mark>.<n>clusters.<minpts>.<minsamps>.annotated_clusters.bed   (R writes .rda)
    clustering/<T>_<C>_<mark>_<n>clusters.<minpts>.<minsamps>.plots_density_scatterplot.{pdf,html}
    differential/   ROTS (2 replicates per condition)   } ChIPbinner functions that
    enrichment/     LOLA region DB                      } snakemake_chipbinner does not run;
<outdir>/exploratory/<mark>_{PCA,correlation}_plot.*    } off by default except PCA/correlation
```

When HDBSCAN finds fewer clusters than `number_of_clusters`, the annotation and density plot are written as empty placeholders, as in snakemake_chipbinner. When it finds more, ChIPbinner keeps the *n* clusters with the lowest mean signal; the port does the same and logs a warning. On the full Cal27 genome, minpts = minsamps = 100 gives 6 clusters, so "3 clusters" keeps just 526 bins and enrichment finds nothing. ChIPbinner's shipped Cal27 clusters cover about 190k bins, so use the larger values of the grid on whole genomes.

Test configs:
- `config/test_synthetic.yaml` exercises every step on synthetic BAMs (chr21 + chr22, made by `tools/make_synthetic_bams.py`) in about 2.5 minutes.
- `config/test_cal27.yaml` starts from ChIPbinner's real Cal27 bigWigs. It covers about 280k bins and is too heavy for the laptop, so run it on the cluster.

## Layout
- `chipbinner/` contains the library:
  - `genomics.py`: I/O and interval overlaps
  - `binning.py`: BAM → bins
  - `normalize.py`: norm_bw, filtering, normalisation factors, merging
  - `cluster.py`
  - `stats.py`: exact ports of `hexbin::hexbin`, `MASS::kde2d` and `base::pretty`
  - `plots.py` and `extra_plots.py`
  - `differential.py`: ROTS
  - `enrichment.py`: LOLA and qvalue
  - `explore.py`: PCA and correlation
  - `pipeline.py`: the workflow
- `tools/make_synthetic_bams.py` writes the synthetic test BAMs. They are labelled synthetic; they are not real data.
- `annotations/` holds ChIPbinner's internal hg38/mm10 genic, intergenic, blacklist and chrom.sizes files, converted from `R/sysdata.rda`.
- `tools/rda_to_bed.py` is the one-off `.rda` → BED/bigWig converter. It uses the `rdata` parser and never runs R code.
- `validation/` holds the checks against ChIPbinner's own outputs.

## Validation against ChIPbinner (downsampled example)
- **norm-bw:** the formula matches R exactly; fitting R's two library-size constants leaves a 1e-14 residual. The shipped R object was built with different library sizes than its shipped input BED gives, so the scores differ by a constant offset of 1.191.
- **pre-clust:** identical to R's shipped matrix (max diff 5e-15) once two things are reproduced: R's file was made with the 1 %/99 % filter that the current R code comments out, and with the columns swapped.
- **generate-clust:** uses scikit-learn's HDBSCAN, or the `hdbscan` package if installed. It agrees with R's labels at ARI 0.988 (305 vs 306 noise bins).

Further checks on the synthetic data, against the real tools and R:
- **Bin counting:** identical to `bedtools makewindows` + `intersect -c`, for both 10 kb and 1.05 kb windows (reads spanning two windows).
- **TMM factors:** match edgeR `calcNormFactors` to 2e-15.
- **`kde2d`:** matches `MASS::kde2d` to a relative 3e-15.
- **`pretty`:** breaks identical to R's.
- **Fisher test:** p-values and conditional odds ratios match R's `fisher.test` to a relative 5e-5.
- **PCA:** variance proportions match R's `prcomp` to 5e-6.
- **Not cross-checked:** ROTS, LOLA and qvalue, because they aren't installed in the local R. They are line-by-line ports of the Bioconductor code.

## Deliberate differences
- `max_x`/`max_y` are optional. They default, like ChIPbinner's `min_*`, to the extent of the hexes holding more than 10 bins.
- In `density-scatter`, the genic panel's opacity range defaults to (0.01, 3), clipped to 1; R uses (0.01, 1), which leaves most hexes nearly black. Use `alpha_range=(0.01, 1)` to get R's look.
- Clusters use a bright palette that stays legible on black (`--cluster-colors alphabet` gives ChIPbinner's palette). The scatter is subsampled to 120k points for display only; centroids and density use every bin. Unclustered bins are shown in grey, and a cluster key with bin counts is added.
- `norm-bw` drops chrM consistently. R drops it from the bin template but not from the score vector.
