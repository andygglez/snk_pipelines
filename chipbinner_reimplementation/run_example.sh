#!/usr/bin/env bash
# Reproduces results_example/ from the ChIPbinner example data (Cal27 H3K36me2, NSD1-KO vs WT).
# Needs python3 with numpy, pandas, plotly, pyBigWig, scikit-learn (kaleido for png/pdf).
set -euo pipefail
cd "$(dirname "$0")"
E=example_data; O=results_example; mkdir -p $O
PY=${PY:-python3}
$PY chipbinner_cli.py genic-scatter \
  --treated $E/Cal27.NSD1_KO.H3K36me2.10kb.norm.bw --wildtype $E/Cal27.WT.H3K36me2.10kb.norm.bw \
  --title "Cal27 H3K36me2 · NSD1-KO vs WT" --xaxis-label WT --yaxis-label NSD1_KO \
  --min-x -5 --min-y -5 --max-x 1 --max-y 1 --formats html,png \
  --out-prefix $O/Cal27.WT_NSD1KO.H3K36me2.10kb.genic_intergenic
$PY chipbinner_cli.py density-scatter \
  --treated $E/Cal27.NSD1_KO.H3K36me2.10kb.norm.bw --wildtype $E/Cal27.WT.H3K36me2.10kb.norm.bw \
  --annotated-clusters $E/Cal27.WT.NSD1_KO.H3K36me2.10kb.annotated_clusters.bed.gz \
  --title H3K36me2 --xaxis-label WT --yaxis-label NSD1_KO \
  --min-x -5 --min-y -5 --max-x 2 --max-y 2 --hexbins 50 --pow 1.1 \
  --include-additional-density-plot true --formats html,png \
  --out-prefix $O/Cal27.WT.NSD1_KO.H3K36me2.10kb_density_scatterplot
