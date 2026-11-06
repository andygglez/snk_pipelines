#!/usr/bin/env python3
"""Run the ChIPbinner workflow from a config (the future `chipbinner:` section of snk-chipseq).

    chipbinner_pipeline.py config/chipbinner.yaml [--force] [--outdir DIR]

Output tree = andygglez/snakemake_chipbinner's (see chipbinner/pipeline.py).
"""
import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from chipbinner.pipeline import Workflow, load_config  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("config")
ap.add_argument("--force", action="store_true", help="recompute outputs that already exist")
ap.add_argument("--outdir", help="override chipbinner.outdir")
a = ap.parse_args()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
for noisy in ("kaleido", "choreographer", "choreo", "asyncio"):
    logging.getLogger(noisy).setLevel(logging.WARNING)
cfg = load_config(a.config)
if a.outdir:
    cfg["outdir"] = a.outdir
Workflow(cfg, force=a.force).run()
