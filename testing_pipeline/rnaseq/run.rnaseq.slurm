#!/bin/bash
#SBATCH --time=0-06:00
#SBATCH --job-name=rna-seq-pipeline
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --account=def-majewski
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=andy.garciagonzalez@mail.mcgill.ca

/lustre09/project/6007495/andygg98/pipelines/.venv/bin/python /lustre09/project/6007495/andygg98/pipelines/pipeline.py \
                                                                run rnaseq -c configs/rnaseq.config.yaml -t 8