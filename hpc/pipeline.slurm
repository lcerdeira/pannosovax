#!/bin/bash
#SBATCH --job-name=pnv-rerun
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/pipeline_%j.out
#SBATCH --error=/home/lshlt19/pannosovax/logs/pipeline_%j.out
#SBATCH --mail-type=END,FAIL

# Re-run completo do PanNosoVax em escala cheia.
bash "$HOME/pannosovax/hpc/run_full_pipeline.sh"
