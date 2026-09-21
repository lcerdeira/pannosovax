#!/bin/bash
#SBATCH --job-name=pnv-surf
#SBATCH --partition=amd
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/surf_%j.out
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
cd ~/pannosovax
for o in kpsc abau spneu; do
  for t in signalp tmhmm; do
    echo "[$(date +%T)] -> $o/$t"
    python scripts/run_surfaceome_tools.py --organism $o --tool $t || echo "  FALHOU $o/$t"
  done
done
for o in kpsc abau spneu; do python scripts/03_surfaceome_filter.py --organism $o || true; done
echo SURFACEOME_DONE
