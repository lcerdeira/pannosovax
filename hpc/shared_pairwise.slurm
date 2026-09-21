#!/bin/bash
#SBATCH --job-name=pnv-pair
#SBATCH --partition=amd
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/pair_%j.out
#SBATCH --mail-type=END,FAIL
# Bloco compartilhado por PAR de patógenos (min-partners=1).
# Exigir os três rendia 2 regiões / 3 epitopos, todos colidindo com os blocos B/MHC,
# e não era falta de estrutura: dobrar a cobertura de A. baumannii não mudou nada.
# Ver o cabeçalho de scripts/04c_select_shared.py.
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
export PATH="$HOME/bin:$PATH"
cd ~/pannosovax
log(){ echo "[$(date '+%F %T')] $*"; }

log "════ 04b · crossmatch, regiões por par ════"
python scripts/04b_structural_crossmatch_v2.py --window 9 --min-plddt 70 --tm 0.5 \
  --min-partners 1 || { log "✗ 04b"; exit 1; }
log "════ 04b2 · validação ════"
python scripts/04b2_validate_shared.py || { log "✗ 04b2"; exit 1; }
log "════ 04d · consistência ════"
python scripts/04d_check_shared_survival.py || log "!! inconsistente"
log "════ 04c · bloco compartilhado ════"
python scripts/04c_select_shared.py || { log "✗ 04c"; exit 1; }
log "════ 08/09/12 · construto ════"
python scripts/08_construct_builder.py || { log "✗ 08"; exit 1; }
python scripts/09_physchem.py || true
python scripts/12_codon_optimize.py || true
python - <<'PY'
import json, os
b="results/08_construct/construct_blocks.json"
if os.path.exists(b):
    d=json.load(open(b)); bl=d.get("blocks",{})
    print(f"  construto: {d.get('length')} aa, {d.get('n_epitopes')} epitopos")
    print("  blocos:", {k:len(v) for k,v in bl.items()})
PY
log PAIR_DONE
