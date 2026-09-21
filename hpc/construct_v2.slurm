#!/bin/bash
#SBATCH --job-name=pnv-v2
#SBATCH --partition=amd
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/v2_%j.out
#SBATCH --mail-type=END,FAIL
# Remonta o construto como v2: o bloco COMPARTILHADO (tese central) precisa do 04c,
# que por sua vez depende de 06/07 já terem rodado. O job anterior monta um construto
# sem esse bloco; aqui refazemos a cauda na ordem correta. Etapas rápidas (minutos).
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
cd ~/pannosovax
log(){ echo "[$(date '+%F %T')] $*"; }

log "→ 04c · seleção do bloco compartilhado"
python scripts/04c_select_shared.py || { log "✗ 04c falhou"; exit 1; }
log "→ 08 · construto (v2, com bloco shared)"
python scripts/08_construct_builder.py || exit 1
log "→ 09 · físico-química"; python scripts/09_physchem.py || true
log "→ 12 · expressão"; python scripts/12_codon_optimize.py || true

python - <<'PY'
import json, os
b="results/08_construct/construct_blocks.json"
if os.path.exists(b):
    d=json.load(open(b)); bl=d.get("blocks",{})
    print(f"  construto: {d.get('length')} aa, {d.get('n_epitopes')} epitopos")
    print("  blocos:", {k:len(v) for k,v in bl.items()})
    print("  SHARED populado:" , len(bl.get("shared",[])) > 0)
PY
echo V2_DONE
