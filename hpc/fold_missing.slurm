#!/bin/bash
#SBATCH --job-name=pnv-fold2
#SBATCH --partition=amd
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=23:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/fold2_%j.out
#SBATCH --mail-type=END,FAIL
# Preenche estruturas faltantes (ESMFold) e refaz o bloco compartilhado sobre a
# cobertura ampliada. Ver scripts/fold_missing_esmfold.py para o porquê.
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
export PATH="$HOME/bin:$PATH"
cd ~/pannosovax
log(){ echo "[$(date '+%F %T')] $*"; }

log "════ ESMFold · estruturas faltantes ════"
python scripts/fold_missing_esmfold.py || log "!! dobramento parcial"

for o in kpsc abau spneu; do
  n=$(ls results/04_structures/$o/*.pdb 2>/dev/null | wc -l)
  c=$(($(wc -l < results/03_surfaceome/${o}_candidates.tsv) - 1))
  printf "  %-6s %s/%s com estrutura\n" "$o" "$n" "$c"
done

log "════ 04b · crossmatch estrutural ════"
python scripts/04b_structural_crossmatch_v2.py --window 9 --min-plddt 70 --tm 0.5 || { log "✗ 04b"; exit 1; }
log "════ 04b2 · validação ════"
python scripts/04b2_validate_shared.py || log "!! nenhuma região validada"
log "════ 04d · consistência ════"
python scripts/04d_check_shared_survival.py || log "!! inconsistente"
log "════ 04c/08/09/12 · construto ════"
python scripts/04c_select_shared.py || log "!! 04c"
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
log FOLD2_DONE
