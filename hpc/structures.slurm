#!/bin/bash
#SBATCH --job-name=pnv-struct
#SBATCH --partition=amd
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=23:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/struct_%j.out
#SBATCH --mail-type=END,FAIL
# Estruturas + bloco compartilhado, recalculados sobre as candidatas CORRIGIDAS.
#
# O bloco compartilhado atual foi derivado do surfaceome piloto: 4 dos 5 epitopos
# que ele põe no construto vêm de proteínas que o filtro rigoroso rejeita. Não dá
# para remendar — o bloco precisa ser redescoberto no conjunto vigente (119/101/34).
#
# Busca no AlphaFold DB é retomável (pula PDB existente); o que não tiver modelo
# fica de fora desta rodada e é contabilizado no resumo.
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
export PATH="$HOME/bin:$PATH"
cd ~/pannosovax
log(){ echo "[$(date '+%F %T')] $*"; }

command -v TMalign >/dev/null || { log "✗ TMalign ausente no PATH"; exit 1; }

log "════ estruturas · AlphaFold DB (todas as candidatas) ════"
python scripts/fetch_structures_afdb.py --per-org 200 --max-len 800 || log "!! busca parcial"

log "════ cobertura estrutural ════"
for o in kpsc abau spneu; do
  n=$(ls results/04_structures/$o/*.pdb 2>/dev/null | wc -l)
  c=$(($(wc -l < results/03_surfaceome/${o}_candidates.tsv) - 1))
  printf "  %-6s %s estruturas de %s candidatas\n" "$o" "$n" "$c"
done

log "════ 04b · crossmatch estrutural (TMalign) ════"
python scripts/04b_structural_crossmatch_v2.py --window 9 --min-plddt 70 --tm 0.5 \
  || { log "✗ 04b falhou"; exit 1; }

log "════ 04b2 · validar regiões contra epitopos seguros ════"
python scripts/04b2_validate_shared.py || { log "✗ 04b2 falhou"; exit 1; }

log "════ 04d · o bloco novo se sustenta? ════"
python scripts/04d_check_shared_survival.py || log "!! ainda inconsistente"

log "════ 04c/08 · remontar o construto com o bloco novo ════"
python scripts/04c_select_shared.py || { log "✗ 04c falhou"; exit 1; }
python scripts/08_construct_builder.py || { log "✗ 08 falhou"; exit 1; }
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

log "STRUCT_DONE"
