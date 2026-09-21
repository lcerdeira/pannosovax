#!/bin/bash
#SBATCH --job-name=pnv-plddt
#SBATCH --partition=amd
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/plddt_%j.out
#SBATCH --mail-type=END,FAIL
# Repara a escala de pLDDT dos modelos ESMFold e refaz o crossmatch.
# Os 57 modelos gravados em 0-1 foram descartados pelo corte --min-plddt 70 sem
# aviso, então a comparação anterior entre pares nunca os enxergou.
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
export PATH="$HOME/bin:$PATH"
cd ~/pannosovax
log(){ echo "[$(date '+%F %T')] $*"; }

log "════ reparo da escala de pLDDT ════"
python - <<'PY'
import sys, pathlib
sys.path.insert(0, "scripts")
from fold_missing_esmfold import normalize_plddt
n = 0
for org in ["kpsc", "abau", "spneu"]:
    for p in sorted(pathlib.Path(f"results/04_structures/{org}").glob("*.pdb")):
        t = p.read_text(); f = normalize_plddt(t)
        if f != t:
            p.write_text(f); n += 1
print(f"  {n} modelos reescalados para 0-100")
PY

log "════ 04b · TRES patogenos (reavaliacao honesta da opcao A) ════"
python scripts/04b_structural_crossmatch_v2.py --window 9 --min-plddt 70 --tm 0.5 \
  --min-partners 2 || { log "✗ 04b(3)"; exit 1; }
cp results/04_shared/shared_regions_v2.tsv results/04_shared/shared_regions_3way.tsv
python scripts/04b2_validate_shared.py 2>&1 | tail -6
cp results/04_shared/shared_validated_v2.tsv results/04_shared/shared_validated_3way.tsv

log "════ 04b · por PAR (opcao B escolhida) ════"
python scripts/04b_structural_crossmatch_v2.py --window 9 --min-plddt 70 --tm 0.5 \
  --min-partners 1 || { log "✗ 04b(1)"; exit 1; }
python scripts/04b2_validate_shared.py || { log "✗ 04b2"; exit 1; }
python scripts/04d_check_shared_survival.py || log "!! inconsistente"
python scripts/04c_select_shared.py || { log "✗ 04c"; exit 1; }
python scripts/08_construct_builder.py || { log "✗ 08"; exit 1; }
python scripts/09_physchem.py 2>&1 | grep -E "length|APROVADO" || true
python scripts/12_codon_optimize.py 2>&1 | grep -E "proteína|CAI|GC global" || true

log "════ 04e · assimetria, agora com todas as estruturas ════"
python scripts/04e_analyze_asymmetry.py || true
log PLDDT_DONE
