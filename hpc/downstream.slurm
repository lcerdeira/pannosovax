#!/bin/bash
#SBATCH --job-name=pnv-down
#SBATCH --partition=amd
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=23:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/down_%j.out
#SBATCH --mail-type=END,FAIL
# A jusante do surfaceome validado: epitopos -> conservação -> segurança ->
# cobertura -> construto. Cada etapa é retomável (cache por proteína no IEDB,
# skip-se-existe nas demais), então resubmeter continua de onde parou.
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
cd ~/pannosovax
T=${SLURM_CPUS_PER_TASK:-8}
log(){ echo "[$(date '+%F %T')] $*"; }
# Contador de falhas. O job precisa sair != 0 quando alguma etapa quebra, senão o
# Slurm marca COMPLETED e o job encadeado com --dependency=afterok dispara sobre
# resultados que não existem. Foi o que aconteceu no run de 18/09: o blastp falhou
# nos três organismos, tudo a jusante saiu vazio, e mesmo assim o job terminou 0.
FAILED=0
run(){ local n="$1"; shift; local t0=$SECONDS
  log "→ $n"; if "$@"; then log "✓ $n ($((SECONDS-t0))s)"; else
    local rc=$?; log "✗ $n FALHOU (rc=$rc)"; FAILED=$((FAILED+1)); return 1; fi; }

log "════ 05 · epitopos (IEDB) — etapa longa, retomável ════"
run epitopes python scripts/run_iedb_epitopes.py

log "════ 04d · consistência do bloco compartilhado (informativo) ════"
# Não conta como falha do pipeline: é uma pendência de manuscrito, não de execução.
python scripts/04d_check_shared_survival.py || log "!! bloco compartilhado inconsistente — ver acima antes de submeter"

log "════ 05c · conservação ════"
for o in kpsc abau spneu; do run "homologs:$o" python scripts/build_homologs_blast.py --organisms $o --threads $T; done
run conservation python scripts/apply_conservation.py

log "════ 06 · segurança ════"
for o in kpsc abau spneu; do for k in mhc1 mhc2 bcell; do
  run "safety:$o:$k" python scripts/06_safety_screen.py --organism $o --class $k; done; done

log "════ 07-08 · cobertura e construto ════"
for k in mhc1 mhc2; do run "coverage:$k" python scripts/07_population_coverage.py --class $k --n 18 --brazil-weight 0.5; done
run construct  python scripts/08_construct_builder.py
run physchem   python scripts/09_physchem.py
run expression python scripts/12_codon_optimize.py

log "════ RESUMO ════"
python - <<'PY'
import pandas as pd, os
def u(f):
    try: d=pd.read_csv(f,sep="\t"); return d["peptide"].nunique() if len(d) else 0
    except: return 0
for o in ["kpsc","abau","spneu"]:
    cand=sum(1 for _ in open(f"results/03_surfaceome/{o}_candidates.tsv"))-1
    safe=sum(u(f"results/06_safety/{o}_{k}_safe.tsv") for k in ["mhc1","mhc2","bcell"])
    print(f"  {o}: candidatas={cand} epitopos_seguros={safe}")
if os.path.exists("results/08_construct/construct.fasta"):
    s="".join(l.strip() for l in open("results/08_construct/construct.fasta") if not l.startswith(">"))
    print(f"  construto: {len(s)} aa")
PY
if [ "$FAILED" -gt 0 ]; then
  log "DOWNSTREAM_FALHOU: $FAILED etapa(s) quebraram — job sai != 0 para barrar o encadeado"
  exit 1
fi
log "DOWNSTREAM_DONE"
