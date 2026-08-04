#!/bin/bash
# PanNosoVax — re-run em escala cheia (encadeado, autônomo).
# Cada estágio é retomável; se o job morrer, resubmeter retoma de onde parou.
# NÃO usa `set -e` global: um estágio que falha loga e o driver segue quando faz
# sentido, mas para nos pré-requisitos duros (core genome, surfaceome).
set -uo pipefail
cd "$HOME/pannosovax"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate pnv-fold

T=${SLURM_CPUS_PER_TASK:-8}
ORGS="kpsc abau spneu"
log(){ echo "[$(date '+%F %T')] $*"; }
stage(){ log "════════ $* ════════"; }

run(){ # run <nome> <cmd...>  — loga, cronometra, NÃO aborta o driver
  local name="$1"; shift
  local t0=$SECONDS
  log "→ $name"
  if "$@"; then log "✓ $name ($((SECONDS-t0))s)"; else log "✗ $name FALHOU (rc=$?)"; return 1; fi
}

stage "01b · download (idempotente — pula existentes)"
run download python scripts/01b_download_proteomes.py --all

stage "02 · core genome (BLAST vs referência)"
for o in $ORGS; do
  if [ -s "results/02_pangenome/${o}_core_proteins.faa" ]; then
    log "✓ core:$o já existe — pulando"
  else
    run "core:$o" python scripts/02_core_genome_blast.py --organism $o --threads $T || exit 2
  fi
done

stage "03 · surfaceome REAL (DeepLocPro + SignalP-6 + DeepTMHMM)"
for o in $ORGS; do
  run "deeploc:$o" python scripts/run_surfaceome_tools.py --organism $o --tool deeploc
  run "signalp:$o" python scripts/run_surfaceome_tools.py --organism $o --tool signalp
  run "tmhmm:$o"   python scripts/run_surfaceome_tools.py --organism $o --tool tmhmm
done
for o in $ORGS; do run "filter:$o" python scripts/03_surfaceome_filter.py --organism $o || exit 3; done

stage "05 · epitopos (IEDB) — o passo LONGO, retomável por cache"
run epitopes python scripts/run_iedb_epitopes.py

stage "05c · conservação (blastp)"
run homologs python scripts/build_homologs_blast.py --threads $T
run conservation python scripts/apply_conservation.py

stage "06 · segurança (self + comensal + mimetismo)"
for o in $ORGS; do for k in mhc1 mhc2 bcell; do
  run "safety:$o:$k" python scripts/06_safety_screen.py --organism $o --class $k
done; done

stage "07 · cobertura HLA (mundo + Brasil)"
for k in mhc1 mhc2; do run "coverage:$k" python scripts/07_population_coverage.py --class $k --n 18 --brazil-weight 0.5; done

stage "08-12 · construto, físico-química, expressão"
run construct  python scripts/08_construct_builder.py
run physchem   python scripts/09_physchem.py
run expression python scripts/12_codon_optimize.py

stage "FIM — resumo"
python - <<'PY'
import pandas as pd, glob, os
def u(f):
    try: d=pd.read_csv(f,sep="\t"); return d["peptide"].nunique() if len(d) else 0
    except: return 0
for o in ["kpsc","abau","spneu"]:
    try:
        pres=pd.read_csv(f"results/02_pangenome/{o}_presence.tsv",sep="\t")
        cand=pd.read_csv(f"results/03_surfaceome/{o}_candidates.tsv",sep="\t")
        safe=sum(u(f"results/06_safety/{o}_{k}_safe.tsv") for k in ["mhc1","mhc2","bcell"])
        print(f"{o}: genomas={int(pres.n_genomes_compared.max())+1} core={int(pres.is_core.sum())} "
              f"candidatas_superficie={len(cand)} epitopos_seguros={safe}")
    except Exception as e: print(f"{o}: incompleto ({e})")
PY
log "PIPELINE_DONE"
