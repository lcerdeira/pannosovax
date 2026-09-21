#!/bin/bash
#SBATCH --job-name=pnv-surf-gp
#SBATCH --partition=amd
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/home/lshlt19/pannosovax/logs/surfgp_%j.out
#SBATCH --mail-type=END,FAIL
# Refaz o surfaceome do pneumococo com o subconjunto AMPLIADO.
#
# O subconjunto submetido ao SignalP-6/DeepTMHMM era {Cellwall, Extracellular} = 63
# proteínas. Isso deixava de fora as 404 de membrana citoplasmática — que em
# Gram-positivo incluem as lipoproteínas EXPOSTAS na superfície, PsaA entre elas. Sem
# submetê-las, não há como o SignalP dizer quais têm lipobox, e o estágio 03 não tem
# como recuperá-las. Agora vão as 467.
#
# As saídas antigas do spneu são movidas para .prev porque as ferramentas pulam se o
# arquivo existe — e o que existe foi calculado sobre o subconjunto estreito.
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh; conda activate pnv-fold
cd ~/pannosovax
log(){ echo "[$(date '+%F %T')] $*"; }

S=results/03_surfaceome
for f in spneu_signalp.tsv spneu_tmhmm.tsv; do
  [ -f "$S/$f" ] && mv -f "$S/$f" "$S/$f.prev" && log "guardado $f -> $f.prev"
done
rm -rf "$S/_cache"/spneu_* 2>/dev/null && log "cache de lotes do spneu limpo"

for t in signalp tmhmm; do
  log "→ spneu/$t (subconjunto ampliado)"
  python scripts/run_surfaceome_tools.py --organism spneu --tool $t || { log "✗ spneu/$t"; exit 1; }
done

# O filtro roda nos três: a exclusão por anotação (ribossomais, chaperonas) vale
# para todos, e kpsc/abau não precisam de nova submissão à nuvem.
for o in kpsc abau spneu; do
  log "→ 03 · filtro $o"
  python scripts/03_surfaceome_filter.py --organism $o || { log "✗ 03/$o"; exit 1; }
done

log "════ CONTROLE POSITIVO ════"
python scripts/03e_check_known_antigens.py || log "!! antígenos conhecidos ainda ausentes"

log "════ RESUMO ════"
for o in kpsc abau spneu; do
  printf "  %-6s candidatas=%s\n" "$o" "$(($(wc -l < $S/${o}_candidates.tsv) - 1))"
done
echo SURFGP_DONE
