#!/bin/bash
# Gera os insumos da Figura 1 do Application Note e compõe a figura.
# Reprodutível: a captura da UI e o grafo do Snakemake são regerados do zero,
# não são imagens guardadas à mão.
#
# Requisitos: graphviz (`dot`), Google Chrome, e o backend do app.
set -euo pipefail
cd "$(dirname "$0")/../.."
WORK=/tmp/figb
PORT=8765
mkdir -p "$WORK"

PY=${PYTHON:-python3}
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

echo "── 1/3 · grafo de regras do Snakemake ──"
$PY -m snakemake -s workflow/Snakefile --rulegraph 2>/dev/null > "$WORK/rulegraph.dot"
$PY - "$WORK" <<'PY'
import re, sys
w = sys.argv[1]
COLOR = {"all":"#6C7B8A","fetch_genomes":"#8A93A5","download_proteomes":"#8A93A5",
         "core_genome":"#3B6EA5","localization":"#2A9D8F","signal_peptide":"#2A9D8F",
         "topology":"#2A9D8F","surfaceome":"#2A9D8F","epitopes":"#E8973A",
         "homologs":"#E8973A","conservation":"#E8973A","safety":"#D1495B",
         "coverage":"#8A6BBF","construct":"#1D3557","physchem":"#1D3557",
         "expression":"#1D3557","report":"#5B6675"}
src = open(f"{w}/rulegraph.dot").read()
def fix(m):
    i, lab = m.group(1), m.group(2)
    return (f'{i}[label = "{lab}", color = "white", fillcolor = "{COLOR.get(lab,"#8A93A5")}", '
            f'style="rounded,filled", fontcolor="white"]')
out = re.sub(r'(\d+)\[label = "([^"]+)", color = "[^"]+", style="rounded"\]', fix, src)
out = out.replace('graph[bgcolor=white, margin=0];',
                  'graph[bgcolor=white, margin=0.2, rankdir=LR, nodesep=0.28, ranksep=0.5];')
out = out.replace('fontsize=10, penwidth=2', 'fontsize=11, penwidth=0')
out = out.replace('edge[penwidth=2, color=grey]',
                  'edge[penwidth=1.6, color="#B9C2CC", arrowsize=0.7]')
open(f"{w}/rulegraph_styled.dot","w").write(out)
PY
dot -Tpng -Gdpi=200 "$WORK/rulegraph_styled.dot" -o "$WORK/dag.png"

echo "── 2/3 · captura da interface ──"
if ! curl -sS -m 5 -o /dev/null "http://127.0.0.1:$PORT/"; then
  echo "   subindo backend..."
  $PY app/backend.py & BACKPID=$!
  trap 'kill $BACKPID 2>/dev/null || true' EXIT
  sleep 6
fi
"$CHROME" --headless --disable-gpu --screenshot="$WORK/ui_exec.png" \
  --window-size=1180,900 --hide-scrollbars --force-device-scale-factor=2 \
  --virtual-time-budget=4000 "http://127.0.0.1:$PORT/#execucao" 2>/dev/null

echo "── 3/3 · compondo a figura ──"
$PY scripts/report/fig_appnote.py
echo "✓ results/report/figures/appnote_fig1.png"
