#!/bin/bash
# Demo do PanNosoVax — roda as etapas rápidas e autossuficientes num conjunto mínimo.
# Saída isolada em results_demo/ (não toca nos results/ reais).
set -euo pipefail
cd "$(dirname "$0")/.."
CFG=demo/config.demo.yaml

echo "── PanNosoVax · demo (3 genomas/organismo) ──"

# 1) coloca as listas de seleção reduzidas onde o pipeline as procura (results_demo/)
mkdir -p results_demo/01_genomes
cp demo/selection/*_selected.tsv results_demo/01_genomes/

# 2) download dos 9 proteomas (idempotente)
echo "→ download (9 proteomas)"
python scripts/01b_download_proteomes.py --all --config "$CFG"

# 3) core genome por BLAST vs referência
echo "→ core genome"
for org in kpsc abau spneu; do
  python scripts/02_core_genome_blast.py --organism "$org" --config "$CFG" --threads 4
done

echo
echo "── resumo ──"
python - <<'PY'
import pandas as pd, glob, os
for org in ["kpsc","abau","spneu"]:
    f=f"results_demo/02_pangenome/{org}_presence.tsv"
    if os.path.exists(f):
        d=pd.read_csv(f,sep="\t")
        print(f"  {org}: {len(d)} proteínas de referência, {int(d.is_core.sum())} core "
              f"(≥95% de {int(d.n_genomes_compared.max())+1} genomas)")
PY
echo "✓ demo concluído — saída em results_demo/"
