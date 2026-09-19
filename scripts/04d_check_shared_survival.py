#!/usr/bin/env python3
"""
Estágio 04d — consistência do bloco compartilhado com o surfaceome em escala cheia.

Por que este script existe. O bloco estruturalmente compartilhado (a tese central do
PanNosoVax) foi derivado no estágio 04b a partir de estruturas preditas para as
candidatas do surfaceome ANTIGO — o de escala piloto, filtrado por palavra-chave de
anotação. O re-run em escala cheia substituiu aquele conjunto por candidatas validadas
de verdade (DeepLocPro + SignalP-6 + DeepTMHMM). Nada garante que as proteínas-fonte do
bloco compartilhado tenham sobrevivido a essa troca.

Se uma delas não for mais candidata de superfície, o epitopo compartilhado que veio dela
não pode ser defendido no manuscrito: estaríamos alegando exposição de superfície para
uma proteína que o nosso próprio filtro, agora rigoroso, rejeita. Isso é falha de
consistência interna, não detalhe de implementação — daí o script sair com código 1.

Este é um teste de consistência, não uma etapa de transformação: não escreve nada em
`shared_structural_epitopes.tsv`. Quem seleciona o bloco é o 04c.

Uso:
    python scripts/04d_check_shared_survival.py
    python scripts/04d_check_shared_survival.py --shared results/04_shared/shared_validated_v2.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, get_logger, write_table

log = get_logger("04d_shared_check")

ORGS = ["kpsc", "abau", "spneu"]


def load_candidates(org: str) -> pd.DataFrame | None:
    """Candidatas de superfície do surfaceome atual, para um organismo."""
    path = ROOT / f"results/03_surfaceome/{org}_candidates.tsv"
    if not path.exists():
        log.error("%s: %s não existe — rode o estágio 03 antes", org, path)
        return None
    return pd.read_csv(path, sep="\t")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shared", default="results/04_shared/shared_validated_v2.tsv",
                    help="tabela de regiões compartilhadas validadas (saída do 04b)")
    ap.add_argument("--out", default="results/04_shared/shared_survival_check.tsv")
    args = ap.parse_args()

    shared_path = ROOT / args.shared
    if not shared_path.exists():
        log.error("%s não existe — rode o estágio 04b antes", shared_path)
        return 2
    shared = pd.read_csv(shared_path, sep="\t")

    # Conjunto de candidatas atual, por organismo. Um organismo sem tabela derruba o
    # teste: seria pior concluir "sobreviveu" por ausência de evidência.
    current: dict[str, set[str]] = {}
    for org in ORGS:
        cand = load_candidates(org)
        if cand is None:
            return 2
        current[org] = set(cand["protein_id"].astype(str))
        log.info("%s: %d candidatas de superfície no surfaceome atual", org, len(current[org]))

    rows = []
    for (org, protein), grp in shared.groupby(["organism", "protein"], sort=False):
        survives = str(protein) in current.get(org, set())
        rows.append({
            "organism": org,
            "protein": protein,
            "n_regioes_compartilhadas": len(grp),
            "n_epitopos_seguros": int(grp["n_epitopos_seguros_na_regiao"].sum()),
            "parceiros": grp["partner_orgs"].iloc[0],
            "ainda_candidata": survives,
        })

    rep = pd.DataFrame(rows).sort_values(
        ["ainda_candidata", "n_epitopos_seguros"], ascending=[True, False])
    write_table(rep, ROOT / args.out, log)

    lost = rep[~rep["ainda_candidata"]]
    kept = rep[rep["ainda_candidata"]]

    log.info("proteínas-fonte do bloco compartilhado: %d", len(rep))
    log.info("  sobreviveram ao surfaceome em escala cheia: %d (%d regiões, %d epitopos)",
             len(kept), int(kept["n_regioes_compartilhadas"].sum()),
             int(kept["n_epitopos_seguros"].sum()))

    if lost.empty:
        log.info("OK — o bloco compartilhado é consistente com o surfaceome atual.")
        return 0

    log.error("  PERDIDAS: %d (%d regiões, %d epitopos ficariam indefensáveis)",
              len(lost), int(lost["n_regioes_compartilhadas"].sum()),
              int(lost["n_epitopos_seguros"].sum()))
    for r in lost.itertuples():
        log.error("    %s %s — %d regiões, %d epitopos seguros",
                  r.organism, r.protein, r.n_regioes_compartilhadas, r.n_epitopos_seguros)
    log.error("Ação: re-rodar 04b sobre as candidatas novas, ou excluir do 04c os "
              "epitopos vindos destas proteínas. NÃO submeter o manuscrito antes.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
