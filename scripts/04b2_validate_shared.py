#!/usr/bin/env python3
"""
Estágio 04b2 — valida as regiões compartilhadas contra os epitopos seguros.

Uma região estruturalmente compartilhada (saída do 04b) só interessa se epitopos que
JÁ passaram pela triagem de segurança caem dentro dela. A região é uma coordenada no
espaço; o epitopo é o que de fato entra no construto.

Este passo existia como manipulação avulsa: `shared_validated_v2.tsv` estava no
repositório sem nenhum script que o produzisse, o que impedia reproduzir o bloco
compartilhado — justamente a tese central do trabalho — e foi como uma tabela derivada
do surfaceome piloto sobreviveu à troca para o conjunto rigoroso sem ninguém notar.

Entrada : results/04_shared/shared_regions_v2.tsv
          results/06_safety/{org}_{classe}_safe.tsv
Saída   : results/04_shared/shared_validated_v2.tsv   (lido pelo 04c e pelo 04d)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, get_logger, write_table

log = get_logger("04b2_validate")

CLASSES = ["mhc1", "mhc2", "bcell"]


def safe_peptides(org: str) -> set[str]:
    peps: set[str] = set()
    for klass in CLASSES:
        path = ROOT / f"results/06_safety/{org}_{klass}_safe.tsv"
        if not path.exists():
            log.warning("%s: %s ausente — rode o estágio 06", org, path.name)
            continue
        df = pd.read_csv(path, sep="\t")
        if len(df):
            peps |= set(df["peptide"].astype(str))
    return peps


def main() -> int:
    src = ROOT / "results/04_shared/shared_regions_v2.tsv"
    if not src.exists():
        log.error("%s ausente — rode o 04b primeiro", src)
        return 2
    regions = pd.read_csv(src, sep="\t")
    log.info("%d regiões compartilhadas na entrada", len(regions))

    safe = {org: safe_peptides(org) for org in regions["organism"].unique()}
    for org, peps in safe.items():
        log.info("  %s: %d peptídeos seguros", org, len(peps))

    rows = []
    for r in regions.itertuples():
        region = str(r.peptide)
        # Um epitopo pertence à região se for substring dela: ambos vêm da mesma
        # proteína, então a contenção em sequência é contenção em posição.
        inside = sorted((p for p in safe.get(r.organism, ()) if p in region), key=len)
        if not inside:
            continue
        rows.append({
            "organism": r.organism,
            "protein": r.protein,
            "region_len": len(region),
            "partner_orgs": r.partner_orgs,
            "n_epitopos_seguros_na_regiao": len(inside),
            "exemplo_epitopo": inside[len(inside) // 2],
        })

    out = pd.DataFrame(rows).sort_values("n_epitopos_seguros_na_regiao", ascending=False)
    write_table(out, ROOT / "results/04_shared/shared_validated_v2.tsv", log)

    if out.empty:
        log.error("nenhuma região compartilhada contém epitopo seguro — o bloco "
                  "compartilhado ficaria vazio")
        return 1
    log.info("%d regiões validadas, de %d proteínas distintas, %d epitopos seguros no total",
             len(out), out["protein"].nunique(), int(out["n_epitopos_seguros_na_regiao"].sum()))
    for org, grp in out.groupby("organism"):
        log.info("  %s: %d regiões (%d proteínas)", org, len(grp), grp["protein"].nunique())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
