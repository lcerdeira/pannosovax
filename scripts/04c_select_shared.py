#!/usr/bin/env python3
"""
Estágio 04c — seleção do bloco de epitopos ESTRUTURALMENTE COMPARTILHADOS.

Este é o bloco que materializa a tese central do PanNosoVax: epitopos que caem em
regiões estruturalmente equivalentes (TM-score alto) entre patógenos distintos,
identificados no estágio 04 por sobreposição estrutural, não por identidade de sequência.

O critério é compartilhamento por PAR de patógenos, não pelos três. Exigir os três
produzia apenas 2 regiões e 3 epitopos, todos já representados nos blocos B/MHC — bloco
vazio. O gargalo é A. baumannii x S. pneumoniae (4 janelas), e ele não é artefato de
cobertura: dobrar as estruturas de A. baumannii (30 -> 76) não mudou o resultado.
Compartilhamento por par rende 314 regiões em 61 proteínas.

O construto v1 foi montado antes deste bloco existir e por isso o ignorou. Aqui
selecionamos de `shared_validated_v2.tsv` os melhores representantes, com dois filtros:

  1. sem colisão de k-mer (>=8) com os epitopos já presentes nos blocos B/MHC-I/MHC-II
     — evita contar o mesmo determinante duas vezes e inflar o construto;
  2. sem redundância interna entre os próprios shared.

A prioridade é para regiões que cruzam a fronteira Gram (ver `crosses_gram`), depois
por profundidade da região e por proteína distinta.

Saída: results/04_shared/shared_structural_epitopes.tsv  (lido pelo estágio 08)
"""
from __future__ import annotations

import argparse

import pandas as pd

from common import GRAM, get_logger, load_config, outpath, write_table

log = get_logger("04c_shared")

KMER = 8


def kmers(seq: str, k: int = KMER) -> set[str]:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else {seq}


def collides(pep: str, pool_kmers: set[str]) -> bool:
    return bool(kmers(pep) & pool_kmers)


def existing_epitopes(cfg: dict) -> list[str]:
    eps = []
    for klass in ("mhc1", "mhc2"):
        p = outpath(cfg, "07_coverage", f"selected_{klass}.tsv")
        if p.exists():
            eps += pd.read_csv(p, sep="\t")["peptide"].tolist()
    for org in cfg["organisms"]:
        p = outpath(cfg, "06_safety", f"{org}_bcell_safe.tsv")
        if p.exists():
            eps += pd.read_csv(p, sep="\t")["peptide"].tolist()
    return eps


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=5, help="quantos epitopos compartilhados")
    ap.add_argument("--min-len", type=int, default=13)
    ap.add_argument("--source", default="shared_validated_v2.tsv",
                    help="tabela de regiões validadas (04b2); use shared_validated_3way.tsv "
                         "para o bloco exigido nos três patógenos")
    ap.add_argument("--out", default="shared_structural_epitopes.tsv")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    src = outpath(cfg, "04_shared", args.source)
    if not src.exists():
        raise SystemExit(f"faltando {src} — rode o estágio 04 (sobreposição estrutural)")

    d = pd.read_csv(src, sep="\t")

    # Prioridade: regiões que atravessam a fronteira Gram-negativo/Gram-positivo.
    #
    # Sobreposição estrutural entre K. pneumoniae e A. baumannii é esperada — são duas
    # Gammaproteobacteria, e é de onde vem a maior parte das 314 regiões. O que sustenta
    # a tese é a equivalência entre um Gram-negativo e o pneumococo, que compartilham
    # nem parede celular nem ancestral próximo: 76 regiões, contra 235 do par fácil.
    # Ordenar por profundidade da região apenas encheria o bloco com o par trivial.
    def crosses_gram(row) -> bool:
        orgs = {row["organism"], *str(row["partner_orgs"]).split("|")}
        return len({GRAM[o] for o in orgs if o in GRAM}) > 1

    d["cruza_gram"] = d.apply(crosses_gram, axis=1)
    d = d.sort_values(["cruza_gram", "n_epitopos_seguros_na_regiao"],
                      ascending=[False, False])
    log.info("%d regiões candidatas, %d delas cruzando a fronteira Gram",
             len(d), int(d["cruza_gram"].sum()))

    # pool de k-mers dos blocos já existentes (evita dupla contagem do mesmo determinante)
    pool = set()
    for e in existing_epitopes(cfg):
        pool |= kmers(e)
    log.info("pool de %d k-mers dos blocos B/MHC existentes", len(pool))

    picked, seen_prot, rows = [], set(), []
    for _, r in d.iterrows():
        # Preferimos o maior peptídeo seguro da região: o representante mediano
        # reprovava no corte de tamanho regiões que continham peptídeo longo.
        ep = str(r.get("maior_epitopo") or r["exemplo_epitopo"])
        if len(ep) < args.min_len:
            continue
        if r["protein"] in seen_prot:
            continue
        if collides(ep, pool):
            log.info("  descartado %s (colide com bloco existente)", ep)
            continue
        picked.append(ep)
        pool |= kmers(ep)                       # também evita colisão entre shared
        seen_prot.add(r["protein"])
        rows.append({"peptide": ep,
                     "organisms": f"{r['organism']}+{r['partner_orgs']}",
                     "anchor_protein": r["protein"],
                     "n_safe_in_region": int(r["n_epitopos_seguros_na_regiao"]),
                     "cruza_gram": bool(r["cruza_gram"]),
                     "note": "regiao estruturalmente compartilhada (TM>=0.5)"})
        if len(picked) >= args.n:
            break

    out = pd.DataFrame(rows)
    write_table(out, outpath(cfg, "04_shared", args.out), log)
    log.info("bloco compartilhado: %d epitopos selecionados", len(out))
    for _, r in out.iterrows():
        log.info("   %-16s  %s%s  (%d seguros na regiao)",
                 r["peptide"], r["organisms"],
                 "  [cruza Gram]" if r["cruza_gram"] else "",
                 r["n_safe_in_region"])


if __name__ == "__main__":
    main()
