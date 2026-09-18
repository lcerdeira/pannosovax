#!/usr/bin/env python3
"""
Estágio 04c — seleção do bloco de epitopos ESTRUTURALMENTE COMPARTILHADOS.

Este é o bloco que materializa a tese central do PanNosoVax: epitopos que caem em
regiões estruturalmente equivalentes (TM-score alto) nas proteínas dos três patógenos,
identificados no estágio 04 (sobreposição estrutural, não identidade de sequência).

O construto v1 foi montado antes deste bloco existir e por isso o ignorou. Aqui
selecionamos de `shared_validated_v2.tsv` os melhores representantes, com dois filtros:

  1. sem colisão de k-mer (>=8) com os epitopos já presentes nos blocos B/MHC-I/MHC-II
     — evita contar o mesmo determinante duas vezes e inflar o construto;
  2. sem redundância interna entre os próprios shared.

Prioriza regiões com maior número de epitopos seguros (mais "profundas") e proteínas
distintas, para maximizar a diversidade estrutural do bloco.

Saída: results/04_shared/shared_structural_epitopes.tsv  (lido pelo estágio 08)
"""
from __future__ import annotations

import argparse

import pandas as pd

from common import get_logger, load_config, outpath, write_table

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
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    src = outpath(cfg, "04_shared", "shared_validated_v2.tsv")
    if not src.exists():
        raise SystemExit(f"faltando {src} — rode o estágio 04 (sobreposição estrutural)")

    d = pd.read_csv(src, sep="\t").sort_values(
        "n_epitopos_seguros_na_regiao", ascending=False)

    # pool de k-mers dos blocos já existentes (evita dupla contagem do mesmo determinante)
    pool = set()
    for e in existing_epitopes(cfg):
        pool |= kmers(e)
    log.info("pool de %d k-mers dos blocos B/MHC existentes", len(pool))

    picked, seen_prot, rows = [], set(), []
    for _, r in d.iterrows():
        ep = str(r["exemplo_epitopo"])
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
                     "note": "regiao estruturalmente compartilhada (TM>=0.5)"})
        if len(picked) >= args.n:
            break

    out = pd.DataFrame(rows)
    write_table(out, outpath(cfg, "04_shared", "shared_structural_epitopes.tsv"), log)
    log.info("bloco compartilhado: %d epitopos selecionados", len(out))
    for _, r in out.iterrows():
        log.info("   %-16s  %s  (%d seguros na regiao)",
                 r["peptide"], r["organisms"], r["n_safe_in_region"])


if __name__ == "__main__":
    main()
