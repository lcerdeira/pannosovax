#!/usr/bin/env python3
"""
Estágio 04c — seleção do bloco de epitopos ESTRUTURALMENTE COMPARTILHADOS.

Este é o bloco que materializa a tese central do PanNosoVax: epitopos que caem em
regiões estruturalmente equivalentes (TM-score alto) entre patógenos distintos,
identificados no estágio 04 por sobreposição estrutural, não por identidade de sequência.

O bloco é montado em DOIS NÍVEIS, porque as evidências de compartilhamento não têm o
mesmo peso e uma alegação uniforme cairia na revisão — basta um epitopo do nível mais
fraco para desmentir a afirmação geral:

  1. `tres_patogenos` — região equivalente nos três. É a alegação mais forte; rende
     poucos epitopos (2), insuficientes para um bloco sozinha;
  2. `par_cruza_gram` — região equivalente entre um Gram-negativo e o pneumococo, que
     não compartilham nem arquitetura de parede nem ancestral próximo.

O par K. pneumoniae x A. baumannii fica FORA: são duas Gammaproteobacteria e respondem
por 2423 das 2682 janelas, então incluí-lo encheria o bloco com o resultado esperado.

O mecanismo por trás do nível 2 é a dobra do transportador ABC de ligação a substrato:
151 das 182 janelas kpsc x spneu e 57 das 77 abau x spneu são ABC-SB x ABC-SB (ver
`04e_analyze_asymmetry.py`). Já kpsc x abau conversa por porina, classe ausente no
Gram-positivo.

Sobre cada candidato aplicamos ainda dois filtros:

  1. sem colisão de k-mer (>=8) com os epitopos já presentes nos blocos B/MHC-I/MHC-II
     — evita contar o mesmo determinante duas vezes e inflar o construto;
  2. sem redundância interna entre os próprios shared, e uma proteína-âncora por vaga.

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
                    help="regiões compartilhadas por PAR de patógenos (04b2)")
    ap.add_argument("--source-3way", default="shared_validated_3way.tsv",
                    help="regiões compartilhadas pelos TRÊS; entram primeiro no bloco")
    ap.add_argument("--out", default="shared_structural_epitopes.tsv")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)

    def crosses_gram(row) -> bool:
        orgs = {row["organism"], *str(row["partner_orgs"]).split("|")}
        return len({GRAM[o] for o in orgs if o in GRAM}) > 1

    # ── bloco em dois níveis ──────────────────────────────────────────────────
    #
    # Nem toda evidência de compartilhamento tem o mesmo peso, e alegar um grau único
    # para o bloco inteiro não sobrevive à revisão: basta um epitopo do nível mais
    # fraco para desmentir a afirmação geral. Então declaramos o nível de cada um.
    #
    #   1. três patógenos — a alegação mais forte, e existe: 33 regiões, 9 validadas.
    #      Sozinha rende apenas 2 epitopos, poucos para um bloco;
    #   2. par cruzando a fronteira Gram — equivalência entre um Gram-negativo e o
    #      pneumococo, que não compartilham parede celular nem ancestral próximo.
    #
    # O par K. pneumoniae x A. baumannii não entra: são duas Gammaproteobacteria e
    # respondem por 2423 das 2682 janelas. Encheria o bloco com o resultado esperado.
    def pair_key(row) -> str:
        # Par não-ordenado (organism, partner) que sustenta a região. Precisa disso
        # separado de "nivel" porque o nível 2 mistura kpsc+spneu e abau+spneu, e
        # ordenar só por profundidade os dois juntos deixa o mais raso desaparecer
        # atrás do mais fundo — no caso, spneu+abau (mais profundo) engolia todas as
        # vagas e kpsc+spneu (182 janelas, o par mais numeroso) não entrava em nenhuma.
        orgs = tuple(sorted({row["organism"], *str(row["partner_orgs"]).split("|")}))
        return "|".join(orgs)

    frames = []
    for fname, tier in [(args.source_3way, "tres_patogenos"),
                        (args.source, "par_cruza_gram")]:
        p = outpath(cfg, "04_shared", fname)
        if not p.exists():
            log.warning("%s ausente — nível '%s' fica de fora", fname, tier)
            continue
        f = pd.read_csv(p, sep="\t")
        if f.empty:
            continue
        f["cruza_gram"] = f.apply(crosses_gram, axis=1)
        if tier == "par_cruza_gram":
            f = f[f["cruza_gram"]]
        f["nivel"] = tier
        f["pair_key"] = f.apply(pair_key, axis=1)
        frames.append(f)
    if not frames:
        raise SystemExit("faltam as tabelas do 04b2 — rode o estágio 04")

    def diversify(f: pd.DataFrame) -> pd.DataFrame:
        """Intercala os pares dentro do nível, cada um ordenado por profundidade.

        Sem isso, o par mais profundo (mais regiões densas) ocupa as vagas todas e um
        par com menos regiões, ainda que real, nunca aparece no construto -- foi o que
        aconteceu com kpsc+spneu (o par mais numeroso, 182 janelas) perante
        spneu+abau. Round-robin por par garante que cada par tenha chance na fila,
        preservando a ordenação por profundidade DENTRO de cada par.
        """
        groups = [g.sort_values("n_epitopos_seguros_na_regiao", ascending=False)
                  for _, g in f.groupby("pair_key", sort=False)]
        if len(groups) <= 1:
            return f.sort_values("n_epitopos_seguros_na_regiao", ascending=False)
        rows, i = [], 0
        while any(len(g) > i for g in groups):
            for g in groups:
                if len(g) > i:
                    rows.append(g.iloc[i])
            i += 1
        return pd.DataFrame(rows)

    frames = [diversify(f) for f in frames]
    d = pd.concat(frames, ignore_index=True)
    for tier, g in d.groupby("nivel", sort=False):
        pares = g["pair_key"].value_counts().to_dict()
        log.info("nível '%s': %d regiões candidatas, pares: %s", tier, len(g), pares)

    # pool de k-mers dos blocos já existentes (evita dupla contagem do mesmo determinante)
    pool = set()
    for e in existing_epitopes(cfg):
        pool |= kmers(e)
    log.info("pool de %d k-mers dos blocos B/MHC existentes", len(pool))

    # A ordem do laço já garante a precedência do nível 1: como `d` vem concatenado com
    # os três patógenos primeiro, um epitopo de par nunca ocupa a vaga de um deles nem
    # bloqueia sua proteína por colisão.
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
                     "nivel": r["nivel"],
                     "pair_key": r["pair_key"],
                     "note": "regiao estruturalmente compartilhada (TM>=0.5)"})
        if len(picked) >= args.n:
            break

    out = pd.DataFrame(rows)
    write_table(out, outpath(cfg, "04_shared", args.out), log)
    log.info("bloco compartilhado: %d epitopos selecionados", len(out))
    if len(out):
        for nivel, g in out.groupby("nivel", sort=False):
            log.info("  nível '%s': %d", nivel, len(g))
        for _, r in out.iterrows():
            log.info("   %-16s  %-22s  %-15s  (%d seguros na regiao)",
                     r["peptide"], r["organisms"], r["nivel"], r["n_safe_in_region"])


if __name__ == "__main__":
    main()
