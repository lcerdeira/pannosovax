#!/usr/bin/env python3
"""
Estágio 07 — seleção do conjunto mínimo de epitopos que maximiza cobertura HLA.

Problema: temos centenas de epitopos seguros e conservados, mas o construto precisa
ser curto (expressão, custo de síntese, imunodominância). Qual subconjunto escolher?

Isso é um **problema de cobertura máxima** (maximum coverage), NP-difícil. A literatura
de vacinas multi-epitopo tipicamente escolhe "os 10 de melhor score", o que é
demonstravelmente subótimo: os 10 melhores costumam ligar aos mesmos alelos comuns e
deixam populações inteiras descobertas.

Aqui usamos guloso com garantia (1 - 1/e) sobre uma função objetivo que otimiza
**cobertura fenotípica ponderada** — a probabilidade de um indivíduo aleatório da
população ter pelo menos um alelo que apresenta pelo menos um epitopo do conjunto —
com um peso combinando frequência mundial e frequência brasileira.

Uso:
    python scripts/07_population_coverage.py --class mhc1 --n 12 --brazil-weight 0.5
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from common import get_logger, load_config, outpath, read_alleles, write_table

log = get_logger("07_coverage")


def phenotypic_coverage(binding: np.ndarray, freqs: np.ndarray) -> float:
    """
    binding: matriz booleana [n_epitopos x n_alelos]
    freqs  : frequência alélica por alelo

    Sob Hardy-Weinberg, a probabilidade de um indivíduo NÃO carregar o alelo i é
    (1 - f_i)^2. A cobertura fenotípica do conjunto é 1 - prod_i (1 - f_i)^2 sobre
    os alelos cobertos por pelo menos um epitopo do conjunto.
    """
    covered = binding.any(axis=0)
    if not covered.any():
        return 0.0
    return float(1.0 - np.prod((1.0 - freqs[covered]) ** 2))


def greedy_select(binding: np.ndarray, freqs: np.ndarray, n: int) -> list[int]:
    """Guloso: a cada passo adiciona o epitopo com maior ganho marginal de cobertura."""
    chosen: list[int] = []
    remaining = set(range(binding.shape[0]))

    for step in range(min(n, binding.shape[0])):
        best_idx, best_gain, best_cov = None, -1.0, 0.0
        current = phenotypic_coverage(binding[chosen], freqs) if chosen else 0.0
        for idx in remaining:
            cov = phenotypic_coverage(binding[chosen + [idx]], freqs)
            gain = cov - current
            if gain > best_gain:
                best_idx, best_gain, best_cov = idx, gain, cov
        if best_idx is None:
            break
        chosen.append(best_idx)
        remaining.discard(best_idx)
        log.info("  passo %2d: epitopo #%d, cobertura=%.4f (+%.4f)",
                 step + 1, best_idx, best_cov, best_gain)
        if best_gain <= 1e-9:
            log.info("  ganho marginal esgotado — parando em %d epitopos", len(chosen))
            break
    return chosen


def build_binding_matrix(df: pd.DataFrame, alleles: list[str]):
    """Matriz epitopo x alelo a partir da tabela longa do IEDB."""
    peptides = df["peptide"].drop_duplicates().tolist()
    pidx = {p: i for i, p in enumerate(peptides)}
    aidx = {a: i for i, a in enumerate(alleles)}
    mat = np.zeros((len(peptides), len(alleles)), dtype=bool)
    for pep, allele in zip(df["peptide"], df["allele"]):
        if allele in aidx:
            mat[pidx[pep], aidx[allele]] = True
    return mat, peptides


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--class", dest="klass", required=True, choices=["mhc1", "mhc2"])
    ap.add_argument("--n", type=int, default=12, help="número máximo de epitopos")
    ap.add_argument("--brazil-weight", type=float, default=0.5,
                    help="0=só frequência mundial, 1=só frequência brasileira")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    klass = args.klass

    frames = []
    for org in cfg["organisms"]:
        p = outpath(cfg, "06_safety", f"{org}_{klass}_safe.tsv")
        if p.exists():
            frames.append(pd.read_csv(p, sep="\t"))
        else:
            log.warning("faltando %s", p)
    if not frames:
        raise SystemExit("nenhum arquivo de epitopos seguros encontrado — rode o estágio 06")
    df = pd.concat(frames, ignore_index=True)

    if cfg.get("construct", {}).get("exclude_free_cysteine"):
        before = df["peptide"].nunique()
        df = df[~df["peptide"].str.contains("C")]
        log.info("exclui epitopos com cisteína: %d -> %d peptídeos únicos",
                 before, df["peptide"].nunique())

    alle = read_alleles(cfg["epitopes"][klass]["alleles_file"])
    w = args.brazil_weight
    alle["freq_mix"] = (1 - w) * alle["freq_world"] + w * alle["freq_brazil"]

    results = []
    for org in df["organism"].unique():
        sub = df[df["organism"] == org]
        mat, peptides = build_binding_matrix(sub, alle["allele"].tolist())
        log.info("%s/%s: %d epitopos únicos, matriz %s", org, klass, len(peptides), mat.shape)

        per_org_n = max(1, args.n // df["organism"].nunique())
        idx = greedy_select(mat, alle["freq_mix"].to_numpy(), per_org_n)

        chosen = [peptides[i] for i in idx]
        cov_world = phenotypic_coverage(mat[idx], alle["freq_world"].to_numpy())
        cov_br = phenotypic_coverage(mat[idx], alle["freq_brazil"].to_numpy())

        for rank, pep in enumerate(chosen, 1):
            row = sub[sub["peptide"] == pep].iloc[0]
            results.append({
                "organism": org, "class": klass, "rank": rank, "peptide": pep,
                "protein_id": row["protein_id"],
                "conservation": row.get("conservation"),
                "n_alleles": int(mat[peptides.index(pep)].sum()),
                "set_coverage_world": cov_world,
                "set_coverage_brazil": cov_br,
            })
        log.info("%s/%s: cobertura mundial=%.1f%%, Brasil=%.1f%% com %d epitopos",
                 org, klass, 100 * cov_world, 100 * cov_br, len(chosen))

    out = pd.DataFrame(results)
    write_table(out, outpath(cfg, "07_coverage", f"selected_{klass}.tsv"), log)


if __name__ == "__main__":
    main()
