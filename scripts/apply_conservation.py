#!/usr/bin/env python3
"""Aplica conservação por k-mer exato aos epitopos brutos e grava os conservados.

conservation(peptide, protein) = fração dos GENOMAS do organismo cujo homólogo da
proteína contém o peptídeo exato (denominador = nº de isolados; ausência conta contra).
Usa results/02_pangenome/{org}_homologs.json (de build_homologs_blast.py).

Processa todos os results/05_epitopes/{org}_{klass}_raw.tsv que existirem; escreve
{org}_{klass}_conserved.tsv (conservation >= config.epitopes.min_epitope_conservation).
Reexecutável conforme novos _raw.tsv forem terminando.
"""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, load_config, write_table, ROOT

log = get_logger("conservation")


def load_homologs(org: str):
    p = ROOT / f"results/02_pangenome/{org}_homologs.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


def main() -> None:
    cfg = load_config()
    thr = cfg["epitopes"]["min_epitope_conservation"]
    raws = sorted((ROOT / "results/05_epitopes").glob("*_raw.tsv"))
    if not raws:
        log.warning("nenhum _raw.tsv encontrado"); return

    hom_cache: dict[str, dict] = {}
    for raw in raws:
        org, klass, _ = raw.stem.split("_", 2)
        if org not in hom_cache:
            hom_cache[org] = load_homologs(org)
        homologs = hom_cache[org]
        if homologs is None:
            log.warning("%s: sem mapa de homólogos — pulando %s", org, raw.name); continue

        df = pd.read_csv(raw, sep="\t")
        if not len(df):
            write_table(df, ROOT / f"results/05_epitopes/{org}_{klass}_conserved.tsv", log); continue

        # cache por (protein_id, peptide) — muitos alelos repetem o mesmo peptídeo
        seen: dict[tuple, float] = {}
        cons = []
        for pid, pep in zip(df["protein_id"], df["peptide"]):
            key = (pid, pep)
            if key not in seen:
                h = homologs.get(pid)
                if not h:
                    seen[key] = float("nan")
                else:
                    n = h["n_genomes"]
                    hits = sum(1 for s in h["homologs"] if pep in s)
                    seen[key] = hits / n
            cons.append(seen[key])
        df["conservation"] = cons
        df["pass_conservation"] = pd.Series(cons) >= thr
        write_table(df, raw, log)  # regrava raw enriquecido
        cons_df = df[df["pass_conservation"].fillna(False)]
        write_table(cons_df, ROOT / f"results/05_epitopes/{org}_{klass}_conserved.tsv", log)
        n_ep = df.drop_duplicates(["protein_id", "peptide"]).shape[0]
        n_cons = cons_df.drop_duplicates(["protein_id", "peptide"]).shape[0]
        log.info("%s/%s: %d epitopos únicos, %d conservados >=%.0f%% (%d linhas com alelo)",
                 org, klass, n_ep, n_cons, 100 * thr, len(cons_df))


if __name__ == "__main__":
    main()
