#!/usr/bin/env python3
"""
Estágio 09 — caracterização físico-química do construto.

Antes de gastar tempo de GPU com predição de estrutura, verificamos se a proteína é
sequer fabricável. Critérios de aprovação (todos com base na literatura de expressão
heteróloga):

  índice de instabilidade < 40      → estável in vitro (Guruprasad et al.)
  GRAVY < 0                         → hidrofílica, solúvel em citoplasma
  meia-vida estimada (E. coli) alta → não degrada durante a expressão
  pI fora da faixa 6.5–7.5          → evita precipitação em pH fisiológico durante purificação
  antigenicidade (VaxiJen) > 0.5    → limiar bacteriano padrão

Uso:
    python scripts/09_physchem.py --fasta results/08_construct/construct.fasta
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from common import get_logger, load_config, outpath, write_table

log = get_logger("09_physchem")

# Meia-vida N-terminal em E. coli (regra do N-end, Tobias et al. 1991)
NEND_ECOLI_STABLE = set("MGASTVP")


def characterize(seq: str) -> dict:
    pa = ProteinAnalysis(seq)
    aa = pa.count_amino_acids()
    return {
        "length": len(seq),
        "molecular_weight_kda": round(pa.molecular_weight() / 1000, 2),
        "theoretical_pi": round(pa.isoelectric_point(), 2),
        "instability_index": round(pa.instability_index(), 2),
        "aliphatic_index": round(
            100 * (aa["A"] + 2.9 * aa["V"] + 3.9 * (aa["I"] + aa["L"])) / len(seq), 2
        ),
        "gravy": round(pa.gravy(), 3),
        "aromaticity": round(pa.aromaticity(), 3),
        "n_negative_DE": aa["D"] + aa["E"],
        "n_positive_KR": aa["K"] + aa["R"],
        "n_cysteine": aa["C"],
        "helix_frac": round(pa.secondary_structure_fraction()[0], 3),
        "turn_frac": round(pa.secondary_structure_fraction()[1], 3),
        "sheet_frac": round(pa.secondary_structure_fraction()[2], 3),
        "ecoli_halflife_class": "estável (>10 h)" if seq[0] in NEND_ECOLI_STABLE
                                else "instável (<2 min)",
    }


def verdict(p: dict) -> dict:
    checks = {
        "estavel": p["instability_index"] < 40,
        "soluvel": p["gravy"] < 0,
        "termoestavel": p["aliphatic_index"] > 70,
        "pi_fora_da_faixa_neutra": not (6.5 <= p["theoretical_pi"] <= 7.5),
        "sem_cisteina_livre_impar": p["n_cysteine"] % 2 == 0,
        "expressavel_ecoli": p["ecoli_halflife_class"].startswith("estável"),
    }
    checks["APROVADO"] = all(checks.values())
    return checks


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fasta = Path(args.fasta) if args.fasta else outpath(cfg, "08_construct", "construct.fasta")
    if not fasta.exists():
        raise SystemExit(f"construto não encontrado: {fasta} (rode o estágio 08)")

    rows = []
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper().replace("*", "")
        props = characterize(seq)
        checks = verdict(props)
        rows.append({"id": rec.id, **props, **checks})

        log.info("── %s ──", rec.id)
        for k, v in props.items():
            log.info("   %-26s %s", k, v)
        for k, v in checks.items():
            log.info("   %-26s %s", k, "OK" if v else "FALHOU")

    df = pd.DataFrame(rows)
    write_table(df, outpath(cfg, "09_physchem", f"{fasta.stem}_properties.tsv"), log)

    if not df["APROVADO"].all():
        log.warning("construto reprovado em pelo menos um critério — "
                    "considere trocar a ordem dos blocos ou o adjuvante e remontar")


if __name__ == "__main__":
    main()
