#!/usr/bin/env python3
"""
Estágio 03e — controle positivo do surfaceome.

Um filtro de superfície só é confiável se recupera os antígenos que a literatura já
estabeleceu para cada patógeno. Este script verifica isso explicitamente, em vez de
deixar a sensibilidade do filtro implícita no número de candidatas.

Foi escrito depois de um falso negativo que passou despercebido: PsaA e PspC estavam
no core proteome do pneumococo, mas o preditor as rotulou `CytoplasmicMembrane` e o
filtro — que herdara a regra de Gram-negativo, onde esse compartimento é a membrana
interna — as descartou. Nenhum número agregado denunciava isso; só a ausência de nomes
conhecidos denuncia.

Os alvos são reconhecidos pela anotação do produto, não pelo accession, porque o
accession varia entre genomas de referência.

Sai 1 se algum antígeno esperado estiver no core mas fora das candidatas.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, get_logger

log = get_logger("03e_antigens")

# Antígenos proteicos com literatura de vacina para cada organismo. Reconhecidos por
# substring na anotação (minúsculas).
KNOWN = {
    "spneu": {
        "PsaA": "adhesin psaa",
        "PspC/CbpA": "pspc domain",
        "pneumolisina": "pneumolysin",
        "LytA": "amidase lyta",
    },
    "kpsc": {
        "OmpA": "ompa",
        "OmpK36": "ompk36",
        "FyuA/siderophore": "siderophore",
        "MrkD/fimbria": "fimbri",
    },
    "abau": {
        "OmpA": "ompa",
        "Omp33-36": "omp33",
        "Bap/adesina": "adhesin",
        "CarO": "caro",
    },
}


def products(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    col = "product" if "product" in df.columns else df.columns[1]
    return df[col].fillna("").astype(str).str.lower()


def main() -> int:
    missing_total = 0
    for org, targets in KNOWN.items():
        cand_path = ROOT / f"results/03_surfaceome/{org}_candidates.tsv"
        full_path = ROOT / f"results/03_surfaceome/{org}_surfaceome_full.tsv"
        if not cand_path.exists() or not full_path.exists():
            log.warning("%s: sem tabelas do estágio 03 — pulando", org)
            continue
        cand, full = products(cand_path), products(full_path)

        hits, absent, not_core = [], [], []
        for name, needle in targets.items():
            if cand.str.contains(needle, regex=False).any():
                hits.append(name)
            elif full.str.contains(needle, regex=False).any():
                absent.append(name)      # está no core, o filtro rejeitou
            else:
                not_core.append(name)    # nem chegou ao core

        log.info("%s: %d/%d antígenos conhecidos recuperados — %s",
                 org, len(hits), len(targets), ", ".join(hits) or "nenhum")
        if not_core:
            log.info("  fora do core (esperado para loci variáveis): %s", ", ".join(not_core))
        if absent:
            missing_total += len(absent)
            log.error("  NO CORE MAS REJEITADOS pelo filtro: %s", ", ".join(absent))

    if missing_total:
        log.error("%d antígeno(s) estabelecido(s) presente(s) no core foram descartados "
                  "pelo surfaceome — o filtro está perdendo superfície real.", missing_total)
        return 1
    log.info("Controle positivo OK.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
