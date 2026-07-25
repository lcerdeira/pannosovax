#!/usr/bin/env python3
"""
Estágio 01b — download dos proteomas anotados.

Atalho deliberado: as montagens RefSeq já vêm com anotação PGAP e arquivo
`protein.faa` pronto. Reanotar tudo com Bakta/Prokka custaria dias de CPU e
mudaria pouco o resultado para proteínas de superfície do core genome, que são
justamente as bem anotadas. Usamos a anotação RefSeq no primeiro passe e
deixamos a reanotação como refinamento opcional.

Isso torna o estágio 02 viável em horas em vez de dias.

Uso:
    python scripts/01b_download_proteomes.py --organism spneu --n 60
    python scripts/01b_download_proteomes.py --all --n 60
"""
from __future__ import annotations

import argparse
import io
import time
import zipfile

import pandas as pd
import requests

from common import get_logger, load_config, outpath

log = get_logger("01b_download")

DL = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
      "{acc}/download?include_annotation_type=PROT_FASTA")


def fetch_one(acc: str, dest, retries: int = 3) -> bool:
    out = dest / f"{acc}.faa"
    if out.exists() and out.stat().st_size > 1000:
        return True
    for attempt in range(retries):
        try:
            r = requests.get(DL.format(acc=acc), timeout=180)
            r.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                names = [n for n in z.namelist() if n.endswith("protein.faa")]
                if not names:
                    log.warning("%s: sem protein.faa no pacote", acc)
                    return False
                data = z.read(names[0])
            out.write_bytes(data)
            return True
        except Exception as exc:                      # noqa: BLE001
            wait = 3 * (2 ** attempt)
            log.warning("%s falhou (%s), retry em %ds", acc, type(exc).__name__, wait)
            time.sleep(wait)
    return False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--n", type=int, default=None,
                    help="sobrescreve o nº de genomas SÓ para piloto/teste. Por padrão usa "
                         "organisms.<org>.n_genomes do config — a config é a fonte da verdade.")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    targets = list(cfg["organisms"]) if args.all else [args.organism]

    for org in targets:
        sel = pd.read_csv(outpath(cfg, "01_genomes", f"{org}_selected.tsv"), sep="\t")
        # O nº de genomas vem do config; --n é override explícito de piloto.
        # Um default silencioso aqui (era 60) fez a análise rodar em escala de piloto
        # enquanto o config declarava 400/300/400 — números do manuscrito não batiam.
        n = args.n if args.n is not None else int(cfg["organisms"][org].get("n_genomes", 0))
        if args.n is not None:
            log.warning("%s: MODO PILOTO — n=%d por --n (config pede %s)", org, args.n,
                        cfg["organisms"][org].get("n_genomes"))
        if n and n < len(sel):
            # amostra equilibrada por região, respeitando o alvo
            sel = (sel.groupby("region", group_keys=False)
                      .apply(lambda g: g.head(max(1, n // sel["region"].nunique())))
                      .head(n))
        log.info("%s: alvo=%d genomas | disponíveis na seleção=%d", org, n or len(sel), len(sel))
        dest = outpath(cfg, "01_genomes", f"{org}_proteomes").parent / f"{org}_proteomes"
        dest.mkdir(parents=True, exist_ok=True)

        log.info("%s: baixando %d proteomas -> %s", org, len(sel), dest)
        ok = 0
        for i, acc in enumerate(sel["accession"], 1):
            if fetch_one(acc, dest):
                ok += 1
            if i % 10 == 0:
                log.info("  %s: %d/%d (%d ok)", org, i, len(sel), ok)
            time.sleep(0.35)          # cortesia com a API
        log.info("%s: %d/%d proteomas obtidos", org, ok, len(sel))


if __name__ == "__main__":
    main()
