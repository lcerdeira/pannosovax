#!/usr/bin/env python3
"""
Preenche as estruturas que o AlphaFold DB não tem, via ESMFold (API ESMAtlas).

A cobertura do AFDB é desigual entre os organismos: 94/119 em K. pneumoniae e 34/34
em S. pneumoniae, mas apenas 30/101 em A. baumannii. Isso não é aleatório — muitas
candidatas de A. baumannii mapeiam para acessos UniProt recentes, para os quais ainda
não há modelo depositado. O efeito é direto sobre a tese: o crossmatch estrutural só
testa pares para os quais existem as DUAS estruturas, então um organismo mal coberto
reduz o bloco compartilhado por um motivo puramente de disponibilidade de banco.

O ESMFold preenche essa lacuna sem GPU local. O pLDDT sai no B-factor, que é o que o
04b lê — mesmo contrato do AFDB.

Limite da API: 400 resíduos. Acima disso a proteína fica sem estrutura nesta rodada e
é contabilizada no resumo, não silenciada.

Uso:
    python scripts/fold_missing_esmfold.py
    python scripts/fold_missing_esmfold.py --organisms abau --max-len 400
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, get_logger

log = get_logger("esmfold")

API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
RETRIES = 3
PAUSE = 1.0          # cortesia com a API pública entre submissões


def normalize_plddt(pdb: str) -> str:
    """Reescreve o B-factor de 0–1 para 0–100, a escala do AlphaFold DB.

    A API do ESMAtlas devolve pLDDT normalizado entre 0 e 1, enquanto os modelos do
    AlphaFold DB usam 0 a 100. Como o 04b mistura as duas origens num mesmo corte de
    confiança, gravar cada uma na sua escala faz o corte rejeitar todo modelo do
    ESMFold sem dizer nada. Gravamos já convertido para que o arquivo em disco tenha
    uma escala só, independente de quem o leia.
    """
    out, bvals = [], []
    for line in pdb.splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            try:
                bvals.append(float(line[60:66]))
            except ValueError:
                pass
    if not bvals or max(bvals) > 1.0:
        return pdb                      # já está em 0–100
    for line in pdb.splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            try:
                b = float(line[60:66]) * 100.0
                line = f"{line[:60]}{b:6.2f}{line[66:]}"
            except ValueError:
                pass
        out.append(line)
    return "".join(out)


def fold(seq: str) -> str | None:
    for attempt in range(RETRIES):
        try:
            r = requests.post(API, data=seq, timeout=300)
            if r.status_code == 200 and r.text.lstrip().startswith(("HEADER", "ATOM")):
                return r.text
            raise RuntimeError(f"HTTP {r.status_code}: {r.text[:120]}")
        except Exception as exc:  # noqa: BLE001
            wait = 5 * (2 ** attempt)
            log.warning("  ESMFold falhou (%s) — retry em %ds", str(exc)[:120], wait)
            time.sleep(wait)
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--organisms", nargs="+", default=["kpsc", "abau", "spneu"])
    ap.add_argument("--max-len", type=int, default=400, help="limite da API ESMAtlas")
    args = ap.parse_args()

    grand = {"ok": 0, "falhou": 0, "longa": 0}
    for org in args.organisms:
        odir = ROOT / f"results/04_structures/{org}"
        odir.mkdir(parents=True, exist_ok=True)
        cand = pd.read_csv(ROOT / f"results/03_surfaceome/{org}_candidates.tsv", sep="\t")
        keep = set(cand["protein_id"].astype(str))
        seqs = {r.id: str(r.seq) for r in
                SeqIO.parse(ROOT / f"results/02_pangenome/{org}_core_proteins.faa", "fasta")
                if r.id in keep}

        missing = [p for p in keep if not (odir / f"{p}.pdb").exists()]
        todo = [p for p in missing if len(seqs.get(p, "")) <= args.max_len]
        too_long = len(missing) - len(todo)
        log.info("%s: %d de %d candidatas sem estrutura; %d dentro do limite de %d aa",
                 org, len(missing), len(keep), len(todo), args.max_len)

        ok = fail = 0
        for i, pid in enumerate(sorted(todo), 1):
            pdb = fold(seqs[pid])
            if pdb is None:
                fail += 1
                log.error("  %s: sem estrutura", pid)
                continue
            (odir / f"{pid}.pdb").write_text(normalize_plddt(pdb))
            ok += 1
            if i % 10 == 0 or i == len(todo):
                log.info("  %s: %d/%d dobradas", org, i, len(todo))
            time.sleep(PAUSE)

        total = len(list(odir.glob("*.pdb")))
        log.info("%s: +%d ESMFold (%d falhas, %d acima de %d aa) -> %d/%d com estrutura",
                 org, ok, fail, too_long, args.max_len, total, len(keep))
        grand["ok"] += ok; grand["falhou"] += fail; grand["longa"] += too_long

    log.info("total: +%d estruturas, %d falhas, %d longas demais para a API",
             grand["ok"], grand["falhou"], grand["longa"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
