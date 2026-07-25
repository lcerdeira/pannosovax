#!/usr/bin/env python3
"""
Estágio 03b — contagem de hélices transmembrana (DeepTMHMM).

Por que isso é um script separado: o DeepTMHMM não tem binário livremente
redistribuível; a via oficial é a API BioLib (pacote `pybiolib`), que exige rede e
uma conta. O restante do pipeline não pode depender disso para rodar.

Ordem de tentativa:
  1. pybiolib -> DTU/DeepTMHMM (predição real, gff3);
  2. binário local `deeptmhmm` ou `tmhmm` no PATH;
  3. nada disso disponível -> TSV só com cabeçalho + WARNING.

No caso 3 o estágio 03 vê n_tm_helices ausente, preenche com 0 e o filtro de
topologia passa a ser inócuo. Isso é deliberado: preferimos declarar que a
topologia ficou **não filtrada** a inventar contagens de hélices. Quem for
publicar precisa rodar a predição de verdade antes.

Entrada : FASTA de proteínas (results/02_pangenome/{org}_core_proteins.faa)
Saída   : TSV com colunas protein_id, n_tm_helices

Uso:
    python scripts/03b_deeptmhmm_client.py \
        --fasta results/02_pangenome/kpsc_core_proteins.faa \
        --out   results/03_surfaceome/kpsc_tmhmm.tsv
    python scripts/03b_deeptmhmm_client.py --fasta x.faa --out y.tsv --batch-size 200
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from common import get_logger, load_config, write_table

log = get_logger("03b_deeptmhmm")

COLUMNS = ["protein_id", "n_tm_helices"]


def parse_gff3(text: str) -> dict[str, int]:
    """DeepTMHMM gff3: linhas '<id>\\tTMhelix\\t<start>\\t<end>'."""
    counts: dict[str, int] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].strip() == "TMhelix":
            counts[parts[0]] = counts.get(parts[0], 0) + 1
        elif "PredHel=" in line:                       # formato TMHMM 2.0 short
            m = re.search(r"PredHel=(\d+)", line)
            if m:
                counts[parts[0]] = int(m.group(1))
    return counts


def run_biolib(records: list, batch_size: int) -> dict[str, int] | None:
    try:
        import biolib  # noqa: F401  (pybiolib)
    except ImportError:
        log.info("pybiolib não instalado — pulando a API BioLib")
        return None

    import biolib

    try:
        app = biolib.load("DTU/DeepTMHMM")
    except Exception as exc:                            # rede/auth/app indisponível
        log.warning("não foi possível carregar DTU/DeepTMHMM via BioLib: %s", exc)
        return None

    counts: dict[str, int] = {}
    for start in range(0, len(records), batch_size):
        chunk = records[start:start + batch_size]
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
            SeqIO.write(chunk, fh, "fasta")
            tmp = fh.name
        log.info("BioLib: lote %d-%d (%d sequências)",
                 start + 1, start + len(chunk), len(chunk))
        try:
            job = app.cli(args=f"--fasta {tmp}", machine="local")
            out = job.get_stdout().decode("utf-8", "replace")
            counts.update(parse_gff3(out))
        except Exception as exc:
            log.warning("lote falhou na BioLib (%s) — abortando essa via", exc)
            return counts or None
        finally:
            Path(tmp).unlink(missing_ok=True)
    return counts


def run_local(records: list, batch_size: int) -> dict[str, int] | None:
    exe = shutil.which("deeptmhmm") or shutil.which("tmhmm")
    if not exe:
        return None
    log.info("usando binário local: %s", exe)

    counts: dict[str, int] = {}
    for start in range(0, len(records), batch_size):
        chunk = records[start:start + batch_size]
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
            SeqIO.write(chunk, fh, "fasta")
            tmp = fh.name
        try:
            proc = subprocess.run([exe, tmp], capture_output=True, text=True, check=False)
            if proc.returncode != 0:
                log.warning("%s retornou %d: %s", exe, proc.returncode,
                            proc.stderr.strip()[:200])
                return counts or None
            counts.update(parse_gff3(proc.stdout))
        finally:
            Path(tmp).unlink(missing_ok=True)
    return counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch-size", type=int, default=500,
                    help="sequências por submissão; lotes grandes estouram o timeout da API")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    load_config(args.config)                            # valida o config cedo
    fasta = Path(args.fasta)
    out = Path(args.out)

    if not fasta.exists():
        log.warning("FASTA de entrada ausente: %s", fasta)
        records = []
    else:
        records = list(SeqIO.parse(fasta, "fasta"))
    log.info("%d sequências em %s", len(records), fasta.name)

    counts = None
    if records:
        counts = run_biolib(records, args.batch_size)
        if not counts:
            counts = run_local(records, args.batch_size)

    if not counts:
        log.warning("DeepTMHMM indisponível (nem pybiolib nem binário local). "
                    "Escrevendo %s só com cabeçalho — o estágio 03 tratará a "
                    "topologia como NÃO FILTRADA (n_tm_helices=0 para todos).", out.name)
        write_table(pd.DataFrame(columns=COLUMNS), out, log)
        return

    df = pd.DataFrame({"protein_id": list(counts), "n_tm_helices": list(counts.values())})
    df = df.sort_values("n_tm_helices", ascending=False)
    write_table(df, out, log)
    log.info("hélices TM: %d proteínas preditas, %d com >1 hélice",
             len(df), int((df["n_tm_helices"] > 1).sum()))


if __name__ == "__main__":
    main()
