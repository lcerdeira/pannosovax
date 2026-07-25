#!/usr/bin/env python3
"""
Estágio 10 — modelo 3D do construto (ColabFold / AlphaFold2).

Por que predizer a estrutura de um construto quimérico é diferente de predizer a de
uma proteína natural: o construto não existe na natureza, então não há MSA profundo.
O AlphaFold vai devolver pLDDT alto nos blocos que correspondem a domínios reais
(adjuvante, epitopos vindos de proteínas conservadas) e pLDDT baixo nos linkers —
e isso é o resultado esperado, não uma falha. Linker flexível deve ter pLDDT baixo;
se um GPGPG aparecer com pLDDT 90 é sinal de que ele está empacotado contra o corpo
da proteína, o que atrapalha o processamento proteolítico do epitopo vizinho.

Por isso o que interessa aqui não é o pLDDT global (número que a literatura reporta
e que diz pouco), e sim o **pLDDT por bloco do construto**, cruzado com o mapa do
estágio 08. Epitopo com pLDDT baixo é epitopo mal apresentado.

Entrada : results/08_construct/construct.fasta
          results/08_construct/construct_map.tsv
Saída   : PDB indicado em --out
          results/10_structure/plddt_by_block.tsv

Uso:
    python scripts/10_structure.py --fasta results/08_construct/construct.fasta \
                                   --out results/10_structure/construct_refined.pdb
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("10_structure")

COLUMNS = ["element", "label", "start", "end", "length", "mean_plddt", "min_plddt"]

NOTEBOOK_HINT = """\
ColabFold não está instalado localmente. Duas opções:

  (a) local:
      pip install "colabfold[alphafold]"
      colabfold_batch --amber --templates --num-recycle 3 {fasta} {outdir}

  (b) Google Colab (gratuito), notebook oficial AlphaFold2.ipynb:
      https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
      Cole a sequência de {fasta}, use num_recycles=3, amber=True, templates=True,
      e salve o PDB de rank 1 em {out}
"""


def run_colabfold(fasta: Path, outdir: Path, out_pdb: Path) -> bool:
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("rodando colabfold_batch — isso leva de minutos a horas")
    proc = subprocess.run(
        ["colabfold_batch", "--amber", "--templates", "--num-recycle", "3",
         str(fasta), str(outdir)],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        log.warning("colabfold_batch falhou: %s", proc.stderr.strip()[-300:])
        return False
    # ColabFold nomeia os modelos com *_rank_001_*.pdb
    ranked = sorted(outdir.glob("*rank_001*.pdb")) or sorted(outdir.glob("*.pdb"))
    if not ranked:
        log.warning("colabfold_batch terminou mas não produziu PDB")
        return False
    shutil.copy2(ranked[0], out_pdb)
    log.info("modelo de melhor rank: %s -> %s", ranked[0].name, out_pdb.name)
    return True


def plddt_per_residue(pdb: Path) -> dict[int, float]:
    """B-factor do CA. Em saída de AlphaFold/ColabFold esse campo carrega o pLDDT."""
    vals: dict[int, float] = {}
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals[int(line[22:26])] = float(line[60:66])
                except ValueError:
                    continue
    return vals


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True, help="caminho do PDB final")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fasta = Path(args.fasta)
    out_pdb = Path(args.out)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    workdir = outpath(cfg, "10_structure", "colabfold", "_").parent

    if out_pdb.exists():
        log.info("%s já existe — pulando a predição", out_pdb.name)
    elif shutil.which("colabfold_batch"):
        run_colabfold(fasta, workdir, out_pdb)
    else:
        hint = NOTEBOOK_HINT.format(fasta=fasta, outdir=workdir, out=out_pdb)
        log.warning("ColabFold ausente — nenhuma estrutura foi predita.\n%s", hint)
        marker = out_pdb.with_suffix(".ABSENT.txt")
        marker.write_text(hint)
        log.warning("marcador escrito em %s", marker)

    plddt_out = outpath(cfg, "10_structure", "plddt_by_block.tsv")
    map_path = outpath(cfg, "08_construct", "construct_map.tsv")

    if not out_pdb.exists() or not map_path.exists():
        falta = "PDB" if not out_pdb.exists() else "construct_map.tsv"
        log.warning("%s ausente — escrevendo %s só com cabeçalho", falta, plddt_out.name)
        write_table(pd.DataFrame(columns=COLUMNS), plddt_out, log)
        return

    plddt = plddt_per_residue(out_pdb)
    if not plddt:
        log.warning("nenhum B-factor legível em %s — o PDB pode não vir do AlphaFold",
                    out_pdb.name)
        write_table(pd.DataFrame(columns=COLUMNS), plddt_out, log)
        return

    cmap = pd.read_csv(map_path, sep="\t")
    rows = []
    for _, blk in cmap.iterrows():
        vals = [plddt[i] for i in range(int(blk["start"]), int(blk["end"]) + 1) if i in plddt]
        rows.append({
            "element": blk["element"], "label": blk["label"],
            "start": int(blk["start"]), "end": int(blk["end"]),
            "length": int(blk["end"]) - int(blk["start"]) + 1,
            "mean_plddt": round(sum(vals) / len(vals), 2) if vals else float("nan"),
            "min_plddt": round(min(vals), 2) if vals else float("nan"),
        })

    df = pd.DataFrame(rows, columns=COLUMNS)
    write_table(df, plddt_out, log)

    eps = df[~df["element"].isin(["linker", "tag"])]
    fracos = eps[eps["mean_plddt"] < 70]
    log.info("pLDDT médio do construto: %.1f", df["mean_plddt"].mean())
    if len(fracos):
        log.warning("%d blocos de epitopo/adjuvante com pLDDT médio <70: %s",
                    len(fracos), ", ".join(fracos["label"].astype(str)))


if __name__ == "__main__":
    main()
