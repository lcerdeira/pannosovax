#!/usr/bin/env python3
"""
Estágio 02 — extração do core genome a partir da saída do Panaroo.

O Panaroo já agrupa genes ortólogos entre centenas de genomas, mas ele não decide
por nós o que é "core". Essa decisão é científica, não técnica:

  * um limiar de 100% (core estrito) é um erro comum na literatura — basta uma
    montagem fragmentada ou um gene truncado por erro de sequenciamento para
    eliminar um antígeno perfeitamente conservado;
  * um limiar frouxo (<90%) traz genes acessórios que existem só em um clado, o
    que derruba a cobertura da vacina justamente nas cepas que não amostramos.

Usamos pangenome.core_threshold (padrão 0.95): presente em >=95% dos genomas.
A presença é contada por genoma, não por cópia — parálogos no mesmo genoma contam
uma vez só, senão genes em famílias expandidas ganhariam vantagem artificial.

A sequência representativa vem do pan_genome_reference.fa (o representante que o
Panaroo escolheu para cada cluster) e é traduzida quando ainda é nucleotídica.
Os alinhamentos por gene são copiados para o estágio 04, que precisa deles para
estimar dN/dS.

Entrada : results/02_pangenome/{org}_panaroo/gene_presence_absence.csv
          results/02_pangenome/{org}_panaroo/pan_genome_reference.fa
          results/02_pangenome/{org}_panaroo/aligned_gene_sequences/
Saída   : results/02_pangenome/{org}_core_proteins.faa
          results/02_pangenome/{org}_gene_presence_absence.csv
          results/02_pangenome/{org}_gene_alignments/

Uso:
    python scripts/02_extract_core.py --organism kpsc
    python scripts/02_extract_core.py --organism spneu --panaroo-dir /outro/caminho
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from common import get_logger, load_config, outpath, write_table

log = get_logger("02_extract_core")

# Colunas de metadados do gene_presence_absence.csv do Panaroo/Roary; tudo o que
# vier depois delas é um genoma.
META_COLS = {
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
    "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
    "Accessory Fragment", "Accessory Order with Fragment", "QC", "Min group size nuc",
    "Max group size nuc", "Avg group size nuc",
}


def presence_matrix(csv_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Lê gene_presence_absence.csv e devolve (tabela, lista de colunas de genoma)."""
    df = pd.read_csv(csv_path, dtype=str, low_memory=False)
    genomes = [c for c in df.columns if c not in META_COLS]
    return df, genomes


def is_nucleotide(seq: str) -> bool:
    letters = set(seq.upper()) - {"-", "N", "X", "*"}
    return bool(letters) and letters <= set("ACGTU")


def to_protein(seq: str) -> str:
    """Traduz se ainda for nucleotídeo; corta para múltiplo de 3 antes de traduzir."""
    if not is_nucleotide(seq):
        return seq.replace("-", "")
    clean = seq.replace("-", "")
    clean = clean[: len(clean) - len(clean) % 3]
    return str(Seq(clean).translate(to_stop=True))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--panaroo-dir", default=None,
                    help="sobrescreve results/02_pangenome/{org}_panaroo")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    org = args.organism
    thr = cfg["pangenome"]["core_threshold"]
    base = outpath(cfg, "02_pangenome", org).parent

    pan = Path(args.panaroo_dir) if args.panaroo_dir else base / f"{org}_panaroo"
    pres_csv = pan / "gene_presence_absence.csv"
    if not pres_csv.exists():
        raise SystemExit(f"saída do Panaroo não encontrada: {pres_csv}. Rode o estágio 02 primeiro.")

    df, genomes = presence_matrix(pres_csv)
    n_genomes = len(genomes)
    if n_genomes == 0:
        raise SystemExit(f"{pres_csv} não tem colunas de genoma reconhecíveis")

    present = df[genomes].notna() & df[genomes].ne("")
    df["n_present"] = present.sum(axis=1)
    df["frac_present"] = df["n_present"] / n_genomes
    df["is_core"] = df["frac_present"] >= thr

    core = df[df["is_core"]]
    log.info("%s: %d genomas, %d clusters, %d core (>=%.0f%%)",
             org, n_genomes, len(df), len(core), 100 * thr)

    write_table(df, base / f"{org}_gene_presence_absence.csv", log)

    # Sequências representativas: o header do pan_genome_reference.fa carrega o
    # nome do cluster, que é a chave de junção com a matriz de presença.
    ref_fa = pan / "pan_genome_reference.fa"
    annot = dict(zip(df["Gene"], df.get("Annotation", pd.Series("", index=df.index)).fillna("")))
    core_names = set(core["Gene"])
    written = 0
    faa = base / f"{org}_core_proteins.faa"

    if ref_fa.exists():
        with open(faa, "w") as out:
            for rec in SeqIO.parse(ref_fa, "fasta"):
                name = rec.description.split(";")[-1].strip() or rec.id
                if name not in core_names:
                    name = rec.id
                if name not in core_names:
                    continue
                prot = to_protein(str(rec.seq))
                if not prot:
                    continue
                out.write(f">{name} {annot.get(name, '')}\n")
                for i in range(0, len(prot), 60):
                    out.write(prot[i:i + 60] + "\n")
                written += 1
        log.info("escrito %s (%d proteínas core)", faa.relative_to(faa.parents[2]), written)
    else:
        log.warning("pan_genome_reference.fa ausente — escrevendo FASTA vazio; "
                    "estágio 03 não terá o que filtrar")
        faa.write_text("")

    # Alinhamentos por gene: só os core interessam ao estágio 04.
    src = pan / "aligned_gene_sequences"
    dst = base / f"{org}_gene_alignments"
    dst.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        copied = 0
        for aln in src.iterdir():
            if aln.suffix in {".fas", ".fa", ".fasta", ".aln"} and aln.stem in core_names:
                shutil.copy2(aln, dst / aln.name)
                copied += 1
        log.info("%s: %d alinhamentos de genes core copiados para %s", org, copied, dst.name)
    else:
        log.warning("aligned_gene_sequences/ ausente (rode panaroo com '-a core') — "
                    "estágio 04 não conseguirá calcular dN/dS")


if __name__ == "__main__":
    main()
