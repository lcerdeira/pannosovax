#!/usr/bin/env python3
"""
Estágio 05 — predição de epitopos B, CD4 (MHC-II) e CD8 (MHC-I).

Usa as APIs públicas e gratuitas do IEDB (nextgen-tools) — sem licença, sem custo.
Se as ferramentas locais (netMHCpan, bepipred) estiverem instaladas, usa-as
preferencialmente porque é muito mais rápido para milhares de sequências.

Regra de ouro do projeto: um epitopo só sobrevive se for **conservado**. Predizemos
sobre a sequência de referência e depois verificamos, no alinhamento do core genome,
em que fração dos genomas aquele k-mer exato aparece. Epitopo predito com nota
excelente mas presente em 60% dos isolados é inútil na prática.

Uso:
    python scripts/05_epitope_predict.py --organism kpsc --class mhc1
    python scripts/05_epitope_predict.py --organism spneu --class bcell
"""
from __future__ import annotations

import argparse
import io
import time

import pandas as pd
import requests
from Bio import SeqIO

from common import get_logger, load_config, outpath, read_alleles, write_table

log = get_logger("05_epitopes")

IEDB_MHC1 = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
IEDB_MHC2 = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
IEDB_BCELL = "http://tools-cluster-interface.iedb.org/tools_api/bcell/"


def _post(url: str, data: dict, retries: int = 4) -> pd.DataFrame:
    for attempt in range(retries):
        try:
            resp = requests.post(url, data=data, timeout=300)
            resp.raise_for_status()
            text = resp.text
            if text.lower().startswith("error"):
                raise RuntimeError(text[:300])
            return pd.read_csv(io.StringIO(text), sep="\t")
        except Exception as exc:                       # noqa: BLE001
            wait = 5 * (2 ** attempt)
            log.warning("IEDB falhou (%s), retry em %ds", exc, wait)
            time.sleep(wait)
    raise RuntimeError(f"IEDB indisponível após {retries} tentativas: {url}")


def predict_mhc(seqs: dict, cfg_cls: dict, klass: str) -> pd.DataFrame:
    alleles = read_alleles(cfg_cls["alleles_file"])["allele"].tolist()
    url = IEDB_MHC1 if klass == "mhc1" else IEDB_MHC2
    lengths = cfg_cls["lengths"]
    frames = []

    for i, (pid, rec) in enumerate(seqs.items(), 1):
        fasta = f">{pid}\n{rec.seq}\n"
        payload = {
            "method": cfg_cls["method"],
            "sequence_text": fasta,
            "allele": ",".join(alleles),
            "length": ",".join(str(x) for x in lengths) if klass == "mhc1" else None,
        }
        payload = {k: v for k, v in payload.items() if v is not None}
        df = _post(url, payload)
        df["protein_id"] = pid
        frames.append(df)
        if i % 25 == 0:
            log.info("  %d/%d proteínas processadas", i, len(seqs))

    out = pd.concat(frames, ignore_index=True)
    rank_col = next((c for c in out.columns if "percentile" in c.lower()), None)
    if rank_col:
        out = out[out[rank_col] <= cfg_cls["percentile_rank_max"]]
        out = out.rename(columns={rank_col: "percentile_rank"})
    return out


def predict_bcell(seqs: dict, cfg_b: dict) -> pd.DataFrame:
    frames = []
    for pid, rec in seqs.items():
        df = _post(IEDB_BCELL, {"method": cfg_b["linear_method"],
                                "sequence_text": str(rec.seq)})
        df["protein_id"] = pid
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    score_col = next((c for c in out.columns if "score" in c.lower()), None)
    if score_col:
        out = out[out[score_col] >= cfg_b["linear_threshold"]]
    return out


def conservation(peptide: str, homologs: list[str]) -> float:
    """Fração de homólogos do core genome que contêm o peptídeo exato."""
    if not homologs:
        return float("nan")
    return sum(1 for h in homologs if peptide in h) / len(homologs)


def load_homologs(cfg: dict, org: str) -> dict[str, list[str]]:
    """Lê os alinhamentos por gene do pangenoma -> {protein_id: [seqs dos isolados]}."""
    aln_dir = outpath(cfg, "02_pangenome", org).parent / f"{org}_gene_alignments"
    if not aln_dir.exists():
        log.warning("alinhamentos ausentes em %s — conservação ficará NaN", aln_dir)
        return {}
    homologs = {}
    for fa in aln_dir.glob("*.fasta"):
        homologs[fa.stem] = [str(r.seq).replace("-", "") for r in SeqIO.parse(fa, "fasta")]
    log.info("%s: alinhamentos carregados para %d genes", org, len(homologs))
    return homologs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--class", dest="klass", required=True,
                    choices=["mhc1", "mhc2", "bcell"])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    org, klass = args.organism, args.klass

    cand = pd.read_csv(outpath(cfg, "03_surfaceome", f"{org}_candidates.tsv"), sep="\t")
    faa = outpath(cfg, "02_pangenome", f"{org}_core_proteins.faa")
    keep = set(cand["protein_id"])
    seqs = {r.id: r for r in SeqIO.parse(faa, "fasta") if r.id in keep}
    log.info("%s/%s: predizendo epitopos em %d proteínas candidatas", org, klass, len(seqs))

    if klass in {"mhc1", "mhc2"}:
        df = predict_mhc(seqs, cfg["epitopes"][klass], klass)
        pep_col = next(c for c in df.columns if c.lower() in {"peptide", "sequence"})
    else:
        df = predict_bcell(seqs, cfg["epitopes"]["bcell"])
        pep_col = next(c for c in df.columns if c.lower() in {"peptide", "sequence"})

    df = df.rename(columns={pep_col: "peptide"})

    homologs = load_homologs(cfg, org)
    df["conservation"] = [
        conservation(p, homologs.get(pid, [])) for p, pid in zip(df["peptide"], df["protein_id"])
    ]
    thr = cfg["epitopes"]["min_epitope_conservation"]
    df["pass_conservation"] = df["conservation"].fillna(0) >= thr
    df["organism"] = org
    df["epitope_class"] = klass

    out = outpath(cfg, "05_epitopes", f"{org}_{klass}_raw.tsv")
    write_table(df, out, log)
    write_table(df[df["pass_conservation"]],
                outpath(cfg, "05_epitopes", f"{org}_{klass}_conserved.tsv"), log)
    log.info("%s/%s: %d epitopos preditos, %d conservados (>=%.0f%%)",
             org, klass, len(df), int(df["pass_conservation"].sum()), 100 * thr)


if __name__ == "__main__":
    main()
