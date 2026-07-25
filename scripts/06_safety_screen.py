#!/usr/bin/env python3
"""
Estágio 06 — triagem negativa de segurança. Este é o filtro que a maioria dos
artigos de vacina multi-epitopo faz pela metade.

Quatro camadas:

  A. Auto-similaridade com o proteoma humano — evita reatividade cruzada óbvia.
     Feito com BLASTp curto (task=blastp-short), que é o modo correto para peptídeos;
     usar blastp padrão em 9-mers é um erro comum e produz falsos negativos.

  B. Mimetismo autoimune por k-mer. Homologia global baixa não garante segurança:
     um 7-mer contíguo idêntico a uma proteína humana já basta para reatividade
     cruzada de células T em vários modelos. Fazemos varredura exata de k-mers.

  C. Similaridade com o **microbioma comensal** respiratório e intestinal.
     Camada ausente em quase toda a literatura. Uma vacina que gera anticorpos
     contra epitopos compartilhados com Neisseria lactamica ou Streptococcus mitis
     pode causar disbiose e abrir nicho — exatamente o que queremos evitar em
     paciente hospitalizado.

  D. Alergenicidade e toxicidade (AllerTOP/AlgPred, ToxinPred) via API.

Saída: results/06_safety/{org}_{class}_safe.tsv + relatório de reprovações por camada.
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("06_safety")


def blast_peptides(peptides: list[str], db: Path, threads: int = 4) -> pd.DataFrame:
    """BLASTp-short dos peptídeos contra uma base de proteínas."""
    if not db.with_suffix(db.suffix + ".phr").exists() and not Path(str(db) + ".phr").exists():
        log.info("criando base BLAST para %s", db.name)
        subprocess.run(["makeblastdb", "-in", str(db), "-dbtype", "prot"], check=True)

    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
        for i, pep in enumerate(peptides):
            fh.write(f">pep{i}\n{pep}\n")
        query = fh.name

    cmd = ["blastp", "-task", "blastp-short", "-query", query, "-db", str(db),
           "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
           "-evalue", "200000", "-max_target_seqs", "5",
           "-num_threads", str(threads)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    Path(query).unlink(missing_ok=True)

    if not proc.stdout.strip():
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length",
                                     "evalue", "bitscore"])
    df = pd.read_csv(pd.io.common.StringIO(proc.stdout), sep="\t", header=None,
                     names=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
    df["pep_index"] = df["qseqid"].str.removeprefix("pep").astype(int)
    return df


def worst_hit(peptides: list[str], db: Path, threads: int) -> pd.Series:
    """Maior identidade **normalizada pela cobertura do epitopo** por peptídeo.

    O blastp-short devolve o melhor HSP local; para peptídeos curtos isso é quase
    sempre uma janela de 5-8 resíduos com pident ~100. Usar pident cru tornaria o
    limiar de 35% sem sentido (todo epitopo o cruza). Medimos, em vez disso, a
    fração do epitopo idêntica a uma janela self/comensal:
        cov_ident(%) = (resíduos idênticos) / (tamanho do epitopo) * 100
                     = pident * align_len / len(peptídeo)
    É isso que dá significado a max_identity_to_self=0.35. 0 se nenhum hit.
    """
    hits = blast_peptides(peptides, db, threads)
    if hits.empty:
        return pd.Series(0.0, index=range(len(peptides)))
    pep_len = {i: len(p) for i, p in enumerate(peptides)}
    hits = hits.copy()
    hits["cov_ident"] = (hits["pident"] * hits["length"]
                         / hits["pep_index"].map(pep_len)).clip(upper=100.0)
    best = hits.groupby("pep_index")["cov_ident"].max()
    return best.reindex(range(len(peptides)), fill_value=0.0)


def build_kmer_set(fasta: Path, k: int) -> set[str]:
    """Todos os k-mers do proteoma humano. ~11M para k=7 — cabe em memória."""
    from Bio import SeqIO
    kmers = set()
    for rec in SeqIO.parse(fasta, "fasta"):
        s = str(rec.seq)
        kmers.update(s[i:i + k] for i in range(len(s) - k + 1))
    log.info("k-mers humanos (k=%d): %d", k, len(kmers))
    return kmers


def max_shared_kmer(peptide: str, human_kmers: set[str], k: int) -> bool:
    """True se o peptídeo compartilha algum k-mer exato com o proteoma humano."""
    return any(peptide[i:i + k] in human_kmers for i in range(len(peptide) - k + 1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--class", dest="klass", required=True,
                    choices=["mhc1", "mhc2", "bcell"])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    saf = cfg["safety"]
    org, klass = args.organism, args.klass
    threads = cfg.get("threads", 4)

    src = outpath(cfg, "05_epitopes", f"{org}_{klass}_conserved.tsv")
    df = pd.read_csv(src, sep="\t").reset_index(drop=True)
    # o mesmo peptídeo aparece em muitas linhas (um por alelo); só BLASTa/varre os únicos
    uniq = sorted(set(df["peptide"].tolist()))
    thr = saf["max_identity_to_self"] * 100
    log.info("%s/%s: triando %d epitopos conservados (%d peptídeos únicos)",
             org, klass, len(df), len(uniq))

    # ── Camada A: proteoma humano ──────────────────────────────────────────────
    human = Path(saf["human_proteome"])
    id_human_u = worst_hit(uniq, human, threads)
    df["id_human"] = df["peptide"].map({p: id_human_u[i] for i, p in enumerate(uniq)})
    df["pass_human"] = df["id_human"] < thr

    # ── Camada B: mimetismo por k-mer exato ────────────────────────────────────
    k = saf["max_shared_kmer"] + 1
    human_kmers = build_kmer_set(human, k)
    shares = {p: max_shared_kmer(p, human_kmers, k) for p in uniq}
    df["shares_human_kmer"] = df["peptide"].map(shares)
    df["pass_mimicry"] = ~df["shares_human_kmer"]

    # ── Camada C: microbioma comensal ──────────────────────────────────────────
    worst_commensal = pd.Series(0.0, index=range(len(uniq)))
    for prot in saf["commensal_proteomes"]:
        p = Path(prot)
        if not p.exists():
            log.warning("proteoma comensal ausente: %s (camada C incompleta)", p)
            continue
        worst_commensal = pd.concat(
            [worst_commensal, worst_hit(uniq, p, threads)], axis=1
        ).max(axis=1)
    df["id_commensal"] = df["peptide"].map({p: worst_commensal[i] for i, p in enumerate(uniq)})
    df["pass_commensal"] = df["id_commensal"] < thr

    # ── Camada D: alergenicidade e toxicidade ──────────────────────────────────
    # Ferramentas web sem API estável; o wrapper grava a fila para submissão em lote
    # e lê o resultado de volta. Ver docs/PROTOCOL.md, seção 6D.
    queue = outpath(cfg, "06_safety", f"{org}_{klass}_allergen_queue.fasta")
    with open(queue, "w") as fh:
        for i, pep in enumerate(uniq):          # fila = peptídeos únicos, na ordem de `uniq`
            fh.write(f">pep{i}\n{pep}\n")
    log.info("fila de alergenicidade/toxicidade escrita: %s (%d únicos)", queue, len(uniq))

    allerg_res = queue.with_name(f"{org}_{klass}_allergen_results.tsv")
    if allerg_res.exists():
        res = pd.read_csv(allerg_res, sep="\t")  # mesma ordem da fila (= uniq)
        alg = {uniq[i]: bool(v) for i, v in enumerate(res["is_allergen"])}
        tox = {uniq[i]: bool(v) for i, v in enumerate(res["is_toxin"])}
        df["is_allergen"] = df["peptide"].map(alg)
        df["is_toxin"] = df["peptide"].map(tox)
        df["pass_allergen"] = ~df["is_allergen"].fillna(False)
        df["pass_toxin"] = ~df["is_toxin"].fillna(False)
    else:
        log.warning("resultados de alergenicidade ainda não disponíveis — "
                    "camada D marcada como PENDENTE")
        df["pass_allergen"] = pd.NA
        df["pass_toxin"] = pd.NA

    gates = ["pass_human", "pass_mimicry", "pass_commensal", "pass_allergen", "pass_toxin"]
    df["safe"] = df[gates].fillna(True).all(axis=1)

    write_table(df, outpath(cfg, "06_safety", f"{org}_{klass}_screened.tsv"), log)
    write_table(df[df["safe"]], outpath(cfg, "06_safety", f"{org}_{klass}_safe.tsv"), log)

    report = pd.DataFrame([{
        "organism": org, "class": klass, "input": len(df),
        **{g.replace("pass_", "failed_"): int((~df[g].fillna(True)).sum()) for g in gates},
        "survivors": int(df["safe"].sum()),
    }])
    write_table(report, outpath(cfg, "06_safety", f"{org}_{klass}_report.tsv"), log)
    log.info("%s/%s: %d/%d epitopos passaram na triagem de segurança",
             org, klass, int(df["safe"].sum()), len(df))


if __name__ == "__main__":
    main()
