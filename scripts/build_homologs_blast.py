#!/usr/bin/env python3
"""Constrói o mapa de homólogos por candidata via blastp (sem panaroo/mafft).

Para cada organismo: concatena os proteomas por genoma (id -> '{genoma}|{orig}'),
faz um BLAST DB, roda blastp das candidatas contra ele e guarda o melhor hit por
genoma (proxy de ortólogo). Saída: results/02_pangenome/{org}_homologs.json
    { protein_id: {"n_genomes": N, "homologs": ["SEQ_do_genoma_1", ...] } }

O denominador de conservação é o nº de genomas do organismo (isolado sem ortólogo
conta contra a conservação — coerente com "presente em X% dos isolados").
"""
from __future__ import annotations
import argparse, json, subprocess, tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, ROOT

log = get_logger("homologs")

PROTEOME_DIR = {org: ROOT / f"results/01_genomes/{org}_proteomes" for org in ("kpsc", "abau", "spneu")}
# gate de ortologia (permissivo — candidatas são core, mas cobre variantes de sequência)
MIN_PIDENT = 60.0
MAX_EVALUE = 1e-10
MIN_QCOV = 60.0


def run(cmd: list[str]) -> None:
    """subprocess com o stderr VISÍVEL em caso de falha.

    Antes isto era `check=True, capture_output=True` puro: quando o blastp falhava, o
    traceback mostrava a linha de comando e engolia a mensagem de erro, que é a única
    parte útil. Custou uma rodada inteira de diagnóstico no HPC.
    """
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        log.error("comando falhou (rc=%d): %s", p.returncode, " ".join(cmd[:3]))
        for line in (p.stderr or p.stdout or "(sem saída)").strip().splitlines()[-10:]:
            log.error("  %s", line)
        raise subprocess.CalledProcessError(p.returncode, cmd, p.stdout, p.stderr)


def write_query(org: str, work: Path) -> Path:
    """FASTA das candidatas ATUAIS, montado a partir do surfaceome + proteoma core.

    Antes este caminho era `results/03b_deeptmhmm/{org}_candidates.faa`, um arquivo que
    sobrou do fluxo piloto e que nenhum estágio escreve — no re-run em escala cheia ele
    simplesmente não existia e o blastp morria nos três organismos. Derivar a consulta
    das mesmas duas fontes que o estágio 05 usa garante que o mapa de homólogos cubra
    exatamente o conjunto de candidatas vigente, em vez de um resquício de outro run.
    """
    cand = pd.read_csv(ROOT / f"results/03_surfaceome/{org}_candidates.tsv", sep="\t")
    keep = set(cand["protein_id"].astype(str))
    faa = ROOT / f"results/02_pangenome/{org}_core_proteins.faa"
    query = work / f"{org}_query.faa"
    n = 0
    with open(query, "w") as out:
        for rec in SeqIO.parse(faa, "fasta"):
            if rec.id in keep:
                out.write(f">{rec.id}\n{rec.seq}\n")
                n += 1
    if n == 0:
        raise SystemExit(f"{org}: nenhuma das {len(keep)} candidatas foi encontrada em {faa.name}")
    if n < len(keep):
        log.warning("%s: %d das %d candidatas sem sequência no proteoma core",
                    org, len(keep) - n, len(keep))
    log.info("%s: consulta com %d candidatas", org, n)
    return query


def build_db(org: str, work: Path) -> tuple[Path, dict[str, str], int]:
    """Concatena proteomas -> FASTA marcado, cria DB, retorna (db, id2seq, n_genomes)."""
    id2seq: dict[str, str] = {}
    genomes = sorted(PROTEOME_DIR[org].glob("*.faa"))
    cat = work / f"{org}_all.faa"
    with open(cat, "w") as out:
        for fa in genomes:
            g = fa.stem
            for rec in SeqIO.parse(fa, "fasta"):
                tag = f"{g}|{rec.id}"
                id2seq[tag] = str(rec.seq)
                out.write(f">{tag}\n{rec.seq}\n")
    run(["makeblastdb", "-in", str(cat), "-dbtype", "prot",
         "-out", str(work / f"{org}_db")])
    log.info("%s: DB com %d proteínas de %d genomas", org, len(id2seq), len(genomes))
    return work / f"{org}_db", id2seq, len(genomes)


def run_blast(query: Path, db: Path, work: Path, threads: int) -> Path:
    out = work / "hits.tsv"
    run([
        "blastp", "-query", str(query), "-db", str(db), "-out", str(out),
        "-outfmt", "6 qseqid sseqid pident evalue bitscore qcovs",
        "-max_target_seqs", "2000", "-evalue", "1e-6", "-num_threads", str(threads),
    ])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--organisms", nargs="+", default=["kpsc", "abau", "spneu"])
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    for org in args.organisms:
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            query = write_query(org, work)
            db, id2seq, n_genomes = build_db(org, work)
            hits = run_blast(query, db, work, args.threads)
            # melhor hit por (query, genoma) por bitscore
            best: dict[str, dict[str, tuple[float, str]]] = {}
            with open(hits) as fh:
                for line in fh:
                    q, s, pid, ev, bits, qcov = line.rstrip("\n").split("\t")
                    if float(pid) < MIN_PIDENT or float(ev) > MAX_EVALUE or float(qcov) < MIN_QCOV:
                        continue
                    genome = s.split("|", 1)[0]
                    bits = float(bits)
                    d = best.setdefault(q, {})
                    if genome not in d or bits > d[genome][0]:
                        d[genome] = (bits, s)
            out = {}
            for q, per_genome in best.items():
                seqs = [id2seq[sid] for _, sid in per_genome.values()]
                out[q] = {"n_genomes": n_genomes, "n_homologs": len(seqs), "homologs": seqs}
            dst = ROOT / f"results/02_pangenome/{org}_homologs.json"
            dst.write_text(json.dumps(out))
            covered = sum(1 for v in out.values() if v["n_homologs"] >= 0.9 * n_genomes)
            log.info("%s: %d candidatas mapeadas; %d com homólogo em >=90%% dos %d genomas -> %s",
                     org, len(out), covered, n_genomes, dst.name)


if __name__ == "__main__":
    main()
