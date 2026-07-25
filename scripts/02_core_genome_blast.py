#!/usr/bin/env python3
"""
Estágio 02 — core genome por ortologia recíproca contra uma referência.

Por que não Panaroo aqui. Panaroo/Roary constroem o pangenoma completo por
clusterização all-vs-all, o que é o correto quando se quer descrever a arquitetura do
pangenoma. Mas nosso objetivo é outro: queremos o conjunto de proteínas **presentes em
quase todos os isolados** para depois filtrar por superfície. Para isso basta uma
abordagem centrada em referência, que é O(N) em vez de O(N²) e roda em minutos.

Método: uma proteína da referência é considerada core se tiver ortólogo
(identidade >= id_min, cobertura >= cov_min) em >= core_threshold dos genomas.

Limitação assumida e reportada: genes ausentes da referência não são descobertos. Para
seleção de antígenos vacinais isso é aceitável — um antígeno que nem sequer existe na
cepa de referência dificilmente seria um bom candidato universal —, mas invalidaria o
uso destes números como descrição do pangenoma. Não os use para isso.

Uso:
    python scripts/02_core_genome_blast.py --organism spneu
    python scripts/02_core_genome_blast.py --all --threads 8
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from common import get_logger, load_config, outpath, write_table

log = get_logger("02_core")

ID_MIN = 90.0        # identidade mínima para considerar ortólogo
COV_MIN = 0.80       # cobertura mínima sobre a proteína de referência


def pick_reference(prot_dir: Path, cfg_ref: str | None) -> Path:
    """Usa a referência do config se ela foi baixada; senão o maior proteoma."""
    if cfg_ref:
        cand = prot_dir / f"{cfg_ref}.faa"
        if cand.exists():
            return cand
    files = sorted(prot_dir.glob("*.faa"), key=lambda p: p.stat().st_size, reverse=True)
    if not files:
        raise SystemExit(f"nenhum proteoma em {prot_dir} — rode o estágio 01b")
    return files[0]


def run_blast(query: Path, db: Path, threads: int) -> pd.DataFrame:
    cmd = ["blastp", "-query", str(query), "-db", str(db),
           "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
           "-evalue", "1e-10", "-max_target_seqs", "1",
           "-num_threads", str(threads)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    if not proc.stdout.strip():
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length",
                                     "qlen", "slen", "evalue", "bitscore"])
    from io import StringIO
    return pd.read_csv(StringIO(proc.stdout), sep="\t", header=None,
                       names=["qseqid", "sseqid", "pident", "length", "qlen",
                              "slen", "evalue", "bitscore"])


def core_for_organism(org: str, cfg: dict, threads: int) -> pd.DataFrame:
    prot_dir = outpath(cfg, "01_genomes", f"{org}_proteomes").parent / f"{org}_proteomes"
    proteomes = sorted(prot_dir.glob("*.faa"))
    ref = pick_reference(prot_dir, cfg["organisms"][org].get("reference"))
    others = [p for p in proteomes if p != ref]
    log.info("%s: referência %s, %d genomas de comparação",
             org, ref.stem, len(others))

    ref_seqs = {r.id: r for r in SeqIO.parse(ref, "fasta")}
    log.info("%s: %d proteínas na referência", org, len(ref_seqs))

    presence: dict[str, int] = defaultdict(int)
    identities: dict[str, list[float]] = defaultdict(list)

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for i, p in enumerate(others, 1):
            db = tmpdir / p.stem
            db.write_bytes(p.read_bytes())
            subprocess.run(["makeblastdb", "-in", str(db), "-dbtype", "prot",
                            "-out", str(db)], check=True, capture_output=True)
            hits = run_blast(ref, db, threads)
            if not hits.empty:
                # Proteínas repetitivas (ex.: adesinas com repetições SSURE) geram
                # múltiplos HSPs contra o mesmo genoma. Sem deduplicar por query, a
                # contagem de presença ultrapassa o número de genomas — observamos
                # frações >1,0 antes desta correção.
                hits = hits.sort_values("bitscore", ascending=False)
                hits = hits.drop_duplicates(subset="qseqid", keep="first")
                hits["cov"] = hits["length"] / hits["qlen"]
                good = hits[(hits["pident"] >= ID_MIN) & (hits["cov"] >= COV_MIN)]
                for qid, pid in zip(good["qseqid"], good["pident"]):
                    presence[qid] += 1
                    identities[qid].append(pid)
            if i % 10 == 0:
                log.info("  %s: %d/%d genomas comparados", org, i, len(others))

    n = len(others)
    thr = cfg["pangenome"]["core_threshold"]
    rows = []
    for pid, rec in ref_seqs.items():
        cnt = presence.get(pid, 0)
        frac = cnt / n if n else 0.0
        ids = identities.get(pid, [])
        rows.append({
            "protein_id": pid,
            "product": rec.description.split(" ", 1)[1] if " " in rec.description else "",
            "length": len(rec.seq),
            "n_genomes_with_ortholog": cnt,
            "n_genomes_compared": n,
            "presence_fraction": round(frac, 4),
            "mean_identity": round(sum(ids) / len(ids), 2) if ids else None,
            "is_core": frac >= thr,
        })
    df = pd.DataFrame(rows).sort_values("presence_fraction", ascending=False)

    core = df[df["is_core"]]
    faa = outpath(cfg, "02_pangenome", f"{org}_core_proteins.faa")
    with open(faa, "w") as fh:
        for pid in core["protein_id"]:
            rec = ref_seqs[pid]
            fh.write(f">{rec.id} {rec.description.split(' ', 1)[-1]}\n{rec.seq}\n")

    write_table(df, outpath(cfg, "02_pangenome", f"{org}_presence.tsv"), log)
    log.info("%s: %d/%d proteínas core (>=%.0f%% dos genomas) -> %s",
             org, len(core), len(df), 100 * thr, faa.name)
    return pd.DataFrame([{
        "organism": org, "reference": ref.stem, "genomes_compared": n,
        "ref_proteins": len(df), "core_proteins": len(core),
        "core_fraction": round(len(core) / len(df), 4) if len(df) else 0,
    }])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    threads = args.threads or cfg.get("threads", 4)
    targets = list(cfg["organisms"]) if args.all else [args.organism]

    summaries = [core_for_organism(org, cfg, threads) for org in targets]
    write_table(pd.concat(summaries, ignore_index=True),
                outpath(cfg, "02_pangenome", "core_summary.tsv"), log)


if __name__ == "__main__":
    main()
