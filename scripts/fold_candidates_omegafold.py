#!/usr/bin/env python3
"""Fola candidatas de superfície via OmegaFold (BioLib cloud) para o estágio 04b.

OmegaFold é single-sequence (sem MSA) — adequado aqui, onde só queremos a DOBRA para
comparar por TM-align. Prioriza classes funcionais cujas dobras podem convergir entre
organismos (porina, TonB, ABC substrate-binding, adesina, lipoproteína).

Saída: results/04_structures/{org}/{protein_id}.pdb   (retomável — pula PDBs existentes)

Uso:
    python scripts/fold_candidates_omegafold.py --per-org 2      # piloto
    python scripts/fold_candidates_omegafold.py --per-org 8 --max-len 600
"""
from __future__ import annotations
import argparse, os, shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, ROOT

log = get_logger("fold")

PRIORITY = ["porin", "tonb", "abc", "substrate-binding", "adhesin", "adesina",
            "siderophore", "lipoprotein", "outer membrane", "pilus", "fimbri"]


def rank_candidates(org: str, per_org: int, max_len: int) -> list[tuple[str, str]]:
    cand = pd.read_csv(ROOT / f"results/03_surfaceome/{org}_candidates.tsv", sep="\t")
    keep = set(cand["protein_id"])
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / f"results/02_pangenome/{org}_core_proteins.faa", "fasta")
            if r.id in keep and len(r.seq) <= max_len}
    prod = dict(zip(cand["protein_id"], cand["product"].astype(str).str.lower()))
    # score de prioridade por palavra-chave funcional; desempata por presença/identidade
    def score(pid):
        p = prod.get(pid, "")
        return sum(1 for kw in PRIORITY if kw in p)
    ordered = sorted(seqs.keys(), key=lambda p: (-score(p), len(seqs[p])))
    return [(p, seqs[p]) for p in ordered[:per_org]]


def fold_one(pid: str, seq: str, out_pdb: Path) -> bool:
    import biolib
    work = out_pdb.parent
    fa = work / f"{pid}.faa"
    fa.write_text(f">{pid}\n{seq}\n")
    prev = os.getcwd()
    os.chdir(work)
    try:
        app = biolib.load("protein-tools/omegafold")
        job = app.cli(args=f"{pid}.faa")
        job.get_stdout()
        tmp = work / f"_out_{pid}"
        if tmp.exists():
            shutil.rmtree(tmp)
        job.save_files(str(tmp))
        pdbs = list(tmp.glob("*.pdb")) or list(tmp.rglob("*.pdb"))
        if not pdbs:
            log.warning("%s: OmegaFold não produziu PDB", pid)
            return False
        shutil.copy2(pdbs[0], out_pdb)
        shutil.rmtree(tmp)
        return True
    except Exception as exc:  # noqa: BLE001
        log.error("%s: folding falhou — %s", pid, str(exc)[:200])
        return False
    finally:
        os.chdir(prev)
        fa.unlink(missing_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-org", type=int, default=2)
    ap.add_argument("--max-len", type=int, default=700)
    ap.add_argument("--organisms", nargs="+", default=["kpsc", "abau", "spneu"])
    args = ap.parse_args()

    for org in args.organisms:
        odir = ROOT / f"results/04_structures/{org}"
        odir.mkdir(parents=True, exist_ok=True)
        picks = rank_candidates(org, args.per_org, args.max_len)
        log.info("%s: folding %d candidatas (<=%d aa)", org, len(picks), args.max_len)
        for pid, seq in picks:
            out_pdb = odir / f"{pid}.pdb"
            if out_pdb.exists() and out_pdb.stat().st_size > 0:
                log.info("  %s já existe — pulando", pid); continue
            ok = fold_one(pid, seq, out_pdb)
            log.info("  %s (%d aa): %s", pid, len(seq), "OK" if ok else "FALHOU")


if __name__ == "__main__":
    main()
