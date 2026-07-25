#!/usr/bin/env python3
"""Busca estruturas prontas do AlphaFold DB (EBI) para as candidatas do estágio 04b.

Em vez de folding local (GPU/BioLib inviável), mapeia RefSeq WP_ -> UniProt e baixa o
modelo AlphaFold pré-computado. O pLDDT vem no B-factor do PDB — exatamente o que o
`04b_structural_crossmatch.py` lê. Retomável (pula PDBs existentes).

Prioriza classes funcionais cujas dobras podem convergir entre organismos.
Saída: results/04_structures/{org}/{protein_id}.pdb

Uso:
    python scripts/fetch_structures_afdb.py --per-org 12 --max-len 800
"""
from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, ROOT

log = get_logger("afdb")
PRIORITY = ["porin", "tonb", "abc", "substrate-binding", "adhesin", "adesina",
            "siderophore", "lipoprotein", "outer membrane", "pilus", "fimbri"]


def uniprot_accs(wp: str, n: int = 5) -> list[str]:
    q = f"xref:refseq-{wp}"
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={requests.utils.quote(q)}"
                     f"&fields=accession&format=tsv&size={n}", timeout=30)
    return [l for l in r.text.strip().split("\n")[1:] if l]


def afdb_pdb_url(acc: str) -> str | None:
    try:
        r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=30)
        if r.status_code == 200 and r.json():
            return r.json()[0].get("pdbUrl")
    except Exception:  # noqa: BLE001
        pass
    return None


def rank_candidates(org: str, per_org: int, max_len: int):
    cand = pd.read_csv(ROOT / f"results/03_surfaceome/{org}_candidates.tsv", sep="\t")
    keep = set(cand["protein_id"])
    lens = {r.id: len(r.seq) for r in SeqIO.parse(ROOT / f"results/02_pangenome/{org}_core_proteins.faa", "fasta")
            if r.id in keep}
    prod = dict(zip(cand["protein_id"], cand["product"].astype(str).str.lower()))
    def score(pid):
        return sum(1 for kw in PRIORITY if kw in prod.get(pid, ""))
    ids = [p for p in cand["protein_id"] if lens.get(p, 1e9) <= max_len]
    ids.sort(key=lambda p: (-score(p), lens.get(p, 0)))
    return ids[:per_org]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-org", type=int, default=12)
    ap.add_argument("--max-len", type=int, default=800)
    ap.add_argument("--organisms", nargs="+", default=["kpsc", "abau", "spneu"])
    args = ap.parse_args()

    for org in args.organisms:
        odir = ROOT / f"results/04_structures/{org}"
        odir.mkdir(parents=True, exist_ok=True)
        picks = rank_candidates(org, args.per_org, args.max_len)
        got = 0
        for wp in picks:
            out = odir / f"{wp}.pdb"
            if out.exists() and out.stat().st_size > 0:
                got += 1; continue
            url = None
            for acc in uniprot_accs(wp):
                url = afdb_pdb_url(acc)
                if url:
                    break
            if not url:
                log.info("  %s: sem modelo AF", wp); continue
            pdb = requests.get(url, timeout=60).text
            out.write_text(pdb)
            got += 1
            log.info("  %s <- %s (%d ATOM)", wp, url.split("/")[-1], pdb.count("\nATOM"))
        log.info("%s: %d/%d candidatas com estrutura", org, got, len(picks))


if __name__ == "__main__":
    main()
