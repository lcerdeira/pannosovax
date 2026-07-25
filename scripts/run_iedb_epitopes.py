#!/usr/bin/env python3
"""Predição de epitopos via IEDB (APIs públicas nextgen) nas candidatas refinadas.

Roda MHC-I, MHC-II e B-cell linear para os 3 organismos, com checkpoint por
proteína (retomável). Escreve na disposição que o pipeline espera:
    results/05_epitopes/{org}_{klass}_raw.tsv

Conservação (k-mer exato em >=95% dos isolados) SÓ é aplicada se existirem os
alinhamentos por gene em results/02_pangenome/{org}_gene_alignments/. Se ausentes,
a coluna `conservation` fica NaN e NÃO se descarta nada silenciosamente — o filtro
de conservação (e a etapa 07) fica pendente até os alinhamentos serem gerados.

Uso:
    python scripts/run_iedb_epitopes.py                 # tudo
    python scripts/run_iedb_epitopes.py --classes mhc2  # só uma classe
    python scripts/run_iedb_epitopes.py --organisms spneu
"""
from __future__ import annotations
import argparse, io, sys, time
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, load_config, read_alleles, write_table, ROOT

log = get_logger("iedb")

URL = {
    "mhc1": "http://tools-cluster-interface.iedb.org/tools_api/mhci/",
    "mhc2": "http://tools-cluster-interface.iedb.org/tools_api/mhcii/",
    "bcell": "http://tools-cluster-interface.iedb.org/tools_api/bcell/",
}


def post(url: str, data: dict, retries: int = 4) -> str:
    for attempt in range(retries):
        try:
            r = requests.post(url, data=data, timeout=300)
            r.raise_for_status()
            t = r.text
            if t.lower().startswith("error") or "<ul class=\"errorlist" in t[:200]:
                raise RuntimeError(t[:200].replace("\n", " "))
            return t
        except Exception as exc:  # noqa: BLE001
            wait = 5 * (2 ** attempt)
            log.warning("IEDB falhou (%s) — retry em %ds", str(exc)[:150], wait)
            time.sleep(wait)
    raise RuntimeError(f"IEDB indisponível após {retries} tentativas: {url}")


def predict_mhc1(pid, seq, alleles, lengths):
    pa, pl = [], []
    for a in alleles:
        for L in lengths:
            pa.append(a); pl.append(str(L))
    txt = post(URL["mhc1"], {"method": "netmhcpan_ba", "sequence_text": f">{pid}\n{seq}\n",
                             "allele": ",".join(pa), "length": ",".join(pl)})
    df = pd.read_csv(io.StringIO(txt), sep="\t")
    rank = next((c for c in df.columns if "percentile" in c.lower() or c.lower() == "rank"), None)
    df = df.rename(columns={rank: "percentile_rank"})
    pep = next(c for c in df.columns if c.lower() == "peptide")
    return df.rename(columns={pep: "peptide"})


def predict_mhc2(pid, seq, alleles):
    txt = post(URL["mhc2"], {"method": "netmhciipan_ba", "sequence_text": f">{pid}\n{seq}\n",
                             "allele": ",".join(alleles)})
    df = pd.read_csv(io.StringIO(txt), sep="\t")
    rank = next((c for c in df.columns if c.lower() in ("rank", "percentile_rank")
                 or "percentile" in c.lower()), None)
    df = df.rename(columns={rank: "percentile_rank"})
    pep = next(c for c in df.columns if c.lower() == "peptide")
    return df.rename(columns={pep: "peptide"})


def predict_bcell(pid, seq, min_len, max_len, thr):
    txt = post(URL["bcell"], {"method": "Bepipred-2.0", "sequence_text": seq})
    df = pd.read_csv(io.StringIO(txt), sep="\t")
    # colunas: Position, Residue, Score, Assignment ('E' = epitopo)
    assign = df["Assignment"].tolist(); res = df["Residue"].tolist(); sc = df["Score"].tolist()
    runs = []
    i = 0
    n = len(assign)
    while i < n:
        if assign[i] == "E":
            j = i
            while j < n and assign[j] == "E":
                j += 1
            peptide = "".join(res[i:j]); score = sum(sc[i:j]) / (j - i)
            if min_len <= len(peptide) <= max_len and score >= thr:
                runs.append({"peptide": peptide, "start": i, "end": j - 1,
                             "length": len(peptide), "score": round(score, 4)})
            i = j
        else:
            i += 1
    return pd.DataFrame(runs)


def load_seqs(org):
    cand = pd.read_csv(ROOT / f"results/03_surfaceome/{org}_candidates.tsv", sep="\t")
    keep = set(cand["protein_id"])
    faa = ROOT / f"results/02_pangenome/{org}_core_proteins.faa"
    return {r.id: str(r.seq) for r in SeqIO.parse(faa, "fasta") if r.id in keep}


def alignments_dir(org):
    d = ROOT / f"results/02_pangenome/{org}_gene_alignments"
    return d if d.exists() else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", nargs="+", default=["mhc2", "mhc1", "bcell"],
                    choices=["mhc1", "mhc2", "bcell"])
    ap.add_argument("--organisms", nargs="+", default=["kpsc", "abau", "spneu"],
                    choices=["kpsc", "abau", "spneu"])
    args = ap.parse_args()

    cfg = load_config()
    ep = cfg["epitopes"]
    al1 = read_alleles(ep["mhc1"]["alleles_file"])["allele"].tolist()
    al2 = read_alleles(ep["mhc2"]["alleles_file"])["allele"].tolist()
    cache = ROOT / "results/05_epitopes/cache"
    cache.mkdir(parents=True, exist_ok=True)

    for org in args.organisms:
        seqs = load_seqs(org)
        for klass in args.classes:
            cdir = cache / f"{org}_{klass}"; cdir.mkdir(exist_ok=True)
            frames = []
            t0 = time.time()
            for i, (pid, seq) in enumerate(seqs.items(), 1):
                cf = cdir / f"{pid}.tsv"
                if cf.exists():
                    d = pd.read_csv(cf, sep="\t") if cf.stat().st_size else pd.DataFrame()
                else:
                    try:
                        if klass == "mhc1":
                            d = predict_mhc1(pid, seq, al1, ep["mhc1"]["lengths"])
                            d = d[d["percentile_rank"] <= ep["mhc1"]["percentile_rank_max"]]
                        elif klass == "mhc2":
                            d = predict_mhc2(pid, seq, al2)
                            d = d[d["percentile_rank"] <= ep["mhc2"]["percentile_rank_max"]]
                        else:
                            d = predict_bcell(pid, seq, ep["bcell"]["min_len"],
                                              ep["bcell"]["max_len"], ep["bcell"]["linear_threshold"])
                    except Exception as exc:  # noqa: BLE001
                        log.error("%s/%s %s: %s", org, klass, pid, str(exc)[:150])
                        continue
                    d = d.copy()
                    d["protein_id"] = pid
                    d.to_csv(cf, sep="\t", index=False)
                if len(d):
                    frames.append(d)
                if i % 20 == 0 or i == len(seqs):
                    rate = (time.time() - t0) / i
                    log.info("%s/%s: %d/%d proteínas (%.1fs/prot, ~%.0fmin restam)",
                             org, klass, i, len(seqs), rate, rate * (len(seqs) - i) / 60)
            out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
            # conservação: só se houver alinhamentos
            adir = alignments_dir(org)
            if adir and len(out):
                homologs = {}
                for fa in adir.glob("*.fasta"):
                    homologs[fa.stem] = [str(r.seq).replace("-", "") for r in SeqIO.parse(fa, "fasta")]
                out["conservation"] = [
                    (sum(1 for h in homologs.get(pid, []) if pep in h) / len(homologs[pid]))
                    if homologs.get(pid) else float("nan")
                    for pep, pid in zip(out["peptide"], out["protein_id"])
                ]
            else:
                if len(out):
                    out["conservation"] = float("nan")
                if not adir:
                    log.warning("%s: sem alinhamentos por gene — conservação NaN "
                                "(filtro de conservação e etapa 07 ficam PENDENTES)", org)
            if len(out):
                out["organism"] = org; out["epitope_class"] = klass
            write_table(out, ROOT / f"results/05_epitopes/{org}_{klass}_raw.tsv", log)
            log.info("%s/%s: %d epitopos brutos gravados", org, klass, len(out))


if __name__ == "__main__":
    main()
