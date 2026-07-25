#!/usr/bin/env python3
"""Roda DeepTMHMM (BioLib cloud) nas 315 candidatas e grava n_tm_helices + sinal.

Uso: python scripts/run_deeptmhmm_candidates.py
Saídas:
  results/03b_deeptmhmm/deeptmhmm_out/   (arquivos brutos: TMRs.gff3, .3line, .md)
  results/03b_deeptmhmm/all_tmhmm.tsv    (protein_id, n_tm_helices, has_signal)
"""
from __future__ import annotations
import os, re, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FASTA = ROOT / "results/03b_deeptmhmm/all_candidates.faa"
WORK = ROOT / "results/03b_deeptmhmm"
OUTDIR = WORK / "deeptmhmm_out"


def parse_gff3(path: Path):
    tm, sig = {}, {}
    ids = []
    for line in path.read_text().splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            m = re.match(r"#\s+(\S+)\s+Length:", line)
            if m:
                ids.append(m.group(1))
                tm.setdefault(m.group(1), 0)
                sig.setdefault(m.group(1), 0)
            continue
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        pid, kind = parts[0], parts[1].strip()
        if kind == "TMhelix":
            tm[pid] = tm.get(pid, 0) + 1
        elif kind == "signal":
            sig[pid] = 1
    return ids, tm, sig


def main():
    import biolib
    OUTDIR.mkdir(parents=True, exist_ok=True)
    # biolib envia arquivos referenciados por caminho relativo existente no cwd
    os.chdir(WORK)
    local_fa = "all_candidates.faa"
    app = biolib.load("DTU/DeepTMHMM")
    print(f"submetendo {local_fa} ...", flush=True)
    job = app.cli(args=f"--fasta {local_fa}")
    out = job.get_stdout().decode("utf-8", "replace")
    sys.stdout.write(out[-600:] + "\n")
    if OUTDIR.exists():
        shutil.rmtree(OUTDIR)
    job.save_files(str(OUTDIR))
    gff = OUTDIR / "TMRs.gff3"
    ids, tm, sig = parse_gff3(gff)
    import pandas as pd
    rows = [{"protein_id": p, "n_tm_helices": tm.get(p, 0), "has_signal": sig.get(p, 0)} for p in ids]
    df = pd.DataFrame(rows).sort_values("n_tm_helices", ascending=False)
    df.to_csv(WORK / "all_tmhmm.tsv", sep="\t", index=False)
    print(f"OK: {len(df)} proteínas; "
          f"{(df.n_tm_helices>1).sum()} com >1 hélice TM; "
          f"{(df.n_tm_helices==0).sum()} com 0 hélices; "
          f"{df.has_signal.sum()} com peptídeo sinal", flush=True)


if __name__ == "__main__":
    main()
