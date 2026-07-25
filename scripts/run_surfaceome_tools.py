#!/usr/bin/env python3
"""Roda as 3 ferramentas de surfaceome no proteoma core e emite os TSVs que o
`03_surfaceome_filter.py` espera (psortb, signalp, tmhmm).

Substitui o passo que nunca rodou (PSORTb/SignalP estavam ausentes). Fontes:
  localização  -> DeepLocPro (CLI local, Gram-aware)         -> {org}_psortb.tsv
  peptídeo sinal -> SignalP-6 (BioLib, DTU)                  -> {org}_signalp.tsv
  topologia    -> DeepTMHMM (BioLib, DTU)                    -> {org}_tmhmm.tsv

DeepLocPro usa vocabulário próprio; mapeamos para o de PSORTb que o config referencia
(gram_negative_ok=[OuterMembrane, Extracellular], gram_positive_ok=[Cellwall, Extracellular]).
Score PSORTb (0-10, limiar 7.5 no config) = probabilidade da classe predita × 10.

Uso:
    python scripts/run_surfaceome_tools.py --organism kpsc --tool deeploc
    python scripts/run_surfaceome_tools.py --organism spneu --tool signalp
    python scripts/run_surfaceome_tools.py --organism kpsc --tool tmhmm
"""
from __future__ import annotations
import argparse, os, shutil, subprocess, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, load_config, ROOT

log = get_logger("surfaceome")

GRAM = {"kpsc": "negative", "abau": "negative", "spneu": "positive"}
# DeepLocPro -> vocabulário PSORTb (o que o config e o 03_surfaceome_filter usam)
LOC_MAP = {
    "Outer Membrane": "OuterMembrane", "Extracellular": "Extracellular",
    "Cell wall & surface": "Cellwall", "Cytoplasmic Membrane": "CytoplasmicMembrane",
    "Cytoplasmic": "Cytoplasmic", "Periplasmic": "Periplasmic",
}


def core_faa(org: str) -> Path:
    return ROOT / f"results/02_pangenome/{org}_core_proteins.faa"


def out_tsv(org: str, name: str) -> Path:
    p = ROOT / f"results/03_surfaceome/{org}_{name}.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


# ── DeepLocPro (CLI) -> formato PSORTb 'long' ────────────────────────────────
def run_deeploc(org: str, fasta: Path):
    import pandas as pd
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["deeplocpro", "-f", str(fasta), "-o", td,
                        "-g", GRAM[org], "-d", "cpu"], check=True)
        csv = next(Path(td).glob("results_*.csv"))
        df = pd.read_csv(csv)
    prob_cols = ["Cell wall & surface", "Extracellular", "Cytoplasmic",
                 "Cytoplasmic Membrane", "Outer Membrane", "Periplasmic"]
    prob_cols = [c for c in prob_cols if c in df.columns]
    dst = out_tsv(org, "psortb")
    with open(dst, "w") as fh:
        for _, r in df.iterrows():
            loc = LOC_MAP.get(str(r["Localization"]), str(r["Localization"]).replace(" ", ""))
            score = float(r[prob_cols].max()) * 10 if prob_cols else float("nan")
            fh.write(f"SeqID: {r['ACC']}\n  Final Prediction:\n  {loc} {score:.2f}\n\n")
    log.info("%s: DeepLocPro -> %s (%d proteínas)", org, dst.name, len(df))


# ── SignalP-6 (BioLib) -> tabular que parse_signalp espera ───────────────────
def run_signalp(org: str, fasta: Path):
    import biolib, pandas as pd
    work = Path(tempfile.mkdtemp())
    local = work / "in.faa"
    shutil.copy(fasta, local)
    prev = os.getcwd(); os.chdir(work)
    try:
        app = biolib.load("DTU/SignalP-6")
        job = app.cli(args="--fastafile in.faa --organism other --format txt --output_dir output")
        job.get_stdout()
        outdir = work / "sp"
        job.save_files(str(outdir))
        pred = next(outdir.rglob("prediction_results.txt"))
    finally:
        os.chdir(prev)
    rows = []
    for line in pred.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        pid = p[0].split()[0]
        prediction = p[1].strip()           # OTHER/SP/LIPO/TAT/TATLIPO/PILIN
        # colunas que parse_signalp lê: protein_id, prediction, p_other, p_sp, p_lipo, p_tat, cs
        p_other = p[2] if len(p) > 2 else "0"
        p_sp = p[3] if len(p) > 3 else "0"
        p_lipo = p[4] if len(p) > 4 else "0"
        p_tat = p[5] if len(p) > 5 else "0"
        cs = p[-1] if "CS pos" in line else ""
        rows.append([pid, prediction, p_other, p_sp, p_lipo, p_tat, cs])
    dst = out_tsv(org, "signalp")
    pd.DataFrame(rows).to_csv(dst, sep="\t", header=False, index=False)
    log.info("%s: SignalP-6 -> %s (%d proteínas)", org, dst.name, len(rows))
    shutil.rmtree(work, ignore_errors=True)


# ── DeepTMHMM (BioLib) -> gff3 (parse_tmhmm conta TMhelix) ───────────────────
def run_tmhmm(org: str, fasta: Path):
    import biolib
    work = Path(tempfile.mkdtemp())
    shutil.copy(fasta, work / "in.faa")
    prev = os.getcwd(); os.chdir(work)
    try:
        app = biolib.load("DTU/DeepTMHMM")
        job = app.cli(args="--fasta in.faa")
        job.get_stdout()
        outdir = work / "tm"
        job.save_files(str(outdir))
        gff = next(outdir.rglob("*.gff3"))
        shutil.copy(gff, out_tsv(org, "tmhmm"))
    finally:
        os.chdir(prev)
    log.info("%s: DeepTMHMM -> %s", org, out_tsv(org, "tmhmm").name)
    shutil.rmtree(work, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--tool", required=True, choices=["deeploc", "signalp", "tmhmm"])
    ap.add_argument("--fasta", default=None, help="default: proteoma core do organismo")
    args = ap.parse_args()
    load_config()
    fasta = Path(args.fasta) if args.fasta else core_faa(args.organism)
    if not fasta.exists():
        raise SystemExit(f"FASTA ausente: {fasta}")
    {"deeploc": run_deeploc, "signalp": run_signalp, "tmhmm": run_tmhmm}[args.tool](args.organism, fasta)


if __name__ == "__main__":
    main()
