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


# Localizações que contam como "acessível a anticorpo" por tipo de parede.
SURFACE_LOC = {"kpsc": {"OuterMembrane", "Extracellular"},
               "abau": {"OuterMembrane", "Extracellular"},
               "spneu": {"Cellwall", "Extracellular"}}


def surface_subset(org: str, fasta: Path) -> Path:
    """FASTA só com as proteínas que o DeepLocPro localizou na superfície.

    SignalP-6 e DeepTMHMM rodam na nuvem (BioLib), onde cada job tem ~15 min de
    overhead e há tempo-limite. Submeter o proteoma core inteiro (~7.800 seqs) é
    inviável — e desnecessário: a localização já elimina ~94%. Rodamos as duas
    ferramentas apenas no subconjunto de superfície, que é o que entra no construto.
    """
    import re
    from Bio import SeqIO
    psortb = out_tsv(org, "psortb")
    if not psortb.exists():
        log.warning("%s: sem %s — rode --tool deeploc primeiro; usando o core inteiro",
                    org, psortb.name)
        return fasta
    txt = psortb.read_text()
    keep = {p for p, loc, _ in
            re.findall(r"SeqID: (\S+)\n  Final Prediction:\n  (\S+) ([\d.]+)", txt)
            if loc in SURFACE_LOC[org]}
    dst = out_tsv(org, "_surface_subset").with_suffix(".faa")
    recs = [r for r in SeqIO.parse(fasta, "fasta") if r.id in keep]
    SeqIO.write(recs, dst, "fasta")
    log.info("%s: subconjunto de superfície = %d de %d proteínas core",
             org, len(recs), sum(1 for _ in SeqIO.parse(fasta, "fasta")))
    return dst


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


# O proteoma core inteiro (~3900 seqs) estoura o tempo-limite da nuvem BioLib
# ("Job exceeded max run time"). Loteamos e retomamos (cache por bloco). Lote de 500
# equilibra: grande o bastante para amortizar o overhead de cold-start de cada job
# BioLib (~10-15 min fixos), pequeno o bastante para caber no tempo-limite.
BATCH = 500
RETRIES = 3


def _batches(fasta: Path, size: int):
    from Bio import SeqIO
    recs = list(SeqIO.parse(fasta, "fasta"))
    for i in range(0, len(recs), size):
        yield i // size, recs[i:i + size]


def _run_biolib_batches(org, tool, fasta, uri, args_tmpl, out_name):
    """Roda `uri` em lotes; devolve a lista de arquivos de saída (um por lote).
    Cache por lote em results/03_surfaceome/_cache/{org}_{tool}/batch_XX.<ext>."""
    import biolib
    from Bio import SeqIO
    cache = out_tsv(org, "_x").parent / "_cache" / f"{org}_{tool}"
    cache.mkdir(parents=True, exist_ok=True)
    app = biolib.load(uri)
    outputs = []
    total = sum(1 for _ in _batches(fasta, BATCH))
    for bi, chunk in _batches(fasta, BATCH):
        cached = cache / f"batch_{bi:03d}_{out_name}"
        if cached.exists() and cached.stat().st_size > 0:
            outputs.append(cached); continue
        # cabeçalho SÓ com o ID: o SignalP-6 nomeia arquivos de plot pelo header inteiro,
        # e descrições longas de produto estouram o limite de 255 chars do FS ao baixar.
        for r in chunk:
            r.description = ""
        last_err = None
        for attempt in range(RETRIES):
            work = Path(tempfile.mkdtemp())
            SeqIO.write(chunk, work / "in.faa", "fasta")
            prev = os.getcwd(); os.chdir(work)
            try:
                job = app.cli(args=args_tmpl)
                job.get_stdout()
                job.save_files(str(work / "out"))
                produced = next((work / "out").rglob(out_name))
                shutil.copy(produced, cached)
                outputs.append(cached)
                log.info("  %s/%s: lote %d/%d ok (%d seqs)", org, tool, bi + 1, total, len(chunk))
                last_err = None
                break
            except Exception as exc:  # noqa: BLE001
                last_err = exc
                log.warning("  %s/%s: lote %d tentativa %d/%d falhou: %s",
                            org, tool, bi + 1, attempt + 1, RETRIES, str(exc)[:150])
            finally:
                os.chdir(prev); shutil.rmtree(work, ignore_errors=True)
        if last_err is not None:
            raise last_err
    return outputs


# ── SignalP-6 (BioLib, loteado) -> tabular que parse_signalp espera ──────────
def run_signalp(org: str, fasta: Path):
    import pandas as pd
    outs = _run_biolib_batches(
        org, "signalp", fasta, "DTU/SignalP-6",
        "--fastafile in.faa --organism other --format txt --output_dir output",
        "prediction_results.txt")
    rows = []
    for pred in outs:
        for line in pred.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            cs = p[-1] if "CS pos" in line else ""
            rows.append([p[0].split()[0], p[1].strip(),
                         p[2] if len(p) > 2 else "0", p[3] if len(p) > 3 else "0",
                         p[4] if len(p) > 4 else "0", p[5] if len(p) > 5 else "0", cs])
    pd.DataFrame(rows).to_csv(out_tsv(org, "signalp"), sep="\t", header=False, index=False)
    log.info("%s: SignalP-6 -> %s (%d proteínas)", org, out_tsv(org, "signalp").name, len(rows))


# ── DeepTMHMM (BioLib, loteado) -> gff3 concatenado ──────────────────────────
def run_tmhmm(org: str, fasta: Path):
    outs = _run_biolib_batches(org, "tmhmm", fasta, "DTU/DeepTMHMM",
                               "--fasta in.faa", "TMRs.gff3")
    dst = out_tsv(org, "tmhmm")
    with open(dst, "w") as fh:
        fh.write("##gff-version 3\n")
        for gff in outs:
            for line in gff.read_text().splitlines():
                if not line.startswith("#") and line.strip():
                    fh.write(line + "\n")
    log.info("%s: DeepTMHMM -> %s (%d lotes)", org, dst.name, len(outs))


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
    # idempotência: se a saída final já existe e tem conteúdo, não refaz (deeploc é caro)
    name = {"deeploc": "psortb", "signalp": "signalp", "tmhmm": "tmhmm"}[args.tool]
    done = out_tsv(args.organism, name)
    if done.exists() and done.stat().st_size > 0:
        log.info("%s/%s: %s já existe — pulando", args.organism, args.tool, done.name)
        return
    # deeploc roda no core inteiro (é local e barato); signalp/tmhmm só na superfície
    if args.tool in ("signalp", "tmhmm") and args.fasta is None:
        fasta = surface_subset(args.organism, fasta)
    {"deeploc": run_deeploc, "signalp": run_signalp, "tmhmm": run_tmhmm}[args.tool](args.organism, fasta)


if __name__ == "__main__":
    main()
