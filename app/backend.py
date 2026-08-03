#!/usr/bin/env python3
"""PanNosoVax Studio — backend local.

Fatia vertical da arquitetura: a UI fala com esta API, que dirige o Snakemake.
Nada de reimplementar pipeline: as regras do workflow/Snakefile continuam sendo
a fonte da verdade (DAG, retomada e provenance vêm de graça).

Endpoints:
    GET  /api/project       -> nome, organismos e config do projeto
    GET  /api/stages        -> etapas, descrição e se estão pendentes (dry-run)
    POST /api/run           -> dispara execução de etapas selecionadas
    GET  /api/status        -> estado do job corrente
    GET  /api/log           -> últimas linhas do log
    POST /api/cancel        -> encerra o job corrente
    GET  /api/results       -> números-chave do run (construto, físico-química)
    GET  /api/figures       -> figuras disponíveis
    GET  /api/figure/{name} -> serve uma figura PNG

Uso (dev):  python app/backend.py    → http://127.0.0.1:8765
"""
from __future__ import annotations

import os
import re
import signal
import subprocess
import sys
import threading
from collections import deque
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

ROOT = Path(__file__).resolve().parent.parent
SNAKEFILE = ROOT / "workflow" / "Snakefile"
PY = sys.executable

# Etapas expostas na v1. Deliberadamente NÃO expomos as que exigem passo manual
# (estrutura 3D, MD, imunossimulação): elas geram "pacote para submissão" e
# entram numa v2, para não prometer ao usuário algo que trava no meio.
STAGES = [
    ("epitopes",  "Predição de epitopos",   "B, MHC-I e MHC-II via IEDB. Lento (horas) — retomável."),
    ("safety",    "Triagem de segurança",   "Self humano, mimetismo 7-mer e microbioma comensal."),
    ("coverage",  "Cobertura populacional", "Conjunto mínimo que maximiza cobertura HLA."),
    ("construct", "Montagem do construto",  "Adjuvante, linkers e blocos de epitopos."),
    ("physchem",  "Físico-química",         "MW, pI, GRAVY, instabilidade, solubilidade."),
    ("expression","Expressão e códon",      "Otimização para E. coli e variante mRNA."),
    ("report",    "Figuras e tabelas",      "Gera as figuras do relatório."),
]

app = FastAPI(title="PanNosoVax Studio")

_job: dict = {"proc": None, "stages": [], "state": "idle", "returncode": None}
_log: deque[str] = deque(maxlen=800)
_lock = threading.Lock()


def _snakemake(*args: str) -> list[str]:
    return [PY, "-m", "snakemake", "-s", str(SNAKEFILE), "-d", str(ROOT), *args]


def _pending() -> set[str]:
    """Etapas que o Snakemake ainda executaria (dry-run)."""
    try:
        out = subprocess.run(_snakemake("-n", "--quiet"), capture_output=True,
                             text=True, timeout=180, cwd=ROOT).stdout
    except Exception:
        return set()
    pend, in_table = set(), False
    for line in out.splitlines():
        if re.match(r"^job\s+count", line):
            in_table = True
            continue
        if in_table:
            m = re.match(r"^(\w[\w_]*)\s+\d+$", line.strip())
            if m:
                pend.add(m.group(1))
    return pend


@app.get("/api/stages")
def stages():
    pend = _pending()
    return {"stages": [{"id": s, "name": n, "help": h, "pending": s in pend} for s, n, h in STAGES]}


class RunReq(BaseModel):
    stages: list[str]
    cores: int = 4


def _pump(proc: subprocess.Popen):
    for line in iter(proc.stdout.readline, ""):
        _log.append(line.rstrip())
    proc.wait()
    with _lock:
        _job["state"] = "done" if proc.returncode == 0 else "error"
        _job["returncode"] = proc.returncode


@app.post("/api/run")
def run(req: RunReq):
    with _lock:
        if _job["state"] == "running":
            return {"ok": False, "error": "Já existe uma execução em andamento."}
        valid = [s for s in req.stages if s in {x[0] for x in STAGES}]
        if not valid:
            return {"ok": False, "error": "Nenhuma etapa válida selecionada."}
        _log.clear()
        proc = subprocess.Popen(
            _snakemake("--cores", str(req.cores), "--rerun-incomplete", *valid),
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, preexec_fn=os.setsid if os.name != "nt" else None)
        _job.update(proc=proc, stages=valid, state="running", returncode=None)
        threading.Thread(target=_pump, args=(proc,), daemon=True).start()
    return {"ok": True, "stages": valid}


@app.get("/api/status")
def status():
    return {"state": _job["state"], "stages": _job["stages"], "returncode": _job["returncode"]}


@app.get("/api/log")
def log(tail: int = 200):
    return {"lines": list(_log)[-tail:]}


@app.post("/api/cancel")
def cancel():
    with _lock:
        proc = _job.get("proc")
        if proc and _job["state"] == "running":
            try:
                if os.name != "nt":
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                else:
                    proc.terminate()
            except Exception as exc:
                return {"ok": False, "error": str(exc)}
            _job["state"] = "cancelled"
            return {"ok": True}
    return {"ok": False, "error": "Nada em execução."}


# ── Projeto: informação de config para a tela inicial ───────────────────────
FIGDIR = ROOT / "results" / "report" / "figures"
CONFIG = ROOT / "config" / "config.yaml"
ORG_LABELS = {"kpsc": "K. pneumoniae (KpSC)", "abau": "A. baumannii",
              "spneu": "S. pneumoniae"}


@app.get("/api/project")
def project():
    import yaml
    info = {"name": "PanNosoVax", "organisms": [], "outdir": "results", "has_results": False}
    try:
        cfg = yaml.safe_load(CONFIG.read_text())
        info["name"] = cfg.get("project", "PanNosoVax")
        info["outdir"] = cfg.get("outdir", "results")
        info["organisms"] = [{"id": o, "label": ORG_LABELS.get(o, o),
                              "n_genomes": cfg["organisms"][o].get("n_genomes")}
                             for o in cfg.get("organisms", {})]
    except Exception as exc:  # noqa: BLE001
        info["error"] = str(exc)
    info["has_results"] = (ROOT / "results" / "08_construct" / "construct.fasta").exists()
    return info


# ── Resultados: resumo e figuras ────────────────────────────────────────────
FIG_CAPTIONS = {
    "F0_graphical_abstract": "Resumo gráfico",
    "F1_funil_epitopos": "Atrição de epitopos",
    "F2_triagem_seguranca": "Triagem de segurança",
    "F3_cobertura_populacional": "Cobertura HLA (mundo vs Brasil)",
    "F4_construto": "Arquitetura do construto",
    "F5_crossmatch_estrutural": "Convergência estrutural",
    "F6_molecula": "Modelo 3D do imunógeno",
    "F8_pan_nosocomial": "Um imunógeno, três patógenos",
}


@app.get("/api/results")
def results_summary():
    """Números-chave do run, se existirem (não inventa nada)."""
    import pandas as pd
    out = {"available": False, "construct": None, "physchem": None}
    fasta = ROOT / "results" / "08_construct" / "construct.fasta"
    if fasta.exists():
        seq = "".join(l for l in fasta.read_text().splitlines() if not l.startswith(">"))
        out["construct"] = {"length": len(seq)}
        out["available"] = True
    pc = ROOT / "results" / "09_physchem" / "construct_properties.tsv"
    if pc.exists():
        try:
            d = pd.read_csv(pc, sep="\t").iloc[0]
            out["physchem"] = {k: (float(d[k]) if k in d else None) for k in
                               ["molecular_weight_kda", "theoretical_pi",
                                "instability_index", "gravy", "n_cysteine"]}
        except Exception:  # noqa: BLE001
            pass
    return out


@app.get("/api/figures")
def figures():
    figs = []
    if FIGDIR.is_dir():
        for p in sorted(FIGDIR.glob("*.png")):
            stem = p.stem
            figs.append({"name": stem,
                         "caption": FIG_CAPTIONS.get(stem, stem.replace("_", " ")),
                         "url": f"/api/figure/{stem}"})
    return {"figures": figs}


@app.get("/api/figure/{name}")
def figure(name: str):
    # sanitiza: só nome-base, sem travessia de diretório
    safe = re.sub(r"[^A-Za-z0-9_.-]", "", name)
    p = FIGDIR / f"{safe}.png"
    if not p.exists():
        return {"error": "figura não encontrada"}
    return FileResponse(p, media_type="image/png")


UI = Path(__file__).resolve().parent / "ui"
if UI.is_dir():
    app.mount("/ui", StaticFiles(directory=UI), name="ui")

    @app.get("/")
    def index():
        return FileResponse(UI / "index.html")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8765, log_level="warning")
