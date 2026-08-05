#!/usr/bin/env python3
"""Figura 1 do Application Note (Paper B).

Painel (a): a interface durante a execução — etapas, estado por etapa, controles.
Painel (b): o grafo de regras do Snakemake que a interface dirige.

A figura existe para sustentar a tese do paper numa imagem só: a camada gráfica
**dirige** um workflow reproduzível, não o substitui.

Entradas (geradas antes por build_appnote_figure.sh):
    /tmp/figb/ui_exec.png   captura da UI (Chrome headless)
    /tmp/figb/dag.png       rulegraph renderizado (graphviz)

Saída: results/report/figures/appnote_fig1.{png,pdf}
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from common import get_logger, ROOT

log = get_logger("fig_appnote")
INK, SUB, LINE = "#16202E", "#5C6A7A", "#DFE5EC"
UI = Path("/tmp/figb/ui_exec.png")
DAG = Path("/tmp/figb/dag.png")


def main():
    for p in (UI, DAG):
        if not p.exists():
            raise SystemExit(f"faltando {p} — rode scripts/report/build_appnote_figure.sh")

    ui, dag = mpimg.imread(UI), mpimg.imread(DAG)
    # Alturas proporcionais ao aspecto real de cada imagem: com a mesma largura, a
    # altura ocupada é largura/aspecto. Sem isso sobra um vão grande entre os painéis.
    a_ui = ui.shape[1] / ui.shape[0]
    a_dag = dag.shape[1] / dag.shape[0]
    h_ui, h_dag = 1 / a_ui, 1 / a_dag

    W = 7.2
    fig = plt.figure(figsize=(W, W * (h_ui + h_dag) + 0.55), facecolor="white")
    gs = fig.add_gridspec(2, 1, height_ratios=[h_ui, h_dag], hspace=0.10,
                          left=0.02, right=0.98, top=0.965, bottom=0.015)

    for i, (img, tag) in enumerate(((ui, "a"), (dag, "b"))):
        ax = fig.add_subplot(gs[i])
        ax.imshow(img)
        ax.axis("off")
        # rótulo do painel fora da imagem, à esquerda — não colide com o conteúdo
        ax.text(-0.012, 1.0, tag, transform=ax.transAxes, ha="right", va="top",
                fontsize=13, fontweight="bold", color=INK)

    out = ROOT / "results/report/figures"
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"appnote_fig1.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    log.info("appnote_fig1.{png,pdf} salvo em %s", out)


if __name__ == "__main__":
    main()
