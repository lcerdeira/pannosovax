#!/usr/bin/env python3
"""Resumo gráfico (graphical abstract) do PanNosoVax numa figura só.

Junta: molécula 3D (hero) + fluxo do pipeline + cobertura + números-chave.
Reaproveita os helpers de fig_molecule.py. Fundo escuro, formato paisagem.

Saída: results/report/figures/F0_graphical_abstract.{png,svg}
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.patches as mpatches
from common import get_logger, ROOT
from report.mol_render import add_molecule

log = get_logger("graph_abs")
BG = "#0E1420"; PANEL = "#161D2B"; INK = "#EAF0F7"; SUB = "#AEB6C2"
ORG = {"kpsc": "#3B6EA5", "abau": "#D1495B", "spneu": "#2A9D8F"}


def render_molecule(ax):
    coords, elem = add_molecule(ax, radius=2.5, nring=12)
    adj = coords[[i for i in range(len(coords)) if elem[i + 1] == "adjuvant"]]
    if len(adj):
        cen = adj.mean(0)
        for s, a in [(760, 0.13), (420, 0.24), (200, 0.6), (90, 1.0)]:
            ax.scatter(*cen, s=s, c="#FFD98A", alpha=a, edgecolors="none",
                       depthshade=False, zorder=6)
    ax.view_init(elev=16, azim=-62)


def _bacterium(ax, cx, cy, color, kind, label):
    if kind == "cocci":
        for dx in (-0.012, 0.012):
            ax.add_patch(mpatches.Circle((cx + dx, cy), 0.017, facecolor=color,
                         edgecolor="white", lw=1.2, transform=ax.transAxes, zorder=4))
    elif kind == "coccobacillus":
        ax.add_patch(mpatches.FancyBboxPatch((cx-0.03, cy-0.017), 0.06, 0.034,
                     boxstyle="round,pad=0,rounding_size=0.017", facecolor=color,
                     edgecolor="white", lw=1.2, transform=ax.transAxes, zorder=4))
    else:
        ax.add_patch(mpatches.FancyBboxPatch((cx-0.04, cy-0.016), 0.08, 0.032,
                     boxstyle="round,pad=0,rounding_size=0.016", facecolor=color,
                     edgecolor="white", lw=1.2, transform=ax.transAxes, zorder=4))
    ax.text(cx, cy - 0.038, label, ha="center", va="top", color=color, fontsize=9.5,
            fontweight="bold", fontstyle="italic", transform=ax.transAxes, zorder=5)


def main():
    fig = plt.figure(figsize=(15, 8.4), facecolor=BG)
    gs = GridSpec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 1],
                  left=0.02, right=0.98, top=0.86, bottom=0.05, wspace=0.06, hspace=0.28)

    # título
    fig.text(0.5, 0.955, "PanNosoVax", ha="center", color=INK, fontsize=30, fontweight="bold")
    fig.text(0.5, 0.905, "vacina multi-epitopo in silico contra pneumonia nosocomial — "
             "KpSC · A. baumannii · S. pneumoniae", ha="center", color=SUB, fontsize=13)

    # molécula (ocupa a coluna esquerda inteira)
    axm = fig.add_subplot(gs[:, 0], projection="3d", facecolor=BG)
    render_molecule(axm)

    # patógenos-alvo ao redor da molécula (overlay 2D sobre a coluna esquerda)
    axo = fig.add_axes([0, 0, 1, 1]); axo.axis("off"); axo.set_xlim(0, 1); axo.set_ylim(0, 1)
    targets = [(0.085, 0.10, "kpsc", "bacillus", "K. pneumoniae"),
               (0.275, 0.065, "spneu", "cocci", "S. pneumoniae"),
               (0.465, 0.10, "abau", "coccobacillus", "A. baumannii")]
    for cx, cy, org, kind, lab in targets:
        axo.annotate("", xy=(cx, cy + 0.035), xytext=(0.275, 0.36),
                     arrowprops=dict(arrowstyle="-|>", color=ORG[org], lw=1.8,
                                     shrinkA=18, shrinkB=4, alpha=0.85), zorder=2)
        _bacterium(axo, cx, cy, ORG[org], kind, lab)
    axo.text(0.275, 0.40, "um imunógeno → três patógenos", ha="center", color=SUB,
             fontsize=9.5, fontstyle="italic", zorder=3)

    # ── painel superior direito: fluxo do pipeline ───────────────────────────
    axf = fig.add_subplot(gs[0, 1]); axf.set_facecolor(BG); axf.axis("off")
    axf.set_xlim(0, 1); axf.set_ylim(0, 1)
    steps = [("3 patógenos", "nicho: pneumonia\nno hospedeiro vulnerável", "#6C7B8A"),
             ("302 antígenos", "core de superfície,\nnunca capsulares", "#2A9D8F"),
             ("epitopos\nconservados", "≥95% dos isolados", "#3B6EA5"),
             ("4 camadas de\nsegurança", "inclui triagem vs\nmicrobioma comensal", "#D1495B"),
             ("cobertura\nHLA-BR", "ponderada p/ o Brasil", "#E8973A")]
    n = len(steps); y = 0.62
    for i, (t, s, c) in enumerate(steps):
        x = 0.5 / n + i / n
        axf.add_patch(mpatches.FancyBboxPatch((x - 0.44 / n, y - 0.16), 0.88 / n, 0.32,
                      boxstyle="round,pad=0.01,rounding_size=0.02", facecolor=PANEL,
                      edgecolor=c, linewidth=1.8, transform=axf.transData))
        axf.text(x, y + 0.05, t, ha="center", va="center", color=c, fontsize=9.5, fontweight="bold")
        axf.text(x, y - 0.09, s, ha="center", va="center", color=SUB, fontsize=6.6)
        if i < n - 1:
            axf.annotate("", xy=(x + 0.62 / n, y), xytext=(x + 0.44 / n, y),
                         arrowprops=dict(arrowstyle="-|>", color=SUB, lw=1.5))
    axf.text(0.5, 0.97, "Do genoma ao imunógeno", ha="center", va="top",
             color=INK, fontsize=12, fontweight="bold")

    # ── painel inferior direito: cobertura + tiles ───────────────────────────
    axb = fig.add_subplot(gs[1, 1]); axb.set_facecolor(BG); axb.axis("off")
    axb.set_xlim(0, 1); axb.set_ylim(0, 1)
    # mini-barras de cobertura (MHC-I mundo/Brasil por organismo)
    sel = ROOT / "results/07_coverage/selected_mhc1.tsv"
    if sel.exists():
        df = pd.read_csv(sel, sep="\t")
        axb.text(0.02, 0.93, "Cobertura populacional (MHC-I)", color=INK, fontsize=11,
                 fontweight="bold", va="top")
        orgs = ["kpsc", "abau", "spneu"]; bw = 0.1
        for i, org in enumerate(orgs):
            g = df[df["organism"] == org]
            if not len(g):
                continue
            cw = g["set_coverage_world"].iloc[0]; cb = g["set_coverage_brazil"].iloc[0]
            xc = 0.16 + i * 0.30
            for j, (val, hatch, lab) in enumerate([(cw, None, "mundo"), (cb, "///", "BR")]):
                axb.add_patch(mpatches.Rectangle((xc + j * bw - 0.02, 0.42), bw * 0.9, 0.34 * val,
                              facecolor=ORG[org], hatch=hatch, edgecolor=BG,
                              alpha=1 if j == 0 else 0.55))
                axb.text(xc + j * bw + 0.025, 0.42 + 0.34 * val + 0.01, f"{100*val:.0f}",
                         ha="center", color=INK, fontsize=7)
            axb.text(xc + bw, 0.37, {"kpsc": "KpSC", "abau": "A. bau.", "spneu": "S. pneu."}[org],
                     ha="center", color=SUB, fontsize=8, fontstyle="italic")
    # tiles de números-chave
    tiles = [("44", "epitopos"), ("~98%", "cobertura HLA"),
             ("3", "patógenos"), ("788 aa", "CAI 0,84 · 0 Cys")]
    for i, (big, small) in enumerate(tiles):
        x = 0.13 + i * 0.25
        axb.text(x, 0.20, big, ha="center", color="#7FD1C2" if i % 2 else "#8FB3E0",
                 fontsize=17, fontweight="bold")
        axb.text(x, 0.07, small, ha="center", color=SUB, fontsize=8)

    for ext in ("png", "svg"):
        fig.savefig(ROOT / f"results/report/figures/F0_graphical_abstract.{ext}",
                    dpi=300, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    log.info("F0_graphical_abstract.{png,svg} salvo")


if __name__ == "__main__":
    main()
