#!/usr/bin/env python3
"""Variações da molécula PanNosoVax a partir do tubo 3D sombreado (mol_render).

Gera:
  F6_molecula.{png,svg}          — versão escura (substitui a antiga)
  F6_molecula_claro.{png,svg}    — fundo claro, para figura de artigo
  F8_pan_nosocomial.{png,svg}    — molécula ao centro + 3 patógenos + setas
  F6_molecula_360.gif            — rotação 360°

Uso: python scripts/report/molecule_variants.py [--only dark light pan gif]
"""
from __future__ import annotations
import argparse, io, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import proj3d
from common import get_logger, ROOT
from report.mol_render import add_molecule, COL

log = get_logger("mol_var")
FIGDIR = ROOT / "results/report/figures"
ORG = {"kpsc": "#3B6EA5", "abau": "#D1495B", "spneu": "#2A9D8F"}
LEG = [("adjuvant", "Adjuvante RS09"), ("helper", "PADRE"), ("bcell", "Epitopos B (×12)"),
       ("mhc2", "Epitopos MHC-II (×15)"), ("mhc1", "Epitopos MHC-I (×17)"),
       ("linker", "Linkers / His-tag")]


def _adjuvant_center(coords, elem):
    idx = [i for i in range(len(coords)) if elem[i + 1] == "adjuvant"]
    return coords[idx].mean(0) if idx else None


def _glow(ax, pt, color="#FFD98A"):
    for s, a in [(900, 0.12), (500, 0.22), (240, 0.55), (110, 1.0)]:
        ax.scatter(*pt, s=s, c=color, alpha=a, edgecolors="none", depthshade=False, zorder=6)


def render(dark=True):
    ink = "#EAF0F7" if dark else "#1A2230"
    sub = "#AEB6C2" if dark else "#5B6675"
    bg = "#0E1420" if dark else "#FFFFFF"
    fig = plt.figure(figsize=(12, 7.6), facecolor=bg)
    ax = fig.add_subplot(111, projection="3d", facecolor=bg)
    ax.set_position([0.0, 0.02, 1.0, 0.82])
    coords, elem = add_molecule(ax, radius=2.5, nring=12)
    ax.view_init(elev=15, azim=-62)
    cen = _adjuvant_center(coords, elem)
    if cen is not None:
        _glow(ax, cen)
    ax.text(*coords[0], "  N", color=ink, fontsize=13, fontweight="bold", zorder=7)
    ax.text(*coords[-1], "C ", color=ink, fontsize=13, fontweight="bold", ha="right", zorder=7)

    fig.text(0.5, 0.965, "PanNosoVax", ha="center", color=ink, fontsize=28, fontweight="bold")
    fig.text(0.5, 0.918, "imunógeno quimérico multi-epitopo pan-nosocomial · 788 aa",
             ha="center", color=sub, fontsize=12.5)
    if cen is not None:
        fig.canvas.draw()
        x2, y2, _ = proj3d.proj_transform(*cen, ax.get_proj())
        fx, fy = fig.transFigure.inverted().transform(ax.transData.transform((x2, y2)))
        fig.text(0.185, 0.75, "Adjuvante RS09", ha="center", color="#F2A93B", fontsize=15, fontweight="bold")
        fig.text(0.185, 0.715, "agonista sintético de TLR4", ha="center", color="#E0AC5E", fontsize=10.5)
        fig.patches.append(mpatches.FancyArrowPatch((0.235, 0.735), (fx, fy),
                           transform=fig.transFigure, arrowstyle="-|>", mutation_scale=18,
                           color="#F2A93B", lw=2))
    handles = [mpatches.Patch(color=COL[e], label=l) for e, l in LEG]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=3,
              frameon=False, fontsize=10, labelcolor=ink, handlelength=1.1, columnspacing=1.6)
    name = "F6_molecula" if dark else "F6_molecula_claro"
    for ext in ("png", "svg"):
        fig.savefig(FIGDIR / f"{name}.{ext}", dpi=300, facecolor=bg, bbox_inches="tight")
    plt.close(fig)
    log.info("%s salvo", name)


def _bacterium(ax, cx, cy, color, kind, label, sublabel):
    """Ícone estilizado do patógeno em coords de eixo (0-1)."""
    if kind == "cocci":       # diplococo (S. pneumoniae)
        for dx in (-0.018, 0.018):
            ax.add_patch(mpatches.Circle((cx + dx, cy), 0.026, facecolor=color,
                         edgecolor="white", lw=1.5, transform=ax.transAxes, zorder=4))
    elif kind == "coccobacillus":  # A. baumannii (curto)
        ax.add_patch(mpatches.FancyBboxPatch((cx-0.045, cy-0.026), 0.09, 0.052,
                     boxstyle="round,pad=0,rounding_size=0.026", facecolor=color,
                     edgecolor="white", lw=1.5, transform=ax.transAxes, zorder=4))
    else:                     # bacilo (KpSC)
        ax.add_patch(mpatches.FancyBboxPatch((cx-0.06, cy-0.024), 0.12, 0.048,
                     boxstyle="round,pad=0,rounding_size=0.024", facecolor=color,
                     edgecolor="white", lw=1.5, transform=ax.transAxes, zorder=4))
    ax.text(cx, cy - 0.06, label, ha="center", va="top", color=color, fontsize=15,
            fontweight="bold", fontstyle="italic", transform=ax.transAxes, zorder=5)
    ax.text(cx, cy - 0.093, sublabel, ha="center", va="top", color="#AEB6C2", fontsize=9.5,
            transform=ax.transAxes, zorder=5)


def render_pan():
    bg = "#0B0F18"; ink = "#EAF0F7"
    fig = plt.figure(figsize=(12, 11), facecolor=bg)
    ax3 = fig.add_axes([0.14, 0.40, 0.72, 0.46], projection="3d", facecolor="none")
    coords, elem = add_molecule(ax3, radius=2.7, nring=12)
    ax3.view_init(elev=14, azim=-62)
    cen = _adjuvant_center(coords, elem)
    if cen is not None:
        _glow(ax3, cen)
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.text(0.5, 0.955, "Um imunógeno, três patógenos", ha="center", color=ink,
             fontsize=27, fontweight="bold")
    fig.text(0.5, 0.912, "PanNosoVax cobre os três principais agentes de pneumonia nosocomial "
             "resistente com um único construto", ha="center", color="#AEB6C2", fontsize=12.5)
    # cápsula central rotulando a molécula
    fig.text(0.5, 0.375, "44 epitopos · adjuvante RS09 · 788 aa", ha="center",
             color="#7F9CC9", fontsize=11, fontweight="bold")
    pts = [(0.17, 0.20, "kpsc", "bacillus", "K. pneumoniae", "species complex · KpSC"),
           (0.5, 0.20, "spneu", "cocci", "S. pneumoniae", "independente de sorotipo"),
           (0.83, 0.20, "abau", "coccobacillus", "A. baumannii", "MDR / XDR")]
    for cx, cy, org, kind, lab, sub in pts:
        ax.annotate("", xy=(cx, cy + 0.052), xytext=(0.5, 0.44),
                    arrowprops=dict(arrowstyle="-|>", color=ORG[org], lw=2.6,
                                    shrinkA=40, shrinkB=6, alpha=0.92,
                                    connectionstyle="arc3,rad=0.05"), zorder=3)
        _bacterium(ax, cx, cy, ORG[org], kind, lab, sub)
    for ext in ("png", "svg"):
        fig.savefig(FIGDIR / f"F8_pan_nosocomial.{ext}", dpi=300, facecolor=bg, bbox_inches="tight")
    plt.close(fig)
    log.info("F8_pan_nosocomial salvo")


def render_gif(frames=36):
    import imageio.v2 as imageio
    bg = "#0E1420"
    fig = plt.figure(figsize=(7, 7), facecolor=bg)
    ax = fig.add_subplot(111, projection="3d", facecolor=bg)
    ax.set_position([0, 0, 1, 1])
    coords, elem = add_molecule(ax, radius=2.6, nring=10)
    cen = _adjuvant_center(coords, elem)
    if cen is not None:
        _glow(ax, cen)
    imgs = []
    for az in np.linspace(0, 360, frames, endpoint=False):
        ax.view_init(elev=12, azim=az)
        buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=85, facecolor=bg)
        buf.seek(0); imgs.append(imageio.imread(buf))
    plt.close(fig)
    imageio.mimsave(FIGDIR / "F6_molecula_360.gif", imgs, duration=0.09, loop=0)
    log.info("F6_molecula_360.gif salvo (%d quadros)", frames)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", choices=["dark", "light", "pan", "gif"])
    args = ap.parse_args()
    todo = args.only or ["dark", "light", "pan", "gif"]
    if "dark" in todo: render(dark=True)
    if "light" in todo: render(dark=False)
    if "pan" in todo: render_pan()
    if "gif" in todo: render_gif()


if __name__ == "__main__":
    main()
