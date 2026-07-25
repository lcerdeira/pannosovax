#!/usr/bin/env python3
"""Figura-molécula (graphical abstract) do PanNosoVax_v1.

Renderiza o backbone 3D do construto (modelo ESMFold, montado por 2 segmentos ≤400 aa
pela limitação da API pública) como um tubo colorido por bloco funcional, com o
adjuvante RS09 em destaque. É ilustrativo/graphical-abstract, não um modelo refinado.

Entrada : results/10_structure/h1_esmfold.pdb, h2_esmfold.pdb
          results/08_construct/construct_map.tsv
Saída   : results/report/figures/F6_molecula.{png,svg}
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
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.patches as mpatches
from common import get_logger, ROOT

log = get_logger("fig_mol")

COL = {"tag": "#9AA0AA", "adjuvant": "#F2A93B", "helper": "#8A6BBF",
       "linker": "#C9CED6", "bcell": "#2A9D8F", "mhc2": "#3B6EA5", "mhc1": "#D1495B"}
BG = "#0E1420"


def ca_coords(pdb: Path):
    xs = []
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            xs.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(xs)


def resi_elements():
    cmap = pd.read_csv(ROOT / "results/08_construct/construct_map.tsv", sep="\t")
    n = int(cmap["end"].max())
    elem = ["linker"] * (n + 1)
    for _, r in cmap.iterrows():
        for i in range(int(r["start"]), int(r["end"]) + 1):
            elem[i] = r["element"]
    return elem  # index 1-based


def catmull(P, k=12):
    """Suaviza a polilinha (Catmull-Rom) para o tubo ficar liso."""
    P = np.asarray(P)
    out = []
    for i in range(len(P) - 1):
        p0 = P[max(i - 1, 0)]; p1 = P[i]; p2 = P[i + 1]; p3 = P[min(i + 2, len(P) - 1)]
        for t in np.linspace(0, 1, k, endpoint=False):
            t2, t3 = t * t, t * t * t
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t
                              + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
                              + (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
    out.append(P[-1])
    return np.array(out)


def main():
    c1 = ca_coords(ROOT / "results/10_structure/h1_esmfold.pdb")
    c2 = ca_coords(ROOT / "results/10_structure/h2_esmfold.pdb")
    # centraliza cada metade e afasta a 2ª ao longo de x -> lê como bi-domínio contínuo
    c1 = c1 - c1.mean(0)
    c2 = c2 - c2.mean(0)
    span = (c1[:, 0].max() - c1[:, 0].min())
    c2 = c2 + np.array([span * 1.15, 0, 0])
    coords = np.vstack([c1, c2])                       # 788 CA
    elem = resi_elements()
    cols = [COL.get(elem[i + 1], "#888") for i in range(len(coords))]

    # suaviza mantendo a cor por resíduo (repete a cor k vezes)
    k = 10
    smooth = catmull(coords, k)
    scol = []
    for i in range(len(coords) - 1):
        scol += [cols[i]] * k
    scol = scol[:len(smooth) - 1]

    segs = np.stack([smooth[:-1], smooth[1:]], axis=1)

    fig = plt.figure(figsize=(12, 7.6), facecolor=BG)
    ax = fig.add_subplot(111, projection="3d", facecolor=BG)
    ax.set_position([0.0, 0.02, 1.0, 0.83])            # ocupa quase toda a figura

    # sombra/halo: tubo grosso escuro por baixo
    halo = Line3DCollection(segs, colors="#05070C", linewidths=15, alpha=0.55,
                            capstyle="round", joinstyle="round", zorder=1)
    ax.add_collection3d(halo)
    # tubo colorido
    tube = Line3DCollection(segs, colors=scol[:len(segs)], linewidths=8.5,
                            capstyle="round", joinstyle="round", zorder=2)
    ax.add_collection3d(tube)

    # destaque do adjuvante: esferas com glow, desenhadas por cima
    adj_idx = [i for i in range(len(coords)) if elem[i + 1] == "adjuvant"]
    cen = None
    if adj_idx:
        ap = coords[adj_idx]
        for s, a in [(900, 0.12), (520, 0.22), (260, 0.55), (110, 1.0)]:
            ax.scatter(ap[:, 0], ap[:, 1], ap[:, 2], s=s, c="#FFD98A",
                       alpha=a, edgecolors="none", depthshade=False, zorder=5)
        cen = ap.mean(0)

    # N e C
    ax.text(*coords[0], "  N", color="white", fontsize=13, fontweight="bold", zorder=6)
    ax.text(*coords[-1], "C  ", color="white", fontsize=13, fontweight="bold", ha="right", zorder=6)

    # estética — limites justos por eixo (preenche o quadro), proporção real
    ax.set_axis_off()
    lo = coords.min(0); hi = coords.max(0); pad = 6
    ax.set_xlim(lo[0] - pad, hi[0] + pad)
    ax.set_ylim(lo[1] - pad, hi[1] + pad)
    ax.set_zlim(lo[2] - pad, hi[2] + pad)
    ax.set_box_aspect(hi - lo + 2 * pad)
    ax.view_init(elev=16, azim=-62)

    fig.text(0.5, 0.965, "PanNosoVax", ha="center", color="white", fontsize=28, fontweight="bold")
    fig.text(0.5, 0.918, "imunógeno quimérico multi-epitopo pan-nosocomial · 788 aa",
             ha="center", color="#AEB6C2", fontsize=12.5)
    # rótulo do adjuvante num canto livre, com seta apontando para a molécula
    if cen is not None:
        from mpl_toolkits.mplot3d import proj3d
        fig.canvas.draw()
        x2, y2, _ = proj3d.proj_transform(cen[0], cen[1], cen[2], ax.get_proj())
        disp = ax.transData.transform((x2, y2))
        fx, fy = fig.transFigure.inverted().transform(disp)
        fig.text(0.185, 0.75, "Adjuvante RS09", ha="center", color="#F2A93B",
                 fontsize=15, fontweight="bold")
        fig.text(0.185, 0.715, "agonista sintético de TLR4", ha="center", color="#E4B15E", fontsize=10.5)
        fig.patches.append(mpatches.FancyArrowPatch(
            (0.235, 0.735), (fx, fy), transform=fig.transFigure,
            arrowstyle="-|>", mutation_scale=18, color="#F2A93B", lw=2))

    order = [("adjuvant", "Adjuvante RS09"), ("helper", "PADRE (T-helper)"),
             ("bcell", "Epitopos B (×12)"), ("mhc2", "Epitopos MHC-II (×15)"),
             ("mhc1", "Epitopos MHC-I (×17)"), ("linker", "Linkers / His-tag")]
    handles = [mpatches.Patch(color=COL[e], label=l) for e, l in order]
    leg = ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
                    ncol=3, frameon=False, fontsize=10, labelcolor="#DCE1E8",
                    handlelength=1.1, columnspacing=1.6)

    for ext in ("png", "svg"):
        fig.savefig(ROOT / f"results/report/figures/F6_molecula.{ext}",
                    dpi=300, facecolor=BG, bbox_inches="tight")
    plt.close(fig)
    log.info("F6_molecula.{png,svg} salvo (%d CA)", len(coords))


if __name__ == "__main__":
    main()
