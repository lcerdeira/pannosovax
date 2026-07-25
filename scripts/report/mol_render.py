#!/usr/bin/env python3
"""Renderizador de 'tubo' 3D sombreado para o backbone do construto (qualidade de figura).

Constrói uma malha tubular real (círculo extrudado ao longo do backbone suavizado, com
quadros de transporte paralelo para não torcer) e sombreia por iluminação difusa. Muito
melhor que plotar uma linha grossa. Reutilizado pelas variações da molécula.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from common import ROOT

COL = {"tag": "#9AA0AA", "adjuvant": "#F2A93B", "helper": "#8A6BBF",
       "linker": "#C9CED6", "bcell": "#2A9D8F", "mhc2": "#3B6EA5", "mhc1": "#D1495B"}


def _hex(c):
    c = c.lstrip("#")
    return np.array([int(c[i:i+2], 16) for i in (0, 2, 4)]) / 255.0


def ca_coords(pdb: Path):
    xs = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            xs.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(xs)


def resi_elements():
    import pandas as pd
    cmap = pd.read_csv(ROOT / "results/08_construct/construct_map.tsv", sep="\t")
    n = int(cmap["end"].max()); elem = ["linker"] * (n + 1)
    for _, r in cmap.iterrows():
        for i in range(int(r["start"]), int(r["end"]) + 1):
            elem[i] = r["element"]
    return elem


def load_backbone():
    """CA dos 2 segmentos ESMFold, montados como bi-domínio contínuo."""
    c1 = ca_coords(ROOT / "results/10_structure/h1_esmfold.pdb")
    c2 = ca_coords(ROOT / "results/10_structure/h2_esmfold.pdb")
    c1 = c1 - c1.mean(0); c2 = c2 - c2.mean(0)
    c2 = c2 + np.array([(c1[:, 0].max() - c1[:, 0].min()) * 1.12, 0, 0])
    return np.vstack([c1, c2])


def _spline(P, pts_per_res=6):
    from scipy.interpolate import splprep, splev
    P = np.asarray(P, float)
    # remove duplicados consecutivos p/ o spline não falhar
    keep = np.r_[True, np.any(np.diff(P, axis=0) != 0, axis=1)]
    P = P[keep]
    m = max(len(P) * pts_per_res, 200)
    tck, _ = splprep(P.T, s=len(P) * 0.5, k=3)
    u = np.linspace(0, 1, m)
    S = np.array(splev(u, tck)).T
    res_of = np.clip((u * (len(keep) - 1)).astype(int), 0, len(keep) - 1)
    # mapeia índice do subconjunto de volta ao índice original
    orig_idx = np.where(keep)[0]
    res_of = orig_idx[np.clip(res_of, 0, len(orig_idx) - 1)]
    return S, res_of


def build_tube(coords, elem, radius=2.4, nring=10, adj_boost=1.7):
    S, res_of = _spline(coords)
    T = np.gradient(S, axis=0)
    T /= (np.linalg.norm(T, axis=1, keepdims=True) + 1e-9)
    # transporte paralelo do quadro
    N = np.zeros_like(S)
    ref = np.array([0, 0, 1.0])
    if abs(T[0] @ ref) > 0.9:
        ref = np.array([0, 1.0, 0])
    N[0] = np.cross(T[0], ref); N[0] /= np.linalg.norm(N[0]) + 1e-9
    for i in range(1, len(S)):
        v = N[i-1] - (N[i-1] @ T[i]) * T[i]
        n = np.linalg.norm(v)
        N[i] = v / n if n > 1e-6 else N[i-1]
    B = np.cross(T, N)
    ang = np.linspace(0, 2*np.pi, nring, endpoint=False)
    # raio por ponto (bump no adjuvante)
    rad = np.array([radius * (adj_boost if elem[res_of[i] + 1] == "adjuvant" else 1.0)
                    for i in range(len(S))])
    rings = (S[:, None, :]
             + (rad[:, None, None] * np.cos(ang)[None, :, None]) * N[:, None, :]
             + (rad[:, None, None] * np.sin(ang)[None, :, None]) * B[:, None, :])
    faces, fres = [], []
    for i in range(len(S) - 1):
        for j in range(nring):
            k = (j + 1) % nring
            faces.append([rings[i, j], rings[i, k], rings[i+1, k], rings[i+1, j]])
            fres.append(res_of[i])
    return np.array(faces), np.array(fres)


def shaded_colors(faces, fres, elem, light=(-0.3, 0.5, 0.8), ambient=0.42):
    L = np.array(light, float); L /= np.linalg.norm(L)
    v0, v1, v2 = faces[:, 0], faces[:, 1], faces[:, 2]
    nrm = np.cross(v1 - v0, v2 - v0)
    nrm /= (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-9)
    diff = np.clip(nrm @ L, 0, 1)
    bright = ambient + (1 - ambient) * diff
    base = np.array([_hex(COL.get(elem[r + 1], "#888")) for r in fres])
    rgb = np.clip(base * bright[:, None], 0, 1)
    return rgb


def add_molecule(ax, radius=2.4, nring=10):
    coords = load_backbone(); elem = resi_elements()
    faces, fres = build_tube(coords, elem, radius=radius, nring=nring)
    rgb = shaded_colors(faces, fres, elem)
    # ordena por profundidade (pintor) — média z após a projeção é aproximada por z do centro
    pc = Poly3DCollection(faces, facecolors=rgb, edgecolors="none",
                          linewidths=0, shade=False)
    pc.set_zsort("average")
    ax.add_collection3d(pc)
    lo, hi = coords.min(0), coords.max(0); pad = radius * 1.6
    ax.set_xlim(lo[0]-pad, hi[0]+pad); ax.set_ylim(lo[1]-pad, hi[1]+pad); ax.set_zlim(lo[2]-pad, hi[2]+pad)
    ax.set_box_aspect(hi - lo + 2*pad)
    ax.set_axis_off()
    return coords, elem
