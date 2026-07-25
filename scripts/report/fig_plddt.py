#!/usr/bin/env python3
"""Perfil de confiança (pLDDT) do construto ESMFold, por resíduo e por bloco.

Usa as 2 metades ESMFold (h1: 1-394, h2: 395-788). O B-factor do PDB ESMFold é o pLDDT.
Escreve results/10_structure/plddt_by_block.tsv e a figura F7.
Nota honesta: pLDDT alto num LINKER flexível seria ruim (deveria ser baixo); o interesse
é o pLDDT dos blocos de EPITOPO — epitopo bem definido é epitopo bem apresentado.

Saída: results/report/figures/F7_plddt.{png,svg}
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
from common import get_logger, ROOT

log = get_logger("fig_plddt")
COL = {"tag": "#9AA0AA", "adjuvant": "#E8973A", "helper": "#8A6BBF",
       "linker": "#C9CED6", "bcell": "#2A9D8F", "mhc2": "#3B6EA5", "mhc1": "#D1495B"}
INK, MUTED, GRID = "#22252A", "#6C7079", "#E4E7EB"


def plddt_series():
    vals = []
    for h in ("h1", "h2"):
        for line in (ROOT / f"results/10_structure/{h}_esmfold.pdb").read_text().splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                vals.append(float(line[60:66]))
    v = np.array(vals)
    return v * 100 if v.max() <= 1.5 else v   # ESMFold grava pLDDT em 0-1


def main():
    plddt = plddt_series()
    cmap = pd.read_csv(ROOT / "results/08_construct/construct_map.tsv", sep="\t")
    # tabela por bloco
    rows = []
    for _, r in cmap.iterrows():
        s, e = int(r["start"]), int(r["end"])
        seg = plddt[s - 1:e]
        if len(seg):
            rows.append({"element": r["element"], "label": r["label"], "start": s, "end": e,
                         "length": e - s + 1, "mean_plddt": round(float(seg.mean()), 1),
                         "min_plddt": round(float(seg.min()), 1)})
    tbl = pd.DataFrame(rows)
    tbl.to_csv(ROOT / "results/10_structure/plddt_by_block.tsv", sep="\t", index=False)

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(12, 6.2), height_ratios=[2, 1.1],
                                  gridspec_kw={"hspace": 0.42})
    # perfil por resíduo, faixas coloridas por bloco
    x = np.arange(1, len(plddt) + 1)
    for _, r in cmap.iterrows():
        ax.axvspan(r["start"] - 0.5, r["end"] + 0.5, color=COL.get(r["element"], "#888"),
                   alpha=0.16, lw=0)
    ax.plot(x, plddt, color=INK, lw=1.0)
    ax.axhline(70, ls="--", color=MUTED, lw=1)
    ax.text(len(plddt), 70, " pLDDT 70", va="center", ha="left", fontsize=8, color=MUTED)
    ax.set_xlim(1, len(plddt)); ax.set_ylim(0, 100)
    ax.set_ylabel("pLDDT (confiança ESMFold)"); ax.set_xlabel("Posição no construto (aa)")
    ax.set_title("Confiança estrutural ao longo do PanNosoVax_v1", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color=GRID, lw=0.8)

    # média por classe de bloco
    order = [("mhc1", "Epitopos MHC-I"), ("mhc2", "Epitopos MHC-II"),
             ("bcell", "Epitopos B"), ("adjuvant", "Adjuvante"),
             ("helper", "PADRE"), ("linker", "Linkers")]
    means = []
    for e, lab in order:
        seg = tbl[tbl["element"] == e]["mean_plddt"]
        means.append((lab, e, seg.mean() if len(seg) else np.nan))
    means = [m for m in means if not np.isnan(m[2])]
    ax2.barh([m[0] for m in means], [m[2] for m in means],
             color=[COL[m[1]] for m in means], edgecolor="white", height=0.68, zorder=3)
    for i, m in enumerate(means):
        ax2.text(m[2] + 1, i, f"{m[2]:.0f}", va="center", fontsize=9, color=INK, fontweight="bold")
    ax2.axvline(70, ls="--", color=MUTED, lw=1)
    ax2.set_xlim(0, 100); ax2.invert_yaxis()
    ax2.set_xlabel("pLDDT médio"); ax2.set_title("Confiança média por tipo de bloco", fontsize=11, fontweight="bold")
    ax2.spines[["top", "right"]].set_visible(False); ax2.grid(axis="x", color=GRID, lw=0.8)

    for ext in ("png", "svg"):
        fig.savefig(ROOT / f"results/report/figures/F7_plddt.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    ep = tbl[tbl["element"].isin(["mhc1", "mhc2", "bcell"])]["mean_plddt"]
    log.info("F7_plddt salvo · pLDDT médio dos epitopos=%.1f · %d/%d blocos de epitopo >=70",
             ep.mean(), int((ep >= 70).sum()), len(ep))


if __name__ == "__main__":
    main()
