#!/usr/bin/env python3
"""Figuras de publicação do PanNosoVax — versão polida, paleta segura p/ daltonismo.

Regra herdada do make_figures.py: figura sem dado não é desenhada. Cada figura lê
TSVs reais do pipeline; nada é preenchido com valor plausível.

Paleta por organismo validada (six-checks do skill dataviz, CVD-safe, contraste ok):
    kpsc  #3B6EA5 (azul)   abau  #D1495B (vermelho)   spneu #2A9D8F (verde-azulado)

Figuras:
    F1  funil de atrição de epitopos (predito -> conservado -> seguro -> selecionado)
    F2  triagem de segurança por camada (a novidade: comensal domina em spneu)
    F3  cobertura populacional HLA (mundo vs Brasil), MHC-I e MHC-II
    F4  arquitetura do construto PanNosoVax_v1
    F5  crossmatch estrutural (TM-score) entre organismos — dobra ABC-SBP pan-patógeno

Uso: python scripts/report/paper_figures.py [--only F1 F3]
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager  # noqa: F401
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from common import get_logger, ROOT

log = get_logger("paper_fig")

ORG = {"kpsc": "#3B6EA5", "abau": "#D1495B", "spneu": "#2A9D8F"}
ORG_LAB = {"kpsc": "KpSC", "abau": "A. baumannii", "spneu": "S. pneumoniae"}
INK, MUTED, GRID = "#22252A", "#6C7079", "#E4E7EB"
FIGDIR = ROOT / "results/report/figures"

plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300,
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "axes.labelcolor": INK, "axes.labelsize": 10,
    "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "legend.frameon": False, "legend.fontsize": 9,
})
CLASSES = ["mhc1", "mhc2", "bcell"]
CLABEL = {"mhc1": "MHC-I (CD8)", "mhc2": "MHC-II (CD4)", "bcell": "B-cell"}


def save(fig, name):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIGDIR / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info("  %s.{png,pdf,svg}", name)


def _uniq(path):
    p = ROOT / path
    if not p.exists():
        return 0
    try:
        d = pd.read_csv(p, sep="\t"); return d["peptide"].nunique() if len(d) else 0
    except Exception:
        return 0


# ── F1 · funil de atrição ────────────────────────────────────────────────────
def f1_funnel():
    stages = ["Preditos", "Conservados\n(≥95%)", "Seguros\n(4 camadas)", "Selecionados\n(cobertura)"]
    sel = {}
    for kl in ("mhc1", "mhc2"):
        p = ROOT / f"results/07_coverage/selected_{kl}.tsv"
        if p.exists():
            for org, g in pd.read_csv(p, sep="\t").groupby("organism"):
                sel[org] = sel.get(org, 0) + len(g)
    data = {}
    for org in ORG:
        pred = sum(_uniq(f"results/05_epitopes/{org}_{k}_raw.tsv") for k in CLASSES)
        cons = sum(_uniq(f"results/05_epitopes/{org}_{k}_conserved.tsv") for k in CLASSES)
        safe = sum(_uniq(f"results/06_safety/{org}_{k}_safe.tsv") for k in CLASSES)
        if pred:
            data[org] = [pred, cons, safe, sel.get(org, 0)]
    if not data:
        return False
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    x = np.arange(len(stages))
    for org, vals in data.items():
        ax.plot(x, vals, "-o", lw=2.4, ms=9, color=ORG[org], label=ORG_LAB[org],
                markeredgecolor="white", markeredgewidth=1.2, zorder=3)
        for xi, v in zip(x, vals):
            dy = 12 if org == "kpsc" else (-16 if org == "spneu" else 12)
            ax.annotate(f"{v:,}".replace(",", "."), (xi, v), textcoords="offset points",
                        xytext=(0, dy), ha="center", fontsize=8.5, color=ORG[org], fontweight="bold")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(stages, fontsize=9.5)
    ax.set_ylabel("Epitopos únicos (escala log)")
    ax.set_title("Atrição de epitopos ao longo do pipeline")
    ax.set_xlim(-0.35, len(stages) - 0.65)
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper right")
    fig.tight_layout(); save(fig, "F1_funil_epitopos"); return True


# ── F2 · triagem de segurança ────────────────────────────────────────────────
def f2_safety():
    # atribuição por prioridade: comensal > mimetismo > self; senão sobrevive
    cats = ["Sobrevive", "Falha: comensal", "Falha: mimetismo 7-mer", "Falha: self humano"]
    ccol = {"Sobrevive": "#2A9D8F", "Falha: comensal": "#D1495B",
            "Falha: mimetismo 7-mer": "#E8973A", "Falha: self humano": "#8A6BBF"}
    counts = {org: dict.fromkeys(cats, 0) for org in ORG}
    for org in ORG:
        seen = {}
        for kl in CLASSES:
            p = ROOT / f"results/06_safety/{org}_{kl}_screened.tsv"
            if not p.exists():
                continue
            d = pd.read_csv(p, sep="\t")
            if not len(d):
                continue
            for _, r in d.drop_duplicates("peptide").iterrows():
                pep = r["peptide"]
                if pep in seen:
                    continue
                if not bool(r.get("pass_commensal", True)):
                    cat = "Falha: comensal"
                elif not bool(r.get("pass_mimicry", True)):
                    cat = "Falha: mimetismo 7-mer"
                elif not bool(r.get("pass_human", True)):
                    cat = "Falha: self humano"
                else:
                    cat = "Sobrevive"
                seen[pep] = cat
        for c in seen.values():
            counts[org][c] += 1
    if not any(sum(v.values()) for v in counts.values()):
        return False
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    orgs = list(ORG)
    y = np.arange(len(orgs))
    for i, org in enumerate(orgs):
        tot = sum(counts[org].values()) or 1
        left = 0
        for c in cats:
            frac = 100 * counts[org][c] / tot
            ax.barh(i, frac, left=left, color=ccol[c], edgecolor="white", linewidth=1.4,
                    height=0.62, zorder=3)
            if frac >= 6:
                ax.text(left + frac / 2, i, f"{frac:.0f}%", ha="center", va="center",
                        color="white", fontsize=8.5, fontweight="bold")
            left += frac
        ax.text(101, i, f"n={tot:,}".replace(",", "."), va="center", ha="left",
                fontsize=8.5, color=MUTED)
    ax.set_yticks(y); ax.set_yticklabels([ORG_LAB[o] for o in orgs], fontstyle="italic")
    ax.set_xlim(0, 108); ax.set_xlabel("Epitopos conservados (%)")
    ax.set_title("Triagem negativa de segurança por camada")
    ax.grid(axis="y", visible=False); ax.invert_yaxis()
    handles = [mpatches.Patch(color=ccol[c], label=c) for c in cats]
    ax.legend(handles=handles, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.55))
    fig.tight_layout(); save(fig, "F2_triagem_seguranca"); return True


# ── F3 · cobertura populacional ──────────────────────────────────────────────
def f3_coverage():
    # cobertura fenotípica FINAL por organismo (mundo vs Brasil) + nº de epitopos.
    # Curva de acumulação é plana (satura no 1º epitopo), então mostramos o alvo atingido.
    frames = {}
    for kl in ("mhc1", "mhc2"):
        p = ROOT / f"results/07_coverage/selected_{kl}.tsv"
        if p.exists():
            frames[kl] = pd.read_csv(p, sep="\t")
    if not frames:
        return False
    orgs = list(ORG)
    fig, axes = plt.subplots(1, len(frames), figsize=(5.4 * len(frames), 4.4), squeeze=False, sharey=True)
    for ax, (kl, df) in zip(axes[0], frames.items()):
        x = np.arange(len(orgs)); bw = 0.36; ticklab = []
        for i, org in enumerate(orgs):
            g = df[df["organism"] == org]
            if not len(g):
                ticklab.append(ORG_LAB[org]); continue
            cw = 100 * g["set_coverage_world"].iloc[0]
            cb = 100 * g["set_coverage_brazil"].iloc[0]
            ticklab.append(f"{ORG_LAB[org]}\n({len(g)} epitopos)")
            ax.bar(x[i] - bw / 2, cw, bw, color=ORG[org], edgecolor="white", zorder=3)
            ax.bar(x[i] + bw / 2, cb, bw, facecolor=ORG[org], edgecolor="white", zorder=3,
                   hatch="////", alpha=0.55)
            ax.text(x[i] - bw / 2, cw + 0.35, f"{cw:.1f}", ha="center", va="bottom", fontsize=8, color=INK)
            ax.text(x[i] + bw / 2, cb + 0.35, f"{cb:.1f}", ha="center", va="bottom", fontsize=8, color=INK)
        ax.axhline(90, ls=":", c=MUTED, lw=1.1, zorder=1)
        ax.set_xticks(x); ax.set_xticklabels(ticklab, fontsize=8.5)
        for t, org in zip(ax.get_xticklabels(), orgs):
            t.set_color(INK)
        ax.set_title(CLABEL[kl]); ax.set_ylim(85, 101); ax.grid(axis="x", visible=False)
    axes[0][0].set_ylabel("Cobertura fenotípica da população (%)")
    axes[0][0].text(-0.45, 90.3, "meta 90%", fontsize=8, color=MUTED, va="bottom")
    handles = [mpatches.Patch(facecolor="#7A7F88", label="Mundo"),
               mpatches.Patch(facecolor="#7A7F88", hatch="////", alpha=0.55, label="Brasil")]
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.04), fontsize=9)
    fig.suptitle("Cobertura HLA do conjunto selecionado — mundo vs Brasil",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout(); save(fig, "F3_cobertura_populacional"); return True


# ── F4 · arquitetura do construto ────────────────────────────────────────────
def f4_construct():
    p = ROOT / "results/08_construct/construct_map.tsv"
    if not p.exists():
        return False
    cmap = pd.read_csv(p, sep="\t")
    col = {"tag": "#9AA0AA", "adjuvant": "#E8973A", "helper": "#8A6BBF",
           "linker": "#DDE1E6", "bcell": "#2A9D8F", "mhc2": "#3B6EA5", "mhc1": "#D1495B"}
    labgroup = {"tag": "His-tag", "adjuvant": "Adjuvante RS09", "helper": "PADRE",
                "bcell": "Epitopos B", "mhc2": "Epitopos MHC-II", "mhc1": "Epitopos MHC-I"}
    total = int(cmap["end"].max())
    fig, ax = plt.subplots(figsize=(11, 3.0))
    for _, r in cmap.iterrows():
        w = r["end"] - r["start"] + 1
        ax.add_patch(mpatches.FancyBboxPatch(
            (r["start"], 0), w, 1, boxstyle="round,pad=0,rounding_size=1.5",
            facecolor=col.get(r["element"], "#888"), edgecolor="white", linewidth=1.0, zorder=3))
    # rótulos dos blocos de epitopos: centrados sobre o span, tier baixo
    for elem in ["bcell", "mhc2", "mhc1"]:
        seg = cmap[cmap["element"] == elem]
        if not len(seg):
            continue
        mid = (seg["start"].min() + seg["end"].max()) / 2
        ax.annotate(f"{labgroup[elem]}\n×{len(seg)}", (mid, 1.12), ha="center", va="bottom",
                    fontsize=9, color=col[elem], fontweight="bold")
    # elementos N-terminais pequenos: rótulos espalhados num tier alto, com linha-guia
    feats = []
    tagN = cmap[cmap["element"] == "tag"].iloc[0]
    feats.append(("His₆ (N)", (tagN["start"] + tagN["end"]) / 2, 4, col["tag"]))
    adj = cmap[cmap["element"] == "adjuvant"]
    if len(adj):
        feats.append(("Adjuvante RS09", (adj["start"].iloc[0] + adj["end"].iloc[0]) / 2, 120, col["adjuvant"]))
    hel = cmap[cmap["element"] == "helper"]
    if len(hel):
        feats.append(("PADRE", (hel["start"].iloc[0] + hel["end"].iloc[0]) / 2, 235, col["helper"]))
    for lab, mid, xlab, c in feats:
        ax.annotate(lab, xy=(mid, 1.02), xytext=(xlab, 1.85), ha="center", va="bottom",
                    fontsize=8.5, color=c, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=c, lw=1, shrinkA=0, shrinkB=2))
    # His6 C-terminal
    tagC = cmap[cmap["element"] == "tag"]
    if len(tagC) > 1:
        ax.annotate("His₆ (C)", (tagC["end"].max(), 1.12), ha="right", va="bottom",
                    fontsize=8, color=col["tag"], fontweight="bold")
    ax.set_xlim(-10, total + 10); ax.set_ylim(-0.95, 2.35)
    ax.axis("off")
    # barra de escala N->C
    ax.annotate("N", (-4, 0.5), ha="right", va="center", fontsize=11, fontweight="bold", color=INK)
    ax.annotate("C", (total + 4, 0.5), ha="left", va="center", fontsize=11, fontweight="bold", color=INK)
    ax.plot([0, total], [-0.55, -0.55], color=MUTED, lw=1)
    for xt in range(0, total + 1, 100):
        ax.plot([xt, xt], [-0.5, -0.6], color=MUTED, lw=1)
        ax.text(xt, -0.75, str(xt), ha="center", va="top", fontsize=7, color=MUTED)
    ax.set_title(f"PanNosoVax_v1 — construto quimérico multi-epitopo ({total} aa)", pad=26)
    fig.tight_layout(); save(fig, "F4_construto"); return True


# ── F5 · crossmatch estrutural ───────────────────────────────────────────────
def f5_structural():
    p = ROOT / "results/04_shared/tm_matrix.tsv"
    if not p.exists():
        return False
    tm = pd.read_csv(p, sep="\t")
    pairs = [("kpsc", "spneu"), ("kpsc", "abau"), ("abau", "spneu")]
    seq = LinearSegmentedColormap.from_list("tm", ["#F3F5F7", "#9AC1D4", "#2A6F8E", "#123B4E"])
    short = lambda p: p.replace("WP_", "")  # encurta rótulo, mantém unicidade
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.0))
    fig.subplots_adjust(wspace=0.75, right=0.9, top=0.82, bottom=0.22)
    for ax, (oa, ob) in zip(axes, pairs):
        sub = tm[(tm.org_a == oa) & (tm.org_b == ob)]
        if not len(sub):
            sub = tm[(tm.org_a == ob) & (tm.org_b == oa)].rename(
                columns={"org_a": "org_b", "org_b": "org_a", "prot_a": "prot_b", "prot_b": "prot_a"})
        ra = list(dict.fromkeys(sub["prot_a"])); rb = list(dict.fromkeys(sub["prot_b"]))
        M = np.zeros((len(ra), len(rb)))
        for _, r in sub.iterrows():
            M[ra.index(r.prot_a), rb.index(r.prot_b)] = r.tm
        im = ax.imshow(M, cmap=seq, vmin=0, vmax=0.9, aspect="auto")
        for i in range(len(ra)):
            for j in range(len(rb)):
                if M[i, j] >= 0.5:
                    ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=7,
                            color="white" if M[i, j] > 0.55 else INK, fontweight="bold")
        ax.set_xticks(range(len(rb))); ax.set_xticklabels([short(p) for p in rb], rotation=90, fontsize=6.5, color=ORG[ob])
        ax.set_yticks(range(len(ra))); ax.set_yticklabels([short(p) for p in ra], fontsize=6.5, color=ORG[oa])
        ax.set_title(f"{ORG_LAB[oa]} × {ORG_LAB[ob]}", fontsize=10, pad=8)
        ax.set_xlabel(ORG_LAB[ob], color=ORG[ob], fontstyle="italic", fontsize=9)
        ax.set_ylabel(ORG_LAB[oa], color=ORG[oa], fontstyle="italic", fontsize=9)
        ax.tick_params(length=0); ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
    cax = fig.add_axes([0.925, 0.22, 0.014, 0.6])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("TM-score (sobreposição de dobra)", fontsize=9)
    cbar.ax.axhline(0.5, color="#D1495B", lw=1.5)
    cbar.ax.text(0.5, 0.5, "0,5", transform=cbar.ax.transAxes, ha="center", va="bottom",
                 fontsize=6.5, color="#D1495B", fontweight="bold")
    fig.suptitle("Convergência estrutural de antígenos de superfície entre patógenos",
                 fontsize=13, fontweight="bold", y=0.96)
    save(fig, "F5_crossmatch_estrutural"); return True


FIGS = {"F1": f1_funnel, "F2": f2_safety, "F3": f3_coverage, "F4": f4_construct, "F5": f5_structural}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", choices=list(FIGS))
    args = ap.parse_args()
    wanted = args.only or list(FIGS)
    ok = 0
    for name in wanted:
        try:
            if FIGS[name]():
                ok += 1
            else:
                log.warning("%s pulada (sem dados)", name)
        except Exception as exc:
            log.warning("%s falhou: %s: %s", name, type(exc).__name__, exc)
    log.info("%d/%d figuras geradas em %s", ok, len(wanted), FIGDIR)


if __name__ == "__main__":
    main()
