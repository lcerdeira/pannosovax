#!/usr/bin/env python3
"""
Relatório — figuras do manuscrito.

Regra que governa este script: **figura sem dado não é desenhada**. Cada figura tem
um TSV de entrada; se o TSV não existe ou está vazio (o pipeline escreve arquivos só
com cabeçalho quando uma ferramenta externa falta), a figura é pulada e o nome dela
entra na lista de faltantes ao final. Nunca há eixo preenchido com valor plausível.

Figuras:
  F1 funil de filtragem      — quantos genes sobram a cada estágio; é a figura que
                               justifica o desenho do pipeline inteiro.
  F2 mapa de blocos          — arquitetura do construto em escala, com pLDDT quando há.
  F3 curva de cobertura      — cobertura acumulada mundo vs Brasil por epitopo
                               adicionado; mostra onde o ganho marginal satura.
  F4 painéis de MD           — RMSD e RMSF das réplicas.
  F5 curso temporal imune    — imunoglobulinas e citocinas do C-ImmSim.

Saída: results/report/figures/{nome}.png (300 dpi) e .pdf

Uso:
    python scripts/report/make_figures.py
    python scripts/report/make_figures.py --only funnel coverage
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from common import get_logger, load_config, outpath  # noqa: E402

log = get_logger("make_figures")

try:                                                   # seaborn é opcional
    import seaborn as sns
    sns.set_theme(style="whitegrid", context="paper")
except Exception:
    plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3,
                         "axes.spines.top": False, "axes.spines.right": False})

PALETTE = {"kpsc": "#1b6ca8", "abau": "#c1452e", "spneu": "#2e8b57"}
skipped: list[str] = []


def save(fig, figdir: Path, name: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(figdir / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("figura %s salva", name)


def load(path: Path) -> pd.DataFrame | None:
    """Lê um TSV; devolve None se ausente ou só com cabeçalho."""
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        log.warning("erro lendo %s: %s", path.name, exc)
        return None
    return df if len(df) else None


def fig_funnel(cfg, figdir) -> bool:
    stages = [
        ("Core", "02_pangenome", "{org}_gene_presence_absence.csv", "is_core"),
        ("Surfaceoma", "03_surfaceome", "{org}_surfaceome_full.tsv", "candidate"),
        ("Seleção", "04_selection", "{org}_dnds.tsv", "pass_selection"),
    ]
    data: dict[str, list[int]] = {}
    for org in cfg["organisms"]:
        counts = []
        for _, sub, tpl, col in stages:
            df = load(outpath(cfg, sub, tpl.format(org=org)))
            counts.append(int(df[col].fillna(False).astype(bool).sum())
                          if df is not None and col in df else 0)
        # Epitopos seguros: última barra do funil.
        n_safe = 0
        for klass in ("mhc1", "mhc2", "bcell"):
            d = load(outpath(cfg, "06_safety", f"{org}_{klass}_safe.tsv"))
            n_safe += len(d) if d is not None else 0
        counts.append(n_safe)
        if any(counts):
            data[org] = counts
    if not data:
        return False

    labels = [s[0] for s in stages] + ["Epitopos seguros"]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    width = 0.8 / len(data)
    for i, (org, counts) in enumerate(data.items()):
        xs = [x + i * width - 0.4 + width / 2 for x in range(len(labels))]
        ax.bar(xs, counts, width, label=org, color=PALETTE.get(org))
        for x, c in zip(xs, counts):
            if c:
                ax.text(x, c, f"{c:,}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yscale("log")
    ax.set_ylabel("Nº de elementos (escala log)")
    ax.set_title("Funil de filtragem por organismo")
    ax.legend(frameon=False)
    save(fig, figdir, "F1_funil_filtragem")
    return True


def fig_blocks(cfg, figdir) -> bool:
    cmap = load(outpath(cfg, "08_construct", "construct_map.tsv"))
    if cmap is None:
        return False
    plddt = load(outpath(cfg, "10_structure", "plddt_by_block.tsv"))

    colors = {"tag": "#9e9e9e", "adjuvant": "#e0932c", "helper": "#8e5aa8",
              "linker": "#d9d9d9", "bcell": "#2e8b57", "mhc2": "#1b6ca8",
              "mhc1": "#c1452e", "shared": "#000000"}
    n = 2 if plddt is not None else 1
    fig, axes = plt.subplots(n, 1, figsize=(9, 2.2 * n), sharex=True, squeeze=False)
    ax = axes[0][0]
    for _, r in cmap.iterrows():
        ax.barh(0, r["end"] - r["start"] + 1, left=r["start"],
                color=colors.get(r["element"], "#777"), edgecolor="white", height=0.6)
        if r["end"] - r["start"] > 8:
            ax.text((r["start"] + r["end"]) / 2, 0, str(r["label"]),
                    ha="center", va="center", fontsize=6, color="white", rotation=90)
    ax.set_yticks([])
    ax.set_title("Arquitetura do construto")

    if plddt is not None:
        ax2 = axes[1][0]
        mid = (plddt["start"] + plddt["end"]) / 2
        ax2.bar(mid, plddt["mean_plddt"], width=plddt["length"],
                color=[colors.get(e, "#777") for e in plddt["element"]])
        ax2.axhline(70, ls="--", c="k", lw=0.8)
        ax2.set_ylabel("pLDDT médio")
        ax2.set_ylim(0, 100)
    axes[-1][0].set_xlabel("Posição no construto (aa)")
    save(fig, figdir, "F2_mapa_construto")
    return True


def fig_coverage(cfg, figdir) -> bool:
    frames = []
    for klass in ("mhc1", "mhc2"):
        d = load(outpath(cfg, "07_coverage", f"selected_{klass}.tsv"))
        if d is not None:
            d = d.copy()
            d["class"] = klass
            frames.append(d)
    if not frames:
        return False
    df = pd.concat(frames, ignore_index=True)

    fig, axes = plt.subplots(1, len(frames), figsize=(5 * len(frames), 3.8), squeeze=False)
    for ax, (klass, sub) in zip(axes[0], df.groupby("class")):
        for org, g in sub.groupby("organism"):
            g = g.sort_values("rank")
            ax.plot(g["rank"], 100 * g["set_coverage_world"], "-o", ms=3,
                    color=PALETTE.get(org), label=f"{org} · mundo")
            ax.plot(g["rank"], 100 * g["set_coverage_brazil"], "--s", ms=3,
                    color=PALETTE.get(org), alpha=0.6, label=f"{org} · Brasil")
        ax.axhline(100 * cfg["coverage"]["target_global"], ls=":", c="k", lw=0.8)
        ax.set_xlabel("Epitopos no conjunto")
        ax.set_ylabel("Cobertura fenotípica (%)")
        ax.set_title(klass.upper())
        ax.legend(frameon=False, fontsize=7)
    save(fig, figdir, "F3_cobertura_populacional")
    return True


def fig_md(cfg, figdir) -> bool:
    summ = load(outpath(cfg, "11_md", "md_summary.tsv"))
    md_dir = outpath(cfg, "11_md", "_").parent
    xvgs = sorted(md_dir.glob("rep*/rmsd.xvg"))
    if summ is None and not xvgs:
        return False

    def read_xvg(p: Path):
        xs, ys = [], []
        for line in p.read_text(errors="replace").splitlines():
            if line.startswith(("@", "#", "&")) or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
                except ValueError:
                    continue
        return xs, ys

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    plotted = False
    for p in xvgs:
        xs, ys = read_xvg(p)
        if xs:
            a1.plot(xs, ys, lw=0.8, label=p.parent.name)
            plotted = True
        f = p.parent / "rmsf.xvg"
        if f.exists():
            fx, fy = read_xvg(f)
            if fx:
                a2.plot(fx, fy, lw=0.8, label=p.parent.name)
                plotted = True
    if not plotted and summ is not None:
        a1.bar(summ["replicate"], summ["rmsd_mean_nm"],
               yerr=summ.get("rmsd_sd_nm"), color="#1b6ca8")
        a1.set_xlabel("Réplica")
        a2.bar(summ["replicate"], summ["rmsf_mean_nm"], color="#c1452e")
        a2.set_xlabel("Réplica")
        plotted = True
    if not plotted:
        plt.close(fig)
        return False

    a1.set_ylabel("RMSD do backbone (nm)")
    a1.set_title("Estabilidade global")
    a2.set_ylabel("RMSF por resíduo (nm)")
    a2.set_title("Flexibilidade local")
    for a in (a1, a2):
        a.legend(frameon=False, fontsize=7)
    save(fig, figdir, "F4_dinamica_molecular")
    return True


def fig_immune(cfg, figdir) -> bool:
    df = load(outpath(cfg, "11_immunosim", "immune_response.tsv"))
    if df is None:
        return False
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True)
    for col, c in [("IgM", "#8e5aa8"), ("IgG1", "#1b6ca8"),
                   ("IgG2", "#2e8b57"), ("Ig_total", "#333333")]:
        if col in df and df[col].notna().any():
            a1.plot(df["day"], df[col], label=col, color=c, lw=1.2)
    a1.set_ylabel("Título (arb.)")
    a1.set_title("Imunoglobulinas")
    for col, c in [("TH_cells", "#e0932c"), ("TC_cells", "#c1452e"),
                   ("B_cells", "#2e8b57"), ("IFNg", "#1b6ca8"), ("IL2", "#8e5aa8")]:
        if col in df and df[col].notna().any():
            a2.plot(df["day"], df[col], label=col, color=c, lw=1.2)
    a2.set_ylabel("Células / citocinas (arb.)")
    a2.set_title("Celular e citocinas")
    for a in (a1, a2):
        a.set_xlabel("Dias")
        a.legend(frameon=False, fontsize=7)
    save(fig, figdir, "F5_imunossimulacao")
    return True


FIGURES = {
    "funnel": fig_funnel,
    "blocks": fig_blocks,
    "coverage": fig_coverage,
    "md": fig_md,
    "immune": fig_immune,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", choices=list(FIGURES), default=None)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    figdir = outpath(cfg, "report", "figures", "_").parent
    figdir.mkdir(parents=True, exist_ok=True)

    wanted = args.only or list(FIGURES)
    for name in wanted:
        try:
            if not FIGURES[name](cfg, figdir):
                skipped.append(name)
        except Exception as exc:
            log.warning("figura '%s' falhou (%s: %s)", name, type(exc).__name__, exc)
            skipped.append(name)

    log.info("%d/%d figuras geradas em %s", len(wanted) - len(skipped), len(wanted), figdir)
    if skipped:
        log.warning("figuras PULADAS por falta de dados: %s", ", ".join(skipped))


if __name__ == "__main__":
    main()
