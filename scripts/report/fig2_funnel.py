#!/usr/bin/env python3
"""Figura 2 e Tabela 2 — funil de filtragem com os números reais do pipeline."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("fig2")

LAB = {"kpsc": "KpSC", "abau": "A. baumannii", "spneu": "S. pneumoniae"}
COL = {"kpsc": "#4C72B0", "abau": "#DD8452", "spneu": "#55A868"}


def main() -> None:
    cfg = load_config()
    rows = []
    for org in ["kpsc", "abau", "spneu"]:
        pres = pd.read_csv(outpath(cfg, "02_pangenome", f"{org}_presence.tsv"), sep="\t")
        ann = pd.read_csv(outpath(cfg, "03_surfaceome", f"{org}_annotation_pass.tsv"),
                          sep="\t")
        meta = pd.read_csv(outpath(cfg, "01_genomes", f"{org}_all_metadata.tsv"), sep="\t")
        rows.append({
            "organismo": LAB[org],
            "genomas_disponiveis": len(meta),
            "genomas_analisados": int(pres["n_genomes_compared"].max()) + 1,
            "proteinas_referencia": len(pres),
            "core_genome": int(pres["is_core"].sum()),
            "excluidas_anotacao": int(ann["excluded_by_annotation"].sum()),
            "candidatas_superficie": int(ann["surface_candidate"].sum()),
        })
    t = pd.DataFrame(rows)
    write_table(t, outpath(cfg, "report", "tabela2_funil.tsv"), log)
    print(t.to_string(index=False))

    steps = ["proteinas_referencia", "core_genome", "candidatas_superficie"]
    names = ["Proteínas da\nreferência", "Core genome\n(≥95% dos genomas)",
             "Candidatas de\nsuperfície"]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    x = range(len(steps))
    for _, r in t.iterrows():
        org = [k for k, v in LAB.items() if v == r["organismo"]][0]
        vals = [r[s] for s in steps]
        ax.plot(x, vals, "o-", lw=2.2, ms=8, color=COL[org],
                label=f"{r['organismo']} (n={r['genomas_analisados']} genomas)")
        for xi, v in zip(x, vals):
            ax.annotate(f"{v:,}".replace(",", "."), (xi, v),
                        textcoords="offset points", xytext=(0, 9),
                        ha="center", fontsize=8.5, color=COL[org])
    ax.set_yscale("log")
    ax.set_xticks(list(x))
    ax.set_xticklabels(names, fontsize=9.5)
    ax.set_ylabel("Número de proteínas (escala log)")
    ax.set_title("Funil de filtragem do PanNosoVax", fontsize=12)
    ax.legend(fontsize=9, frameon=False)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    out = outpath(cfg, "report", "figures", "fig2_funil.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    log.info("figura 2 escrita")


if __name__ == "__main__":
    main()
