#!/usr/bin/env python3
"""
Figura S1 — composição geográfica real dos genomas disponíveis.

Esta figura carrega um argumento próprio do artigo: o viés de sequenciamento nos
genomas **completos** tem forma diferente do viés já descrito para o conjunto que
inclui montagens rascunho. Reportá-lo é uma contribuição, não um detalhe de métodos.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from common import get_logger, load_config, outpath

log = get_logger("fig_s1")

LABELS = {"kpsc": "KpSC", "abau": "$A.\\ baumannii$", "spneu": "$S.\\ pneumoniae$"}
REGION_PT = {
    "asia": "Ásia", "europe": "Europa", "north_america": "América do Norte",
    "latin_america": "América Latina", "africa": "África", "other": "Outros/NI",
}
ORDER = ["asia", "europe", "north_america", "latin_america", "africa", "other"]
COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#B0B0B0"]


def main() -> None:
    cfg = load_config()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))

    summary = []
    for ax, org in zip(axes, ["kpsc", "abau", "spneu"]):
        path = outpath(cfg, "01_genomes", f"{org}_all_metadata.tsv")
        if not path.exists():
            log.warning("faltando %s — pulando", path)
            continue
        d = pd.read_csv(path, sep="\t")
        n = len(d)
        counts = d["region"].value_counts()
        vals = [counts.get(r, 0) for r in ORDER]
        pct = [100 * v / n for v in vals]

        bars = ax.barh([REGION_PT[r] for r in ORDER], pct, color=COLORS)
        ax.invert_yaxis()
        ax.set_xlabel("% dos genomas")
        ax.set_title(f"{LABELS[org]}\nn = {n:,} genomas completos".replace(",", "."),
                     fontsize=11)
        ax.set_xlim(0, max(pct) * 1.25)
        for bar, v, p in zip(bars, vals, pct):
            ax.text(bar.get_width() + max(pct) * 0.02, bar.get_y() + bar.get_height() / 2,
                    f"{v}", va="center", fontsize=8.5)
        ax.spines[["top", "right"]].set_visible(False)

        npais = d["country"].dropna().str.split(":").str[0].str.strip().nunique()
        for r in ORDER:
            summary.append({"organismo": LABELS[org].replace("$", "").replace("\\ ", " "),
                            "regiao": REGION_PT[r], "n": counts.get(r, 0),
                            "pct": round(100 * counts.get(r, 0) / n, 2),
                            "total_genomas": n, "paises_distintos": npais})

    fig.suptitle("Composição geográfica dos genomas completos disponíveis no RefSeq "
                 "(julho de 2026)", fontsize=12, y=1.02)
    fig.tight_layout()

    out = outpath(cfg, "report", "figures", "figS1_geografia.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    pd.DataFrame(summary).to_csv(
        outpath(cfg, "report", "tabelaS1_geografia.tsv"), sep="\t", index=False)
    log.info("figura S1 e tabela S1 escritas em %s", out.parent)


if __name__ == "__main__":
    main()
