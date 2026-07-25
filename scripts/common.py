"""Utilitários compartilhados pelo pipeline PanNosoVax."""
from __future__ import annotations

import logging
import os
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config" / "config.yaml"


def load_config(path: str | os.PathLike | None = None) -> dict:
    with open(path or CONFIG_PATH) as fh:
        return yaml.safe_load(fh)


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=os.environ.get("PANNOSOVAX_LOGLEVEL", "INFO"),
        format="%(asctime)s %(levelname)-7s [%(name)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger(name)


def outpath(cfg: dict, *parts: str) -> Path:
    p = ROOT / cfg["outdir"]
    for part in parts:
        p = p / part
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def write_table(df: pd.DataFrame, path: Path, log: logging.Logger | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    if log:
        log.info("escrito %s (%d linhas)", path.relative_to(ROOT), len(df))


def read_alleles(path: str | os.PathLike) -> pd.DataFrame:
    """Lê config/alleles_mhc*.txt -> DataFrame[allele, freq_world, freq_brazil]."""
    rows = []
    with open(ROOT / path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            rows.append(
                {
                    "allele": parts[0],
                    "freq_world": float(parts[1]) if len(parts) > 1 else float("nan"),
                    "freq_brazil": float(parts[2]) if len(parts) > 2 else float("nan"),
                }
            )
    return pd.DataFrame(rows)


GRAM = {"kpsc": "negative", "abau": "negative", "spneu": "positive"}
