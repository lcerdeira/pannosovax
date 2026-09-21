#!/usr/bin/env python3
"""
Estágio 04e — por que K. pneumoniae compartilha estrutura com o pneumococo e
A. baumannii não?

O crossmatch rendeu 182 janelas em kpsc x spneu e apenas 4 em abau x spneu, com
cobertura estrutural comparável (105 e 76 modelos). A diferença é grande demais para
passar sem explicação num manuscrito, e há duas famílias de causa possíveis.

Artefato do método, que invalidaria a comparação:
  A. qualidade das estruturas — abau recebeu 46 modelos de ESMFold contra 11 de kpsc, e
     ESMFold prediz a partir de sequência única, sem MSA. Se os modelos de abau tiverem
     pLDDT sistematicamente menor, o filtro `--min-plddt 70` do 04b poda as janelas de
     abau antes de qualquer biologia entrar em jogo;
  B. tamanho e número de proteínas comparáveis.

Biologia, que seria resultado:
  C. composição funcional — se as janelas kpsc x spneu se concentram numa classe (a
     hipótese é transportador ABC de ligação a substrato, cuja dobra é das mais
     conservadas entre bactérias distantes) e abau tiver poucas proteínas dessa classe
     entre as candidatas, a assimetria é real.

O script não decide sozinho: imprime as três evidências lado a lado.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, get_logger

log = get_logger("04e_assim")

ORGS = ["kpsc", "abau", "spneu"]
CLASSES = {
    "ABC substrate-binding": ["abc transporter substrate", "substrate-binding"],
    "porina/OM": ["porin", "outer membrane"],
    "TonB/sideróforo": ["tonb", "siderophore"],
    "adesina/pilus": ["adhesin", "pilus", "fimbri"],
    "lipoproteína": ["lipoprotein"],
    "peptidase/hidrolase": ["peptidase", "hydrolase", "amidase"],
}


def plddt_stats(pdb: Path) -> tuple[float, str]:
    """pLDDT médio (B-factor dos CA) e origem do modelo."""
    vals, source = [], "AFDB"
    try:
        for line in pdb.read_text().splitlines():
            if line.startswith("TITLE") and "ESMFOLD" in line.upper():
                source = "ESMFold"
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    except Exception:  # noqa: BLE001
        return float("nan"), "?"
    return (sum(vals) / len(vals) if vals else float("nan")), source


def classify(product: str) -> str:
    p = (product or "").lower()
    for name, kws in CLASSES.items():
        if any(k in p for k in kws):
            return name
    return "outra"


def main() -> int:
    windows = ROOT / "results/04_shared/shared_structural_epitopes_v2.tsv"
    if not windows.exists():
        log.error("%s ausente — rode o 04b", windows)
        return 2
    w = pd.read_csv(windows, sep="\t")

    # ── A. qualidade das estruturas ───────────────────────────────────────────
    print("\n=== A. qualidade dos modelos (o confundidor a descartar) ===")
    rows = []
    for org in ORGS:
        for pdb in sorted((ROOT / f"results/04_structures/{org}").glob("*.pdb")):
            m, src = plddt_stats(pdb)
            rows.append({"organism": org, "protein": pdb.stem, "plddt": m, "fonte": src})
    q = pd.DataFrame(rows)
    if q.empty:
        log.error("nenhuma estrutura encontrada")
        return 2
    print(q.groupby(["organism", "fonte"])["plddt"]
           .agg(n="size", media="mean", mediana="median",
                abaixo_de_70=lambda s: int((s < 70).sum())).round(1).to_string())
    print("\n  pLDDT médio por organismo:")
    print(q.groupby("organism")["plddt"].agg(n="size", media="mean").round(1).to_string())

    # ── B. o que entra em cada par ────────────────────────────────────────────
    print("\n=== B. janelas e proteínas por par ===")
    pares = w.groupby(["org_a", "org_b"]).agg(
        janelas=("tm_score", "size"), tm_medio=("tm_score", "mean"),
        tm_max=("tm_score", "max"), prot_a=("prot_a", "nunique"),
        prot_b=("prot_b", "nunique")).round(3)
    print(pares.to_string())

    # ── C. composição funcional ───────────────────────────────────────────────
    print("\n=== C. composição funcional das candidatas COM estrutura ===")
    comp = {}
    for org in ORGS:
        cand = pd.read_csv(ROOT / f"results/03_surfaceome/{org}_candidates.tsv", sep="\t")
        have = set(q[q["organism"] == org]["protein"])
        cand = cand[cand["protein_id"].astype(str).isin(have)]
        comp[org] = cand["product"].map(classify).value_counts()
    print(pd.DataFrame(comp).fillna(0).astype(int).to_string())

    print("\n=== C2. classe funcional das proteínas que DERAM match, por par ===")
    prod = {}
    for org in ORGS:
        c = pd.read_csv(ROOT / f"results/03_surfaceome/{org}_candidates.tsv", sep="\t")
        prod.update({f"{org}:{p}": classify(s)
                     for p, s in zip(c["protein_id"].astype(str), c["product"])})
    w = w.copy()
    w["classe_a"] = [prod.get(f"{a}:{p}", "?") for a, p in zip(w["org_a"], w["prot_a"])]
    w["classe_b"] = [prod.get(f"{b}:{p}", "?") for b, p in zip(w["org_b"], w["prot_b"])]
    for (a, b), g in w.groupby(["org_a", "org_b"]):
        print(f"\n  {a} x {b}  ({len(g)} janelas)")
        top = (g.groupby(["classe_a", "classe_b"]).size()
                .sort_values(ascending=False).head(5))
        print(top.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
