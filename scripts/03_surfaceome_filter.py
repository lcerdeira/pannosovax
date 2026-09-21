#!/usr/bin/env python3
"""
Estágio 03 — do core genome ao "surfaceome vacinável".

Um antígeno só é útil se o anticorpo consegue alcançá-lo. Filtramos por:
  1. localização subcelular (PSORTb) compatível com a parede celular do organismo;
  2. presença de peptídeo-sinal (SignalP 6) OU ancoragem de superfície reconhecida;
  3. no máximo N hélices transmembrana (DeepTMHMM) — proteína politópica não expressa bem;
  4. faixa de tamanho tratável;
  5. bônus funcional para adesinas, porinas e receptores de sideróforo, que são as
     classes com maior taxa histórica de sucesso em vacinas bacterianas.

Entrada : results/02_pangenome/{org}_core_proteins.faa
          results/03_surfaceome/{org}_psortb.tsv, _signalp.tsv, _tmhmm.tsv
Saída   : results/03_surfaceome/{org}_candidates.tsv
"""
from __future__ import annotations

import argparse
import re

import pandas as pd
from Bio import SeqIO

from common import GRAM, get_logger, load_config, outpath, write_table

log = get_logger("03_surfaceome")


def parse_psortb(path) -> pd.DataFrame:
    """PSORTb saída 'long' -> DataFrame[protein_id, localization, psortb_score]."""
    records, cur = [], {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("SeqID:"):
                if cur:
                    records.append(cur)
                cur = {"protein_id": line.split(":", 1)[1].strip().split()[0]}
            elif re.match(r"^\s*Final Prediction:", line):
                cur["_await_final"] = True
            elif cur.get("_await_final") and line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    cur["localization"] = parts[0]
                    try:
                        cur["psortb_score"] = float(parts[-1])
                    except ValueError:
                        cur["psortb_score"] = float("nan")
                    cur.pop("_await_final", None)
    if cur:
        records.append(cur)
    return pd.DataFrame(records).drop(columns=["_await_final"], errors="ignore")


def parse_signalp(path) -> pd.DataFrame:
    """SignalP 6 summary tabular."""
    df = pd.read_csv(path, sep="\t", comment="#", header=None,
                     names=["protein_id", "prediction", "p_other", "p_sp",
                            "p_lipo", "p_tat", "cs_position"], engine="python")
    df["has_signal"] = df["prediction"].ne("OTHER")
    return df[["protein_id", "prediction", "has_signal", "cs_position"]]


def parse_tmhmm(path) -> pd.DataFrame:
    """DeepTMHMM gff3 ou TMHMM short -> contagem de hélices TM."""
    counts: dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 3 and parts[2].strip() in {"TMhelix", "TMhelix_out", "TMhelix_in"}:
                counts[parts[0]] = counts.get(parts[0], 0) + 1
            elif "PredHel=" in line:                       # formato TMHMM short
                pid = parts[0]
                m = re.search(r"PredHel=(\d+)", line)
                if m:
                    counts[pid] = int(m.group(1))
    return pd.DataFrame({"protein_id": list(counts), "n_tm_helices": list(counts.values())})


def functional_score(product: str, keywords: list[str]) -> int:
    product = (product or "").lower()
    return sum(1 for kw in keywords if kw.lower() in product)


def excluded_by_annotation(product: str, keywords: list[str]) -> bool:
    """Proteína que a anotação identifica como citoplasmática abundante.

    O DeepLocPro rotulou 12 proteínas ribossomais do pneumococo como
    "Extracellular" — 67% das candidatas daquele organismo. São falsos-positivos
    conhecidos: ribossomos, chaperonas e fatores de alongamento são abundantes e
    confundem preditores treinados em composição de aminoácidos. Nenhuma é
    alcançável por anticorpo, e cada uma ocupa uma vaga na seleção por cobertura.
    """
    product = (product or "").lower()
    return any(kw.lower() in product for kw in keywords)


LPXTG = re.compile(r"LP.TG")


def surface_anchored_gram_positive(seq: str, signalp_pred: str) -> bool:
    """Evidência de ancoragem exposta na superfície de um Gram-positivo.

    Duas formas reconhecidas: lipobox (SignalP-6 chama LIPO — via Sec/SPII, que
    ancora a proteína no folheto EXTERNO da membrana, como PsaA) e motivo LPXTG
    reconhecido por sortase, que liga covalentemente à parede (como PspC). Sem
    membrana externa para encobri-las, ambas ficam acessíveis ao anticorpo.
    """
    if str(signalp_pred).upper() == "LIPO":
        return True
    # o motivo de sortase fica no extremo C-terminal, antes da cauda hidrofóbica
    return bool(LPXTG.search(str(seq)[-60:]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    sf = cfg["surfaceome"]
    org = args.organism
    base = outpath(cfg, "03_surfaceome", org).parent

    faa = outpath(cfg, "02_pangenome", f"{org}_core_proteins.faa")
    seqs = {r.id: r for r in SeqIO.parse(faa, "fasta")}
    log.info("%s: %d proteínas core na entrada", org, len(seqs))

    df = pd.DataFrame({
        "protein_id": list(seqs),
        "length": [len(r.seq) for r in seqs.values()],
        "product": [r.description.split(" ", 1)[1] if " " in r.description else ""
                    for r in seqs.values()],
    })

    for parser, fname, label in [
        (parse_psortb, f"{org}_psortb.tsv", "psortb"),
        (parse_signalp, f"{org}_signalp.tsv", "signalp"),
        (parse_tmhmm, f"{org}_tmhmm.tsv", "tmhmm"),
    ]:
        path = base / fname
        if path.exists():
            df = df.merge(parser(path), on="protein_id", how="left")
        else:
            log.warning("%s ausente — estágio 03 depende dele (%s)", fname, label)

    df["n_tm_helices"] = df.get("n_tm_helices", pd.Series(0, index=df.index)).fillna(0)
    allowed = sf["gram_negative_ok"] if GRAM[org] == "negative" else sf["gram_positive_ok"]

    score_ok = df.get("psortb_score", 0) >= sf["psortb_min_score"]
    df["pass_localization"] = df.get("localization", "").isin(allowed) & score_ok

    # Gram-positivo: recuperar as ancoradas à superfície que o preditor mandou para
    # a membrana citoplasmática. Em Gram-negativo esse rótulo significa membrana
    # interna e continua sendo motivo de exclusão.
    df["anchored_surface"] = False
    if GRAM[org] == "positive" and sf.get("gram_positive_membrane_if_anchored"):
        at_membrane = df.get("localization", "").eq(sf["gram_positive_membrane_loc"]) & score_ok
        df["anchored_surface"] = at_membrane & df.apply(
            lambda r: surface_anchored_gram_positive(
                seqs[r["protein_id"]].seq, r.get("prediction", "")), axis=1)
        df["pass_localization"] = df["pass_localization"] | df["anchored_surface"]
        log.info("%s: %d proteínas de membrana recuperadas por ancoragem de superfície "
                 "(lipobox ou LPXTG)", org, int(df["anchored_surface"].sum()))

    df["excluded_by_annotation"] = df["product"].map(
        lambda p: excluded_by_annotation(p, sf.get("exclude_annotation_keywords", []))
    )
    n_excl = int((df["excluded_by_annotation"] & df["pass_localization"]).sum())
    if n_excl:
        log.info("%s: %d proteínas descartadas por anotação citoplasmática apesar da "
                 "localização predita", org, n_excl)
    df["pass_topology"] = df["n_tm_helices"] <= sf["max_tm_helices"]
    df["pass_length"] = df["length"].between(sf["min_length"], sf["max_length"])
    df["pass_export"] = df.get("has_signal", False).fillna(False) | df["pass_localization"]
    df["functional_score"] = df["product"].map(
        lambda p: functional_score(p, sf["functional_boost_keywords"])
    )

    df["candidate"] = (df["pass_localization"] & df["pass_topology"]
                       & df["pass_length"] & df["pass_export"]
                       & ~df["excluded_by_annotation"])
    df = df.sort_values(["candidate", "functional_score", "psortb_score"],
                        ascending=False)

    write_table(df, base / f"{org}_surfaceome_full.tsv", log)
    write_table(df[df["candidate"]], base / f"{org}_candidates.tsv", log)

    kept = int(df["candidate"].sum())
    log.info("%s: %d/%d proteínas aprovadas como candidatas de superfície (%.1f%%)",
             org, kept, len(df), 100 * kept / max(len(df), 1))


if __name__ == "__main__":
    main()
