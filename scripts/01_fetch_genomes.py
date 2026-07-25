#!/usr/bin/env python3
"""
Estágio 01 — aquisição de genomas com amostragem estratificada.

A maioria dos estudos de vacinologia reversa baixa "os N primeiros genomas" do NCBI,
o que enviesa fortemente para isolados de Europa/EUA e para clones de surto
super-sequenciados (ST258 em KpSC, GC2 em A. baumannii). Um antígeno que parece
conservado nesse conjunto pode não ser conservado no mundo real.

Aqui fazemos amostragem estratificada por região geográfica e por sequence type,
com teto por ST para não deixar um clone dominar.

Uso:
    python scripts/01_fetch_genomes.py --organism kpsc
    python scripts/01_fetch_genomes.py --all --dry-run
"""
from __future__ import annotations

import argparse
import json
import time
import subprocess
from collections import defaultdict

import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("01_fetch")

REGION_MAP = {
    "latin_america": {"Brazil", "Argentina", "Chile", "Colombia", "Mexico", "Peru",
                      "Uruguay", "Paraguay", "Bolivia", "Ecuador", "Venezuela",
                      "Costa Rica", "Cuba", "Guatemala", "Panama"},
    "africa": {"Nigeria", "South Africa", "Kenya", "Ghana", "Egypt", "Ethiopia",
               "Tanzania", "Uganda", "Malawi", "Senegal", "Morocco", "Tunisia",
               "Zambia", "Mozambique", "Gambia"},
    "asia": {"China", "India", "Japan", "South Korea", "Thailand", "Vietnam",
             "Pakistan", "Bangladesh", "Indonesia", "Malaysia", "Singapore",
             "Philippines", "Nepal", "Taiwan", "Israel", "Turkey", "Iran", "Saudi Arabia"},
    "europe": {"United Kingdom", "France", "Germany", "Italy", "Spain", "Netherlands",
               "Portugal", "Poland", "Sweden", "Norway", "Denmark", "Belgium",
               "Switzerland", "Austria", "Greece", "Czech Republic", "Ireland", "Russia"},
    "north_america": {"USA", "United States", "Canada"},
}


def region_of(country: str | None) -> str:
    if not country:
        return "other"
    country = country.split(":")[0].strip()
    for region, members in REGION_MAP.items():
        if country in members:
            return region
    return "other"


NCBI_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

# A API usa COMPLETE_GENOME/CHROMOSOME; o config usa a forma curta.
LEVEL_MAP = {"complete": "complete_genome", "chromosome": "chromosome",
             "scaffold": "scaffold", "contig": "contig"}


def fetch_metadata(taxid: int, assembly_levels: list[str]) -> pd.DataFrame:
    """
    Consulta a REST API do NCBI Datasets v2 com paginação.

    Optamos pela API em vez do binário `datasets` para eliminar uma dependência de
    instalação e tornar o estágio reprodutível em qualquer máquina com rede.
    """
    import requests

    levels = [LEVEL_MAP.get(x, x) for x in assembly_levels]
    rows, token, page = [], None, 0

    while True:
        # a API espera o parâmetro repetido, não uma lista separada por vírgula
        params = {
            "filters.assembly_level": levels,
            "filters.assembly_source": "refseq",
            "page_size": 1000,
        }
        if token:
            params["page_token"] = token

        resp = requests.get(f"{NCBI_API}/genome/taxon/{taxid}/dataset_report",
                            params=params, timeout=180)
        resp.raise_for_status()
        data = resp.json()

        reports = data.get("reports", [])
        if page == 0:
            log.info("taxid=%s: %s genomas disponíveis no NCBI",
                     taxid, data.get("total_count", "?"))

        for rec in reports:
            ai = rec.get("assembly_info", {})
            # nem todo atributo de BioSample traz 'value' — registros antigos omitem
            attrs = {a.get("name"): a.get("value")
                     for a in ai.get("biosample", {}).get("attributes", [])
                     if a.get("name")}
            rows.append({
                "accession": rec.get("accession"),
                "organism": rec.get("organism", {}).get("organism_name"),
                "strain": rec.get("organism", {}).get("infraspecific_names", {}).get("strain"),
                "level": ai.get("assembly_level"),
                "country": attrs.get("geo_loc_name") or attrs.get("country"),
                "collection_date": attrs.get("collection_date"),
                "host": attrs.get("host"),
                "isolation_source": attrs.get("isolation_source"),
                "release_date": ai.get("release_date"),
                "contig_n50": rec.get("assembly_stats", {}).get("contig_n50"),
                "total_length": rec.get("assembly_stats", {}).get("total_sequence_length"),
            })

        token = data.get("next_page_token")
        page += 1
        if not token:
            break
        time.sleep(0.4)          # cortesia com a API pública

    df = pd.DataFrame(rows)
    if df.empty:
        log.warning("taxid=%s não retornou genomas", taxid)
        return df
    df["region"] = df["country"].map(region_of)
    return df


def stratified_sample(df: pd.DataFrame, n: int, strata: dict[str, float], seed: int = 42):
    """Amostra n genomas respeitando as proporções por região, com fallback."""
    rng = pd.Series(range(len(df))).sample(frac=1, random_state=seed).index
    df = df.iloc[rng].reset_index(drop=True)

    quotas = {region: int(round(n * frac)) for region, frac in strata.items()}
    picked, counts = [], defaultdict(int)

    for region, quota in quotas.items():
        pool = df[df["region"] == region]
        take = pool.head(quota)
        picked.append(take)
        counts[region] = len(take)
        if len(take) < quota:
            log.warning("região %s: só %d/%d genomas disponíveis", region, len(take), quota)

    chosen = pd.concat(picked) if picked else df.head(0)
    deficit = n - len(chosen)
    if deficit > 0:
        rest = df[~df["accession"].isin(chosen["accession"])].head(deficit)
        chosen = pd.concat([chosen, rest])
        log.info("preenchendo déficit de %d genomas com amostra irrestrita", deficit)

    return chosen.reset_index(drop=True), dict(counts)


def download(accessions: list[str], dest) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    zip_path = dest / "genomes.zip"
    cmd = ["datasets", "download", "genome", "accession", *accessions,
           "--include", "genome,protein,gff3", "--filename", str(zip_path)]
    log.info("baixando %d genomas -> %s", len(accessions), zip_path)
    subprocess.run(cmd, check=True)
    subprocess.run(["unzip", "-q", "-o", str(zip_path), "-d", str(dest)], check=True)
    zip_path.unlink()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="só metadados, sem download")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    targets = list(cfg["organisms"]) if args.all else [args.organism]
    if not targets or targets == [None]:
        ap.error("informe --organism ou --all")

    summary = []
    for org in targets:
        spec = cfg["organisms"][org]
        meta = pd.concat(
            [fetch_metadata(t, spec["assembly_level"]) for t in spec["taxids"]],
            ignore_index=True,
        ).drop_duplicates("accession")
        log.info("%s: %d genomas candidatos", org, len(meta))

        chosen, counts = stratified_sample(
            meta, spec["n_genomes"], cfg["geographic_strata"]
        )
        write_table(meta, outpath(cfg, "01_genomes", f"{org}_all_metadata.tsv"), log)
        write_table(chosen, outpath(cfg, "01_genomes", f"{org}_selected.tsv"), log)

        summary.append({"organism": org, "candidates": len(meta),
                        "selected": len(chosen), **counts})

        if not args.dry_run:
            download(chosen["accession"].tolist(),
                     outpath(cfg, "01_genomes", org).parent / org)

    write_table(pd.DataFrame(summary),
                outpath(cfg, "01_genomes", "sampling_summary.tsv"), log)
    log.info("estágio 01 concluído")


if __name__ == "__main__":
    main()
