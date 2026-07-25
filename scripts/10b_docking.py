#!/usr/bin/env python3
"""
Estágio 10b — docking do construto contra receptores da imunidade inata.

Se o adjuvante embutido (RS09, agonista de TLR4) não engajar o receptor, o construto
vira um peptídeo inerte. Docking não prova ativação, mas descarta o caso em que a
fusão enterra o adjuvante numa interface inacessível — que é uma falha de desenho
comum e barata de detectar.

Usamos HADDOCK 2.4 em vez de docking rígido porque HADDOCK é *guiado por informação*:
declaramos como "resíduos ativos" o adjuvante (no ligante) e o bolso de ligação
conhecido do receptor (no alvo), e o programa amostra ao redor disso. Docking cego
entre duas proteínas grandes gera poses sem sentido biológico com scores bonitos —
outro erro recorrente na literatura de vacinas in silico.

Aviso de interpretação: o HADDOCK score é uma soma ponderada empírica
(1.0*vdW + 0.2*elec + 1.0*desolv + 0.1*AIR); ele **não** é energia livre de ligação
e não deve ser reportado em kcal/mol. Comparações só valem entre poses do mesmo par.

Entrada : results/10_structure/construct_refined.pdb
          config docking.receptors (códigos PDB, baixados via RCSB REST se faltarem)
Saída   : results/10_docking/{receptor}/  (arquivos de job)
          results/10_docking/docking_scores.tsv

Uso:
    python scripts/10b_docking.py --pdb results/10_structure/construct_refined.pdb
    python scripts/10b_docking.py --pdb x.pdb --parse-only
"""
from __future__ import annotations

import argparse
import re
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("10b_docking")

COLUMNS = ["receptor", "haddock_score", "vdw", "elec", "desolv", "buried_sasa", "cluster_size"]
RCSB = "https://files.rcsb.org/download/{code}.pdb"

# Resíduos de interface conhecidos por receptor — é isso que transforma o docking
# cego em docking guiado. Numeração da cadeia A de cada entrada do PDB.
KNOWN_INTERFACE = {
    "TLR4": [263, 264, 265, 289, 290, 328, 329, 349, 351, 375, 376, 386],   # bolso MD-2/LPS
    "TLR2": [318, 319, 320, 321, 347, 349, 376, 377, 396, 397, 398],        # bolso de lipopeptídeo
    "MHCII_DRB1_0101": [9, 11, 13, 26, 28, 30, 47, 52, 58, 61, 65, 68, 71, 74],  # sulco peptídico
}


def fetch_pdb(code: str, dst: Path) -> bool:
    if dst.exists() and dst.stat().st_size > 0:
        return True
    url = RCSB.format(code=code.upper())
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            dst.write_bytes(resp.read())
        log.info("baixado %s -> %s", code, dst.name)
        return True
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        log.warning("falha ao baixar %s do RCSB (%s)", code, exc)
        return False


def construct_active_residues(pdb: Path, cfg: dict) -> list[int]:
    """Resíduos ativos do ligante = bloco do adjuvante, lido do mapa do estágio 08."""
    cmap = outpath(cfg, "08_construct", "construct_map.tsv")
    if cmap.exists():
        df = pd.read_csv(cmap, sep="\t")
        adj = df[df["element"] == "adjuvant"]
        if len(adj):
            r = adj.iloc[0]
            return list(range(int(r["start"]), int(r["end"]) + 1))
    log.warning("construct_map.tsv sem bloco 'adjuvant' — usando os 30 primeiros resíduos")
    return list(range(1, 31))


def write_air(path: Path, act_a: list[int], act_b: list[int]) -> None:
    """Restrições ambíguas (AIR) no formato CNS que o HADDOCK consome."""
    lines = []
    for lst, seg, other, oseg in ((act_a, "A", act_b, "B"), (act_b, "B", act_a, "A")):
        for res in lst:
            block = " or\n".join(
                f"        (resid {o} and segid {oseg})" for o in other
            )
            lines.append(
                f"assign (resid {res} and segid {seg})\n"
                f"       (\n{block}\n       )  2.0 2.0 0.0\n"
            )
    path.write_text("\n".join(lines))


def write_run_param(path: Path, rec_pdb: Path, lig_pdb: Path, air: Path, name: str) -> None:
    path.write_text(
        "\n".join([
            f"AMBIG_TBL={air.resolve()}",
            "HADDOCK_DIR=/opt/haddock2.4",
            "N_COMP=2",
            f"PDB_FILE1={rec_pdb.resolve()}",
            f"PDB_FILE2={lig_pdb.resolve()}",
            "PROJECT_DIR=./",
            f"PROT_SEGID_1=A",
            f"PROT_SEGID_2=B",
            f"RUN_NUMBER=1",
            "# HADDOCK 2.4: gere o run com  haddock2.4  neste diretório,",
            f"# depois edite run1/run.cns (structures_0=1000) e rode: cd run1 && csh run.cns",
            f"# projeto: {name}",
            "",
        ])
    )


def parse_haddock(stats: Path) -> dict | None:
    """
    Lê o structures/it1/water/*.stats ou o file.list do HADDOCK.

    Formato típico do .stats: cabeçalho com nomes de coluna e uma linha por cluster.
    Aceitamos também o padrão 'HADDOCK score ... -123.4 +/- 2.1' do cluster summary.
    """
    text = stats.read_text(errors="replace")
    out: dict[str, float] = {}
    patterns = {
        "haddock_score": r"HADDOCK score[^-\d]*(-?[\d.]+)",
        "vdw": r"Van der Waals energy[^-\d]*(-?[\d.]+)",
        "elec": r"Electrostatic energy[^-\d]*(-?[\d.]+)",
        "desolv": r"Desolvation energy[^-\d]*(-?[\d.]+)",
        "buried_sasa": r"Buried Surface Area[^-\d]*(-?[\d.]+)",
        "cluster_size": r"Cluster size[^-\d]*(-?[\d.]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text, re.IGNORECASE)
        if m:
            out[key] = float(m.group(1))
    return out or None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", required=True, help="modelo do construto (estágio 10)")
    ap.add_argument("--parse-only", action="store_true",
                    help="não prepara jobs, apenas coleta resultados existentes")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    lig = Path(args.pdb)
    receptors: dict[str, str] = cfg["docking"]["receptors"]
    base = outpath(cfg, "10_docking", "_").parent
    out_tsv = outpath(cfg, "10_docking", "docking_scores.tsv")

    if not lig.exists():
        log.warning("modelo do construto ausente (%s) — rode o estágio 10", lig)
        write_table(pd.DataFrame(columns=COLUMNS), out_tsv, log)
        return

    active_lig = construct_active_residues(lig, cfg)
    rows = []

    for name, code in receptors.items():
        jobdir = base / name
        jobdir.mkdir(parents=True, exist_ok=True)
        rec_pdb = jobdir / f"{code}.pdb"

        if not args.parse_only:
            if fetch_pdb(code, rec_pdb):
                air = jobdir / "ambig.tbl"
                write_air(air, KNOWN_INTERFACE.get(name, list(range(1, 21))), active_lig)
                write_run_param(jobdir / "run.param", rec_pdb, lig, air, name)
                log.info("job preparado em %s (receptor %s / PDB %s)", jobdir, name, code)
            else:
                log.warning("sem estrutura do receptor %s — job não preparado", name)

        stats = next(
            (p for p in [
                jobdir / "run1" / "structures" / "it1" / "water" / "file.list",
                jobdir / "run1" / "clusters_haddock-sorted.stat",
                jobdir / "haddock_results.txt",
            ] if p.exists()),
            None,
        )
        if stats:
            parsed = parse_haddock(stats)
            if parsed:
                rows.append({"receptor": name, **{c: parsed.get(c) for c in COLUMNS[1:]}})
                log.info("%s: HADDOCK score %s", name, parsed.get("haddock_score"))
                continue
        log.warning("sem resultados de HADDOCK para %s — rode o job em %s", name, jobdir)

    if not rows:
        log.warning("nenhum resultado de docking disponível; %s sai só com cabeçalho. "
                    "Prepare-se para rodar o HADDOCK 2.4 nos diretórios acima ou pelo "
                    "servidor web (https://wenmr.science.uu.nl/haddock2.4/).", out_tsv.name)
    write_table(pd.DataFrame(rows, columns=COLUMNS), out_tsv, log)


if __name__ == "__main__":
    main()
