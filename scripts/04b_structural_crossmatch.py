#!/usr/bin/env python3
"""
Estágio 04b — epitopos estruturalmente compartilhados entre os três organismos.

*** ESTÁGIO EXPLORATÓRIO. ***  Não é um resultado consolidado do pipeline; é uma
hipótese testável que o construto carrega em bloco separado, para que possa ser
avaliada (ou descartada) sem contaminar o resto do desenho.

A ideia: K. pneumoniae, A. baumannii e S. pneumoniae não compartilham identidade de
sequência apreciável nos antígenos de superfície — buscar epitopos pan-patógeno por
BLAST não dá nada, e é por isso que a literatura de vacinas multi-patógeno
essencialmente não existe. Mas proteínas de superfície convergem em **dobra**:
barris beta de porina, domínios de ligação a substrato tipo-ABC, receptores
TonB-dependentes. Se duas proteínas de organismos diferentes se superpõem
estruturalmente (TM-score alto), existem alças de superfície ocupando a mesma
posição espacial mesmo com sequências distintas.

Um anticorpo reconhece superfície, não sequência. A hipótese é que um peptídeo
desenhado sobre essa região comum possa induzir anticorpos com reatividade cruzada.
Ela pode estar errada — reatividade cruzada estrutural entre epítopos lineares
curtos é rara. Por isso o rótulo EXPLORATÓRIO é registrado na própria saída.

Critério: TM-score >= structural_tm_threshold (padrão 0.5, o ponto acima do qual
duas estruturas têm quase certamente a mesma dobra, Xu & Zhang 2010).

Entrada : results/04_structures/{org}/*.pdb   (modelos preditos)
          results/04_selection/{org}_dnds.tsv (ranking dos candidatos)
Saída   : results/04_shared/shared_structural_epitopes.tsv

Uso:
    python scripts/04b_structural_crossmatch.py
    python scripts/04b_structural_crossmatch.py --top 10 --window 15
"""
from __future__ import annotations

import argparse
import itertools
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("04b_crossmatch")

COLUMNS = ["peptide", "organisms", "tm_score", "note"]
DEFAULT_TM = 0.5

# Resíduo -> código de uma letra, para extrair o peptídeo direto do PDB.
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def top_candidates(cfg: dict, org: str, n: int) -> list[str]:
    p = outpath(cfg, "04_selection", f"{org}_dnds.tsv")
    if not p.exists():
        log.warning("%s ausente — usando todos os modelos disponíveis de %s", p.name, org)
        return []
    df = pd.read_csv(p, sep="\t")
    if "pass_selection" in df.columns:
        df = df[df["pass_selection"].fillna(False).astype(bool)]
    if "frac_purifying" in df.columns:
        df = df.sort_values("frac_purifying", ascending=False)
    return df["protein_id"].astype(str).head(n).tolist()


def models_for(cfg: dict, org: str, keep: list[str]) -> list[Path]:
    d = outpath(cfg, "04_structures", org, "_").parent
    if not d.is_dir():
        return []
    pdbs = sorted(d.glob("*.pdb"))
    if keep:
        wanted = set(keep)
        sel = [p for p in pdbs if p.stem in wanted]
        if sel:
            return sel
    return pdbs


def run_tmalign(a: Path, b: Path) -> float | None:
    """TM-align: pegamos o TM-score normalizado pela cadeia mais curta (o maior dos dois)."""
    proc = subprocess.run(["TMalign", str(a), str(b)],
                          capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        return None
    scores = [float(m) for m in re.findall(r"TM-score=\s*([0-9.]+)", proc.stdout)]
    return max(scores) if scores else None


def chain_residues(pdb: Path) -> list[tuple[int, str, float]]:
    """CA de cada resíduo: (numero, aa1, b-factor). O B-factor de modelo ColabFold é o pLDDT."""
    out = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                aa = AA3.get(line[17:20].strip().upper())
                if not aa:
                    continue
                try:
                    out.append((int(line[22:26]), aa, float(line[60:66])))
                except ValueError:
                    continue
    return out


def surface_windows(pdb: Path, window: int, min_plddt: float) -> list[str]:
    """
    Janelas contíguas de alta confiança usadas como proxy de alça exposta.

    Proxy deliberado e grosseiro: sem DSSP/SASA aqui, usamos pLDDT alto como
    indicador de região bem definida. Um cálculo real de SASA é o próximo passo se
    a hipótese sobreviver.
    """
    res = chain_residues(pdb)
    peps = []
    for i in range(len(res) - window + 1):
        block = res[i:i + window]
        if min(b for _, _, b in block) >= min_plddt:
            peps.append("".join(aa for _, aa, _ in block))
    return peps


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top", type=int, default=15, help="candidatos por organismo")
    ap.add_argument("--window", type=int, default=12, help="tamanho do peptídeo emitido")
    ap.add_argument("--min-plddt", type=float, default=70.0)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    thr = cfg.get("structure", {}).get("structural_tm_threshold",
                                       cfg.get("structural_tm_threshold", DEFAULT_TM))
    out_path = outpath(cfg, "04_shared", "shared_structural_epitopes.tsv")

    log.warning("=== ESTÁGIO EXPLORATÓRIO === a hipótese de reatividade cruzada por "
                "sobreposição estrutural NÃO está validada; os peptídeos emitidos aqui "
                "entram no construto em bloco separado e rotulado.")
    log.info("limiar TM-score = %.2f", thr)

    orgs = list(cfg["organisms"])
    models = {org: models_for(cfg, org, top_candidates(cfg, org, args.top)) for org in orgs}
    for org in orgs:
        log.info("%s: %d modelos PDB encontrados", org, len(models[org]))

    if not shutil.which("TMalign") or not all(models.values()):
        motivo = "TMalign não está no PATH" if not shutil.which("TMalign") \
            else "faltam modelos preditos para pelo menos um organismo"
        log.warning("%s — escrevendo %s só com cabeçalho", motivo, out_path.name)
        write_table(pd.DataFrame(columns=COLUMNS), out_path, log)
        return

    # Pares entre organismos diferentes apenas; pares intra-organismo não informam nada.
    pair_hits: dict[tuple[str, str], list[tuple[Path, Path, float]]] = {}
    for oa, ob in itertools.combinations(orgs, 2):
        hits = []
        for pa in models[oa]:
            for pb in models[ob]:
                tm = run_tmalign(pa, pb)
                if tm is not None and tm >= thr:
                    hits.append((pa, pb, tm))
        pair_hits[(oa, ob)] = sorted(hits, key=lambda h: -h[2])
        log.info("%s x %s: %d pares com TM-score >= %.2f", oa, ob, len(hits), thr)

    rows = []
    for (oa, ob), hits in pair_hits.items():
        for pa, pb, tm in hits[:20]:
            common = set(surface_windows(pa, args.window, args.min_plddt)) & \
                set(surface_windows(pb, args.window, args.min_plddt))
            for pep in sorted(common):
                rows.append({
                    "peptide": pep,
                    "organisms": f"{oa}|{ob}",
                    "tm_score": round(tm, 3),
                    "note": f"EXPLORATORIO: {pa.stem} vs {pb.stem}, janela {args.window} aa, "
                            f"pLDDT>={args.min_plddt}",
                })

    df = pd.DataFrame(rows, columns=COLUMNS)
    if not df.empty:
        # Um peptídeo visto em mais de um par cobre mais de dois organismos: prioridade.
        df = df.drop_duplicates(["peptide", "organisms"]).sort_values("tm_score", ascending=False)
    else:
        log.warning("nenhuma janela de superfície comum encontrada — a hipótese não se "
                    "sustentou com estes modelos e limiares")
    write_table(df, out_path, log)


if __name__ == "__main__":
    main()
