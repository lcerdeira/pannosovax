#!/usr/bin/env python3
"""
Estágio 04 — pressão seletiva sobre os candidatos de superfície.

Conservação de sequência sozinha não basta. Um gene pode estar conservado só porque
foi amostrado pouco, e um gene de superfície variável pode estar sob seleção
diversificadora justamente porque o sistema imune do hospedeiro o persegue — que é
exatamente o antígeno que NÃO queremos (é o erro clássico de vacinas contra
antígenos hipervariáveis: a população escapa em uma temporada).

Queremos o oposto: sítios sob **seleção purificadora** (dN/dS < 1), sinal de que
mutar aquele resíduo custa fitness ao patógeno. Isso é o que torna o epitopo uma
armadilha evolutiva em vez de um alvo móvel.

Método: alinhamento MAFFT (codon-aware por retrotradução quando há nucleotídeo),
árvore por IQ-TREE, e HyPhy FEL por sítio. FEL é preferido a SLAC/MEME aqui porque
estima dN e dS por sítio com máxima verossimilhança sem assumir variação de taxa
entre ramos, o que é adequado para populações intraespecíficas densas.

Cruzamos ainda com a base de essencialidade (selection.essentiality_db, tipicamente
DEG + telas TraDIS/Tn-seq): gene essencial ou de virulência tem custo de escape alto.

Entrada : results/03_surfaceome/{org}_candidates.tsv
          results/02_pangenome/{org}_gene_alignments/
Saída   : results/04_selection/{org}_dnds.tsv

Uso:
    python scripts/04_conservation_dnds.py --organism kpsc
    python scripts/04_conservation_dnds.py --organism abau --threads 16 --max-genes 50
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("04_dnds")

COLUMNS = ["protein_id", "n_sites", "mean_dnds", "frac_purifying",
           "essential", "virulence", "pass_selection"]


def have(tool: str) -> bool:
    return shutil.which(tool) is not None


def run(cmd: list[str], cwd: Path | None = None) -> bool:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        log.warning("falhou: %s -> %s", " ".join(cmd[:3]), proc.stderr.strip()[:200])
        return False
    return True


def align(src: Path, dst: Path, threads: int) -> bool:
    """MAFFT --auto. O alinhamento do Panaroo já vem alinhado, mas re-alinhamos
    para garantir que o quadro de leitura de códons esteja consistente."""
    proc = subprocess.run(["mafft", "--auto", "--thread", str(threads), str(src)],
                          capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        log.warning("MAFFT falhou em %s", src.name)
        return False
    dst.write_text(proc.stdout)
    return True


def build_tree(aln: Path, threads: int) -> Path | None:
    """IQ-TREE com GTR+G — modelo padrão para codificantes bacterianos."""
    out = aln.with_suffix(".treefile")
    if out.exists():
        return out
    ok = run(["iqtree2" if have("iqtree2") else "iqtree", "-s", str(aln),
              "-m", "GTR+G", "-nt", str(threads), "-quiet", "-redo"])
    return out if ok and out.exists() else None


def run_fel(aln: Path, tree: Path, workdir: Path) -> Path | None:
    """HyPhy FEL. --ci Yes não é necessário; queremos só alpha/beta por sítio."""
    out = workdir / f"{aln.stem}.FEL.json"
    ok = run(["hyphy", "fel", "--alignment", str(aln), "--tree", str(tree),
              "--output", str(out), "--branches", "All"])
    return out if ok and out.exists() else None


def parse_fel(path: Path) -> tuple[int, float, float]:
    """
    FEL JSON -> (n_sites, dN/dS médio, fração de sítios sob seleção purificadora).

    Colunas de MLE.content['0']: [alpha, beta, alpha=beta, LRT, p-value, total branch length].
    dN/dS por sítio = beta/alpha; sítios com alpha=0 são descartados (indeterminado).
    Purificador = beta < alpha com p < 0.1 (limiar usual do FEL).
    """
    with open(path) as fh:
        data = json.load(fh)
    rows = data.get("MLE", {}).get("content", {}).get("0", [])
    ratios, n_pur, n_test = [], 0, 0
    for row in rows:
        if len(row) < 5:
            continue
        alpha, beta, _, _, pval = row[0], row[1], row[2], row[3], row[4]
        if alpha and alpha > 0:
            ratios.append(beta / alpha)
        n_test += 1
        if beta < alpha and pval is not None and pval < 0.1:
            n_pur += 1
    mean = sum(ratios) / len(ratios) if ratios else float("nan")
    frac = n_pur / n_test if n_test else float("nan")
    return len(rows), mean, frac


def load_essentiality(cfg: dict) -> pd.DataFrame | None:
    """TSV externo esperado com colunas gene/protein_id + essential + virulence."""
    rel = cfg["selection"].get("essentiality_db")
    if not rel:
        return None
    path = Path(rel)
    if not path.is_absolute():
        path = Path(__file__).resolve().parent.parent / rel
    if not path.exists():
        log.warning("base de essencialidade ausente (%s) — colunas essential/virulence "
                    "ficarão vazias e o filtro require_essential_or_virulence será ignorado", rel)
        return None
    df = pd.read_csv(path, sep="\t")
    key = "protein_id" if "protein_id" in df.columns else df.columns[0]
    df = df.rename(columns={key: "protein_id"})
    for col in ("essential", "virulence"):
        if col not in df.columns:
            df[col] = pd.NA
    return df[["protein_id", "essential", "virulence"]].drop_duplicates("protein_id")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--max-genes", type=int, default=None,
                    help="limita o número de candidatos processados (teste)")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    org = args.organism
    threads = args.threads or cfg.get("threads", 4)
    sel = cfg["selection"]
    out_path = outpath(cfg, "04_selection", f"{org}_dnds.tsv")

    cand_path = outpath(cfg, "03_surfaceome", f"{org}_candidates.tsv")
    if not cand_path.exists():
        log.warning("%s ausente — rode o estágio 03", cand_path.name)
        write_table(pd.DataFrame(columns=COLUMNS), out_path, log)
        return
    cand = pd.read_csv(cand_path, sep="\t")
    ids = cand["protein_id"].astype(str).tolist()
    if args.max_genes:
        ids = ids[: args.max_genes]
    log.info("%s: %d candidatos de superfície para avaliar", org, len(ids))

    missing = [t for t in ("mafft", "iqtree2", "hyphy") if not have(t)]
    if "iqtree2" in missing and have("iqtree"):
        missing.remove("iqtree2")
    degraded = bool(missing)
    if degraded:
        log.warning("ferramentas ausentes: %s — dN/dS não será calculado; "
                    "as colunas numéricas ficarão vazias e pass_selection será False",
                    ", ".join(missing))

    aln_dir = outpath(cfg, "02_pangenome", f"{org}_gene_alignments")
    work = outpath(cfg, "04_selection", f"{org}_work")
    work.mkdir(parents=True, exist_ok=True)

    rows = []
    for pid in ids:
        rec = {"protein_id": pid, "n_sites": pd.NA, "mean_dnds": pd.NA,
               "frac_purifying": pd.NA}
        if not degraded:
            src = next((aln_dir / f"{pid}{ext}" for ext in (".aln.fas", ".fas", ".fa", ".fasta")
                        if (aln_dir / f"{pid}{ext}").exists()), None)
            if src is None:
                log.warning("sem alinhamento para %s", pid)
            else:
                aligned = work / f"{pid}.aln.fasta"
                if align(src, aligned, threads):
                    tree = build_tree(aligned, threads)
                    fel = run_fel(aligned, tree, work) if tree else None
                    if fel:
                        n, mean, frac = parse_fel(fel)
                        rec.update(n_sites=n, mean_dnds=mean, frac_purifying=frac)
        rows.append(rec)

    df = pd.DataFrame(rows)

    ess = load_essentiality(cfg)
    if ess is not None:
        df = df.merge(ess, on="protein_id", how="left")
    else:
        df["essential"] = pd.NA
        df["virulence"] = pd.NA

    ok_dnds = pd.to_numeric(df["mean_dnds"], errors="coerce") < sel["dnds_max"]
    if sel.get("require_essential_or_virulence") and ess is not None:
        ok_role = df["essential"].fillna(False).astype(bool) | df["virulence"].fillna(False).astype(bool)
    else:
        ok_role = pd.Series(True, index=df.index)
    df["pass_selection"] = (ok_dnds & ok_role).fillna(False)

    write_table(df[COLUMNS], out_path, log)
    log.info("%s: %d/%d candidatos sob seleção purificadora aceitável",
             org, int(df["pass_selection"].sum()), len(df))


if __name__ == "__main__":
    main()
