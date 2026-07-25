#!/usr/bin/env python3
"""
Relatório — preenchimento dos marcadores ⟨PENDENTE:chave⟩ do manuscrito.

Regra inegociável: **nunca inventar um valor**. Se o TSV que alimenta um marcador
não existe, ou existe só com cabeçalho, o marcador fica exatamente como está e o
script imprime qual arquivo o resolveria. Um manuscrito com marcador visível é
constrangedor; um manuscrito com número fabricado é fraude.

Marcadores textuais (conclusao, discussao_dados, financiamento, conflitos,
contribuicoes, referencias, repo_url) são de autoria humana — este script não os
toca, apenas os lista como pendências abertas.

Atenção ao caractere: os delimitadores são U+27E8 ⟨ e U+27E9 ⟩ (MATHEMATICAL LEFT/
RIGHT ANGLE BRACKET), não '<' '>' nem ‹›. Buscar pelos errados é o motivo usual de
o script "não achar nada".

Entrada : manuscript/manuscript.md + results/**/*.tsv
Saída   : manuscript/manuscript_filled.md

Uso:
    python scripts/report/fill_manuscript.py
    python scripts/report/fill_manuscript.py --input manuscript/manuscript.md --strict
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from common import ROOT, get_logger, load_config, outpath  # noqa: E402

log = get_logger("fill_manuscript")

MARKER = re.compile("⟨PENDENTE:([A-Za-z0-9_]+)⟩")

# Marcadores que dependem de texto humano, não de dados.
AUTHORED = {"conclusao", "conclusion_summary", "discussao_dados", "financiamento",
            "conflitos", "contribuicoes", "referencias", "repo_url", "xx"}

# Onde procurar cada marcador resolvível por dados.
SOURCES = {
    "n_genomes_total": "results/01_genomes/{org}_selected.tsv",
    "md_ns": "config/config.yaml (md.ns)",
    "results_summary": "results/07_coverage/selected_*.tsv + results/09_physchem/*",
    "secao_resultados": "results/0*/*.tsv (contagens dos estágios 02-08)",
}


def n_rows(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    return len(df) or None


def resolve_n_genomes_total(cfg) -> str | None:
    total = 0
    for org in cfg["organisms"]:
        n = n_rows(outpath(cfg, "01_genomes", f"{org}_selected.tsv"))
        if n is None:
            log.info("  n_genomes_total: falta %s_selected.tsv", org)
            return None
        total += n
    return f"{total:,}".replace(",", ".")


def resolve_md_ns(cfg) -> str | None:
    # Vem do config, mas só se a MD tiver de fato rodado — declarar 100 ns sem
    # trajetória seria afirmar um resultado inexistente.
    summ = outpath(cfg, "11_md", "md_summary.tsv")
    if n_rows(summ) is None:
        log.info("  md_ns: md_summary.tsv ausente/vazio — a simulação não rodou")
        return None
    return str(cfg["md"]["ns"])


def stage_counts(cfg) -> dict[str, int]:
    counts: dict[str, int] = {}
    for org in cfg["organisms"]:
        for key, sub, fname in [
            (f"core_{org}", "02_pangenome", f"{org}_gene_presence_absence.csv"),
            (f"cand_{org}", "03_surfaceome", f"{org}_candidates.tsv"),
            (f"sel_{org}", "04_selection", f"{org}_dnds.tsv"),
        ]:
            n = n_rows(outpath(cfg, sub, fname))
            if n:
                counts[key] = n
    for klass in ("mhc1", "mhc2"):
        n = n_rows(outpath(cfg, "07_coverage", f"selected_{klass}.tsv"))
        if n:
            counts[klass] = n
    return counts


def resolve_results_summary(cfg) -> str | None:
    c = stage_counts(cfg)
    if not c:
        return None
    parts = []
    for org in cfg["organisms"]:
        if f"cand_{org}" in c:
            parts.append(f"{org}: {c[f'cand_{org}']} candidatos de superfície"
                         + (f", {c[f'sel_{org}']} sob seleção purificadora"
                            if f"sel_{org}" in c else ""))
    n_ep = c.get("mhc1", 0) + c.get("mhc2", 0)
    if n_ep:
        parts.append(f"{n_ep} epitopos MHC selecionados por cobertura populacional")
    prop = outpath(cfg, "09_physchem", "construct_properties.tsv")
    if prop.exists():
        try:
            df = pd.read_csv(prop, sep="\t")
            row = df.iloc[0].to_dict() if len(df) else {}
            if "length" in row:
                parts.append(f"construto final de {int(row['length'])} aa")
        except Exception:
            pass
    return "; ".join(parts) + "." if parts else None


def resolve_secao_resultados(cfg) -> str | None:
    c = stage_counts(cfg)
    if not c:
        return None
    lines = ["| Organismo | Genes core | Candidatos de superfície | Aprovados na seleção |",
             "|---|---|---|---|"]
    for org in cfg["organisms"]:
        if f"cand_{org}" not in c:
            continue
        lines.append(f"| {cfg['organisms'][org]['label']} | {c.get(f'core_{org}', '—')} "
                     f"| {c[f'cand_{org}']} | {c.get(f'sel_{org}', '—')} |")
    return "\n".join(lines) if len(lines) > 2 else None


RESOLVERS = {
    "n_genomes_total": resolve_n_genomes_total,
    "md_ns": resolve_md_ns,
    "results_summary": resolve_results_summary,
    "secao_resultados": resolve_secao_resultados,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default=str(ROOT / "manuscript" / "manuscript.md"))
    ap.add_argument("--output", default=str(ROOT / "manuscript" / "manuscript_filled.md"))
    ap.add_argument("--strict", action="store_true",
                    help="sai com código 1 se restar algum marcador de dados não resolvido")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    src = Path(args.input)
    if not src.exists():
        raise SystemExit(f"manuscrito não encontrado: {src}")
    text = src.read_text()

    found = MARKER.findall(text)
    log.info("%d marcadores no manuscrito (%d chaves distintas)", len(found), len(set(found)))

    resolved: dict[str, str] = {}
    for key in sorted(set(found)):
        if key in AUTHORED:
            continue
        fn = RESOLVERS.get(key)
        if fn is None:
            continue
        try:
            val = fn(cfg)
        except Exception as exc:
            log.warning("resolvedor de '%s' falhou: %s", key, exc)
            val = None
        if val is not None:
            resolved[key] = val

    out = MARKER.sub(lambda m: resolved.get(m.group(1), m.group(0)), text)
    dst = Path(args.output)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(out)

    remaining = sorted(set(MARKER.findall(out)))
    log.info("escrito %s", dst)

    print("\n── Marcadores resolvidos ─────────────────────────────")
    for k, v in resolved.items():
        preview = v if len(v) <= 90 else v[:87] + "..."
        print(f"  ✓ {k}: {preview}")
    if not resolved:
        print("  (nenhum)")

    print("\n── Marcadores PENDENTES ──────────────────────────────")
    data_pending = []
    for k in remaining:
        if k in AUTHORED:
            print(f"  ✎ {k}: texto de autoria humana — escreva à mão")
        else:
            data_pending.append(k)
            print(f"  ✗ {k}: seria preenchido por {SOURCES.get(k, 'fonte não mapeada')}")
    if not remaining:
        print("  (nenhum)")
    print()

    if args.strict and data_pending:
        raise SystemExit(f"{len(data_pending)} marcadores de dados sem fonte disponível")


if __name__ == "__main__":
    main()
