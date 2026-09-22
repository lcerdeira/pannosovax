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
    "bloco_compartilhado_resumo": "results/04_shared/shared_structural_epitopes.tsv (04c)",
    "mecanismo_estrutural": "results/04_shared/asymmetry_summary.tsv (04e)",
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
    # Separador de milhar em inglês (vírgula), não o brasileiro (ponto) — mesmo
    # ponto que resolve_n_surface_total tinha, escondido aqui: "1.056" no
    # manuscrito em inglês lê como "um vírgula zero cinco seis", não "mil e
    # cinquenta e seis". n_core_total (abaixo) já usava vírgula; este não.
    return f"{total:,}"


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
        # presence.tsv NÃO é só genes core: é o proteoma de referência inteiro, com
        # uma coluna booleana `is_core`. Contar linhas (n_rows) em vez de somar
        # `is_core` inflou "core_{org}" para o total de genes testados — para kpsc,
        # 5865 em vez dos 3803 que de fato são core. Isso produzia uma tabela de
        # resultados (secao_resultados) cuja soma (11975) contradizia o próprio
        # n_core_total (7841, correto) duas linhas abaixo no manuscrito — o tipo de
        # inconsistência que salta aos olhos de qualquer revisor.
        p = outpath(cfg, "02_pangenome", f"{org}_presence.tsv")
        if p.exists():
            try:
                n_core = int(pd.read_csv(p, sep="\t")["is_core"].fillna(False).astype(bool).sum())
                if n_core:
                    counts[f"core_{org}"] = n_core
            except Exception:
                pass
        for key, sub, fname in [
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
    # Os dois manuscritos ativos (Paper A e Paper B, em manuscript/npj-vaccines/ e
    # manuscript/bioinformatics-appnote/) estão em inglês; manuscript_pt_base.md é um
    # rascunho legado sem commits desde o inicial. Antes esta função escrevia em
    # português direto no marcador do manuscrito em inglês — passou despercebido
    # porque ninguém tinha rodado o preenchimento de ponta a ponta ainda.
    c = stage_counts(cfg)
    if not c:
        return None
    parts = []
    for org in cfg["organisms"]:
        if f"cand_{org}" in c:
            label = cfg["organisms"][org]["label"]
            parts.append(f"*{label}*: {c[f'cand_{org}']} surface candidates"
                         + (f", {c[f'sel_{org}']} under purifying selection"
                            if f"sel_{org}" in c else ""))
    n_ep = c.get("mhc1", 0) + c.get("mhc2", 0)
    if n_ep:
        parts.append(f"{n_ep} MHC epitopes selected by population coverage")
    prop = outpath(cfg, "09_physchem", "construct_properties.tsv")
    if prop.exists():
        try:
            df = pd.read_csv(prop, sep="\t")
            row = df.iloc[0].to_dict() if len(df) else {}
            if "length" in row:
                parts.append(f"a final construct of {int(row['length'])} aa")
        except Exception:
            pass
    return "; ".join(parts) + "." if parts else None


def resolve_secao_resultados(cfg) -> str | None:
    c = stage_counts(cfg)
    if not c:
        return None
    lines = ["| Organism | Core genes | Surface candidates | Passed selection |",
             "|---|---|---|---|"]
    for org in cfg["organisms"]:
        if f"cand_{org}" not in c:
            continue
        lines.append(f"| *{cfg['organisms'][org]['label']}* | {c.get(f'core_{org}', '—')} "
                     f"| {c[f'cand_{org}']} | {c.get(f'sel_{org}', '—')} |")
    return "\n".join(lines) if len(lines) > 2 else None


def resolve_n_core_total(cfg) -> str | None:
    """Total de proteínas core somando os organismos (coluna is_core do presence.tsv)."""
    total = 0
    for org in cfg["organisms"]:
        p = outpath(cfg, "02_pangenome", f"{org}_presence.tsv")
        if not p.exists():
            return None
        try:
            total += int(pd.read_csv(p, sep="\t")["is_core"].fillna(False).astype(bool).sum())
        except Exception:
            return None
    return f"{total:,}" if total else None


def resolve_n_surface_total(cfg) -> str | None:
    """Proteínas em compartimento acessível a anticorpo, somando os organismos.

    Lê `{org}_candidates.tsv` — a saída final e autoritativa do estágio 03 (a coluna
    `candidate` do filtro completo, ver scripts/03_surfaceome_filter.py) — em vez de
    reimplementar a regra de localização aqui. A versão anterior reimplementava:
    aceitava só as localizações "puras" do config (`gram_negative_ok`/
    `gram_positive_ok`) e ignorava tanto a recuperação de Gram-positivo ancorado
    quanto a exclusão por anotação citoplasmática (estágio 03e). Isso a deixou
    desincronizada da correção do surfaceome: devolvia 475 (número do surfaceome
    ANTIGO, keyword-based, citado como limitação no §2.2) quando o valor correto,
    coerente com `results_summary` e `secao_resultados`, é 254. Publicar 475 aqui
    teria reintroduzido silenciosamente, no próprio manuscrito, o mesmo erro que o
    pipeline foi corrigido para não cometer.
    """
    c = stage_counts(cfg)
    total = sum(c.get(f"cand_{org}", 0) for org in cfg["organisms"])
    return str(total) if total else None


def resolve_cobertura_resumo(cfg) -> str | None:
    """Cobertura fenotípica mundo/Brasil do conjunto selecionado, por classe."""
    parts = []
    for klass, label in (("mhc1", "MHC-I"), ("mhc2", "MHC-II")):
        p = outpath(cfg, "07_coverage", f"selected_{klass}.tsv")
        if not p.exists():
            continue
        try:
            d = pd.read_csv(p, sep="\t")
        except Exception:
            continue
        if not len(d):
            continue
        w = 100 * d["set_coverage_world"].max()
        b = 100 * d["set_coverage_brazil"].max()
        parts.append(f"{label}: {len(d)} epitopes, {w:.1f}% coverage (world) "
                     f"and {b:.1f}% (Brazil)")
    return "; ".join(parts) + "." if parts else None


def resolve_construto_resumo(cfg) -> str | None:
    """Tamanho, composição e propriedades do construto final."""
    import json
    parts = []
    blocks = outpath(cfg, "08_construct", "construct_blocks.json")
    if blocks.exists():
        try:
            b = json.loads(blocks.read_text())
            n = b.get("n_epitopes")
            bl = b.get("blocks", {})
            if n:
                parts.append(f"The final construct carries {n} epitopes "
                             f"({len(bl.get('bcell', []))} B-cell, "
                             f"{len(bl.get('mhc2', []))} MHC-II, "
                             f"{len(bl.get('mhc1', []))} MHC-I"
                             + (f", {len(bl['shared'])} structurally shared"
                                if bl.get("shared") else "") + ")")
        except Exception:
            pass
    prop = outpath(cfg, "09_physchem", "construct_properties.tsv")
    if prop.exists():
        try:
            r = pd.read_csv(prop, sep="\t").iloc[0]
            parts.append(f"totalling {int(r['length'])} aa and {r['molecular_weight_kda']:.1f} kDa, "
                         f"with pI {r['theoretical_pi']:.2f}, instability index "
                         f"{r['instability_index']:.1f} and {int(r['n_cysteine'])} free cysteines")
        except Exception:
            pass
    return ", ".join(parts) + "." if parts else None


def resolve_bloco_compartilhado_resumo(cfg) -> str | None:
    """Descreve o bloco compartilhado por nível, a partir da tabela que o 04c grava.

    Não hardcoda contagem nem proteínas: lê `nivel` e `pair_key` de
    `shared_structural_epitopes.tsv`, que é a mesma tabela que o construto usa.
    """
    p = outpath(cfg, "04_shared", "shared_structural_epitopes.tsv")
    if not p.exists():
        return None
    try:
        d = pd.read_csv(p, sep="\t")
    except Exception:
        return None
    if not len(d) or "nivel" not in d.columns:
        return None

    NAMES = {"kpsc": "*K. pneumoniae*", "abau": "*A. baumannii*", "spneu": "*S. pneumoniae*"}
    def pair_label(pk: str) -> str:
        return " and ".join(NAMES.get(o, o) for o in str(pk).split("|"))

    parts = []
    three = d[d["nivel"] == "tres_patogenos"]
    if len(three):
        parts.append(f"{len(three)} epitope(s) fall in a region structurally shared by "
                     f"all three pathogens")
    pair = d[d["nivel"] == "par_cruza_gram"]
    if len(pair):
        by_pair = pair["pair_key"].value_counts()
        detail = "; ".join(f"{n} anchored in {pair_label(pk)}" for pk, n in by_pair.items())
        parts.append(f"{len(pair)} additional epitope(s) are shared between exactly one "
                     f"Gram-negative pathogen and *S. pneumoniae* ({detail})")
    if not parts:
        return None
    return "; ".join(parts) + f"; none are shared between the two Gram-negatives alone. " \
        f"All {len(d)} cross a Gram-positive/Gram-negative boundary."


def resolve_mecanismo_estrutural(cfg) -> str | None:
    """Fração de janelas estruturais dominadas por uma classe funcional, por par.

    Fonte: results/04_shared/asymmetry_summary.tsv, escrito por
    scripts/04e_analyze_asymmetry.py — não recalculado nem digitado aqui.
    """
    p = outpath(cfg, "04_shared", "asymmetry_summary.tsv")
    if not p.exists():
        return None
    try:
        d = pd.read_csv(p, sep="\t")
    except Exception:
        return None
    if not len(d):
        return None
    NAMES = {"kpsc": "*K. pneumoniae*", "abau": "*A. baumannii*", "spneu": "*S. pneumoniae*"}
    # scripts/04e_analyze_asymmetry.py classifies in Portuguese (matches the rest of the
    # pipeline's comments/logs); the manuscript is in English, so translate here rather
    # than in the analysis script.
    CLASS_EN = {"ABC substrate-binding": "ABC-transporter substrate-binding proteins",
               "porina/OM": "outer-membrane porins",
               "TonB/sideróforo": "TonB-dependent siderophore receptors",
               "adesina/pilus": "adhesins/pilus proteins",
               "lipoproteína": "lipoproteins",
               "peptidase/hidrolase": "peptidases/hydrolases",
               "outra": "other folds", "nenhuma classe repetida": "no repeated class"}
    parts = []
    for r in d.itertuples():
        cls = CLASS_EN.get(r.classe_dominante, r.classe_dominante)
        parts.append(f"{NAMES.get(r.org_a, r.org_a)}–{NAMES.get(r.org_b, r.org_b)}: "
                     f"{r.n_classe_dominante}/{r.n_janelas} windows ({100*r.frac_classe_dominante:.0f}%) "
                     f"are {cls}")
    return "; ".join(parts) + "."


RESOLVERS = {
    "n_genomes_total": resolve_n_genomes_total,
    "md_ns": resolve_md_ns,
    "results_summary": resolve_results_summary,
    "secao_resultados": resolve_secao_resultados,
    "n_core_total": resolve_n_core_total,
    "n_surface_total": resolve_n_surface_total,
    "cobertura_resumo": resolve_cobertura_resumo,
    "construto_resumo": resolve_construto_resumo,
    "bloco_compartilhado_resumo": resolve_bloco_compartilhado_resumo,
    "mecanismo_estrutural": resolve_mecanismo_estrutural,
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
