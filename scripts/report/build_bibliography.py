#!/usr/bin/env python3
"""Monta o .bib do manuscrito buscando metadados REAIS no CrossRef.

Por que não escrever as entradas à mão: metadados de citação (volume, páginas, ano,
DOI) são exatamente o tipo de detalhe que se lembra errado. Uma referência com DOI
trocado é um erro de integridade, não de formatação. Aqui cada entrada vem da API do
CrossRef e carrega o DOI que a própria API devolveu.

O que este script NÃO faz: decidir se o artigo encontrado é o certo. A busca é por
texto e pode trazer o trabalho errado (por exemplo, uma revisão que cita a ferramenta
em vez do artigo original). Por isso toda entrada sai marcada com `verify = {...}` e
o relatório final lista o que precisa de conferência humana.

Uso:
    python scripts/report/build_bibliography.py            # gera o .bib
    python scripts/report/build_bibliography.py --check    # só relata, não escreve
"""
from __future__ import annotations
import argparse, json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "manuscript" / "npj-vaccines" / "references.bib"
API = "https://api.crossref.org/works"
MAILTO = "louise.cerdeira@lshtm.ac.uk"   # "polite pool" do CrossRef

# (chave bibtex, consulta, o que a referência sustenta no manuscrito)
REFS = [
    ("murray2022amr", "Global burden of bacterial antimicrobial resistance in 2019 systematic analysis",
     "carga global de AMR (Introdução)"),
    ("who2024bppl", "WHO bacterial priority pathogens list 2024 antibiotic resistant bacteria",
     "lista de patógenos prioritários da OMS"),
    ("weinberger2011serotype", "Serotype replacement in disease after pneumococcal vaccination",
     "substituição de sorotipo pós-PCV"),
    ("hansen2022deeptmhmm", "DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks",
     "topologia de membrana (Métodos 4.3)"),
    ("teufel2022signalp6", "SignalP 6.0 predicts all five types of signal peptides using protein language models",
     "peptídeo sinal (Métodos 4.3)"),
    ("moreno2024deeplocpro", "DeepLocPro prediction of subcellular localization of prokaryotic proteins",
     "localização subcelular procariótica (Métodos 4.3)"),
    ("vita2019iedb", "The Immune Epitope Database IEDB 2019 update",
     "predição de epitopos (Métodos 4.4)"),
    ("reynisson2020netmhcpan", "NetMHCpan-4.1 and NetMHCIIpan-4.0 improved predictions of MHC antigen presentation",
     "predição de ligação a MHC"),
    ("jespersen2017bepipred", "BepiPred-2.0 sequence-based B-cell epitope prediction using conformational epitopes",
     "predição de epitopos de célula B"),
    ("bui2006popcoverage", "Predicting population coverage of T-cell epitope-based diagnostics and vaccines",
     "cobertura populacional (Métodos 4.6)"),
    ("zhang2005tmalign", "TM-align a protein structure alignment algorithm based on the TM-score",
     "sobreposição estrutural (Resultados 2.6)"),
    ("xu2010tmscore", "How significant is a protein structure similarity with TM-score of 0.5",
     "limiar de TM-score para mesma dobra"),
    ("varadi2024alphafolddb", "AlphaFold Protein Structure Database 2024 providing structure coverage",
     "modelos estruturais pré-computados"),
    ("camacho2009blast", "BLAST+ architecture and applications",
     "ortologia e triagem de homologia"),
    ("molder2021snakemake", "Sustainable data analysis with Snakemake",
     "reprodutibilidade do pipeline"),
    ("guruprasad1990instability", "Correlation between stability of a protein and its dipeptide composition instability index",
     "índice de instabilidade (Métodos 4.7)"),
    ("sharp1987cai", "The codon adaptation index a measure of directional synonymous codon usage bias",
     "CAI (otimização de códons)"),
    ("tobias1991nend", "The N-end rule in bacteria",
     "regra do N-terminal (desenho do construto)"),
    ("alexander1994padre", "Development of high potency universal DR-restricted helper epitopes PADRE",
     "epitopo T-helper promíscuo PADRE"),
    ("shanmugam2012rs09", "Synthetic Toll-like receptor 4 agonist RS09 adjuvant",
     "adjuvante molecular RS09"),
]


def crossref(query: str, rows: int = 3) -> list[dict]:
    url = (f"{API}?query.bibliographic={urllib.parse.quote(query)}&rows={rows}"
           f"&select=DOI,title,author,container-title,issued,volume,page,type&mailto={MAILTO}")
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return json.load(r)["message"]["items"]
    except Exception as exc:  # noqa: BLE001
        print(f"    ! CrossRef falhou: {str(exc)[:80]}", file=sys.stderr)
        return []


# Itens que o CrossRef devolve e que NÃO são o artigo original: relatórios de
# revisão, recomendações, erratas, notícias. Eram a maior fonte de erro na v1.
JUNK = ("peer review report", "faculty opinions", "recommendation of", "news brief",
        "correction to", "erratum", "corrigendum", "comment on", "reply to",
        "editorial", "in this issue", "author response")


def score(it: dict, query: str, key: str = "") -> tuple[float, list[str]]:
    """Pontua um candidato do CrossRef e devolve (score, motivos de suspeita).

    A chave bibtex codifica o ano esperado (ex.: murray2022amr). Divergência de ano
    é o sinal mais confiável de que veio o artigo errado — tipicamente um comentário
    ou resumo publicado depois, que tem título quase idêntico ao original.
    """
    title = (it.get("title") or [""])[0]
    tl = title.lower()
    journal = (it.get("container-title") or [""])[0]
    year = it.get("issued", {}).get("date-parts", [[None]])[0][0]
    qwords = {w.lower() for w in re.findall(r"[A-Za-z]{4,}", query)}
    twords = {w.lower() for w in re.findall(r"[A-Za-z]{4,}", title)}
    overlap = len(qwords & twords) / max(len(qwords), 1)

    s, flags = overlap, []
    if any(j in tl for j in JUNK):
        s -= 1.0
        flags.append("parece revisão/errata/notícia, não o artigo original")
    if it.get("type") == "posted-content" or "biorxiv" in journal.lower():
        s -= 0.25
        flags.append("preprint")
    if not journal:
        s -= 0.2
        flags.append("sem revista")
    if year is None:
        s -= 0.3
        flags.append("sem ano")
    if overlap < 0.35:
        flags.append(f"baixa sobreposição de título ({overlap:.0%})")
    m = re.search(r"(19|20)\d{2}", key)
    if m and year:
        expected = int(m.group(0))
        if abs(year - expected) > 1:
            s -= 0.5
            flags.append(f"ano diverge: esperado ~{expected}, veio {year}")
    return s, flags


def pick_best(items: list[dict], query: str, key: str = "") -> tuple[dict, list[str]]:
    """Escolhe o candidato de maior score; devolve seus motivos de suspeita."""
    ranked = sorted(((score(it, query, key), it) for it in items),
                    key=lambda x: x[0][0], reverse=True)
    (best_score, flags), best = ranked[0]
    if best_score < 0.35:
        flags = list(dict.fromkeys(flags + ["score baixo — provavelmente o artigo errado"]))
    return best, flags


def to_bibtex(key: str, it: dict, purpose: str, query: str) -> str:
    authors = " and ".join(
        f"{a.get('family','')}, {a.get('given','')}".strip(", ")
        for a in it.get("author", [])[:12] if a.get("family"))
    title = (it.get("title") or [""])[0].replace("{", "").replace("}", "")
    journal = (it.get("container-title") or [""])[0]
    year = it.get("issued", {}).get("date-parts", [[None]])[0][0]
    entry = [f"@article{{{key},"]
    if authors:
        entry.append(f"  author  = {{{authors}}},")
    entry.append(f"  title   = {{{{{title}}}}},")
    if journal:
        entry.append(f"  journal = {{{journal}}},")
    if year:
        entry.append(f"  year    = {{{year}}},")
    if it.get("volume"):
        entry.append(f"  volume  = {{{it['volume']}}},")
    if it.get("page"):
        entry.append(f"  pages   = {{{it['page']}}},")
    entry.append(f"  doi     = {{{it.get('DOI','')}}},")
    entry.append(f"  note    = {{Sustenta: {purpose}. VERIFICAR: entrada obtida por busca "
                 f"automática no CrossRef (consulta: \"{query}\") — conferir se é o artigo "
                 f"original e não uma revisão/preprint superado.}}")
    entry.append("}")
    return "\n".join(entry)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="só relata, não escreve o .bib")
    args = ap.parse_args()

    print(f"Buscando {len(REFS)} referências no CrossRef...\n")
    entries, suspect = [], []
    for key, query, purpose in REFS:
        items = crossref(query)
        if not items:
            suspect.append((key, "NENHUM resultado"))
            print(f"  ✗ {key}: nada encontrado")
            continue
        it, flags = pick_best(items, query, key)
        title = (it.get("title") or ["?"])[0]
        year = it.get("issued", {}).get("date-parts", [[None]])[0][0]
        journal = (it.get("container-title") or ["(sem revista)"])[0]
        if flags:
            suspect.append((key, "; ".join(flags)))
        mark = "⚠" if flags else "✓"
        print(f"  {mark} {key}: {title[:60]} ({year}, {journal[:28]})")
        entries.append(to_bibtex(key, it, purpose, query))
        time.sleep(0.4)  # cortesia com a API

    if args.check:
        print("\n(--check: nada escrito)")
    else:
        header = (
            "% Bibliografia do PanNosoVax — Paper A (npj Vaccines)\n"
            "%\n"
            "% Entradas obtidas automaticamente do CrossRef por\n"
            "% scripts/report/build_bibliography.py. Os metadados vêm da API, mas a\n"
            "% ESCOLHA do artigo é por busca textual e PODE ESTAR ERRADA.\n"
            "% Toda entrada tem nota 'VERIFICAR'. Conferir antes de submeter.\n\n")
        existing = OUT.read_text() if OUT.exists() else ""
        keep = "\n".join(b for b in re.findall(r"@\w+\{[^@]*?\n\}", existing)
                         if "david2025vaccineswatch" in b)
        OUT.write_text(header + (keep + "\n\n" if keep else "") + "\n\n".join(entries) + "\n")
        print(f"\n✓ escrito {OUT.relative_to(ROOT)} ({len(entries)} entradas)")

    if suspect:
        print("\n── PRECISA DE CONFERÊNCIA HUMANA ────────────────────")
        for k, why in suspect:
            print(f"  ⚠ {k}: {why}")


if __name__ == "__main__":
    main()
