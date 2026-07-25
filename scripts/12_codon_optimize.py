#!/usr/bin/env python3
"""
Estágio 12 — otimização de códons, verificação de clonagem e desenho da variante mRNA.

Faz três coisas:

  1. Retrotradução com uso de códons de E. coli K-12, maximizando CAI mas evitando o
     erro clássico de usar sempre o códon mais frequente: isso cria trechos de
     tRNA-depleção e favorece estruturas secundárias no mRNA. Usamos amostragem
     proporcional à frequência, com o códon ótimo forçado apenas nos primeiros
     ~15 códons (onde a eficiência de iniciação realmente importa).

  2. Verifica sítios de restrição do MCS do pET-28a que precisam estar ausentes do
     inserto, além de GC local, repetições e sítios internos de Shine-Dalgarno.

  3. Gera a variante mRNA-LNP: mesma proteína, mas com UTRs otimizados, cauda poli-A,
     e substituição de uridina por N1-metilpseudouridina anotada. É a segunda
     plataforma do artigo — permite comparar custo/velocidade sem redesenhar o antígeno.

Uso:
    python scripts/12_codon_optimize.py --fasta results/08_construct/construct.fasta
"""
from __future__ import annotations

import argparse
import random
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from common import get_logger, load_config, outpath, write_table

log = get_logger("12_codon")

# Uso de códons de E. coli K-12 (frequência relativa dentro de cada família de sinônimos)
ECOLI_USAGE: dict[str, dict[str, float]] = {
    "A": {"GCG": 0.36, "GCC": 0.27, "GCA": 0.21, "GCT": 0.16},
    "R": {"CGC": 0.40, "CGT": 0.38, "CGG": 0.10, "CGA": 0.06, "AGA": 0.04, "AGG": 0.02},
    "N": {"AAC": 0.55, "AAT": 0.45},
    "D": {"GAT": 0.63, "GAC": 0.37},
    "C": {"TGC": 0.55, "TGT": 0.45},
    "Q": {"CAG": 0.65, "CAA": 0.35},
    "E": {"GAA": 0.68, "GAG": 0.32},
    "G": {"GGC": 0.40, "GGT": 0.34, "GGG": 0.15, "GGA": 0.11},
    "H": {"CAT": 0.57, "CAC": 0.43},
    "I": {"ATT": 0.51, "ATC": 0.42, "ATA": 0.07},
    "L": {"CTG": 0.50, "TTA": 0.13, "TTG": 0.13, "CTT": 0.10, "CTC": 0.10, "CTA": 0.04},
    "K": {"AAA": 0.76, "AAG": 0.24},
    "M": {"ATG": 1.00},
    "F": {"TTT": 0.57, "TTC": 0.43},
    "P": {"CCG": 0.52, "CCA": 0.19, "CCT": 0.16, "CCC": 0.13},
    "S": {"AGC": 0.28, "TCT": 0.15, "TCC": 0.15, "AGT": 0.15, "TCG": 0.15, "TCA": 0.12},
    "T": {"ACC": 0.44, "ACG": 0.27, "ACT": 0.17, "ACA": 0.13},
    "W": {"TGG": 1.00},
    "Y": {"TAT": 0.57, "TAC": 0.43},
    "V": {"GTG": 0.37, "GTT": 0.26, "GTC": 0.22, "GTA": 0.15},
    "*": {"TAA": 0.64, "TGA": 0.29, "TAG": 0.07},
}

RESTRICTION = {
    "NdeI": "CATATG", "XhoI": "CTCGAG", "BamHI": "GGATCC", "EcoRI": "GAATTC",
    "HindIII": "AAGCTT", "NcoI": "CCATGG", "SacI": "GAGCTC", "NotI": "GCGGCCGC",
    "SalI": "GTCGAC", "XbaI": "TCTAGA",
}

SHINE_DALGARNO = "AGGAGG"


def optimize(protein: str, seed: int = 7, head_optimal: int = 15,
             alpha: float = 1.0) -> str:
    """
    alpha controla o compromisso entre CAI e diversidade de códons:
      alpha = 0  → uniforme (CAI baixo, sem viés)
      alpha = 1  → proporcional ao uso natural de E. coli
      alpha → ∞  → sempre o códon ótimo (CAI máximo, mas cria tRNA-depleção
                   e estrutura secundária no mRNA)
    """
    rng = random.Random(seed)
    dna = []
    for i, aa in enumerate(protein):
        table = ECOLI_USAGE.get(aa)
        if table is None:
            raise ValueError(f"aminoácido não reconhecido na posição {i + 1}: {aa!r}")
        if i < head_optimal:
            codon = max(table, key=table.get)      # iniciação: usa o códon ótimo
        else:
            codons = list(table)
            weights = [table[c] ** alpha for c in codons]
            codon = rng.choices(codons, weights=weights, k=1)[0]
        dna.append(codon)
    return "".join(dna) + "TAA"


def optimize_to_target(protein: str, target_cai: float, seed: int = 7) -> tuple[str, float]:
    """Aumenta alpha até bater o CAI alvo, mantendo o menor alpha que funciona —
    ou seja, o máximo de diversidade de códons compatível com a meta de expressão."""
    best = None
    for alpha in [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0]:
        dna = optimize(protein, seed=seed, alpha=alpha)
        score = cai(dna[:-3], protein)
        best = (dna, score, alpha)
        if score >= target_cai:
            log.info("   CAI alvo atingido com alpha=%.1f", alpha)
            return dna, score
    log.warning("   CAI alvo %.2f não atingido nem com alpha=10 (melhor: %.4f)",
                target_cai, best[1])
    return best[0], best[1]


def cai(dna: str, protein: str) -> float:
    """Codon Adaptation Index (Sharp & Li 1987), média geométrica dos pesos relativos."""
    import math
    logs = []
    for i, aa in enumerate(protein):
        codon = dna[i * 3:(i + 1) * 3]
        table = ECOLI_USAGE[aa]
        w = table.get(codon, 1e-6) / max(table.values())
        logs.append(math.log(w))
    return round(math.exp(sum(logs) / len(logs)), 4)


def gc_content(s: str) -> float:
    return round((s.count("G") + s.count("C")) / len(s), 4)


def sliding_gc(s: str, window: int = 60) -> list[float]:
    return [gc_content(s[i:i + window]) for i in range(0, len(s) - window + 1, 10)]


def audit(dna: str, avoid: list[str], gc_range: list[float]) -> dict:
    issues = defaultdict(list)
    for name, site in RESTRICTION.items():
        for m in re.finditer(site, dna):
            key = "restricao_proibida" if name in avoid else "restricao_presente"
            issues[key].append(f"{name}@{m.start() + 1}")
    for m in re.finditer(SHINE_DALGARNO, dna):
        issues["shine_dalgarno_interno"].append(f"@{m.start() + 1}")
    for m in re.finditer(r"(.)\1{7,}", dna):
        issues["homopolimero"].append(f"{m.group(1)}x{len(m.group(0))}@{m.start() + 1}")

    win = sliding_gc(dna)
    lo, hi = gc_range
    bad = [i for i, g in enumerate(win) if not (lo <= g <= hi)]
    if bad:
        issues["gc_local_fora_da_faixa"] = [f"janela {i}" for i in bad[:10]]
    return dict(issues)


def repair(dna: str, protein: str, avoid: list[str], max_rounds: int = 40) -> str:
    """Reamostra códons nas regiões problemáticas até limpar os sítios proibidos."""
    for round_i in range(max_rounds):
        problems = audit(dna, avoid, [0.0, 1.0]).get("restricao_proibida", [])
        if not problems:
            return dna
        pos = int(problems[0].split("@")[1]) - 1
        codon_i = pos // 3
        rng = random.Random(1000 + round_i * 31 + codon_i)
        chars = list(dna)
        for ci in range(max(0, codon_i - 1), min(len(protein), codon_i + 3)):
            aa = protein[ci]
            table = ECOLI_USAGE[aa]
            if len(table) == 1:
                continue
            codons, weights = zip(*table.items())
            new = rng.choices(codons, weights=weights, k=1)[0]
            chars[ci * 3:(ci + 1) * 3] = list(new)
        dna = "".join(chars)
    log.warning("não foi possível remover todos os sítios proibidos em %d rodadas", max_rounds)
    return dna


def design_mrna(cds: str) -> dict:
    """Variante mRNA-LNP. UTRs de alta expressão (alfa-globina humana / AES-mtRNR1)."""
    utr5 = "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC"     # + Kozak
    utr3 = ("GCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCA"
            "CCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC")
    polya = "A" * 120
    return {
        "utr5": utr5,
        "cds": cds,
        "utr3": utr3,
        "polya": polya,
        "full": utr5 + cds + utr3 + polya,
        "cap": "Cap1 (m7G ppp Am) — co-transcricional com CleanCap AG",
        "modificacao": "substituição total de uridina por N1-metilpseudouridina (m1Ψ)",
        "nota": "reduz ativação de TLR7/8 e RIG-I, aumenta meia-vida e tradução",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    exp = cfg["expression"]
    fasta = Path(args.fasta) if args.fasta else outpath(cfg, "08_construct", "construct.fasta")
    if not fasta.exists():
        raise SystemExit(f"construto não encontrado: {fasta}")

    rows = []
    for rec in SeqIO.parse(fasta, "fasta"):
        protein = str(rec.seq).upper().replace("*", "")
        dna, _ = optimize_to_target(protein, exp["target_cai_min"])
        dna = repair(dna, protein, exp["restriction_avoid"])
        cds_no_stop = dna[:-3]

        score = cai(cds_no_stop, protein)   # recalculado após o reparo dos sítios
        issues = audit(dna, exp["restriction_avoid"], exp["gc_range"])

        log.info("── %s ──", rec.id)
        log.info("   proteína           %d aa", len(protein))
        log.info("   CDS                %d nt (com códon de parada)", len(dna))
        log.info("   CAI                %.4f (alvo >= %.2f)", score, exp["target_cai_min"])
        log.info("   GC global          %.1f%%", 100 * gc_content(dna))
        if issues:
            for k, v in issues.items():
                log.info("   %-22s %s", k, ", ".join(v[:6]))
        else:
            log.info("   auditoria          sem problemas")

        mrna = design_mrna(cds_no_stop + "TGA")

        with open(outpath(cfg, "12_expression", f"{rec.id}_cds_ecoli.fasta"), "w") as fh:
            fh.write(f">{rec.id}_CDS_ecoli_optimized CAI={score}\n")
            for i in range(0, len(dna), 60):
                fh.write(dna[i:i + 60] + "\n")

        with open(outpath(cfg, "12_expression", f"{rec.id}_mrna.fasta"), "w") as fh:
            fh.write(f">{rec.id}_mRNA 5UTR+CDS+3UTR+polyA len={len(mrna['full'])}\n")
            for i in range(0, len(mrna["full"]), 60):
                fh.write(mrna["full"][i:i + 60] + "\n")

        rows.append({
            "id": rec.id, "protein_aa": len(protein), "cds_nt": len(dna),
            "cai": score, "gc": gc_content(dna),
            "cai_ok": score >= exp["target_cai_min"],
            "forbidden_sites": len(issues.get("restricao_proibida", [])),
            "mrna_nt": len(mrna["full"]),
            "vector": exp["vector"], "host": exp["host"],
        })

    write_table(pd.DataFrame(rows), outpath(cfg, "12_expression", "expression_summary.tsv"), log)


if __name__ == "__main__":
    main()
