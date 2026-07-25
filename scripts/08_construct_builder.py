#!/usr/bin/env python3
"""
Estágio 08 — montagem do construto quimérico multi-epitopo pan-nosocomial.

Arquitetura (N → C):

  [6xHis] - [Adjuvante RS09] -EAAAK- [PADRE] -GPGPG- [Epitopos B ...KK...]
          -GPGPG- [Epitopos MHC-II ...GPGPG...] - [Epitopos MHC-I ...AAY...]
          - [Bloco compartilhado] - [6xHis]

Escolha dos linkers, que não é arbitrária:
  EAAAK  — hélice alfa rígida; isola o adjuvante do resto para que ele mantenha a
           conformação que engaja TLR4 sem interferência estérica dos epitopos.
  GPGPG  — quebra junções imunogênicas espúrias e favorece processamento por
           MHC-II; é o linker canônico para epitopos de célula T helper.
  AAY    — sítio preferencial de clivagem do proteassomo em mamíferos; garante que
           cada epitopo MHC-I seja liberado como unidade correta.
  KK     — sítio de clivagem de catepsina B; separa epitopos B preservando a
           acessibilidade do anticorpo.

O ponto novo do projeto está no **bloco compartilhado**: epitopos que ocorrem em
região estruturalmente equivalente nos três organismos (identificados no estágio 04
por sobreposição estrutural, não por identidade de sequência). Se a resposta a esse
bloco for protetora, um único imunógeno cobre os três patógenos com uma fração da
carga peptídica.

Uso:
    python scripts/08_construct_builder.py
    python scripts/08_construct_builder.py --demo   # roda com epitopos de exemplo
"""
from __future__ import annotations

import argparse
import json

import pandas as pd

from common import get_logger, load_config, outpath, write_table

log = get_logger("08_construct")

# Adjuvantes moleculares disponíveis (sequências publicadas, domínio público)
ADJUVANTS = {
    "RS09": ("APPHALS", "agonista sintético de TLR4, 7 aa — mínimo custo de síntese"),
    "PADRE": ("AKFVAAWTLKAAA", "epitopo T-helper promíscuo pan-HLA-DR"),
    "FliC_D1": (
        "MAQVINTNSLSLLTQNNLNKSQSALGTAIERLSSGLRINSAKDDAAGQAIANRFTANIKGLTQASRNANDGISIAQTTEGALNEINNNLQRVRELAVQSAN",
        "domínio D1 da flagelina de S. typhimurium, agonista de TLR5",
    ),
    "HBD2": ("GIGDPVTCLKSGAICHPVFCPRRYKQIGTCGLPGTKCCKKP",
             "beta-defensina 2 humana, quimioatraente de células dendríticas via CCR6"),
}


def assemble(blocks: dict[str, list[str]], cfg_c: dict) -> tuple[str, list[dict]]:
    """Monta a sequência e devolve o mapa de coordenadas de cada elemento."""
    linkers = cfg_c["linkers"]
    parts: list[tuple[str, str, str]] = []      # (tipo, rótulo, sequência)

    if cfg_c.get("his_tag"):
        # Cauda N-terminal derivada do pET-28a: MGSSHHHHHHSSGLVPRGSH.
        # Começar a proteína diretamente em His seria um erro — pela regra do N-terminal
        # (Tobias et al. 1991) a histidina é desestabilizante em E. coli e derruba a
        # meia-vida para <2 min. A metionina iniciadora do vetor resolve isso e ainda
        # traz o sítio de trombina (LVPRGS) para remover a cauda após a purificação.
        parts.append(("tag", "pET28a-N-His6", "MGSSHHHHHHSSGLVPRGSH"))

    adj_name = cfg_c["adjuvant"]
    adj_seq, adj_note = ADJUVANTS[adj_name]
    parts.append(("adjuvant", adj_name, adj_seq))
    parts.append(("linker", "EAAAK", cfg_c["adjuvant_linker"]))

    if cfg_c.get("include_padre"):
        parts.append(("helper", "PADRE", ADJUVANTS["PADRE"][0]))
        parts.append(("linker", "GPGPG", linkers["mhc2"]))

    block_linker = {"bcell": linkers["bcell"], "mhc2": linkers["mhc2"],
                    "mhc1": linkers["mhc1"], "shared": linkers["mhc2"]}

    for block in cfg_c["block_order"]:
        if block in {"adjuvant", "padre"}:
            continue
        eps = blocks.get(block, [])
        if not eps:
            log.warning("bloco '%s' vazio — pulando", block)
            continue
        lk = block_linker[block]
        for i, ep in enumerate(eps):
            parts.append((block, f"{block}_{i + 1}", ep))
            if i < len(eps) - 1:
                parts.append(("linker", lk, lk))
        parts.append(("linker", lk, lk))         # separador entre blocos

    if parts and parts[-1][0] == "linker":
        parts.pop()
    if cfg_c.get("his_tag"):
        parts.append(("linker", "GPGPG", linkers["mhc2"]))
        parts.append(("tag", "His6-C", "HHHHHH"))

    seq, mapping, pos = "", [], 1
    for kind, label, s in parts:
        seq += s
        mapping.append({"element": kind, "label": label, "sequence": s,
                        "start": pos, "end": pos + len(s) - 1})
        pos += len(s)

    return seq, mapping


def load_blocks(cfg: dict) -> dict[str, list[str]]:
    blocks: dict[str, list[str]] = {"bcell": [], "mhc2": [], "mhc1": [], "shared": []}

    for klass in ("mhc1", "mhc2"):
        p = outpath(cfg, "07_coverage", f"selected_{klass}.tsv")
        if p.exists():
            blocks[klass] = pd.read_csv(p, sep="\t").sort_values("rank")["peptide"].tolist()

    no_cys = cfg.get("construct", {}).get("exclude_free_cysteine")
    bcell = []
    for org in cfg["organisms"]:
        p = outpath(cfg, "06_safety", f"{org}_bcell_safe.tsv")
        if p.exists():
            d = pd.read_csv(p, sep="\t")
            if no_cys:
                d = d[~d["peptide"].str.contains("C")]
            sort_col = "conservation" if "conservation" in d else d.columns[0]
            bcell += d.sort_values(sort_col, ascending=False).head(4)["peptide"].tolist()
    blocks["bcell"] = bcell

    p = outpath(cfg, "04_shared", "shared_structural_epitopes.tsv")
    if p.exists():
        blocks["shared"] = pd.read_csv(p, sep="\t")["peptide"].tolist()

    return blocks


DEMO_BLOCKS = {
    # Epitopos ilustrativos para validar a mecânica de montagem. NÃO são candidatos
    # reais — o pipeline substitui por preditos assim que os estágios 05-07 rodam.
    "bcell": ["DKQYNSTAQNVDKAI", "GNTSDNQKGSNSQVA"],
    "mhc2": ["FVKLLNAKFSYIALS", "IYRLLGTQSKAIVLA"],
    "mhc1": ["KLNDLVEKL", "YLQPRTFLL", "SVFTDLNAA"],
    "shared": ["GTYSAKLNA"],
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--demo", action="store_true",
                    help="monta com epitopos de exemplo (teste da mecânica)")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    blocks = DEMO_BLOCKS if args.demo else load_blocks(cfg)
    if args.demo:
        log.warning("MODO DEMO — epitopos ilustrativos, não candidatos reais")

    total = sum(len(v) for v in blocks.values())
    if total == 0:
        raise SystemExit("nenhum epitopo disponível. Rode os estágios 05-07, "
                         "ou use --demo para testar a montagem.")

    seq, mapping = assemble(blocks, cfg["construct"])

    tag = "demo_" if args.demo else ""
    fasta = outpath(cfg, "08_construct", f"{tag}construct.fasta")
    with open(fasta, "w") as fh:
        fh.write(f">PanNosoVax_v1 len={len(seq)} n_epitopes={total}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")

    write_table(pd.DataFrame(mapping), outpath(cfg, "08_construct", f"{tag}construct_map.tsv"), log)
    with open(outpath(cfg, "08_construct", f"{tag}construct_blocks.json"), "w") as fh:
        json.dump({"blocks": blocks, "length": len(seq), "n_epitopes": total},
                  fh, indent=2, ensure_ascii=False)

    log.info("construto montado: %d aa, %d epitopos, %d elementos",
             len(seq), total, len(mapping))
    log.info("adjuvante: %s (%s)", cfg["construct"]["adjuvant"],
             ADJUVANTS[cfg["construct"]["adjuvant"]][1])
    print(f"\n>PanNosoVax_v1 ({len(seq)} aa)")
    for i in range(0, len(seq), 60):
        print(seq[i:i + 60])


if __name__ == "__main__":
    main()
