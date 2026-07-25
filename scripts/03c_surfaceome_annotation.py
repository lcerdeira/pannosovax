#!/usr/bin/env python3
"""
Estágio 03c — primeiro passe do surfaceome por anotação funcional.

**Este passe NÃO substitui PSORTb/SignalP/DeepTMHMM.** É uma triagem preliminar para
quando essas ferramentas ainda não estão instaladas, e serve para estimar a ordem de
grandeza do funil e priorizar quais proteínas submeter à predição rigorosa.

Método: classificação por família funcional a partir da anotação PGAP do RefSeq, com
listas separadas de inclusão e de exclusão. A exclusão importa tanto quanto a inclusão —
proteínas ribossomais, chaperonas e enzimas metabólicas citoplasmáticas são o principal
ruído das triagens baseadas em texto, e várias delas aparecem em superfície na
literatura como "moonlighting", o que não as torna bons alvos de anticorpo.

Toda proteína aprovada aqui deve passar por PSORTb + SignalP + DeepTMHMM antes de entrar
no construto. O campo `needs_rigorous_check` marca isso explicitamente.

Uso:
    python scripts/03c_surfaceome_annotation.py --all
"""
from __future__ import annotations

import argparse
import re

import pandas as pd

from common import GRAM, get_logger, load_config, outpath, write_table

log = get_logger("03c_surface")

# Famílias com localização de superfície bem estabelecida
INCLUDE = {
    # "porin" sozinho captura aquaporina/aquagliceroporina, que são canais de membrana
    # interna e não porinas de membrana externa — daí a exclusão explícita abaixo.
    "outer_membrane": r"outer membrane|\bOmpA\b|\bOmpC\b|\bOmpF\b|\bOmpK[0-9]*\b|"
                      r"(?<!aqua)(?<!aquaglycero)porin\b|\bOmp[0-9]+\b",
    "tonb_receptor": r"TonB-dependent|siderophore receptor|ferric.*receptor|heme receptor",
    "adhesin": r"adhesin|invasin|intimin|autotransporter|hemagglutinin",
    "pilus_fimbria": r"pilus|pili|fimbri|curli|flagell",
    "lipoprotein": r"lipoprotein",
    "abc_substrate_binding": r"ABC transporter substrate-binding|periplasmic binding|"
                             r"solute-binding",
    "choline_binding": r"choline-binding|PspC|CbpA",
    "sortase_lpxtg": r"sortase|LPXTG|cell wall anchor|cell wall surface anchor",
    "secreted_toxin": r"cytolysin|hemolysin|leukocidin|pneumolysin",
    "capsule_surface_enzyme": r"neuraminidase|sialidase|hyaluronidase|peptidoglycan hydrolase",
    "efflux_outer": r"TolC|efflux.*outer membrane",
    "secretion_system": r"type (IV|VI|II|III) secretion",
}

# Ruído clássico de triagem por texto — citoplasmáticas, ribossomais, metabólicas
EXCLUDE = (
    r"ribosomal|ribosome|tRNA|rRNA|DNA (polymerase|gyrase|ligase|helicase|topoisomerase)|"
    r"RNA polymerase|transcription|translation initiation|elongation factor|"
    r"chaperon|GroEL|DnaK|HtpG|heat shock|cold shock|"
    r"dehydrogenase|synthetase|synthase|reductase|kinase|phosphatase|isomerase|"
    r"transferase|hydratase|carboxylase|aldolase|racemase|deaminase|"
    r"transposase|integrase|recombinase|restriction|methyltransferase|"
    r"hypothetical protein|DUF[0-9]+|"
    r"cytoplasmic|intracellular|"
    r"inner membrane|cytochrome|ATP synthase|NADH|"
    r"aquaporin|aquaglyceroporin|MIP/aquaporin"
)

# Gram-positivos não possuem membrana externa: a família é incoerente para eles e
# qualquer atribuição desse tipo é falso positivo de anotação.
GRAM_POSITIVE_FORBIDDEN = {"outer_membrane", "tonb_receptor", "efflux_outer"}


def classify(product: str) -> tuple[str | None, bool]:
    """Devolve (família, excluido)."""
    p = (product or "")
    if re.search(EXCLUDE, p, flags=re.I):
        return None, True
    for fam, pat in INCLUDE.items():
        if re.search(pat, p, flags=re.I):
            return fam, False
    return None, False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=["kpsc", "abau", "spneu"])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    sf = cfg["surfaceome"]
    targets = list(cfg["organisms"]) if args.all else [args.organism]

    summary = []
    for org in targets:
        d = pd.read_csv(outpath(cfg, "02_pangenome", f"{org}_presence.tsv"), sep="\t")
        core = d[d["is_core"]].copy()

        fam_excl = core["product"].map(classify)
        core["family"] = [f for f, _ in fam_excl]
        core["excluded_by_annotation"] = [e for _, e in fam_excl]
        core["pass_length"] = core["length"].between(sf["min_length"], sf["max_length"])

        # coerência com a parede celular do organismo
        if GRAM[org] == "positive":
            incoerente = core["family"].isin(GRAM_POSITIVE_FORBIDDEN)
            if incoerente.any():
                log.warning("%s (Gram-positivo): %d proteínas descartadas por família "
                            "incompatível com ausência de membrana externa",
                            org, int(incoerente.sum()))
            core.loc[incoerente, "family"] = None

        core["surface_candidate"] = (core["family"].notna()
                                     & ~core["excluded_by_annotation"]
                                     & core["pass_length"])
        core["needs_rigorous_check"] = core["surface_candidate"]   # todas precisam
        core["gram"] = GRAM[org]

        cand = core[core["surface_candidate"]].sort_values(
            ["presence_fraction", "mean_identity"], ascending=False)

        base = outpath(cfg, "03_surfaceome", f"{org}_annotation_pass.tsv")
        write_table(core, base, log)
        write_table(cand, base.with_name(f"{org}_candidates_annotation.tsv"), log)

        byfam = cand["family"].value_counts().to_dict()
        log.info("%s: %d core -> %d candidatos de superfície (%.1f%%)",
                 org, len(core), len(cand), 100 * len(cand) / max(len(core), 1))
        for fam, n in sorted(byfam.items(), key=lambda x: -x[1]):
            log.info("     %-26s %d", fam, n)

        summary.append({"organism": org, "core": len(core),
                        "excluded": int(core["excluded_by_annotation"].sum()),
                        "surface_candidates": len(cand), **byfam})

    write_table(pd.DataFrame(summary).fillna(0),
                outpath(cfg, "03_surfaceome", "annotation_pass_summary.tsv"), log)
    log.warning("passe preliminar — submeta os candidatos a PSORTb/SignalP/DeepTMHMM "
                "antes de qualquer conclusão")


if __name__ == "__main__":
    main()
