#!/usr/bin/env python3
"""
Estágio 04b v2 — epitopos compartilhados por CORRESPONDÊNCIA ESTRUTURAL (não sequência).

*** EXPLORATÓRIO (mais ambicioso que o 04b v1). ***

A v1 extraía janelas por interseção de sequência idêntica — o que é autocontraditório:
proteínas com a mesma dobra e sequência divergente (Seq_ID ~20%) não têm k-mers idênticos,
então a v1 sempre devolvia vazio. A v2 usa a superposição do TMalign: em cada par
inter-organismo com TM-score >= limiar, lê o alinhamento estrutural (linha de match onde
":" = par de resíduos a d < 5 Å) e procura JANELAS CONTÍGUAS onde ambos os resíduos estão
(i) estruturalmente sobrepostos e (ii) com pLDDT alto (proxy de região bem definida/exposta).

Para cada janela emite os DOIS peptídeos (um por organismo): sequências diferentes ocupando
a mesma posição no espaço. A hipótese é que um anticorpo contra essa superfície comum possa
ter reatividade cruzada. É especulativa e vai para bloco SEPARADO e ROTULADO do construto.

Entrada : results/04_structures/{org}/*.pdb  (modelos AlphaFold; ver fetch_structures_afdb.py)
Saída   : results/04_shared/shared_structural_epitopes_v2.tsv
          results/04_shared/shared_regions_v2.tsv   (regiões que cobrem >2 organismos)

Uso:
    python scripts/04b_structural_crossmatch_v2.py --window 9 --min-plddt 70 --tm 0.5
"""
from __future__ import annotations
import argparse, glob, itertools, os, re, subprocess
from pathlib import Path

import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, load_config, ROOT

log = get_logger("04b_v2")

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def ca_plddt(pdb: Path) -> list[float]:
    """pLDDT (B-factor do CA) por resíduo, na ordem da cadeia."""
    vals = []
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    vals.append(0.0)
    return vals


def tmalign(a: Path, b: Path):
    """Roda TMalign; devolve (tm_ref, seq1, match, seq2) ou None."""
    p = subprocess.run(["TMalign", str(a), str(b)], capture_output=True, text=True)
    if p.returncode != 0:
        return None
    tms = [float(x) for x in re.findall(r"TM-score=\s*([0-9.]+)", p.stdout)]
    tm = max(tms) if tms else 0.0
    lines = p.stdout.splitlines()
    idx = next((i for i, l in enumerate(lines) if l.startswith('(":"')), None)
    if idx is None or idx + 3 >= len(lines):
        return None
    seq1, match, seq2 = lines[idx + 1], lines[idx + 2], lines[idx + 3]
    n = min(len(seq1), len(match), len(seq2))
    return tm, seq1[:n], match[:n], seq2[:n]


def shared_windows(aln, plddt1, plddt2, window, min_plddt):
    """Janelas de colunas ':' contíguas com pLDDT>=min nos dois. Devolve lista de dicts."""
    tm, seq1, match, seq2 = aln
    # mapeia coluna -> índice de resíduo (0-based) em cada cadeia
    i1 = i2 = -1
    cols = []  # (col, res1, res2, close_and_confident)
    for c in range(len(match)):
        c1 = seq1[c] != "-"
        c2 = seq2[c] != "-"
        if c1: i1 += 1
        if c2: i2 += 1
        good = (match[c] == ":") and c1 and c2 \
            and i1 < len(plddt1) and i2 < len(plddt2) \
            and plddt1[i1] >= min_plddt and plddt2[i2] >= min_plddt
        cols.append((seq1[c], seq2[c], i1, i2, good))
    # corridas de good
    out = []
    k = 0
    while k < len(cols):
        if cols[k][4]:
            j = k
            while j < len(cols) and cols[j][4]:
                j += 1
            run = cols[k:j]
            if len(run) >= window:
                pep1 = "".join(r[0] for r in run)
                pep2 = "".join(r[1] for r in run)
                seqid = sum(1 for r in run if r[0] == r[1]) / len(run)
                mp = (sum(plddt1[r[2]] for r in run) + sum(plddt2[r[3]] for r in run)) / (2 * len(run))
                out.append({"len": len(run), "pep1": pep1, "pep2": pep2,
                            "res1_start": run[0][2] + 1, "res2_start": run[0][3] + 1,
                            "window_seqid": round(seqid, 2), "mean_plddt": round(mp, 1)})
            k = j
        else:
            k += 1
    return out


def org_of(pdb_path: str) -> str:
    return Path(pdb_path).parent.name


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=9)
    ap.add_argument("--min-plddt", type=float, default=70.0)
    ap.add_argument("--tm", type=float, default=0.5)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    load_config(args.config)

    orgs = ["kpsc", "abau", "spneu"]
    models = {o: sorted(glob.glob(str(ROOT / f"results/04_structures/{o}/*.pdb"))) for o in orgs}
    plddt = {p: ca_plddt(Path(p)) for o in orgs for p in models[o]}
    prod = {}
    for o in orgs:
        d = pd.read_csv(ROOT / f"results/03_surfaceome/{o}_candidates.tsv", sep="\t")
        prod.update(dict(zip(d["protein_id"], d["product"].astype(str))))

    log.warning("=== EXPLORATÓRIO v2 === epitopos por sobreposição estrutural; bloco separado "
                "e rotulado. Reatividade cruzada de epitopo linear NÃO está validada.")

    rows = []
    for oa, ob in itertools.combinations(orgs, 2):
        npairs = 0
        for a in models[oa]:
            for b in models[ob]:
                aln = tmalign(Path(a), Path(b))
                if not aln or aln[0] < args.tm:
                    continue
                wins = shared_windows(aln, plddt[a], plddt[b], args.window, args.min_plddt)
                if wins:
                    npairs += 1
                pa, pb = Path(a).stem, Path(b).stem
                for w in wins:
                    rows.append({
                        "tm_score": round(aln[0], 3),
                        "org_a": oa, "prot_a": pa, "pep_a": w["pep1"], "res_a_start": w["res1_start"],
                        "org_b": ob, "prot_b": pb, "pep_b": w["pep2"], "res_b_start": w["res2_start"],
                        "window_len": w["len"], "window_seqid": w["window_seqid"],
                        "mean_plddt": w["mean_plddt"],
                        "product_a": prod.get(pa, "")[:50], "product_b": prod.get(pb, "")[:50],
                        "note": "EXPLORATORIO_v2: janela estruturalmente sobreposta (d<5A, pLDDT>=%d)" % int(args.min_plddt),
                    })
        log.info("%s x %s: %d pares com janela estrutural compartilhada", oa, ob, npairs)

    df = pd.DataFrame(rows)
    out = ROOT / "results/04_shared/shared_structural_epitopes_v2.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    log.info("escrito %s (%d janelas compartilhadas)", out.name, len(df))

    # regiões que cobrem >2 organismos: mesma proteína/janela aparecendo com parceiros de 2 orgs
    if not df.empty:
        cover = {}
        for _, r in df.iterrows():
            for org, prot, pep, start in [(r.org_a, r.prot_a, r.pep_a, r.res_a_start),
                                          (r.org_b, r.prot_b, r.pep_b, r.res_b_start)]:
                key = (prot, int(start))
                d = cover.setdefault(key, {"organism": org, "protein": prot, "peptide": pep,
                                           "res_start": int(start), "partners": set()})
                d["partners"].add(r.org_b if org == r.org_a else r.org_a)
        multi = [{"organism": v["organism"], "protein": v["protein"], "peptide": v["peptide"],
                  "res_start": v["res_start"], "n_partner_orgs": len(v["partners"]),
                  "partner_orgs": "|".join(sorted(v["partners"]))}
                 for v in cover.values() if len(v["partners"]) >= 2]
        mdf = pd.DataFrame(multi).sort_values("n_partner_orgs", ascending=False) if multi else pd.DataFrame()
        mout = ROOT / "results/04_shared/shared_regions_v2.tsv"
        mdf.to_csv(mout, sep="\t", index=False)
        log.info("regiões cobrindo >=3 organismos (2 parceiros): %d -> %s",
                 len(mdf), mout.name)


if __name__ == "__main__":
    main()
