#!/usr/bin/env python3
"""
Estágio 11b — imunossimulação do esquema prime-boost (C-ImmSim).

O que a simulação agrega: o construto pode ter epitopos excelentes e ainda assim
falhar no esquema — resposta primária sem memória, ou IgM que nunca troca de classe
para IgG. C-ImmSim modela a dinâmica de células B/T num autômato celular e mostra se
o segundo e o terceiro reforço realmente elevam o platô de IgG e a população de
células B de memória.

Restrição prática: C-ImmSim é serviço web (Kraken/UniMoRe) e **não tem API estável**.
Fingir uma chamada programática aqui produziria números inventados, que é o pior
resultado possível num artigo. Então este script faz duas coisas honestas:

  1. prepara o arquivo de entrada e escreve as instruções exatas de submissão;
  2. se o usuário já baixou o resultado, converte para TSV padronizado.

Parâmetros do esquema vêm de config immunosim (doses, interval_days, simulation_steps).
Nota de calibração: 1 passo de simulação ≈ 8 horas reais; 28 dias entre doses são
84 passos.

Entrada : results/08_construct/construct.fasta
          results/11_immunosim/cimmsim_raw.txt  (baixado manualmente, opcional)
Saída   : results/11_immunosim/cimmsim_input.txt
          results/11_immunosim/SUBMISSAO.md
          results/11_immunosim/immune_response.tsv

Uso:
    python scripts/11b_immune_sim.py --fasta results/08_construct/construct.fasta
    python scripts/11b_immune_sim.py --fasta x.fasta --raw meus_resultados.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from common import get_logger, load_config, outpath, write_table

log = get_logger("11b_immunosim")

COLUMNS = ["day", "IgM", "IgG1", "IgG2", "Ig_total", "B_cells",
           "TH_cells", "TC_cells", "IFNg", "IL2"]

STEP_HOURS = 8          # 1 passo do C-ImmSim ≈ 8 h

# Sinônimos que aparecem nas várias exportações do C-ImmSim.
ALIASES = {
    "day": ["day", "days", "time", "t"],
    "IgM": ["igm", "ig m"],
    "IgG1": ["igg1", "ig g1"],
    "IgG2": ["igg2", "ig g2"],
    "Ig_total": ["ig_total", "igtotal", "ig total", "total ig", "igm+igg"],
    "B_cells": ["b cells", "b_cells", "b-cell", "blym", "b lymphocytes"],
    "TH_cells": ["th cells", "th_cells", "t helper", "thlym", "th"],
    "TC_cells": ["tc cells", "tc_cells", "t cytotoxic", "tclym", "tc", "ctl"],
    "IFNg": ["ifn-g", "ifng", "ifn g", "interferon"],
    "IL2": ["il-2", "il2", "il 2"],
}


def write_input(fasta: Path, dst: Path) -> str:
    recs = list(SeqIO.parse(fasta, "fasta"))
    if not recs:
        raise SystemExit(f"FASTA vazio: {fasta}")
    seq = str(recs[0].seq)
    with open(dst, "w") as fh:
        fh.write(f">{recs[0].id}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return seq


def write_instructions(dst: Path, cfg_i: dict, seq_len: int, inp: Path) -> list[int]:
    doses = cfg_i.get("doses", 3)
    interval = cfg_i.get("interval_days", 28)
    steps = cfg_i.get("simulation_steps", 1050)
    inject_steps = [1 + i * int(round(interval * 24 / STEP_HOURS)) for i in range(doses)]

    dst.write_text(f"""\
# Submissão ao C-ImmSim — instruções

C-ImmSim é serviço web sem API estável. Este passo é **manual e obrigatório**;
nenhum número é gerado automaticamente aqui.

1. Abra: https://kraken.iac.rm.cnr.it/C-IMMSIM/
2. Escolha "Protein antigen" e cole o conteúdo de:
       {inp}
   (construto de {seq_len} aa)
3. Parâmetros (de config/config.yaml, seção immunosim):
       Random seed .................. 12345
       Simulation volume ............ 10 µL
       Simulation steps ............. {steps}   (1 passo ≈ {STEP_HOURS} h ≈ {steps * STEP_HOURS / 24:.0f} dias)
       Number of injections ......... {doses}
       Injection at time steps ...... {', '.join(map(str, inject_steps))}
                                      (intervalo de {interval} dias entre doses)
       Adjuvant ..................... none
         → o adjuvante já está fundido no construto; marcar adjuvante externo
           duplicaria o estímulo e inflaria a resposta artificialmente.
       LPS .......................... no
       HLA ..........................  deixe o padrão (A0101/A0201/B0702/B0801,
                                       DRB1_0101/DRB1_0401)
4. Baixe o resultado tabular e salve como:
       {inp.parent / 'cimmsim_raw.txt'}
5. Rode de novo:  python scripts/11b_immune_sim.py --fasta <construto>
   → gera immune_response.tsv com as colunas padronizadas do pipeline.

O que olhar no resultado:
  * IgG1 deve subir a cada dose e atingir platô mais alto que o da dose anterior;
    se o platô não sobe, não há memória e o esquema prime-boost não está funcionando.
  * IgM alto com IgG baixo = falha de troca de classe (falta ajuda de T CD4).
  * IFN-γ e IL-2 sustentados indicam polarização Th1, desejável contra bactérias
    intracelulares facultativas e para opsonização eficiente.
""")
    return inject_steps


def normalize(raw: Path) -> pd.DataFrame | None:
    """Converte a exportação do C-ImmSim para as colunas padronizadas do pipeline."""
    for sep in ("\t", ",", r"\s+"):
        try:
            df = pd.read_csv(raw, sep=sep, comment="#", engine="python")
            if df.shape[1] > 1:
                break
        except Exception:
            continue
    else:
        return None
    if df.shape[1] < 2:
        return None

    lower = {str(c).strip().lower(): c for c in df.columns}
    out = pd.DataFrame()
    for target, names in ALIASES.items():
        src = next((lower[n] for n in names if n in lower), None)
        if src is None:
            src = next((orig for low, orig in lower.items()
                        if any(low.startswith(n) for n in names)), None)
        out[target] = pd.to_numeric(df[src], errors="coerce") if src else pd.NA
        if src is None:
            log.warning("coluna '%s' não encontrada na exportação — ficará vazia", target)

    if out["day"].isna().all():
        # Sem coluna de tempo: o índice da linha é o passo de simulação.
        out["day"] = [round(i * STEP_HOURS / 24, 2) for i in range(len(out))]
    return out[COLUMNS]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--raw", default=None,
                    help="arquivo baixado do C-ImmSim (padrão: results/11_immunosim/cimmsim_raw.txt)")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fasta = Path(args.fasta)
    if not fasta.exists():
        raise SystemExit(f"construto ausente: {fasta}. Rode o estágio 08.")

    inp = outpath(cfg, "11_immunosim", "cimmsim_input.txt")
    seq = write_input(fasta, inp)
    steps = write_instructions(outpath(cfg, "11_immunosim", "SUBMISSAO.md"),
                               cfg["immunosim"], len(seq), inp)
    log.info("entrada preparada (%d aa); doses nos passos %s", len(seq), steps)

    out_tsv = outpath(cfg, "11_immunosim", "immune_response.tsv")
    raw = Path(args.raw) if args.raw else outpath(cfg, "11_immunosim", "cimmsim_raw.txt")

    if not raw.exists():
        log.warning("resultado do C-ImmSim não encontrado (%s). C-ImmSim é serviço WEB "
                    "sem API — a submissão é manual. Siga results/11_immunosim/SUBMISSAO.md "
                    "e rode este script de novo. Escrevendo %s só com cabeçalho.",
                    raw.name, out_tsv.name)
        write_table(pd.DataFrame(columns=COLUMNS), out_tsv, log)
        return

    df = normalize(raw)
    if df is None or df.empty:
        log.warning("não consegui interpretar %s como tabela do C-ImmSim — "
                    "escrevendo só o cabeçalho", raw.name)
        write_table(pd.DataFrame(columns=COLUMNS), out_tsv, log)
        return

    write_table(df, out_tsv, log)
    igg = pd.to_numeric(df["IgG1"], errors="coerce")
    if igg.notna().any():
        log.info("IgG1 máximo = %.1f no dia %.1f", igg.max(), df.loc[igg.idxmax(), "day"])


if __name__ == "__main__":
    main()
