#!/usr/bin/env python3
"""Estágio 10 — dobra o construto COMPLETO (788 aa) com ESMFold local (HuggingFace).

Por que local: a API pública do ESMAtlas recusa sequências >400 aa (HTTP 413), e o
ColabFold exige Colab/GPU. O `EsmForProteinFolding` do transformers é PyTorch puro
(não precisa do openfold/CUDA), então roda em CPU no macOS.

Memória é o gargalo real: a atenção do trunk escala ~O(L²) e L=788. Usamos
`set_chunk_size` para trocar velocidade por memória. Se ainda assim estourar,
reduza --chunk.

Saída: results/10_structure/construct_refined.pdb  (pLDDT no B-factor, escala 0-100)

Uso:
    python scripts/fold_construct_esmfold_local.py --chunk 32
"""
from __future__ import annotations
import argparse, time, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import get_logger, ROOT

log = get_logger("esmfold")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=32, help="chunk do trunk (menor = menos RAM)")
    ap.add_argument("--fasta", default=str(ROOT / "results/08_construct/construct.fasta"))
    ap.add_argument("--out", default=str(ROOT / "results/10_structure/construct_refined.pdb"))
    args = ap.parse_args()

    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding
    from Bio import SeqIO

    rec = next(SeqIO.parse(args.fasta, "fasta"))
    seq = str(rec.seq).upper().replace("*", "")
    log.info("sequência: %s (%d aa)", rec.id, len(seq))

    torch.set_grad_enabled(False)
    log.info("carregando facebook/esmfold_v1 (~2.6 GB no primeiro uso)...")
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model.eval()
    model.trunk.set_chunk_size(args.chunk)
    log.info("modelo pronto; chunk=%d; dobrando em CPU (pode levar dezenas de minutos)...", args.chunk)

    t0 = time.time()
    inputs = tok([seq], return_tensors="pt", add_special_tokens=False)
    out = model(**inputs)
    dt = time.time() - t0
    log.info("fold concluído em %.1f min", dt / 60)

    pdb = model.output_to_pdb(out)[0]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(pdb)

    plddt = out["plddt"][0, :, 1].tolist()      # pLDDT por resíduo (CA)
    mean = sum(plddt) / len(plddt)
    log.info("escrito %s · pLDDT médio=%.1f · >=70 em %d/%d resíduos",
             args.out, mean, sum(1 for p in plddt if p >= 70), len(plddt))


if __name__ == "__main__":
    main()
