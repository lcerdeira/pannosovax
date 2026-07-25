# Protocolo de execução — do zero ao manuscrito

Ordem de execução, dependências externas e pontos onde é preciso intervenção manual.

## Passo 0 — Ambiente e dados externos

```bash
mamba env create -f environment.yml
conda activate pannosovax
mkdir -p data/external
```

Baixe uma vez (os únicos downloads pesados do projeto):

| Arquivo | Origem | Tamanho |
|---------|--------|---------|
| `UP000005640_human.fasta` | UniProt, proteoma humano de referência | ~35 MB |
| `bakta_db` | `bakta_db download --output data/external/bakta_db --type full` | ~65 GB (ou `--type light`, ~2 GB) |
| `commensal_respiratory.fasta` | UniProt: *N. lactamica*, *S. mitis*, *S. oralis*, *Moraxella catarrhalis*, *Corynebacterium* spp. | ~50 MB |
| `commensal_gut.fasta` | HMP / UHGG, representantes de espécie | ~200 MB |
| `deg_essential.tsv` | DEG (Database of Essential Genes) + telas TraDIS publicadas | ~5 MB |

Se o disco for limitado, `bakta --type light` e Prokka funcionam; a anotação fica um pouco
menos rica, mas o pipeline não quebra.

## Passo 1 — Rodar o pipeline

```bash
snakemake -c8 --use-conda -n     # dry-run: confira o grafo antes
snakemake -c8 --use-conda
```

Tempo estimado num notebook comum: **18–30 h**, dominado pela anotação (passo 02) e pelas
chamadas ao IEDB (passo 05). Ambos são retomáveis — o Snakemake não refaz o que já existe.

## Passo 2 — Intervenções manuais inevitáveis

Quatro serviços não têm API estável e exigem submissão pelo navegador. O pipeline para,
escreve o arquivo de entrada e diz exatamente o que fazer.

### 2A · Alergenicidade e toxicidade (estágio 06, camada D)

O estágio 06 escreve `results/06_safety/{org}_{classe}_allergen_queue.fasta`.

1. Submeta em AllerTOP v2.0 (`ddg-pharmfac.net/AllerTOP`) e ToxinPred 3.0.
2. Salve os resultados como `{org}_{classe}_allergen_results.tsv`, no mesmo diretório,
   com as colunas `is_allergen` e `is_toxin` (booleanas), uma linha por peptídeo, **na
   mesma ordem da fila**.
3. Rode de novo: `snakemake -c8 --use-conda --forcerun safety`.

Enquanto os resultados não existem, a camada D fica marcada `PENDENTE` e os epitopos
passam pelas camadas A–C. Não confunda isso com aprovação.

### 2B · Estrutura 3D (estágio 10)

Se `colabfold_batch` não estiver instalado localmente, o script escreve um arquivo
`*.ABSENT.txt` com o comando pronto para colar no ColabFold no Google Colab. Rode lá
(gratuito, ~15 min para 200 aa), baixe o PDB de melhor ranking e coloque em
`results/10_structure/construct_refined.pdb`.

Depois: refino em GalaxyRefine, validação em SAVES v6 (ERRAT, Verify3D, PROCHECK).
Guarde os relatórios em `results/10_structure/validation/`.

### 2C · Docking (estágio 10b)

O script baixa os PDBs dos receptores e gera os arquivos de restrição e a configuração
do HADDOCK 2.4. Submeta no servidor web (conta acadêmica gratuita), baixe o resultado e
coloque em `results/10_docking/haddock_run/`. Rode de novo para a análise.

Alternativa sem espera de fila: ClusPro, também gratuito.

### 2D · Imunossimulação (estágio 11b)

C-ImmSim é exclusivamente web. O script escreve `results/11_immunosim/SUBMISSAO.md` com
os parâmetros já calculados (timesteps das três doses, dado o intervalo de 28 dias do
config). Submeta, baixe o CSV e coloque no diretório indicado.

## Passo 3 — Figuras e manuscrito

```bash
python scripts/report/make_figures.py
python scripts/report/fill_manuscript.py
```

`fill_manuscript.py` resolve os marcadores `⟨PENDENTE:xx⟩` a partir dos TSVs em
`results/` e escreve `manuscript/manuscript_filled.md`. Ele **nunca inventa valor**: se o
TSV de origem não existir, o marcador permanece e o script relata qual arquivo o
preencheria. Ao final, revise a lista de marcadores restantes — ela é a sua lista de
pendências para submissão.

## Passo 4 — Antes de submeter

- [ ] Nenhum `⟨PENDENTE:⟩` restante no manuscrito
- [ ] Camada D da segurança concluída (não `PENDENTE`)
- [ ] Três réplicas de MD convergidas — verifique platô de RMSD, não só a média
- [ ] Tabela S1 com todos os acessos de genoma
- [ ] Código público com DOI (Zenodo, gratuito)
- [ ] Declaração de uso de IA conforme a política do periódico
- [ ] Limitações revisadas: elas devem estar mais fortes do que o instinto pede

## Ordem de importância, se o tempo for curto

Se você precisa de um manuscrito submetível rápido, os estágios 01–09 e 12 já constituem
um artigo completo. Docking e MD (10–11) fortalecem, mas são complementares — vários
periódicos aceitam o desenho sem eles, desde que as limitações digam isso claramente.
O que **não** pode faltar: o filtro de conservação (05), as quatro camadas de segurança
(06) e a otimização de cobertura (07). São eles que sustentam a novidade do trabalho.
