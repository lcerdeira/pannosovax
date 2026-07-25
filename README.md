# PanNosoVax — POC de vacina multi-epitopo pan-nosocomial *in silico*

Desenho, do zero e com custo próximo de zero, de um imunógeno quimérico contra os três
principais agentes de pneumonia nosocomial e adquirida na comunidade com resistência
crescente: **Klebsiella pneumoniae species complex (KpSC)**, **Acinetobacter baumannii**
e **Streptococcus pneumoniae**.

## A tese (o "fora da caixa")

A vacinologia reversa clássica produz um candidato por patógeno, cada um preso ao
sorotipo capsular. Isso é caro, lento e — no caso do pneumococo — já falhou uma vez por
substituição de sorotipo. Este projeto inverte três premissas:

1. **Um construto, três patógenos.** Em vez de tratar KpSC, *A. baumannii* e
   *S. pneumoniae* como três programas, tratamos o **nicho** (pulmão + trato
   respiratório em paciente hospitalizado) como o alvo. O construto é uma proteína
   quimérica única com blocos de epitopos dos três organismos, mais um bloco de
   epitopos **compartilhados por homologia estrutural** (não por identidade de
   sequência) — a parte realmente nova.
2. **Independência de cápsula.** Só entram antígenos proteicos de superfície do
   *core genome*, presentes em ≥95% dos genomas e independentes de sorotipo K/O/capsular.
   Isso é a resposta direta ao problema de substituição de sorotipo das PCVs.
3. **Escape caro por construção.** Antígenos são priorizados por seleção purificadora
   (dN/dS < 1) e por acoplamento a fitness/virulência. Um mutante de escape paga um
   preço fisiológico — o que não é verdade para a maioria dos alvos de vacinologia reversa.

Dois diferenciais adicionais de segurança e equidade, ambos gratuitos:

- **Triagem negativa ampliada:** além do proteoma humano, filtramos contra o
  **microbioma comensal respiratório e intestinal** (evita disbiose vacinal) e contra
  painéis de alérgenos e de mimetismo autoimune.
- **Cobertura populacional ponderada por HLA brasileiro/latino-americano**, além da
  cobertura global. A maioria dos artigos de multi-epitopo otimiza implicitamente para
  HLA europeu.

## Estágios do pipeline

| # | Estágio | Ferramenta principal | Custo |
|---|---------|----------------------|-------|
| 01 | Aquisição de genomas (NCBI Datasets, BV-BRC) | `datasets`, `ncbi-genome-download` | 0 |
| 02 | Core genome por ortologia recíproca vs referência (BLAST, O(N)) | BLASTp | 0 |
| 03 | Surfaceome: localização subcelular, peptídeo-sinal, hélices TM | **DeepLocPro**, SignalP-6, DeepTMHMM | 0 |
| 04 | Conservação, dN/dS e ligação a fitness | MAFFT + IQ-TREE + HyPhy | 0 |
| 05 | Predição de epitopos B, CD4 (MHC-II) e CD8 (MHC-I) | IEDB API, NetMHCpan, BepiPred, DiscoTope | 0 |
| 06 | Triagem de segurança (humano, microbioma, alérgenos, toxinas) | BLASTp, AllerTOP, ToxinPred | 0 |
| 07 | Cobertura populacional HLA (global + BR/LatAm) | IEDB Population Coverage | 0 |
| 08 | Montagem do construto quimérico + adjuvante | script próprio | 0 |
| 09 | Físico-química, estrutura 3D, refino | ProtParam, ColabFold, GalaxyRefine | 0 |
| 10 | Docking com TLR2/TLR4 e dinâmica molecular | HADDOCK/ClusPro + GROMACS | 0 |
| 11 | Simulação imunológica (prime-boost) | C-ImmSim | 0 |
| 12 | Otimização de códons e clonagem *in silico* | script próprio + SnapGene Viewer | 0 |

## Documentação

| Documento | Conteúdo |
|-----------|----------|
| [docs/STRATEGY.md](docs/STRATEGY.md) | **Leia primeiro.** Sequência POC → preprint → grant, proteção de prioridade, e a decisão sobre o formato da plataforma |
| [docs/POSITIONING.md](docs/POSITIONING.md) | Relação com vaccines.watch: o que eles cobrem, o que deixaram aberto, o que citar |
| [docs/PLATFORM.md](docs/PLATFORM.md) | Especificação do índice de antígenos e do site estático |
| [docs/PROTOCOL.md](docs/PROTOCOL.md) | Execução passo a passo, incluindo as intervenções manuais inevitáveis |
| [docs/COST.md](docs/COST.md) | Custo real, incluindo o que só aparece na validação experimental |

Manuscrito em [manuscript/manuscript.md](manuscript/manuscript.md); versão Word gerada em
`manuscript/PanNosoVax_manuscrito.docx` via:

```bash
pandoc manuscript/_submission.md -o manuscript/PanNosoVax_manuscrito.docx
```

## Como rodar

```bash
mamba env create -f environment.yml
conda activate pannosovax
```

Cada script em `scripts/` roda isolado (`python scripts/03_surfaceome_filter.py --help`).
O re-run completo em escala é orquestrado por `hpc/run_full_pipeline.sh` (Slurm).

> O `workflow/Snakefile` está sendo **reconciliado** com os scripts que efetivamente
> produzem os resultados — ver "Estado atual". Até lá, o caminho canônico é o driver em `hpc/`.

## Dois entregáveis / dois papers

Ver [manuscript/SPLIT.md](manuscript/SPLIT.md). O desenho da vacina vai para *npj Vaccines*;
o app (**PanNosoVax Studio**, em `app/`) para uma *Application Note* em *Bioinformatics*.

## Estado atual (desenvolvimento ativo)

- **Re-run em escala cheia (400/300/356 genomas) em andamento** no HPC. Os números e figuras
  mudam até isso assentar — a análise inicial rodou em escala de piloto (173 genomas).
- **Surfaceome sendo validado de verdade** (DeepLocPro + SignalP-6 + DeepTMHMM); a versão
  anterior dependia de anotação e estava marcada `needs_rigorous_check`.
- **Snakefile em reconciliação** com os scripts reais.

Os números do manuscrito são preenchidos por `scripts/report/fill_manuscript.py` a partir dos
TSVs em `results/`. **Nenhum número é inventado** — arquivos sem dado não viram figura.

## Licença

GNU General Public License v3.0 — ver [LICENSE](LICENSE).

## Aviso

Este é um estudo computacional. Nada aqui substitui validação *in vitro*/*in vivo*.
O manuscrito declara isso explicitamente na seção de limitações, como deve.
