# Dataset de demonstração — PanNosoVax

Um conjunto **mínimo e real** para verificar que o pipeline roda numa máquina limpa,
em poucos minutos, sem os dias de compute do run em escala cheia.

## O que é

- **3 genomas reais por organismo** (9 no total), acessos do NCBI listados em `selection/`
- **`config.demo.yaml`** — mesmos limiares do pipeline, mas `n_genomes: 3` e saída isolada
  em `results_demo/` (não toca nos `results/` reais)

Acessos escolhidos (ver `selection/{org}_selected.tsv`):

| organismo | genomas |
|---|---|
| KpSC | GCF_020525545.1, GCF_051122745.1, GCF_041217655.1 |
| *A. baumannii* | GCF_051294685.1, GCF_024139075.1, GCF_052216905.1 |
| *S. pneumoniae* | GCF_900692555.1, GCF_900692575.1, GCF_900693085.1 |

## Como rodar

```bash
conda activate pannosovax
bash demo/run_demo.sh
```

Isso executa, em `results_demo/`, as etapas **autossuficientes e rápidas**:
download dos 9 proteomas → core genome (BLAST). Leva ~2–3 min e prova que a aquisição e a
definição de core genome funcionam ponta a ponta.

## Escopo (honesto)

As etapas seguintes (surfaceome, IEDB, segurança, cobertura, construto) ainda usam caminhos
fixos em `results/`; o **demo end-to-end completo** é liberado junto com a reconciliação do
`workflow/Snakefile` (ver README principal, "Estado atual"), que unifica a saída por `outdir`.
Até lá, este demo cobre a aquisição + core genome, que é o suficiente para o CI e para uma
verificação rápida de instalação.

O app **PanNosoVax Studio** (`app/`) usará este mesmo dataset como "projeto de exemplo".
