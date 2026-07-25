# Plataforma companheira: rastreamento de diversidade de antígenos proteicos

## A lacuna, em uma frase

vaccines.watch responde "qual formulação polissacarídica cobre a população circulante hoje?".
**Ninguém responde** "quais antígenos proteicos de superfície permanecem conservados, sob
seleção purificadora, e imunologicamente acessíveis — e como isso varia por região e ao longo
do tempo?"

Essa segunda pergunta é a que importa para qualquer vacina que pretenda escapar da esteira de
substituição de sorotipo. E ela é continuamente móvel: um antígeno conservado em 2015 pode não
ser em 2030. Um pipeline de execução única (que é o que nossa POC é hoje) não captura isso.
Uma plataforma captura.

## Nome

Evite o sufixo `.watch`. É a convenção de marca do Centre for Genomic Pathogen Surveillance
— pathogen.watch, vaccines.watch, microreact — e adotá-la leria como carona de marca,
justamente o que o [posicionamento](POSITIONING.md) recomenda evitar.

Sugestões que sinalizam complementaridade sem ambiguidade de origem: **CoreAntigen Atlas**,
**ProteoVax**, **PanAntigen**. Preferência: *CoreAntigen Atlas* — descreve exatamente o que é.

## Arquitetura mínima viável

O ponto crítico de custo: **não reimplementar a ingestão de genomas.** Eles já fazem isso a
cada 4 horas e disponibilizam tudo via Pathogenwatch e ENA/SRA. Consumir em vez de duplicar
reduz o custo de operação em uma ordem de grandeza e reforça a complementaridade em vez da
competição.

```
  Pathogenwatch / ENA  ──▶  ingestão incremental (só metadados + assemblies)
           │
           ▼
  ┌──────────────────────────────────────────────────┐
  │  núcleo: nosso pipeline 02→07 já implementado    │
  │  pangenoma → surfaceome → dN/dS → epitopo →      │
  │  segurança → cobertura HLA                       │
  └──────────────────────────────────────────────────┘
           │
           ▼
  índice de antígenos versionado (Parquet + DuckDB)
           │
           ▼
  API REST  ──▶  front-end estático (sem servidor de aplicação)
```

Custo de operação estimado: **US$ 5–20/mês** (armazenamento de objetos + página estática).
DuckDB sobre Parquet elimina a necessidade de banco gerenciado; um recomputo semanal em
runner gratuito de CI dá conta do volume.

## As quatro telas que justificam a existência da plataforma

Cada uma responde a algo que vaccines.watch estruturalmente não pode responder, porque não
carrega dados proteicos.

**1 · Índice de conservação por antígeno.** Para cada proteína core de superfície: fração de
isolados que a contêm, identidade média, dN/dS por sítio. Ordenável, filtrável por espécie.
A tela que um desenvolvedor de vacina abre primeiro.

**2 · Deriva temporal.** Conservação do antígeno por ano de coleta. Um antígeno cuja curva
está caindo é um alvo que está sendo perdido — sinal de alerta precoce, invisível numa análise
de instantâneo único. É o análogo proteico do monitoramento pós-rollout deles.

**3 · Estratificação geográfica com honestidade sobre o viés.** Conservação por região, com o
denominador sempre visível. Se a Nigéria contribui com 4 genomas, a tela mostra "4", não uma
porcentagem tranquilizadora. O preprint deles mostra que 79,2% dos genomas vêm de países de
alta renda e 1,5% de países de baixa renda; qualquer plataforma que exiba médias globais sem
esse denominador está mentindo por omissão.

**4 · Cobertura HLA por população.** A curva de cobertura acumulada do nosso estágio 07,
interativa: escolha a população, veja quantos epitopos são necessários para atingir 90%.
Esta é a tela sem equivalente em lugar nenhum, e é o argumento de equidade tornado tangível.

## Faseamento realista

| Fase | Entrega | Esforço | Depende de |
|------|---------|---------|-----------|
| 0 | Artigo do PanNosoVax (POC, três patógenos, execução única) | ~2 meses | pipeline atual, já pronto |
| 1 | Índice estático publicado: TSV + Parquet + DOI no Zenodo | ~2 semanas | fase 0 |
| 2 | Front-end estático com as telas 1 e 4 | ~1 mês | fase 1 |
| 3 | Recomputo automático e tela 2 (deriva temporal) | ~2 meses | fase 2 |
| 4 | Artigo de recurso da plataforma | ~1 mês | fase 3 |

A fase 1 sozinha já é citável e já marca presença — e custa duas semanas, não dois anos.
Comece por ela.

## Por que isto é defensável a longo prazo

Uma plataforma vive de manutenção contínua, e manutenção depende de financiamento. O nosso
diferencial não é ter mais genomas — nunca teremos — mas ter **a camada analítica que eles
escolheram não construir**, ancorada num pipeline de desenho de imunógeno de ponta a ponta.

Isso também é o que torna o conjunto financiável: CARB-X, GARDP, FAPESP e Wellcome financiam
especificamente alvos de AMR e equidade em saúde global. Um recurso que documenta antígenos
proteicos independentes de sorotipo, com cobertura HLA explícita para populações
sub-representadas, encaixa nesses editais com precisão incomum.
