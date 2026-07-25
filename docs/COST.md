# Custo real do projeto

## Fase computacional (este repositório): **US$ 0 – 120**

| Item | Opção gratuita | Custo |
|------|----------------|-------|
| Genomas | NCBI RefSeq, BV-BRC, PATRIC | 0 |
| Anotação | Bakta / Prokka (local) | 0 |
| Pangenoma | Panaroo | 0 |
| Localização subcelular | PSORTb 3.0 (local), SignalP 6 (licença acadêmica gratuita) | 0 |
| Topologia TM | DeepTMHMM (API pública) | 0 |
| Epitopos | IEDB Analysis Resource (API aberta, sem chave) | 0 |
| Estrutura 3D | ColabFold no Google Colab (tier gratuito) | 0 |
| Refino/validação | GalaxyRefine, SAVES v6, PROCHECK | 0 |
| Docking | HADDOCK 2.4 (acadêmico), ClusPro | 0 |
| Dinâmica molecular | GROMACS em CPU local ou Colab | 0 |
| Imunossimulação | C-ImmSim | 0 |
| Compute | notebook comum resolve; Colab Pro opcional | 0–120 |

O único gasto real é energia elétrica e, opcionalmente, ~US$ 10/mês de Colab Pro se
você quiser rodar as dinâmicas de 100 ns em GPU em vez de esperar dias em CPU.

## O que o dinheiro compra depois (validação úmida — fora deste repositório)

Aqui é onde a honestidade importa. Um artigo *in silico* é publicável sozinho, mas ele
é uma **hipótese**, não uma vacina. Escala de custo para levar adiante:

| Etapa | Custo aproximado (USD) | Prazo |
|-------|------------------------|-------|
| Síntese do gene + clonagem em pET-28a | 300 – 800 | 3–4 semanas |
| Expressão e purificação (IMAC, 10 mg) | 1.500 – 4.000 | 4–6 semanas |
| ELISA com soro de convalescentes | 2.000 – 5.000 | 4 semanas |
| Imunização em camundongo (n=30, 3 grupos) | 8.000 – 20.000 | 3 meses |
| Ensaio de opsonofagocitose (OPKA) | 5.000 – 15.000 | 6 semanas |
| Desafio letal em modelo murino de pneumonia | 20.000 – 60.000 | 4 meses |
| **Subtotal pré-clínico exploratório** | **~40.000 – 100.000** | **~12 meses** |
| Toxicologia GLP + lote GMP fase I | 2 – 8 milhões | 2–3 anos |

**Estratégia de custo baixo para a validação inicial:** peptídeos sintéticos dos
epitopos individuais (US$ 50–150 cada) testados por ELISPOT de IFN-γ contra PBMC de
doadores saudáveis validam a imunogenicidade dos epitopos T por **menos de US$ 3.000**,
sem precisar expressar a proteína inteira. É o teste de maior retorno por dólar e o
próximo passo lógico depois deste artigo.

## Fontes de financiamento compatíveis (Brasil)

- FAPESP Auxílio Regular / Jovem Pesquisador — cobre a fase pré-clínica exploratória
- CNPq Universal — cobre síntese e ELISPOT
- Programa de bolsas PIBIC/mestrado — mão de obra do pipeline computacional
- Editais Fiocruz/Butantan de parceria em vacinas bacterianas
- CARB-X e GARDP — financiam especificamente alvos de resistência antimicrobiana
