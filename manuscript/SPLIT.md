# PanNosoVax — separação em dois manuscritos

Decisão (2026-07-23): o trabalho gera **duas contribuições distintas** que devem ser
publicadas separadamente. Enfiar a ferramenta no paper da vacina dilui o argumento
biológico e arrisca desk reject por escopo em revista de vacinas.

| | Paper A — a vacina | Paper B — a ferramenta |
|---|---|---|
| **Alvo** | **npj Vaccines** | **Bioinformatics (OUP) — Application Note** |
| **Pasta** | `manuscript/npj-vaccines/` | `manuscript/bioinformatics-appnote/` |
| **Pergunta** | É possível um imunógeno único, sorotipo-independente e seguro contra a flora, cobrindo os três patógenos? | Como a comunidade roda esse desenho sem linha de comando? |
| **Unidade** | achado biológico + construto | software maduro, testado, disponível |
| **Limite** | ~sem limite rígido (research article) | ~1300–2000 palavras + 1 figura |
| **Revisão** | meses | meses; exige o software **já no ar** |
| **Fator de impacto** | npj Vaccines ~7–9 (confirmar ano) | Bioinformatics ~4–6 (confirmar ano) |

**Sequenciamento (consequência real):** o npj Vaccines pode ser submetido assim que
traduzido/revisado. O Application Note **só pode ser submetido depois que o app estiver
empacotado, testado e público** (binários mac/win + repositório + DOI Zenodo). Não são
paralelos na submissão, embora o desenvolvimento seja.

**Citação cruzada:** o Paper A cita o Paper B na seção de disponibilidade de código/dados
("o pipeline está disponível como PanNosoVax Studio, ref B"); o Paper B cita o Paper A como
o caso de uso que o motivou e valida.

## Alocação de conteúdo (do manuscrito atual)

Do `manuscript/manuscript.md` (rascunho único, em português):

**→ Paper A (npj Vaccines):**
- Introdução inteira (AMR, o trio, falha das estratégias capsulares, substituição de sorotipo)
- Métodos 2.1–2.8 (amostragem, pangenoma, surfaceome, filtro evolutivo, epitopos,
  **as 4 camadas de segurança**, cobertura máxima, construto) — descritos como MÉTODO
  científico, não como features de software
- Métodos 2.9–2.10 (estrutura/MD, imunossimulação/manufatura) — como validação in silico
- Resultados 3.1–3.4 e as figuras F1–F8
- As DUAS teses de novidade: triagem vs microbioma comensal + cobertura ponderada pelo Brasil

**→ Paper B (Bioinformatics AppNote):**
- A arquitetura GUI-sobre-Snakemake (novo — não está no manuscrito atual)
- Empacotamento multiplataforma, retomada, o problema do tempo de execução e como a UI lida
- Disponibilidade, dataset de demonstração, instalação em um clique
- NÃO repete a biologia: cita o Paper A

**Regra de ouro para não sobrepor:** o Paper A pode mencionar "implementamos um pipeline
reproduzível (Snakemake)" em uma frase, mas todo o *como usar / instalar / a interface* é
exclusivo do Paper B. Nenhuma figura é compartilhada.

## Estado
- [x] Estrutura de pastas criada
- [x] Esqueletos dos dois papers (outline + seções)
- [x] **Paper A escrito em inglês** (`npj-vaccines/manuscript_en.md`, ~2000 palavras)
- [x] Métodos reescritos descrevendo o pipeline real (sem enquadramento de "ferramenta")
- [x] **Bibliografia via CrossRef** (21 entradas, `scripts/report/build_bibliography.py`);
      3 entradas sinalizadas para conferência humana
- [x] **Paper B escrito** (`bioinformatics-appnote/manuscript_en.md`, 931 palavras de corpo)
- [ ] Preencher números do Paper A (automático, quando o re-run terminar)
- [ ] Paper A: autoria/afiliações, Reporting Summary do portfólio Nature
- [ ] Paper B: release com binários, DOI Zenodo, Figura 1, demo end-to-end
- [ ] Verificar as 3 referências sinalizadas (snakemake, RS09, DeepTMHMM preprint)
