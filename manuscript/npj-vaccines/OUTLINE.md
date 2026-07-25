# Paper A — npj Vaccines

**Título de trabalho:** *A pan-nosocomial multi-epitope immunogen against three WHO-priority
respiratory pathogens, designed with commensal-microbiome safety screening and Brazil-weighted
HLA coverage*

**Tipo:** Article (npj Vaccines). **Idioma:** inglês.

## O argumento em uma frase
Um único imunógeno quimérico pode, in silico, cobrir *K. pneumoniae*, *A. baumannii* e
*S. pneumoniae* usando apenas antígenos core de superfície (não capsulares → independente de
sorotipo), sob seleção purificadora, filtrado contra o proteoma humano **e o microbioma
comensal**, e otimizado para cobertura HLA **ponderada pela população brasileira**.

## As duas teses de novidade (o que faz o paper ser publicável)
1. **Triagem negativa contra o microbioma comensal** — ausente na literatura de vacinas
   multi-epitopo. Evita disbiose; empiricamente, derruba 77% dos epitopos de *S. pneumoniae*
   por homologia com estreptococos comensais (Fig. 2).
2. **Cobertura de epitopos ponderada por frequências HLA brasileiras**, não só mundiais.

## Estrutura e figuras
| Seção | Conteúdo | Figuras |
|---|---|---|
| Abstract | 150 palavras estruturado | — |
| Introduction | AMR; o trio; falência das estratégias capsulares; substituição de sorotipo; a lacuna | — |
| Results | funil; validação por controles; segurança por camada; cobertura HLA; construto; físico-química; convergência estrutural exploratória | F1–F5, F8 |
| Methods | 2.1–2.10 (ver SPLIT.md) como método científico | — |
| Discussion | limitações (in silico; falta validação úmida); o achado ABC-SBP; próximos passos | — |
| Data/Code availability | genomas (acessos), pipeline como PanNosoVax Studio (→ Paper B) | — |

**Figuras já prontas** (`results/report/figures/`): F1 funil, F2 segurança, F3 cobertura,
F4 construto, F5 convergência estrutural, F8 pan-nosocomial (candidata a Graphical Abstract).

## Pendências específicas desta revista
- [ ] Tradução PT→EN
- [ ] Abstract estruturado (Background/Methods/Results/Conclusions) ≤ 150 palavras
- [ ] Reporting summary / checklists do portfólio Nature
- [ ] Confirmar que nenhuma figura de software entra aqui
- [ ] Declaração de disponibilidade citando o Paper B (mesmo que "in preparation")
- [ ] Números finais dos TSVs (não inventar — usar fill_manuscript.py)
