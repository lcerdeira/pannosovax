# Paper B — Bioinformatics (OUP), Application Note

**Título de trabalho:** *PanNosoVax Studio: a desktop application for reproducible
multi-epitope vaccine design without the command line*

**Tipo:** Application Note. **Limite:** ~2 páginas (≈1300–2000 palavras), **1 figura**,
referências enxutas.

## O argumento em uma frase
Desenho de vacina reversa é acessível a quem programa; este app entrega o mesmo pipeline
reproduzível (Snakemake) a pesquisadores clínicos que não usam linha de comando — sem abrir
mão de retomada, proveniência e reprodutibilidade.

## Por que isso é publicável (e não só "mais uma GUI")
- **A lacuna é real**: ferramentas de vacinologia reversa existentes ou são servidores web
  com limites de submissão, ou são pipelines CLI. O público que mais precisa (clínicos,
  epidemiologistas, saúde pública) é justamente o que menos usa terminal.
- **GUI sobre Snakemake**, não em vez dele: a interface escreve a config e dirige o DAG, então
  reprodutibilidade e proveniência ficam intactas — é o que diferencia de um wrapper qualquer.
- **Honestidade sobre tempo**: o desenho leva horas (IEDB). O app trata isso como requisito de
  produto (checkpoint, fechar e voltar, estimativa na tela), não como detalhe.

## Estrutura (formato da revista)
| Seção | Conteúdo |
|---|---|
| **Abstract** | *Summary* (o que é, para quem) · *Availability and implementation* (URL, licença, binários mac/win, DOI) · *Contact* · *Supplementary information* |
| **1 Introduction** | a lacuna: vacinologia reversa exige CLI; público-alvo clínico; o que existe e por que não serve |
| **2 Implementation** | arquitetura em 3 camadas (UI → Snakemake → scripts); etapas cobertas; retomada; empacotamento mac/win via CI; ambiente fixado |
| **3 Results / Use case** | o caso PanNosoVax (→ cita Paper A): 3 patógenos, 302 antígenos, construto de 788 aa; dataset demo roda em ~10 min |
| **4 Conclusion** | limitações (etapas manuais: ColabFold, C-ImmSim, AllerTOP) e roadmap |

**Figura única (a que temos direito):** painel com (a) captura da interface mostrando as
etapas e o progresso, e (b) o DAG do Snakemake — comunica "fácil de usar" e "reprodutível"
na mesma imagem. *Não* reutilizar nenhuma figura do Paper A.

## Requisitos DUROS da revista (bloqueiam submissão)
- [ ] Software **público e funcional** — repositório aberto + licença OSI
- [ ] **Binários** macOS e Windows (via GitHub Actions)
- [ ] Documentação de instalação e uso
- [ ] **Dataset de demonstração** que roda rápido
- [ ] DOI arquivado (Zenodo)
- [ ] Testes automatizados mínimos
- [ ] Suporte declarado a plataformas

## Estado do software (2026-07-23)
- [x] Backbone Snakemake (`workflow/Snakefile`, 19 regras) — validado, DAG constrói
- [x] Backend FastAPI (`app/backend.py`) dirigindo o Snakemake — fatia vertical funcionando
- [x] UI web (`app/ui/index.html`) lendo estado real das etapas
- [ ] Empacotamento (pywebview + PyInstaller) e CI mac/win
- [ ] Ambiente fixado / dataset demo / testes / documentação

**Consequência:** este manuscrito só é submetível quando os itens acima estiverem fechados.
Escrever o texto antes disso é desperdício — a revista rejeita software indisponível.
