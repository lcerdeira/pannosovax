---
title: "Um imunógeno quimérico pan-nosocomial desenhado por vacinologia reversa
        estrutural contra Klebsiella pneumoniae species complex, Acinetobacter
        baumannii e Streptococcus pneumoniae: desenho in silico independente de
        sorotipo com otimização de cobertura HLA"
short_title: "Vacina multi-epitopo pan-nosocomial in silico"
keywords: [vacinologia reversa, multi-epitopo, resistência antimicrobiana, ESKAPE,
           imunoinformática, cobertura populacional HLA, independente de sorotipo]
target_journals:
  - Frontiers in Immunology (Vaccines and Molecular Therapeutics)
  - Scientific Reports
  - npj Vaccines
  - Vaccine
  - PLOS ONE
---

> **Estado do manuscrito.** Estrutura, métodos, racional e discussão estão completos.
> Marcadores no formato PENDENTE (entre colchetes angulares) indicam números que serão preenchidos automaticamente por
> `scripts/report/fill_manuscript.py` quando o pipeline rodar de ponta a ponta.
> Nenhum resultado foi inventado ou estimado.

---

## Resumo

**Contexto.** *Klebsiella pneumoniae* species complex (KpSC), *Acinetobacter baumannii*
e *Streptococcus pneumoniae* concentram grande parte da carga global de pneumonia
bacteriana grave e figuram entre as prioridades máximas da lista de patógenos da
Organização Mundial da Saúde para pesquisa de novos antimicrobianos. Nenhum dos dois
primeiros possui vacina licenciada. O terceiro possui, mas suas vacinas conjugadas são
dependentes de sorotipo capsular e vêm sofrendo erosão de eficácia por substituição de
sorotipo. O denominador comum — pneumonia em hospedeiro hospitalizado ou vulnerável —
sugere que a unidade natural de intervenção talvez não seja a espécie, e sim o nicho.

**Objetivo.** Desenhar, inteiramente *in silico* e a custo desprezível, um imunógeno
quimérico multi-epitopo único, independente de sorotipo, dirigido simultaneamente aos
três patógenos, e avaliar sua plausibilidade estrutural, imunológica e de manufatura.

**Métodos.** Partimos de 5.789 genomas completos ou em nível de cromossomo disponíveis no RefSeq
(4.290 KpSC, 1.143 *A. baumannii*, 356 *S. pneumoniae*), dos quais 1.056 foram
selecionados por amostragem estratificada por região geográfica para reduzir o viés de
super-representação de clones de surto e de países de sequenciamento intensivo. Definimos o core genome de cada organismo por
Panaroo (presença ≥95%) e restringimos ao surfaceome por PSORTb, SignalP 6 e DeepTMHMM.
Priorizamos proteínas sob seleção purificadora (HyPhy FEL, dN/dS < 1) e associadas a
essencialidade ou virulência. Epitopos de células B lineares, T CD4 (MHC-II) e T CD8
(MHC-I) foram preditos via IEDB e filtrados por conservação ≥95% entre isolados. Aplicamos
uma triagem negativa de quatro camadas — proteoma humano, mimetismo autoimune por 7-mer
exato, **proteomas de comensais respiratórios e intestinais**, e alergenicidade/toxicidade.
A seleção final do conjunto de epitopos foi tratada como problema de cobertura máxima e
resolvida por algoritmo guloso com garantia (1−1/e) sobre cobertura fenotípica ponderada
por frequências alélicas mundiais **e brasileiras**. O construto foi montado com adjuvante
molecular RS09 e linkers específicos por classe (EAAAK, GPGPG, AAY, KK), modelado por
ColabFold, refinado, acoplado a TLR4/TLR2 e submetido a dinâmica molecular de
⟨PENDENTE:md_ns⟩ ns em triplicata, além de imunossimulação de esquema prime-boost.

**Resultados.** ⟨PENDENTE:results_summary⟩

**Conclusões.** ⟨PENDENTE:conclusion_summary⟩ O desenho é apresentado explicitamente como
hipótese testável; propomos um caminho de validação experimental escalonado cujo primeiro
passo custa menos de US$ 3.000.

---

## 1. Introdução

A resistência antimicrobiana deixou de ser uma projeção. Estimativas de carga global
atribuem milhões de mortes anuais a infecções bacterianas resistentes, e três organismos
aparecem de forma persistente no topo das listas de prioridade: *Klebsiella pneumoniae*,
*Acinetobacter baumannii* e *Streptococcus pneumoniae*. Os dois primeiros são
Gram-negativos nosocomiais com repertório crescente de carbapenemases; o terceiro é um
Gram-positivo de circulação comunitária cuja resistência a betalactâmicos e macrolídeos
avança de forma consistente.

A resposta vacinal a esse trio é desigual e, em todos os casos, insatisfatória.

Para **KpSC**, não há vacina licenciada. Os candidatos mais avançados são bioconjugados
de polissacarídeo O, cuja cobertura esbarra na diversidade de sorotipos O e K e na
capacidade do organismo de trocar loci capsulares por recombinação.

Para ***A. baumannii***, também não há vacina licenciada, apesar de mais de uma década
de candidatos promissores em modelo murino — OmpA, Ata, Bap, vesículas de membrana
externa. A translação esbarra na heterogeneidade entre linhagens e na dificuldade de
padronizar correlatos de proteção.

Para ***S. pneumoniae***, existem vacinas conjugadas eficazes, e é justamente aí que
está a lição mais instrutiva. As PCV7, PCV13, PCV15 e PCV20 protegem contra os sorotipos
que contêm — e a cada expansão de valência, o nicho ecológico liberado é ocupado por
sorotipos não vacinais. O fenômeno de substituição de sorotipo é a demonstração empírica
de que uma estratégia baseada em cápsula é uma corrida armamentista contra um genoma
naturalmente transformável.

A magnitude desse problema foi recentemente quantificada em escala global. Ao analisar
mais de 100 mil genomas públicos de alta qualidade dos três patógenos, David e
colaboradores [@david2025vaccineswatch] demonstraram, por meio da plataforma
vaccines.watch, que os sorotipos incluídos na PCV13 correspondem a apenas 36,2%
(11.907/32.918) dos genomas públicos de pneumococo, subindo a 87,4% na PCV21 — ao custo
de uma escalada de valência que já atingiu 21 sorotipos. Os mesmos autores documentam a
ausência quase completa de pipeline vacinal para os dois Gram-negativos: um único
candidato baseado em antígeno O para KpSC chegou a ensaio de fase 1/2 (Kleb4V,
NCT04959344) e nenhum candidato está em desenvolvimento clínico ativo para
*A. baumannii*.

Esse trabalho oferece, portanto, o mapa quantitativo do problema que motiva o presente
estudo — e, de forma explícita, restringe-se a **alvos polissacarídicos** (sorotipos
capsulares do pneumococo, antígenos K e O de KpSC, K e OC de *A. baumannii*). O espaço
de antígenos **proteicos** conservados, que é onde uma solução independente de sorotipo
necessariamente reside, permanece sistematicamente inexplorado nessa escala. É esse o
espaço que ocupamos aqui, e entendemos o presente trabalho como complementar — e não
concorrente — àquele esforço de vigilância.

### 1.1 Três premissas que este trabalho inverte

**Premissa 1 — um patógeno, um programa de vacina.** Os três organismos compartilham um
nicho anatômico (via aérea inferior), uma população de risco (paciente hospitalizado,
idoso, imunocomprometido) e uma janela de intervenção (pré-internação ou admissão). Do
ponto de vista de saúde pública, o alvo relevante é *pneumonia bacteriana grave em
hospedeiro vulnerável*, não cada espécie isoladamente. Propomos tratar o nicho como
unidade de desenho e construir um imunógeno único.

**Premissa 2 — o antígeno protetor está na cápsula.** Restringimos deliberadamente o
espaço de busca a proteínas de superfície do core genome, presentes em ≥95% dos isolados
e independentes de tipagem capsular. Isso sacrifica a imunogenicidade conhecida dos
polissacarídeos em troca de algo que as PCVs não têm: invariância ao sorotipo.

**Premissa 3 — basta que o antígeno seja conservado.** Conservação observada é um
instantâneo; ela não diz se o patógeno *pode* variar aquele sítio sem custo. Por isso
adicionamos um filtro evolutivo: só entram proteínas sob seleção purificadora
demonstrável (dN/dS < 1 por sítio, HyPhy FEL) e ligadas a essencialidade ou virulência.
A intenção é que o escape imunológico exija ao patógeno um custo de aptidão — tornando o
antígeno **evolutivamente ancorado**, não apenas conservado.

### 1.2 Dois problemas metodológicos silenciosos que corrigimos

A literatura de vacinas multi-epitopo *in silico* cresceu muito rápido, e com ela dois
hábitos que consideramos defeitos, ambos corrigíveis sem custo adicional:

**(a) Triagem de segurança incompleta.** O padrão é fazer BLAST contra o proteoma humano
e declarar segurança. Isso ignora duas coisas. Primeiro, homologia global baixa não exclui
mimetismo: um 7-mer contíguo idêntico a uma proteína humana já é suficiente para
reatividade cruzada de células T em diversos modelos experimentais, e passa despercebido
por BLAST em peptídeos curtos se o modo `blastp-short` não for usado. Segundo — e isso é
praticamente ausente da literatura — ninguém verifica homologia com o **microbioma
comensal**. Uma vacina respiratória que induza anticorpos contra epitopos compartilhados
com *Neisseria lactamica*, *Streptococcus mitis* ou *Moraxella* comensais pode perturbar
uma flora protetora e abrir nicho para o próprio patógeno que se quer combater. Incluímos
essa camada.

**(b) Cobertura populacional eurocêntrica e seleção gulosa por score.** A escolha dos
epitopos costuma ser "os N de melhor score", o que é demonstravelmente subótimo: os
melhores ligantes tendem a se concentrar nos mesmos alelos HLA comuns, deixando
populações inteiras sem cobertura. Além disso, os painéis de referência usados refletem
frequências alélicas predominantemente europeias. Tratamos a seleção como problema formal
de cobertura máxima e otimizamos cobertura fenotípica ponderada por frequências alélicas
mundiais **e da população brasileira miscigenada**, reportando ambas separadamente.

### 1.3 O elemento exploratório: epitopos compartilhados por equivalência estrutural

A parte mais especulativa e mais interessante do desenho. KpSC e *A. baumannii* são
Gram-negativos e compartilham famílias de proteínas de membrana externa com dobramento
de barril-β conservado (OmpA-like, receptores TonB-dependentes). *S. pneumoniae* é
Gram-positivo, e homologia de sequência com os outros dois é essencialmente nula. Porém,
proteínas de superfície funcionalmente análogas — adesinas, proteínas ligadoras de
substrato de transportadores ABC — podem apresentar **epitopos conformacionais em posição
estruturalmente equivalente** mesmo sem ancestralidade detectável por sequência.

Buscamos esses casos por sobreposição estrutural (TM-align sobre modelos preditos) em vez
de alinhamento de sequência, e reservamos um bloco do construto para eles. Se a hipótese
se sustentar experimentalmente, uma fração pequena da carga peptídica cobriria os três
organismos simultaneamente. Tratamos isso como hipótese exploratória e a reportamos com
esse rótulo — não como resultado estabelecido.

---

## 2. Materiais e Métodos

Todo o código, configuração e ambiente estão disponíveis em
`https://github.com/⟨PENDENTE:repo_url⟩`, com versões fixadas em `environment.yml` e
orquestração por Snakemake para reprodutibilidade completa a partir de um comando.

### 2.1 Amostragem de genomas

Genomas completos e de nível cromossomo foram obtidos do NCBI RefSeq via `datasets` CLI
para os táxons listados na Tabela 1. A amostragem foi **estratificada por região
geográfica** (América Latina 25%, Ásia 25%, Europa 20%, África 15%, América do Norte 15%),
com preenchimento irrestrito do déficit quando a cota regional não pôde ser atendida.

O motivo é substantivo, não cosmético, e a magnitude do viés já foi medida de forma
independente. Na revisão de mais de 100 mil genomas públicos desses mesmos três patógenos
conduzida por David e colaboradores [@david2025vaccineswatch], 79,2% (79.515/100.381) dos
genomas provêm de países de alta renda e apenas 1,5% (1.470/100.381) de países de baixa
renda; 47,1% dos genomas de pneumococo vêm dos Estados Unidos, e sete países respondem por
79,7% do total. Mais da metade dos países do mundo não contribui com nenhum genoma, fração
que sobe a 65,9% (164/249) para *A. baumannii* e 77,9% (194/249) para *S. pneumoniae*.

Um antígeno que pareça conservado num conjunto com essa composição pode estar apenas
refletindo a clonalidade da amostra e a geografia de quem sequencia. A estratificação é uma
correção parcial e imperfeita desse viés — ela redistribui o que existe, e não pode criar
representatividade onde não há dado. Reportamos a composição resultante explicitamente
(Tabela S1) para que o leitor julgue o alcance da correção.

### 2.2 Anotação e pangenoma

Anotação com Bakta. Pangenoma com Panaroo em `--clean-mode strict`, identidade 0.95.
Genes core definidos como presentes em ≥95% dos genomas. Alinhamentos por gene gerados
com MAFFT para uso posterior em conservação e dN/dS.

### 2.3 Definição do surfaceome vacinável

Filtros sequenciais sobre as proteínas core:

1. **Localização subcelular** (PSORTb 3.0, escore ≥7.5): membrana externa ou
   extracelular para os Gram-negativos; parede celular ou extracelular para o pneumococo.
2. **Exportação**: peptídeo-sinal por SignalP 6 ou localização de superfície confirmada.
3. **Topologia**: no máximo 1 hélice transmembrana (DeepTMHMM). Proteínas politópicas
   são notoriamente difíceis de expressar e purificar em forma nativa.
4. **Tamanho**: 100–1200 aminoácidos.
5. **Bônus funcional**: adesinas, porinas, receptores de sideróforo TonB-dependentes,
   pili/fímbrias, proteínas ligadoras de colina e substratos de ABC receberam prioridade,
   por serem as classes com maior taxa histórica de sucesso em vacinas bacterianas.

### 2.4 Filtro evolutivo

Para cada gene candidato, dN/dS por sítio via HyPhy FEL sobre o alinhamento de codons e
árvore IQ-TREE. Retivemos genes com evidência de seleção purificadora predominante
(dN/dS < 1). Cruzamos com essencialidade (DEG e telas TraDIS/Tn-seq publicadas) e com
anotação de virulência, exigindo pelo menos um dos dois.

### 2.5 Predição de epitopos

- **MHC-I**: NetMHCpan-4.1 BA via IEDB, 9- e 10-mers, painel de 27 alelos de referência,
  percentil ≤1.0.
- **MHC-II**: NetMHCIIpan via IEDB, 15-mers, 27 alelos, percentil ≤2.0.
- **Células B lineares**: BepiPred 2.0, limiar 0.55, comprimento 12–22.
- **Células B conformacionais**: DiscoTope 3.0 sobre os modelos estruturais.
- **Antigenicidade**: VaxiJen v2.0, modelo bacteriano, limiar 0.5.

**Filtro de conservação de epitopo.** Cada peptídeo predito foi buscado como
correspondência exata nas sequências homólogas de todos os isolados do respectivo
organismo; retivemos apenas os presentes em ≥95%. Este filtro é aplicado *antes* de
qualquer consideração de score, e é o mais agressivo do pipeline.

### 2.6 Triagem negativa de segurança (quatro camadas)

- **A — Proteoma humano**: BLASTp em modo `blastp-short` contra UP000005640; descarte se
  identidade local ≥35%.
- **B — Mimetismo por k-mer**: descarte se qualquer 7-mer contíguo do peptídeo ocorrer
  exatamente no proteoma humano.
- **C — Microbioma comensal**: mesmo critério da camada A contra proteomas de comensais
  respiratórios (*Neisseria lactamica*, *Streptococcus mitis*, *S. oralis*, *Moraxella*
  spp., *Corynebacterium* spp.) e de representantes intestinais do Human Microbiome Project.
- **D — Alergenicidade e toxicidade**: AllerTOP v2.0 / AlgPred 2.0 e ToxinPred 3.0.

### 2.7 Seleção do conjunto por cobertura máxima

Formalizamos a escolha do subconjunto de epitopos como problema de cobertura máxima.
Sob equilíbrio de Hardy-Weinberg, a probabilidade de um indivíduo carregar ao menos um
alelo do conjunto coberto *S* é

  Cobertura(S) = 1 − ∏(1 − f_i)²  para todo alelo i coberto por algum epitopo de S,

onde f_i é a frequência do alelo i na população de interesse. Maximizamos essa função por
algoritmo guloso, que para funções submodulares monótonas garante ao menos (1 − 1/e) ≈
63% do ótimo. Usamos uma frequência combinada f = (1−w)·f_mundial + w·f_Brasil, com
w = 0.5, e reportamos as coberturas mundial e brasileira separadamente.

### 2.8 Montagem do construto

Arquitetura N→C: cauda His6 derivada do pET-28a (MGSSHHHHHHSSGLVPRGSH) — adjuvante RS09
(APPHALS, agonista sintético de TLR4) — EAAAK — PADRE — GPGPG — bloco de epitopos B
(linker KK) — bloco MHC-II (GPGPG) — bloco MHC-I (AAY) — bloco compartilhado — His6.

A escolha de cada linker é funcional: EAAAK forma hélice rígida e isola conformacionalmente
o adjuvante; GPGPG interrompe epitopos juncionais espúrios e favorece processamento MHC-II;
AAY é sítio preferencial de clivagem do imunoproteassomo, liberando cada epitopo MHC-I
como unidade correta; KK é sítio de catepsina B.

Registramos uma correção de desenho detectada pelo próprio pipeline: a versão inicial
iniciava a proteína diretamente na cauda de histidinas, o que, pela regra do N-terminal,
classifica a proteína como desestabilizada em *E. coli* (meia-vida < 2 min). A adoção da
cauda N-terminal completa do pET-28a, iniciada por metionina, corrige o problema e ainda
incorpora o sítio de trombina para remoção da cauda pós-purificação.

### 2.9 Estrutura, docking e dinâmica molecular

Modelagem por ColabFold (AlphaFold2 com MMseqs2), refino por GalaxyRefine, validação por
gráfico de Ramachandran, ERRAT e Verify3D. Docking proteína–proteína com HADDOCK 2.4
contra TLR4/MD-2 (PDB 4G8A), TLR2 (2Z7X) e HLA-DR1 (1DLH). Dinâmica molecular em GROMACS,
campo de força CHARMM36, água TIP3P, ⟨PENDENTE:md_ns⟩ ns, três réplicas independentes com
sementes distintas, analisando RMSD, RMSF, raio de giro, ligações de hidrogênio e
componentes principais.

### 2.10 Imunossimulação e manufatura

C-ImmSim com três doses em intervalos de 28 dias. Otimização de códons para *E. coli* K-12
com calibração automática do parâmetro de compromisso entre CAI e diversidade de códons,
evitando o erro de usar sempre o códon ótimo (que gera depleção de tRNA e estrutura
secundária). Auditoria de sítios de restrição do MCS do pET-28a, sequências
Shine-Dalgarno internas, homopolímeros e GC em janela deslizante. Desenhamos em paralelo
uma variante **mRNA-LNP** com a mesma sequência proteica, UTRs de alta expressão, cauda
poli-A de 120 nt, Cap1 e substituição integral de uridina por N1-metilpseudouridina.

---

## 3. Resultados

### 3.1 Disponibilidade e composição geográfica dos genomas completos

A consulta ao RefSeq (julho de 2026, níveis *complete* e *chromosome*) recuperou 4.290
genomas para KpSC — somando *K. pneumoniae* (3.870), *K. quasipneumoniae* (217),
*K. variicola* (203) e *K. quasivariicola* (216) —, 1.143 para *A. baumannii* e apenas
**356** para *S. pneumoniae* (Tabela 1, Figura S1).

Essa assimetria é o primeiro resultado relevante. O pneumococo, apesar de ser o mais
estudado dos três e o único com vacinas licenciadas, é o que menos dispõe de montagens
fechadas: 356 genomas completos contra 5.307 em nível de *scaffold*. A genômica
pneumocócica foi construída sobre sequenciamento de leituras curtas voltado a tipagem
capsular, e montagens fechadas nunca foram prioridade.

A composição geográfica desses genomas completos revela um viés cuja **forma difere** do
já descrito para o conjunto que inclui montagens rascunho. David e colaboradores
(David et al., 2025) reportaram, sobre >100 mil genomas majoritariamente rascunho,
predomínio de países de alta renda e, no caso do pneumococo, dos Estados Unidos (47,1%).
No subconjunto de genomas completos observamos um padrão distinto:

- **KpSC**: 55,6% dos genomas vêm da Ásia, com a China isoladamente respondendo por
  **32,9%** (1.411/4.290) do total mundial. África contribui com 1,0% (45) e América
  Latina com 3,1% (131). 87 países representados.
- ***A. baumannii***: 44,5% da Ásia (China 21,3%), 16,7% da América do Norte. África 1,7%
  (19) e América Latina 6,5% (74). 60 países representados.
- ***S. pneumoniae***: EUA 22,8%, Japão 14,0%. **Apenas 3 genomas (0,8%) de toda a América
  Latina** e 35 países representados no mundo inteiro.

Ou seja: para os dois Gram-negativos, o viés dos genomas completos não é do Norte Global
mas da Ásia Oriental, e é ainda mais concentrado em um único país do que o viés descrito
para os rascunhos. Para o pneumococo, o problema não é a forma do viés mas a escassez
absoluta — 356 genomas de 35 países é uma base estreita para qualquer afirmação sobre
conservação global.

Reportamos isso não como detalhe metodológico mas como achado: **qualquer alegação de
"conservação global" derivada de genomas completos herda uma geografia muito específica**,
e a direção do viés depende de qual nível de montagem se exige.

### 3.2 Validação do pipeline por controles positivos e negativos conhecidos

Antes de aplicar o pipeline a antígenos ainda não descritos, verificamos se ele recupera
os antígenos proteicos de *S. pneumoniae* já estabelecidos na literatura. O core genome
pneumocócico obtido (1.392 de 2.133 proteínas da referência, 65,3%; 52 genomas comparados)
contém:

| Antígeno conhecido | Recuperado no core | Anotação encontrada |
|---|---|---|
| PsaA | sim | metal ABC transporter substrate-binding protein |
| Pneumolisina | sim | cholesterol-dependent cytolysin pneumolysin |
| PhtE (família histidina-tríade) | sim | pneumococcal histidine triad protein PhtE |
| NanB (neuraminidase) | sim | neuraminidase NanB |
| PcsB | sim | peptidoglycan hydrolase PcsB |
| PspC / proteínas ligadoras de colina | sim (2) | PspC domain-containing protein |
| Enolase de superfície | sim | surface-displayed alpha-enolase |
| **PspA** | **não — classificado como acessório** | pneumococcal surface protein A |

O caso do **PspA é o controle mais informativo**, e é negativo por acerto. Aplicando o
critério de ortologia de ≥90% de identidade aminoacídica, o PspA da referência não
encontrou ortólogo em **nenhum** dos 52 genomas comparados (fração de presença = 0,00).
Isso não é falha do método: reflete a biologia conhecida do PspA, cuja diversidade em
clados com identidade inter-clado frequentemente abaixo de 70% é precisamente a razão
histórica pela qual vacinas baseadas nele tiveram desempenho inconsistente entre
populações.

Um pipeline que se propõe a selecionar antígenos invariantes deveria recuperar os
antígenos conservados e rejeitar os variáveis. É o que se observa. Consideramos essa
concordância uma verificação necessária — ainda que não suficiente — da validade do
critério de conservação adotado.

### 3.3 Funil de filtragem: do core genome ao surfaceome

Aplicamos o pipeline a um subconjunto-piloto estratificado por região (60 genomas para
KpSC e *A. baumannii*, 53 para *S. pneumoniae*), usando a anotação PGAP do RefSeq e
ortologia recíproca contra genoma de referência (≥90% de identidade, ≥80% de cobertura).

**Tabela 2 — Funil de filtragem.**

| Etapa | KpSC | *A. baumannii* | *S. pneumoniae* |
|---|---:|---:|---:|
| Genomas completos disponíveis (RefSeq) | 4.290 | 1.143 | 356 |
| Genomas analisados no piloto | 60 | 60 | 53 |
| Proteínas na referência | 5.766 | 4.034 | 2.133 |
| Core genome (≥95% dos genomas) | 3.684 | 2.526 | 1.391 |
| Excluídas por anotação citoplasmática | 1.686 | 1.318 | 688 |
| **Candidatas de superfície** | **188** | **88** | **39** |

A fração de core genome ficou entre 62,6% e 65,3% das proteínas da referência nos três
organismos — coerente com o esperado para essas espécies e um indício de que o critério
de ortologia adotado não é nem excessivamente permissivo nem restritivo demais.

O surfaceome resultante corresponde a 5,1% (KpSC), 3,5% (*A. baumannii*) e 2,8%
(*S. pneumoniae*) do respectivo core genome, totalizando **315 proteínas candidatas**. A
composição por família é dominada, nos Gram-negativos, por proteínas de membrana externa
e receptores TonB-dependentes; no pneumococo, por proteínas ligadoras de substrato de
transportadores ABC e proteínas ancoradas por sortase — exatamente a distribuição que a
diferença de arquitetura de parede celular prevê.

Duas correções foram introduzidas durante a execução, ambas detectadas pelos próprios
dados e registradas aqui por transparência metodológica. Primeiro, proteínas com domínios
repetitivos — notadamente a adesina PavB, com repetições SSURE — geravam múltiplos
alinhamentos locais contra o mesmo genoma, inflando a contagem de presença acima do
número de genomas comparados (frações >1,0); passou-se a deduplicar por proteína-consulta
antes da contagem. Segundo, a triagem por anotação classificava aquaporinas e
aquagliceroporinas como porinas de membrana externa, atribuindo membrana externa a um
organismo Gram-positivo; incluiu-se verificação de coerência entre família predita e
arquitetura de parede celular.

### 3.4 Demais resultados

⟨PENDENTE:secao_resultados⟩

Estrutura prevista da seção, com as tabelas e figuras já definidas:

- **Tabela 1** — Amostragem de genomas por organismo e região.
- **Tabela 2** — Funil de filtragem: core → surfaceome → seleção purificadora → epitopos
  → seguros → selecionados. Uma linha por estágio, uma coluna por organismo.
- **Tabela 3** — Epitopos finais: sequência, proteína de origem, classe, conservação,
  número de alelos ligados, antigenicidade.
- **Tabela 4** — Propriedades físico-químicas do construto.
- **Tabela 5** — Cobertura populacional por região, mundial vs. Brasil.
- **Figura 1** — Fluxograma do pipeline.
- **Figura 2** — Funil de filtragem (diagrama de Sankey).
- **Figura 3** — Mapa do construto com blocos e linkers anotados.
- **Figura 4** — Estrutura 3D prevista, colorida por bloco, com validação de Ramachandran.
- **Figura 5** — Complexo com TLR4/MD-2 e detalhe das interações de interface.
- **Figura 6** — RMSD, RMSF e raio de giro ao longo da dinâmica, três réplicas.
- **Figura 7** — Imunossimulação: IgM/IgG, populações de células B e T, citocinas.
- **Figura 8** — Curva de cobertura populacional acumulada por número de epitopos,
  mundial vs. Brasil — a figura que sustenta o argumento da seleção por cobertura máxima.
- **Figura S1** — Distribuição geográfica e por sequence type da amostra de genomas.

---

## 4. Discussão

⟨PENDENTE:discussao_dados⟩ Os parágrafos abaixo independem dos números finais.

### 4.1 O que este desenho tenta resolver

A contribuição não está em aplicar vacinologia reversa — isso é rotina há duas décadas —
mas em três decisões de escopo. Tratar o nicho clínico, e não a espécie, como unidade de
desenho. Excluir deliberadamente antígenos capsulares, aceitando perda de imunogenicidade
conhecida em troca de invariância ao sorotipo. E exigir que o antígeno esteja
evolutivamente ancorado, de modo que o escape custe aptidão ao patógeno.

A terceira decisão merece ênfase. A história das PCVs mostra que um alvo variável leva a
uma corrida sem linha de chegada. Um alvo sob seleção purificadora forte e acoplado a
função essencial oferece, ao menos em princípio, uma trégua mais duradoura. Não temos
como demonstrar isso *in silico* — é uma predição falsificável, e é assim que a
apresentamos.

### 4.2 A camada do microbioma

Consideramos a triagem contra proteomas comensais a contribuição metodológica mais
imediatamente transferível deste trabalho. Ela é gratuita, leva minutos e é praticamente
ausente da literatura de vacinas multi-epitopo. O raciocínio é direto: a vacina será
administrada a pacientes que estão prestes a receber antibióticos de amplo espectro, ou
seja, a um microbioma já fragilizado. Induzir imunidade adaptativa cruzada contra
comensais protetores da via aérea nessa população é um risco de segurança que sequer
aparece nos desenhos convencionais, porque ninguém mede.

### 4.3 Equidade em vacinologia computacional

Painéis de referência HLA e conjuntos de treinamento de preditores refletem
majoritariamente populações europeias. Um construto otimizado nesse referencial pode ter
cobertura sistematicamente inferior em populações africanas, ameríndias e miscigenadas —
exatamente as que concentram a maior carga de pneumonia bacteriana grave. Reportar
cobertura mundial e brasileira lado a lado, e otimizar explicitamente para ambas, é uma
correção barata e que deveria ser padrão.

### 4.4 Complementaridade com a vigilância de alvos polissacarídicos

Vale explicitar como este trabalho se relaciona com o esforço de vigilância genômica de
alvos vacinais representado por vaccines.watch [@david2025vaccineswatch]. As duas
abordagens tratam dos mesmos três patógenos e são, em nosso entender, complementares em
sentido estrito.

Aquela plataforma responde a uma pergunta de formulação e monitoramento: dada a população
de patógenos circulante, qual composição polissacarídica maximiza cobertura, e como essa
cobertura se desloca após a introdução de uma vacina. É a pergunta certa para o paradigma
conjugado, e os dados que ela produz são a melhor descrição disponível do problema.

O presente trabalho responde a uma pergunta anterior e distinta: existe um conjunto de
alvos cuja eficácia **não** dependa de formulação, porque não varia com o sorotipo. Uma
resposta afirmativa exige sair do espaço polissacarídico — precisamente o espaço que
aquela plataforma cobre — e entrar no espaço proteico, que ela não cobre por escolha de
escopo declarada.

A relação natural entre as duas linhas é de retroalimentação. A vigilância de sorotipos
identifica onde a estratégia atual está perdendo terreno e com que velocidade; o desenho
de antígenos proteicos conservados propõe alvos que, por hipótese, não estão sujeitos a
essa erosão. Nenhuma das duas substitui a outra, e um recurso que rastreasse continuamente
a diversidade de antígenos proteicos — hoje inexistente — seria o análogo proteico daquele
esforço.

### 4.5 Limitações

Enunciadas sem atenuação, porque um artigo *in silico* que subestima suas limitações é
pior do que inútil:

1. **Nada aqui foi validado experimentalmente.** Este é um trabalho gerador de hipótese.
   Predição de ligação a MHC não é imunogenicidade; imunogenicidade não é proteção.
2. **Preditores de epitopo de célula B linear têm desempenho modesto.** BepiPred 2.0
   opera com AUC em torno de 0.6 em conjuntos independentes. Os epitopos B do construto
   são, portanto, a parte menos confiável do desenho.
3. **Epitopos juncionais.** A concatenação pode criar epitopos artefatuais nas junções.
   Os linkers mitigam, mas não eliminam; a verificação exige ressíntese e teste.
4. **Simulação imunológica é ilustrativa.** C-ImmSim é um modelo baseado em agentes com
   parametrização genérica; suas curvas não são preditivas de titulação real.
5. **Modelos estruturais de construtos quiméricos têm confiança baixa nos linkers.**
   Regiões de linker são intrinsecamente desordenadas e o pLDDT nelas deve ser lido como
   indicação de flexibilidade, não como erro corrigível.
6. **A hipótese de epitopos estruturalmente compartilhados é especulativa.** Equivalência
   estrutural não implica reatividade cruzada de anticorpos.
7. **O viés de amostragem foi mitigado, não eliminado.** Regiões com menor capacidade de
   sequenciamento continuam sub-representadas nos repositórios públicos.
8. **O adjuvante RS09 tem menos evidência clínica** do que alternativas como flagelina ou
   agonistas de TLR9; a escolha foi guiada por custo de síntese e o desenho prevê troca.

### 4.6 Caminho de validação proposto

Escalonado por custo crescente, de modo que cada etapa possa matar a hipótese antes da
próxima:

1. **Peptídeos sintéticos + ELISPOT de IFN-γ** contra PBMC de doadores saudáveis com HLA
   tipado (< US$ 3.000). Testa diretamente a imunogenicidade dos epitopos T e é o teste de
   maior retorno por dólar.
2. **ELISA com soro de convalescentes** de infecções pelos três organismos, para verificar
   se os epitopos B são reconhecidos por resposta natural.
3. **Expressão em *E. coli*, purificação e caracterização** (dicroísmo circular para
   confirmar o conteúdo de hélice previsto).
4. **Imunização em camundongo** com leitura de IgG total, subclasses (perfil Th1/Th2) e
   proliferação de esplenócitos.
5. **Ensaio de opsonofagocitose** — o correlato mais próximo de proteção para os três.
6. **Desafio letal** em modelo murino de pneumonia, um organismo por vez.

---

## 5. Conclusão

⟨PENDENTE:conclusao⟩

Apresentamos um desenho computacional completo e reprodutível de um imunógeno quimérico
único dirigido a três patógenos respiratórios prioritários, construído sobre três
inversões deliberadas de premissa — nicho em vez de espécie, proteína core em vez de
cápsula, ancoragem evolutiva em vez de conservação observada — e sobre duas correções
metodológicas de baixo custo: triagem de segurança contra o microbioma comensal e
otimização formal de cobertura HLA com peso para populações sub-representadas.

O construto é uma hipótese testável, não uma vacina. Sua utilidade depende inteiramente da
validação experimental, cujo primeiro passo decisivo custa menos que uma diária de UTI.

---

## Declarações

**Disponibilidade de dados e código.** Todo o código, configuração e ambiente estão em
`https://github.com/⟨PENDENTE:repo_url⟩` sob licença MIT. Os genomas são públicos
(NCBI RefSeq); os acessos exatos estão na Tabela S1. O pipeline reproduz integralmente
os resultados com `snakemake -c8 --use-conda`.

**Financiamento.** ⟨PENDENTE:financiamento⟩

**Conflitos de interesse.** ⟨PENDENTE:conflitos⟩

**Contribuições dos autores.** ⟨PENDENTE:contribuicoes⟩

**Uso de IA.** As análises computacionais foram implementadas com auxílio de assistente de
IA na escrita de código e estruturação do texto. Todo o desenho experimental, as decisões
metodológicas, a interpretação dos resultados e a redação final são de responsabilidade
dos autores. Declaração conforme as diretrizes do periódico.

---

## Referências

⟨PENDENTE:referencias⟩ — ver `manuscript/references.bib`. As citações a inserir cobrem,
no mínimo: carga global de resistência antimicrobiana (GBD/Murray et al.); lista de
patógenos prioritários da OMS; substituição de sorotipo pós-PCV; PSORTb; SignalP 6;
DeepTMHMM; Panaroo; Bakta; IQ-TREE; HyPhy FEL; NetMHCpan-4.1; NetMHCIIpan; BepiPred 2.0;
DiscoTope 3.0; VaxiJen; AllerTOP; ToxinPred; IEDB Population Coverage; regra do N-terminal
(Bachmair/Tobias); índice de instabilidade (Guruprasad); CAI (Sharp & Li); ColabFold;
HADDOCK; GROMACS/CHARMM36; C-ImmSim; RS09; PADRE; e a literatura de antígenos proteicos
de KpSC, *A. baumannii* e pneumococo.
