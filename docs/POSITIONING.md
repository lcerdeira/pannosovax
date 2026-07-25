# Posicionamento em relação a vaccines.watch

## O que a fonte realmente é

**Preprint:** David S, Abudahab K, Couto N, Daningrat WOD, Yeats C, Molloy A, Ashton PM,
Alikhan N-F, Aanensen DM, e a NIHR Global Health Research Unit on Genomics and enabling
data for Surveillance of AMR. *Monitoring of vaccine targets and interventions using global
genome data: vaccines.watch.* bioRxiv 2025.06.13.659488, postado em **18 de junho de 2025**.
Licença CC-BY-NC. Financiamento: NIHR (NIHR133307) e Gates Foundation (INV-025280).

É o grupo do **Centre for Genomic Pathogen Surveillance (Oxford)** — a mesma equipe do
Pathogenwatch. Não é um grupo pequeno nem sem recursos.

## Correção factual importante

O campo `published: NA` na API do bioRxiv **não significa que o artigo foi rejeitado**.
Significa apenas que o bioRxiv ainda não vinculou um DOI de periódico ao preprint. Para um
artigo de recurso/plataforma vindo de um consórcio financiado por NIHR e Gates, treze meses
entre preprint e versão final é inteiramente normal — e é bem possível que já esteja
publicado, ou em revisão avançada, sem que o vínculo tenha sido registrado.

**Verificação feita em julho de 2026:** não localizei versão publicada em periódico — o
preprint segue como preprint. Mas isso não indica projeto abandonado nem rejeitado. O mesmo
grupo depositou um artigo do **Pathogenwatch** no medRxiv em **março de 2026**
(10.64898/2026.03.18.26348693), e a plataforma vaccines.watch continua no ar ingerindo
genomas a cada 4 horas. O grupo está ativo e produzindo. Vale reconferir na véspera da
submissão, porque a citação precisa ser atualizada se a versão de periódico sair.

## O que eles cobrem — e o que deliberadamente não cobrem

Isto é o achado que define nossa estratégia.

vaccines.watch é **exclusivamente sobre alvos polissacarídicos**:

| Patógeno | Alvos rastreados |
|----------|------------------|
| *S. pneumoniae* | sorotipos capsulares (base de todas as PCVs licenciadas) |
| KpSC | cápsula (K) e antígeno O do LPS |
| *A. baumannii* | cápsula (K) e outer core do LOS (OC) |

Proteínas aparecem no artigo **uma única vez**, na introdução, como contexto histórico
("vacinologia reversa... MenB"). Não há nenhum antígeno proteico na plataforma. E a seção de
trabalhos futuros deles é explícita: pretendem estender para **mais alvos polissacarídicos**
de outros patógenos, dados de coleções personalizadas e mais metadados.

**Ou seja: eles mapearam o problema com precisão e escolheram não atacar essa parte dele.**
O espaço de antígenos proteicos conservados está aberto, e não por acidente — está fora do
escopo declarado deles.

## Por que isso é uma boa notícia, e não uma disputa

Nossa POC é anticapsular por construção. A tese central do PanNosoVax — só proteínas core de
superfície, independência de sorotipo, ancoragem evolutiva — é **exatamente a resposta ao
problema que os dados deles documentam**. Os números do preprint são a melhor seção de
motivação que nosso artigo poderia ter:

- PCV13 cobre apenas **36,2%** (11.907/32.918) dos genomas públicos globais de pneumococo;
  PCV21 chega a 87,4% — mas ao custo de uma corrida de valência crescente que já vai a 21
  sorotipos.
- **79,2%** (79.515/100.381) de todos os genomas vêm de países de alta renda. Apenas
  **1,5%** (1.470/100.381) de países de baixa renda.
- 47,1% dos genomas de pneumococo vêm dos EUA. Sete países respondem por 79,7% do total.
- Mais da metade dos países não contribui com **nenhum** genoma — subindo para 65,9%
  (164/249) em *A. baumannii* e 77,9% (194/249) em *S. pneumoniae*.
- Só existe um candidato de vacina para KpSC em ensaio (Kleb4V, fase 1/2, NCT04959344) e
  **nenhum** em desenvolvimento clínico ativo para *A. baumannii*.

Aquele 1,5% de países de baixa renda é a justificativa quantitativa, medida por terceiros,
da nossa amostragem estratificada por região e da ponderação por frequências HLA brasileiras.
Não precisamos argumentar que o viés existe — eles já mediram.

## A estratégia recomendada: complementar e citar, não competir

**Não competir na plataforma deles.** Eles têm >100 mil genomas, ingestão automática a cada
4 horas, a infraestrutura do Pathogenwatch por baixo, financiamento Gates/NIHR e treze meses
de vantagem. Essa corrida está perdida antes de começar, e a comunidade de genômica de
pneumococo e *Klebsiella* é pequena o suficiente para que eles sejam, muito provavelmente,
os nossos revisores.

**Ocupar o nicho adjacente**, que está genuinamente vago:

| | vaccines.watch | PanNosoVax (nosso) |
|---|---|---|
| Alvo | polissacarídeo (K, O, OC, sorotipo) | proteína core de superfície |
| Pergunta | qual formulação cobre a população atual? | qual antígeno não depende de formulação? |
| Horizonte | vigilância e monitoramento pós-rollout | desenho de imunógeno independente de sorotipo |
| Falha que endereça | cobertura incompleta | substituição de sorotipo |
| Saída | painel de diversidade de alvos | construto quimérico + validação in silico |

São perguntas complementares sobre os mesmos três patógenos. Um artigo que se posiciona como
a peça proteica ao lado da peça polissacarídica deles é **mais publicável**, não menos: ganha
uma motivação quantitativa robusta de graça, e não convida a uma briga de prioridade.

## Sobre "marcar território"

Vale dizer isto com todas as letras, porque é a diferença entre um bom artigo e um problema
de carreira.

O preprint deles, postado em junho de 2025, **já é a marcação de território**. É exatamente
para isso que preprints existem. Pegar os pontos fortes de um preprint ainda não publicado e
correr para publicar antes é *scooping* — e a licença CC-BY-NC permite reuso **com
atribuição**, não apropriação de ideias.

Numa comunidade deste tamanho, isso seria reconhecido de imediato, pelas pessoas de quem
Louise mais precisa como colaboradoras e revisoras. O custo reputacional é desproporcional a
qualquer ganho.

A boa notícia é que o caminho honesto é também o mais forte: citá-los como motivação,
declarar a complementaridade de forma explícita, e ocupar o espaço proteico que eles
deixaram aberto por escolha. Isso é ciência normal e saudável — construir sobre um preprint
com citação não requer permissão de ninguém.

## Sobre contatar o grupo — sequência corrigida

> **Nota de revisão.** Uma versão anterior deste documento recomendava escrever ao grupo
> antes de publicar. **Essa recomendação está retirada.** Confirmou-se que não há vínculo
> com o consórcio NIHR, e a posição é de pesquisadora independente, sem financiamento
> estabelecido, diante de um grupo bem financiado atuando nos mesmos três patógenos.

A regra passa a ser: **prioridade pública primeiro, conversa depois — se houver.**

Não é preciso atribuir má-fé a ninguém para justificar isso. A proteção é barata e cobre
igualmente o cenário muito mais comum de **descoberta paralela honesta**, em que alguém
chega ao mesmo resultado de forma independente e a prioridade se perde sem que ninguém
tenha feito nada errado.

Depois que o preprint estiver no ar com data, contatar deixa de ter custo e passa a ter
vantagem — inclusive a de citação recíproca. Antes disso, não.

Ver [STRATEGY.md](STRATEGY.md) para a sequência completa de proteção de prioridade e para
a decisão sobre o formato da plataforma.
