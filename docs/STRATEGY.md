# Estratégia: POC → preprint → grant, e o que a plataforma deve ser

## Correção da recomendação anterior

Eu havia sugerido contatar o grupo do CGPS antes de publicar. **Para o seu caso, a
sequência estava invertida.** Quem não tem financiamento nem consórcio não deve revelar
uma linha de trabalho não publicada a um grupo mais bem financiado que atua nos mesmos
três patógenos. Não é preciso atribuir má-fé a ninguém para chegar a essa conclusão: a
proteção estrutural é barata, e protege igualmente contra o cenário muito mais comum e
igualmente danoso de **descoberta paralela honesta** — alguém chegar ao mesmo lugar
sozinho, seis meses antes, sem nunca ter ouvido falar de você.

A sequência correta é: **estabelecer prioridade pública primeiro, conversar depois (ou
nunca).** Depois que o preprint está no ar com data, revelar não custa mais nada — e
inclusive passa a ser vantajoso.

## Proteção de prioridade: o que fazer, em ordem

Um preprint datado é a proteção real. Ele torna a apropriação ineficaz e visível, porque
a prioridade passa a ser pública e verificável por qualquer pessoa. É exatamente para
isso que preprints existem.

1. **Registre o domínio agora.** Custa ~US$ 15/ano e nomes são por ordem de chegada. Faça
   isso hoje, muito antes de precisar dele.
2. **Preprint no bioRxiv assim que o POC estiver sólido.** Não espere o artigo estar
   perfeito; espere ele estar defensável. Preprint é carimbo de data, não versão final.
3. **DOI no Zenodo para código e dados no mesmo dia do preprint.** Vincule os dois. Isso
   estabelece prioridade sobre o *método*, não só sobre o texto.
4. **Caderno de laboratório com datas** — commits no git já servem, desde que o repositório
   seja tornado público na data do preprint.
5. Só então, se quiser, converse com quem for.

Ordem inversa a essa é que cria exposição.

## A pergunta da plataforma: nenhuma das duas opções

Você colocou como "plataforma completa **ou** site de consulta". Minha sugestão é que
a divisão certa é outra, e ela muda a decisão:

> **O ativo científico é o índice de antígenos versionado e citável. O site é só a
> vitrine.**

Construa o conjunto de dados primeiro. O site vem de graça em cima dele.

### Por que não a plataforma viva

Uma plataforma com ingestão contínua é **passivo de manutenção**, e é precisamente onde o
CGPS é mais forte: eles já têm >100 mil genomas, ingestão a cada 4 horas e a infraestrutura
do Pathogenwatch por baixo, financiada por NIHR e Gates. Competir ali é escolher o único
terreno onde a assimetria de recursos é máxima.

Pior para o seu objetivo específico: **um serviço no ar é um compromisso operacional que o
avaliador de grant enxerga como custo recorrente**, não como capacidade. E o modo de falha
é cruel — uma plataforma que cai durante a avaliação do edital é pior do que plataforma
nenhuma.

### Por que não só um site de consulta

Um site de consulta sem um conjunto de dados versionado por trás não é citável, não gera
DOI, e não constitui evidência de capacidade. É marketing sem ativo.

### O que fazer: índice versionado + site estático

```
   dados: Parquet versionado, DOI Zenodo a cada release  ← o ativo citável
             │
   site: HTML estático + DuckDB-WASM (consulta roda no navegador do usuário)
             │
   custo: ~US$ 0–5/mês, sem servidor, sem banco, sem plantão
```

Propriedades que interessam à sua estratégia:

- **Citável.** Cada release (v1.0, v1.1…) tem DOI próprio. Citações acumulam.
- **Não apodrece.** Sem backend, não há o que cair. Um site estático de 2026 ainda estará
  no ar em 2030 sem manutenção.
- **Não promete o que não pode cumprir.** "Índice v1.0, dados de julho de 2026" é uma
  afirmação honesta e permanentemente verdadeira. "Atualizado continuamente" é uma dívida.
- **Compete onde você é forte.** A camada analítica — conservação, dN/dS, epitopo,
  cobertura HLA — é o que eles não têm e declararam não pretender construir.
- **Evidência de capacidade para o grant.** Dado + código + site funcionando = você
  demonstrou que consegue executar. É isso que o edital compra.

### Faseamento

| Fase | Entrega | Esforço | Quando |
|------|---------|---------|--------|
| 0 | POC roda de ponta a ponta nos 3 patógenos | ~1 mês | agora |
| 1 | **Preprint + Zenodo no mesmo dia** | ~2 semanas | imediatamente após a fase 0 |
| 2 | Índice v1.0 em Parquet com DOI | ~2 semanas | após fase 1 |
| 3 | Site estático com consulta e curva de cobertura HLA | ~1 mês | após fase 2 |
| 4 | Submissão a periódico + submissão a edital | ~1 mês | após fase 3 |

Ingestão contínua **não aparece nesta lista**. Ela entra na proposta de grant como
*trabalho a financiar*, não como coisa já feita — que é exatamente onde ela tem valor para
você.

## Sobre a redação da proposta de grant

Não posicione o projeto como concorrente do vaccines.watch. Posicione como a peça que
falta, citando-os. Editais gostam de ver que o proponente conhece o cenário e ocupa uma
lacuna real em vez de duplicar esforço — e "eles cobrem polissacarídeo, nós cobrimos
proteína, os dados deles quantificam a lacuna que atacamos" é um parágrafo muito forte
de justificativa.

Editais compatíveis: FAPESP (Auxílio Regular, Jovem Pesquisador), CNPq Universal, CARB-X,
GARDP, Wellcome Trust, e chamadas de AMR da Fiocruz/Butantan. CARB-X e GARDP financiam
especificamente alvos de resistência antimicrobiana; o argumento de equidade — cobertura
HLA para populações sub-representadas, ancorado no dado de que 1,5% dos genomas públicos
vêm de países de baixa renda — encaixa com precisão incomum no escopo declarado deles.

## Resumo em uma linha

Publique o preprint antes de falar com qualquer pessoa; construa o **conjunto de dados
citável**, não o serviço; e trate a ingestão contínua como aquilo que você pede dinheiro
para fazer, não como aquilo que você faz de graça antes de pedir.
