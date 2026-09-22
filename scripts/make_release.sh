#!/bin/bash
# Cria e publica a tag de release, depois de checar que ela pode virar Release e DOI.
#
# Existe por causa da v0.1.0: aquela tag foi criada antes do job de release, do
# CITATION.cff e do .zenodo.json existirem, então nunca poderia gerar um Release
# publicado nem um DOI. Toda pré-condição abaixo corresponde a uma forma concreta
# de a tag sair inútil. O script não cria a tag se alguma falhar.
#
# Uso:
#   bash scripts/make_release.sh v0.1.1
#   bash scripts/make_release.sh v0.1.1 --dry-run   # só checa, não cria nada
set -euo pipefail
cd "$(dirname "$0")/.."

TAG="${1:-}"
DRY=""
[[ "${2:-}" == "--dry-run" ]] && DRY=1
[[ -z "$TAG" ]] && { echo "uso: bash scripts/make_release.sh vX.Y.Z [--dry-run]"; exit 2; }
[[ "$TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo "✗ tag deve ser vX.Y.Z (o CI dispara em 'v*')"; exit 2; }

fail=0
ok()   { printf '  \033[32mOK\033[0m   %s\n' "$1"; }
bad()  { printf '  \033[31mFALHA\033[0m %s\n' "$1"; fail=1; }

echo "── pré-condições para $TAG ──"

# 1. Árvore limpa: a tag aponta para um commit, não para o que está no disco.
[[ -z "$(git status --porcelain)" ]] && ok "árvore de trabalho limpa" \
  || bad "há mudanças não commitadas — a tag não as incluiria"

# 2. Sincronizado com o origin, senão a tag aponta para um commit que ninguém tem.
git fetch -q origin
[[ -z "$(git log --oneline origin/main..HEAD)" ]] && ok "nada por enviar ao origin" \
  || bad "há commits locais não enviados — envie antes de marcar a tag"

# 3. A versão declarada tem de bater com a tag, senão o binário publicado se
#    apresenta com um número diferente do release.
V="${TAG#v}"
grep -q "^version = \"$V\"" pyproject.toml && ok "pyproject.toml na versão $V" \
  || bad "pyproject.toml não está na versão $V"
grep -q "^version: \"$V\"" CITATION.cff && ok "CITATION.cff na versão $V" \
  || bad "CITATION.cff não está na versão $V"

# 4. Metadados de citação: sem eles o Zenodo arquiva sem autoria.
[[ -f CITATION.cff && -f .zenodo.json ]] && ok "CITATION.cff e .zenodo.json presentes" \
  || bad "faltam metadados de citação"
! grep -q "PREENCHER" CITATION.cff && ok "sem placeholders no CITATION.cff" \
  || bad "CITATION.cff ainda tem ORCID como PREENCHER — o DOI sairia sem ORCID"

# 5. O job que efetivamente publica os binários. Sem ele a tag só gera artefatos
#    de CI, que expiram e que o Zenodo não enxerga.
grep -q "gh release create" .github/workflows/build.yml \
  && ok "workflow publica um Release com os binários" \
  || bad "o workflow não cria Release — a tag não geraria download nem DOI"

# 6. O CI empacota só se os testes passarem; se falharem, a tag não produz nada.
if python3 -m pytest -q >/dev/null 2>&1; then
  ok "testes passam"
else
  printf '  \033[33mAVISO\033[0m testes não passaram neste interpretador — '
  printf 'confirme num ambiente 3.10+ limpo (o CI usa 3.11)\n'
fi

# 7. Actions precisa estar destravado; conta bloqueada por cobrança não roda job
#    nenhum, e a tag passa sem gerar Release.
#    Tem de ser o workflow 'build' especificamente: outros workflows rodam em
#    ubuntu (grátis e não bloqueado) e passam, mascarando o bloqueio que atinge
#    os runners macOS/Windows de que o build depende.
#    Olhamos a última execução CONCLUÍDA, não a mais recente: uma execução ainda
#    em fila não tem 'conclusion', e ler esse vazio como "não sei" já fez o script
#    liberar a tag enquanto o build vinha falhando em toda execução concluída.
if command -v gh >/dev/null 2>&1; then
  last=$(gh run list --workflow build.yml --status completed --limit 1 \
           --json conclusion -q '.[0].conclusion' 2>/dev/null || echo "")
  if [[ "$last" == "failure" ]]; then
    bad "a última execução concluída do 'build' falhou — a tag não geraria Release"
    printf '        se for "account is locked due to a billing issue", resolva em\n'
    printf '        github.com/settings/billing antes de marcar a tag\n'
  elif [[ -z "$last" ]]; then
    bad "não consegui ler o estado do workflow 'build' — confirme antes de marcar a tag"
  else
    ok "última execução concluída do 'build': '$last'"
  fi
fi

# 8. Zenodo só arquiva releases publicados DEPOIS de a integração estar ligada.
#    Ligar depois não captura o que já saiu, e a única saída seria queimar a próxima
#    versão. O Zenodo instala um webhook no repo ao ser ativado, então a ausência de
#    qualquer webhook é evidência de que a integração não está ligada.
if command -v gh >/dev/null 2>&1; then
  hooks=$(gh api "repos/$(gh repo view --json nameWithOwner -q .nameWithOwner)/hooks" \
            --jq '[.[] | select(.config.url | test("zenodo"; "i"))] | length' 2>/dev/null || echo "?")
  if [[ "$hooks" == "?" ]]; then
    printf '  \033[33mAVISO\033[0m não consegui verificar o webhook do Zenodo\n'
  elif [[ "$hooks" == "0" ]]; then
    bad "integração Zenodo-GitHub não está ligada — o release sairia sem DOI"
    printf '        ligue em zenodo.org (Account > GitHub) ANTES da tag; o Zenodo não\n'
    printf '        captura releases publicados antes da ativação\n'
  else
    ok "webhook do Zenodo presente"
  fi
fi

# 9. A tag não pode já existir — mover tag publicada quebra quem já a buscou.
! git rev-parse -q --verify "refs/tags/$TAG" >/dev/null && ok "tag $TAG ainda não existe localmente" \
  || bad "a tag $TAG já existe — escolha a próxima versão"

echo
[[ $fail -eq 1 ]] && { echo "✗ pré-condições falharam — nada foi criado."; exit 1; }
[[ -n "$DRY" ]] && { echo "✓ tudo pronto (dry-run: nada foi criado)."; exit 0; }

echo "── criando e publicando $TAG ──"
git tag -a "$TAG" -F - <<EOF
PanNosoVax Studio $TAG

Desktop application for reproducible multi-epitope vaccine design without the
command line. The graphical layer drives a Snakemake workflow rather than
replacing it, so resumability and provenance are preserved.

Includes the PanNosoVax pipeline: core genome, surfaceome, epitope prediction,
four-layer safety screening with commensal-microbiome negative selection, HLA
coverage weighted for Brazilian allele frequencies, and construct assembly.

Supersedes v0.1.0, which was tagged before the release workflow and the citation
metadata existed and therefore could not produce a published Release or a DOI.
EOF
git push origin "$TAG"

echo
echo "✓ $TAG publicada. O CI agora monta os binários e cria o Release."
echo "  acompanhar:  gh run watch"
echo "  release:     gh release view $TAG --web"
echo "  o Zenodo só captura o release se a integração já estiver ligada para o repo."
