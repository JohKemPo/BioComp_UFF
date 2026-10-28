#!/usr/bin/env bash
# Gera as figuras do dataflow do AkôFlow (docs/akoflow/03-dataflow.md):
#   dataflow-geral.pdf/.png e dataflow-zika21.pdf/.png
#
# O que é preciso ter instalado (fora do ambiente conda do projeto):
#   - Node.js com npm (testado com node 22.22.3)
#   - @mermaid-js/mermaid-cli 12.0.0, instalado em $MMDC_PREFIX
#     (padrão /tmp/mmdc-ptm; o script instala sozinho se não existir)
#   - ghostscript (gs) e poppler-utils (pdfinfo, pdffonts, pdftotext)
#
# Layout: a versão LR do Mermaid sai com 4864 px de largura (geral) e
# 6739 px (Zika-21 não-LR). Para 170 mm usamos a mesma fonte com
# `flowchart TD` (única troca: a primeira linha do bloco) e espaçamento
# menor (cfg). Nós e arestas não mudam: o script confere a contagem no SVG.
#
# Limite conhecido: NENHUMA das duas figuras atinge 7 pt no tamanho final
# de 170 mm. O script mede o texto no PDF final e reporta a faixa. Ver
# docs/akoflow/03-dataflow.md, seção "Figuras".
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
DOC="$ROOT/docs/akoflow/03-dataflow.md"
OUT="$ROOT/docs/akoflow/figuras"
MMDC_PREFIX="${MMDC_PREFIX:-/tmp/mmdc-ptm}"
MMDC="$MMDC_PREFIX/node_modules/.bin/mmdc"
LARGURA_PT=481.89   # 170 mm em pontos
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

if [ ! -x "$MMDC" ]; then
  npm install --prefix "$MMDC_PREFIX" @mermaid-js/mermaid-cli@12.0.0 >/dev/null 2>&1
fi

cat > "$TMP/cfg.json" <<'EOF'
{"flowchart":{"nodeSpacing":18,"rankSpacing":28,"padding":6,"wrappingWidth":150},"fontSize":14}
EOF

# Extrai o k-ésimo bloco ```mermaid``` do documento (1 = geral, 2 = Zika-21).
extrai() {
  awk -v k="$1" '/^```mermaid/{n++; inb=(n==k); next} /^```/{inb=0} inb' "$DOC"
}

# Gera uma figura: $1 = índice do bloco, $2 = nome do arquivo de saída.
gera() {
  local bloco="$1" nome="$2" mmd="$TMP/$2.mmd"
  extrai "$bloco" > "$mmd"
  grep -q '^flowchart LR$' "$mmd" || { echo "bloco $bloco não começa em 'flowchart LR'" >&2; exit 1; }
  sed -i '1s/^flowchart LR$/flowchart TD/' "$mmd"

  "$MMDC" -q -c "$TMP/cfg.json" -i "$mmd" -o "$TMP/$nome-raw.pdf" -b white
  "$MMDC" -q -c "$TMP/cfg.json" -i "$mmd" -o "$TMP/$nome.svg" -b white

  # Contagem no SVG renderizado (confere com a origem).
  python3 - "$TMP/$nome.svg" "$mmd" <<'PY'
import re, sys
svg = open(sys.argv[1], encoding="utf-8").read()
src = open(sys.argv[2], encoding="utf-8").read()
arestas_src = len(re.findall(r'-->|-\.->|==>', src))
nos_svg = len(re.findall(r'class="node ', svg))
arestas_svg = len(re.findall(r'class="[^"]*flowchart-link', svg))
print(f"  nós: SVG={nos_svg}  arestas: SVG={arestas_svg}  arestas na origem={arestas_src}")
PY

  # Ajusta a página à largura de 170 mm, mantendo o vetor.
  local H
  H=$(pdfinfo "$TMP/$nome-raw.pdf" | LC_ALL=C awk -v w="$LARGURA_PT" '/Page size/{printf "%.2f", $5*w/$3}')
  gs -q -o "$OUT/$nome.pdf" -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS="$LARGURA_PT" \
     -dDEVICEHEIGHTPOINTS="$H" -dFIXEDMEDIA -dPDFFitPage -dCompatibilityLevel=1.6 \
     "$TMP/$nome-raw.pdf"

  # PNG 300 dpi, fundo branco (gs já compõe sobre branco).
  gs -q -sDEVICE=png16m -r300 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
     -o "$OUT/$nome.png" "$OUT/$nome.pdf"

  # Fontes: nenhuma pode estar sem embutir.
  if pdffonts "$OUT/$nome.pdf" | tail -n +3 | awk '$(NF-4)=="no"{f=1} END{exit !f}'; then
    echo "ERRO: fonte não embutida em $nome.pdf" >&2; exit 1
  fi

  # Tamanho final do texto (altura das palavras, em pt, no PDF de 170 mm).
  pdftotext -bbox "$OUT/$nome.pdf" - | grep '<word' \
    | LC_ALL=C sed -E 's/.*yMin="([0-9.]+)".*yMax="([0-9.]+)".*/\2 \1/' \
    | LC_ALL=C awk '{printf "%.2f\n", $1-$2}' | sort -g \
    | LC_ALL=C awk '{a[NR]=$1} END{printf "  texto no PDF de 170 mm: min %.1f  mediana %.1f  max %.1f pt  (n=%d palavras)\n", a[1], a[int(NR/2)+1], a[NR], NR}'
  echo "  gerado: $OUT/$nome.pdf e $OUT/$nome.png"
}

mkdir -p "$OUT"
echo "dataflow-geral:"
gera 1 dataflow-geral
echo "dataflow-zika21:"
gera 2 dataflow-zika21

# Figura de impressão (170 mm, 7 pt ou mais): Graphviz a partir de dag-artigo.dot.
# Fonte separada do Mermaid; o texto tem 9,5 pt na origem, em escala 0,856 (563 pt -> 482 pt).
echo "dag-artigo:"
dot -Tpdf "$OUT/dag-artigo.dot" -o "$TMP/dag-raw.pdf"
H=$(pdfinfo "$TMP/dag-raw.pdf" | LC_ALL=C awk -v w="$LARGURA_PT" '/Page size/{printf "%.2f", $5*w/$3}')
gs -q -o "$OUT/dag-artigo.pdf" -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS="$LARGURA_PT" \
   -dDEVICEHEIGHTPOINTS="$H" -dFIXEDMEDIA -dPDFFitPage -dCompatibilityLevel=1.6 "$TMP/dag-raw.pdf"
gs -q -sDEVICE=png16m -r300 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -o "$OUT/dag-artigo.png" "$OUT/dag-artigo.pdf"
echo "  gerado: $OUT/dag-artigo.pdf e $OUT/dag-artigo.png"
