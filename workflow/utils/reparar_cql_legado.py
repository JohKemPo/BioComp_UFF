"""Repara aspas simples não escapadas dentro do blob `value: '{json}'` em
`.cql` legados — DEC-052.

Um gerador antigo (não presente neste repositório; anterior ao
`neo4jProcessing.py` atual, que já usa `json.dumps` por campo com aspas
duplas) serializava metadados como um blob JSON inteiro dentro de uma string
Cypher de aspas simples, sem escapar as aspas simples do próprio texto
(ex.: "Cote d'Ivoire", presente em cepas históricas de Zika/Spondweni). O
resultado é Cypher sintaticamente inválido: a primeira aspa interna fecha a
string do Cypher antes da hora, e tudo que vem depois — inclusive o `;` real
e a instrução seguinte — vira lixo dentro do "bloco" seguinte, até a PRÓXIMA
aspa desencontrada devolver o estado certo por acidente. Medido em
`Zika_Virus_Singapura_Medium_11seq`: 32 pares de instruções `MERGE
(m:Metadata ...)` fundidas em um bloco só (ver DEC-052 e
`Backend/tests/unit/test_parse_cql_blocks.py`).

Este script varre cada linha `MERGE (m:Metadata {value: '...'})`, extrai o
conteúdo entre a primeira aspa após `value:` e o `'})` de fechamento real
(reconhecido por âncora de sufixo, não por contagem de aspas), escapa toda
aspa simples interna que ainda não estiver escapada, e confere que o JSON
resultante (desescapando `\\'` → `'`) decodifica sem erro antes de gravar.

Uso:
    python workflow/utils/reparar_cql_legado.py <arquivo.cql> [--apply]

Sem `--apply`, roda em modo de conferência (mostra quantos blobs seriam
corrigidos e falha se algum não decodificar). Com `--apply`, sobrescreve o
arquivo (faz backup em `<arquivo>.bak-preDEC052` antes, se ainda não existir).
"""
import argparse
import json
import re
import shutil
import sys
from pathlib import Path

PREFIX = "MERGE (m:Metadata {value: '"
SUFFIX_RE = re.compile(r"'\}\)(\s*;?\s*)$")


def escapar_aspas_internas(conteudo: str) -> str:
    """Escapa toda aspa simples não precedida de '\\' dentro de `conteudo`."""
    saida = []
    i = 0
    while i < len(conteudo):
        c = conteudo[i]
        if c == "\\" and i + 1 < len(conteudo):
            saida.append(conteudo[i : i + 2])
            i += 2
            continue
        if c == "'":
            saida.append("\\'")
            i += 1
            continue
        saida.append(c)
        i += 1
    return "".join(saida)


def reparar_linha(linha: str, stats: dict) -> str:
    idx = linha.find(PREFIX)
    if idx == -1:
        return linha

    prefixo = linha[: idx + len(PREFIX)]
    resto = linha[idx + len(PREFIX) :]
    m = SUFFIX_RE.search(resto)
    if not m:
        stats["sem_fechamento"] += 1
        return linha

    conteudo = resto[: m.start()]
    sufixo = resto[m.start() :]

    stats["blobs"] += 1

    corrigido = escapar_aspas_internas(conteudo)
    if corrigido == conteudo:
        return linha  # já estava bem escapado

    # Confere que o JSON decodifica depois de desescapar — nunca grava um
    # blob que não seja um JSON válido.
    try:
        json.loads(corrigido.replace("\\'", "'"))
    except json.JSONDecodeError as e:
        stats["falhas_validacao"] += 1
        raise ValueError(f"JSON inválido após reparo: {e}") from e

    stats["corrigidos"] += 1
    return prefixo + corrigido + sufixo


def reparar_arquivo(caminho: Path, aplicar: bool) -> dict:
    stats = {"blobs": 0, "corrigidos": 0, "sem_fechamento": 0, "falhas_validacao": 0}
    linhas_saida = []
    with caminho.open("r", encoding="utf-8") as f:
        for numero, linha in enumerate(f, 1):
            try:
                linhas_saida.append(reparar_linha(linha, stats))
            except ValueError as e:
                raise ValueError(f"{caminho}:{numero}: {e}") from e

    if aplicar and stats["corrigidos"] > 0:
        backup = caminho.with_suffix(caminho.suffix + ".bak-preDEC052")
        if not backup.exists():
            shutil.copy2(caminho, backup)
        caminho.write_text("".join(linhas_saida), encoding="utf-8")

    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arquivo", type=Path)
    parser.add_argument("--apply", action="store_true", dest="aplicar")
    args = parser.parse_args()

    resultado = reparar_arquivo(args.arquivo, args.aplicar)
    modo = "APLICADO" if args.aplicar else "conferência (sem gravar)"
    print(f"{args.arquivo} [{modo}]")
    print(f"  blobs Metadata encontrados: {resultado['blobs']}")
    print(f"  corrigidos (aspa não escapada): {resultado['corrigidos']}")
    print(f"  sem padrão de fechamento reconhecido: {resultado['sem_fechamento']}")
    if resultado["sem_fechamento"]:
        sys.exit(1)
