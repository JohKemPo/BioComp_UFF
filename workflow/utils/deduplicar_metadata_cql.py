"""Deduplica `Metadata`/`Feature`/`Qualifier` em `.cql` já gerados — G1.

`docs/science/05-grafo-neo4j.md §G1`: até a correção de `neo4jProcessing.py`,
cada ocorrência de um terminal em cada subárvore recriava com `CREATE` o nó
`Metadata` e todo o registro GenBank abaixo dele. Os `.cql` das reexecuções
de 2026-09 têm milhões de `Qualifier` por isso (5,29 M em VARV-49).

Regenerar o `.cql` exige reexecutar a mineração inteira: o gerador consome o
resultado em memória de `SubtreeMiner.find_exact_subsets`, que não é gravado
em disco. O `.cql` legado, porém, contém toda a informação necessária — cada
bloco de metadado traz a subárvore, o `terminal_hash` e o registro completo —,
então este script o reescreve no formato do gerador corrigido, usando o
MESMO `DeduplicadorMetadata`: o resultado é byte a byte o que o gerador novo
produziria (verificado em `test_neo4j_processing.py` contra o golden
pré-correção).

Blocos que não são de metadado (User, Tree, Subtree, Support) passam
inalterados, com o espaço em branco original.

O corte em blocos segue a regra de `Backend/src/services/cql_batch_service.py
::parse_cql_blocks` (C-5e): `;` só termina bloco fora de string simples ou
dupla, com `\\` como escape. Diferença deliberada: aqui NÃO se removem
comentários `//` — em `Zika_ZIKV480_reexec_20260901` há 500 ocorrências de
`//` dentro de descrições do GenBank ("ZIKV/Homo_sapiens/USA//2016/...").

Uso:
    python -m workflow.utils.deduplicar_metadata_cql <arquivo.cql> [--saida X]

Grava `<arquivo>.reparado.cql` (nunca sobrescreve o original), verifica o
resultado contra o original e imprime as contagens antes/depois.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import re
import resource
import sys
import time
from typing import Iterator

from workflow.utils.neo4jProcessing import (
    USER_PLACEHOLDER,
    DeduplicadorMetadata,
    linha_match_subtree,
    project_name_from_output_path,
)

_ESPECIAIS = re.compile(r"[\\'\";]")
_ABRE_LEGADO = "CREATE (m:Metadata {"
_LIGACAO_LEGADA = "CREATE (child)-[:HAS_METADATA]->(m)"
_RE_MATCH_SUBTREE = re.compile(
    r"^MATCH \(child:Subtree \{name: '([^']*)', uid: '" + re.escape(USER_PLACEHOLDER) + r"'\}\)$")
_RE_HASH_FINAL = re.compile(r"(?:^|, )terminal_hash: (-?\d+)$")
_RE_FEATURE = re.compile(r"^CREATE \(m\)-\[:HAS_FEATURE\]->\(f_\d+:Feature \{")
_RE_QUALIFIER = re.compile(r"^CREATE \(f_\d+\)-\[:HAS_QUALIFIER\]->\(q_\d+_\d+:Qualifier \{")
_RE_PROJETO = re.compile(r"CREATE \(t:Tree \{name: '[^']*', uid: '[^']*', project: '([^']*)'\}\)")


def segmentos(texto: str) -> Iterator[str]:
    """Fatia `texto` em segmentos que terminam em `;` fora de string.

    Cada segmento inclui o espaço em branco que o precede, de modo que
    ``"".join(segmentos(t)) == t``. O último segmento pode não terminar em `;`
    (sobra final do arquivo).

    A regra de aspas/escape é a de `parse_cql_blocks` (Backend), percorrendo
    só os caracteres que mudam de estado — equivalente ao laço caractere a
    caractere, mas em velocidade de `re`.
    """
    inicio = 0
    aspa_simples = aspa_dupla = False
    pular_ate = -1
    for achado in _ESPECIAIS.finditer(texto):
        i = achado.start()
        if i < pular_ate:
            continue
        c = texto[i]
        if c == "\\":
            pular_ate = i + 2
            continue
        if c == "'" and not aspa_dupla:
            aspa_simples = not aspa_simples
        elif c == '"' and not aspa_simples:
            aspa_dupla = not aspa_dupla
        elif c == ";" and not aspa_simples and not aspa_dupla:
            yield texto[inicio:i + 1]
            inicio = i + 1
    if inicio < len(texto):
        yield texto[inicio:]


def blocos(texto: str) -> Iterator[str]:
    """Blocos sem o espaço em branco das bordas, como `parse_cql_blocks`."""
    for seg in segmentos(texto):
        b = seg.strip()
        if b:
            yield b


def eh_bloco_legado(bloco: str) -> bool:
    linhas = bloco.split("\n", 2)
    return len(linhas) > 1 and linhas[1].startswith(_ABRE_LEGADO)


def ler_bloco_legado(bloco: str) -> tuple:
    """Decompõe um bloco de metadado do gerador antigo.

    Returns
    -------
    tuple
        ``(subtree_name, terminal_hash, prop_body, linhas_feature)`` com
        `terminal_hash` como texto e `linhas_feature` literal.

    Raises
    ------
    ValueError
        Se o bloco não tiver exatamente a forma que o gerador antigo emitia.
    """
    linhas = bloco.split("\n")
    m = _RE_MATCH_SUBTREE.match(linhas[0])
    if (len(linhas) < 4 or not m or linhas[0] != linha_match_subtree(m.group(1))
            or not linhas[1].startswith(_ABRE_LEGADO) or not linhas[1].endswith("})")
            or linhas[2] != _LIGACAO_LEGADA or linhas[-1] != ";"):
        raise ValueError(f"bloco de Metadata fora do formato esperado: {bloco[:200]!r}")
    prop_body = linhas[1][len(_ABRE_LEGADO):-2]
    h = _RE_HASH_FINAL.search(prop_body)
    if not h:
        raise ValueError(f"Metadata sem terminal_hash: {linhas[1][:200]!r}")
    linhas_feature = linhas[3:-1]
    for linha in linhas_feature:
        if not (_RE_FEATURE.match(linha) or _RE_QUALIFIER.match(linha)):
            raise ValueError(f"linha inesperada num bloco de Metadata: {linha[:200]!r}")
    return m.group(1), h.group(1), prop_body, linhas_feature


def projeto_do_texto(texto: str) -> str:
    """Projeto gravado nos nós `Tree` do arquivo; exige um único valor."""
    projetos = set(_RE_PROJETO.findall(texto))
    if len(projetos) != 1:
        raise ValueError(
            f"esperado exatamente um `project` nos nós Tree, achado {sorted(projetos)} — "
            "arquivo anterior à rotulagem? ver rotular_projeto_cql_legado.py")
    return projetos.pop()


def reparar_texto(texto: str, project_name: str) -> Iterator[str]:
    """Reescreve os blocos de metadado no formato do gerador corrigido.

    Gera pedaços de texto; ``"".join(...)`` é o arquivo reparado. Rodar sobre
    um arquivo já reparado não muda nada (não há mais bloco legado).

    Raises
    ------
    ValueError
        Se duas ocorrências completas do mesmo terminal divergirem — nesse
        caso deduplicar perderia informação, e o reparo não decide sozinho.
    """
    dedup = DeduplicadorMetadata(project_name)
    for seg in segmentos(texto):
        bloco = seg.strip()
        if not bloco or not eh_bloco_legado(bloco):
            yield seg
            continue
        prefixo = seg[:len(seg) - len(seg.lstrip())]
        subtree_name, h, prop_body, linhas_feature = ler_bloco_legado(bloco)
        yield prefixo + "\n".join(
            dedup.instrucoes(subtree_name, h, prop_body, linhas_feature))
        if dedup.conflitos:
            raise ValueError(f"terminal_hash {dedup.conflitos[0]} com registros divergentes")


# ---------------------------------------------------------------------------
# Verificação: lê original e reparado de forma independente do deduplicador
# e confere que só a duplicação mudou.
# ---------------------------------------------------------------------------

_RE_MERGE_META = re.compile(
    r"^MERGE \(m:Metadata \{uid: '" + re.escape(USER_PLACEHOLDER)
    + r"', project: (\"(?:[^\"\\]|\\.)*\"), terminal_hash: (-?\d+)\}\)$")


def contar_original(texto: str) -> dict:
    """Contagens e conteúdo por terminal de um `.cql` no formato antigo."""
    outros, ocorrencias = [], []
    completo = {}
    n_feat = n_qual = 0
    for b in blocos(texto):
        if not eh_bloco_legado(b):
            outros.append(b)
            continue
        sub, h, prop_body, feats = ler_bloco_legado(b)
        ocorrencias.append((sub, h))
        n_feat += sum(1 for x in feats if _RE_FEATURE.match(x))
        n_qual += sum(1 for x in feats if _RE_QUALIFIER.match(x))
        if DeduplicadorMetadata.eh_completo(prop_body, feats):
            conteudo = (prop_body, tuple(feats))
            if completo.setdefault(h, conteudo) != conteudo:
                raise ValueError(f"terminal_hash {h}: ocorrências completas divergentes")
    return {"outros": outros, "ocorrencias": ocorrencias, "completo": completo,
            "Metadata": len(ocorrencias), "Feature": n_feat, "Qualifier": n_qual,
            "HAS_METADATA": len(ocorrencias), "terminais": len({h for _, h in ocorrencias})}


def contar_reparado(texto: str, project_name: str) -> dict:
    """Contagens de um `.cql` no formato novo, lido sem o deduplicador."""
    outros, ligacoes = [], []
    criacoes = collections.defaultdict(list)
    projeto_json = json.dumps(project_name)
    n_feat = n_qual = 0
    for b in blocos(texto):
        linhas = b.split("\n")
        if len(linhas) >= 2 and _RE_MERGE_META.match(linhas[1]) and linhas[0].startswith("MATCH (child:Subtree"):
            m = _RE_MERGE_META.match(linhas[1])
            sub = _RE_MATCH_SUBTREE.match(linhas[0]).group(1)
            assert linhas[2:] == ["CREATE (child)-[:HAS_METADATA]->(m)", ";"], b[:300]
            assert m.group(1) == projeto_json, m.group(1)
            ligacoes.append((sub, m.group(2)))
        elif _RE_MERGE_META.match(linhas[0]):
            m = _RE_MERGE_META.match(linhas[0])
            assert m.group(1) == projeto_json, m.group(1)
            assert linhas[1].startswith("SET m += {") and linhas[1].endswith("}"), linhas[1][:200]
            assert linhas[-1] == ";"
            feats = linhas[2:-1]
            if feats:
                assert feats[0] == "WITH m WHERE NOT EXISTS { (m)-[:HAS_FEATURE]->() }", feats[0]
                feats = feats[1:]
            n_feat += sum(1 for x in feats if _RE_FEATURE.match(x))
            n_qual += sum(1 for x in feats if _RE_QUALIFIER.match(x))
            criacoes[m.group(2)].append((linhas[1][len("SET m += {"):-1], tuple(feats)))
        else:
            assert "Metadata" not in b.split("\n", 1)[0], b[:200]
            outros.append(b)
    return {"outros": outros, "ligacoes": ligacoes, "criacoes": criacoes,
            "Metadata": len(criacoes), "Feature": n_feat, "Qualifier": n_qual,
            "HAS_METADATA": len(ligacoes), "terminais": len({h for _, h in ligacoes})}


def verificar(original: dict, reparado: dict) -> list:
    """Lista de divergências (vazia = reparo correto).

    Invariantes: (1) blocos que não são de metadado idênticos e na mesma
    ordem; (2) a sequência de ligações (subárvore, terminal) é a mesma das
    ocorrências originais — uma aresta HAS_METADATA por ocorrência; (3) um
    `Metadata` por terminal; (4) o registro completo de cada terminal aparece
    numa única criação, idêntico ao original, e as demais criações do mesmo
    terminal só têm identidade e nenhuma feature.
    """
    erros = []
    if original["outros"] != reparado["outros"]:
        erros.append("blocos não-Metadata divergem")
    if original["ocorrencias"] != reparado["ligacoes"]:
        erros.append("sequência de ligações HAS_METADATA diverge das ocorrências originais")
    if set(reparado["criacoes"]) != {h for _, h in original["ocorrencias"]}:
        erros.append("conjunto de terminais com Metadata diverge")
    for h, criacoes in reparado["criacoes"].items():
        esperado = original["completo"].get(h)
        completas = [c for c in criacoes if DeduplicadorMetadata.eh_completo(c[0], list(c[1]))]
        vazias = [c for c in criacoes if c not in completas]
        if esperado is None:
            if completas:
                erros.append(f"{h}: criação completa sem ocorrência completa no original")
        elif completas != [esperado]:
            erros.append(f"{h}: registro completo ausente, duplicado ou alterado")
        if any(c[1] for c in vazias):
            erros.append(f"{h}: criação vazia com features")
    return erros


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("cql")
    ap.add_argument("--saida", help="padrão: <arquivo sem .cql>.reparado.cql")
    args = ap.parse_args(argv)

    origem = os.path.abspath(args.cql)
    saida = args.saida or re.sub(r"\.cql$", "", origem) + ".reparado.cql"
    if os.path.abspath(saida) == origem:
        ap.error("a saída não pode sobrescrever o original")

    t0 = time.time()
    with open(origem, encoding="utf-8") as f:
        texto = f.read()
    projeto = projeto_do_texto(texto)
    projeto_caminho = project_name_from_output_path(os.path.dirname(origem))
    if projeto != projeto_caminho:
        print(f"AVISO: project no arquivo ({projeto}) difere do diretório ({projeto_caminho})")

    temporario = saida + ".parcial"
    with open(temporario, "w", encoding="utf-8") as f:
        for pedaco in reparar_texto(texto, projeto):
            f.write(pedaco)

    with open(temporario, encoding="utf-8") as f:
        reparado_txt = f.read()
    antes = contar_original(texto)
    depois = contar_reparado(reparado_txt, projeto)
    erros = verificar(antes, depois)
    if erros:
        print("VERIFICAÇÃO FALHOU — saída mantida em", temporario)
        for e in erros[:20]:
            print("  -", e)
        return 1
    os.replace(temporario, saida)

    print(f"projeto: {projeto}")
    print(f"{'':14}{'antes':>12}{'depois':>12}")
    for k in ("Metadata", "Feature", "Qualifier", "HAS_METADATA", "terminais"):
        print(f"{k:14}{antes[k]:>12}{depois[k]:>12}")
    print(f"bytes         {os.path.getsize(origem):>12}{os.path.getsize(saida):>12}")
    print(f"blocos não-Metadata inalterados: {len(antes['outros'])}")
    print(f"saída: {saida}")
    print(f"tempo {time.time() - t0:.1f}s, pico de memória "
          f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
