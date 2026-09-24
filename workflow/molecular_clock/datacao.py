"""
Datação por mínimos quadrados (LSD2, embutido no IQ-TREE) sobre árvore já inferida.

Porte para produção da lógica validada em E12.1-E12.7
(`docs/science/15-datacao-leve-e12.md`). O script de análise
`docs/science/scripts/datacao_leve_lsd2.py` **não é importado**: ele mora no
repositório pai (`docs/`), e o submódulo `BioComp_UFF` tem ciclo próprio — um
`import` atravessando essa fronteira quebraria quem clonasse só o submódulo. O
script continua existindo como **referência**: a tabela de diff de
`docs/science/16-relogio-molecular-producao.md` roda os dois sobre a mesma
entrada e confere que dão o mesmo número.

O que mudou em relação ao script, e por quê (cada item está medido no
documento 16):

- **`-T 1 -seed S` por padrão.** Com ``-T 4`` a reotimização de comprimento de
  ramo muda o **ponto** estimado (tMRCA 1685,1 → 1682,4 no Zika-21) e, mesmo
  com semente fixa, duas corridas diferem; com ``-T 1 -seed S`` são idênticas
  até o intervalo. É o D17 do RAxML-NG em outro binário.
- **Data de clado pela bipartição canônica** (`canonical_item_id`), não pelo
  conjunto de pontas cru nem pelo rótulo interno do Biopython (D3, DEC-123).
- **Mínimo de datas contado no grupo interno**: com ``--date-no-outgroup`` o
  grupo externo sai da árvore de tempo, e data de ponta removida não calibra.
- **Sinal temporal é portão**, não só diagnóstico: inclinação raiz-ponta ≤ 0
  no grupo interno recusa a datação (E12 §4).
"""

from __future__ import annotations

import datetime as _dt
import os
import re
import subprocess
from typing import Dict, FrozenSet, Iterable, List, Optional, Set, Tuple

from workflow.stability.clade_identity import (
    canonical_bipartition,
    canonical_digest,
    canonical_item_id,
    strip_accession_version,
)

__all__ = [
    "MESES",
    "extrair_collection_date",
    "converter_data",
    "data_decimal",
    "data_anotada_para_decimal",
    "ler_rotulos_do_alinhamento",
    "ler_arvore_de_origem",
    "verificar_grupo_externo",
    "escrever_entrada",
    "escrever_arquivo_datas",
    "sinal_temporal",
    "minimo_de_pontas_datadas",
    "rodar_lsd2",
    "interpretar_linha_de_resultado",
    "datas_por_clado",
    "biparticoes_do_grupo_interno",
]

MESES = {
    "JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
    "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12,
}


# --------------------------------------------------------------------------- #
# Datas de coleta
# --------------------------------------------------------------------------- #
def extrair_collection_date(gb_path: str) -> Dict[str, Optional[str]]:
    """
    Lê `collection_date` de cada registro de um GenBank flat file.

    Parameters
    ----------
    gb_path : str
        `out/outputs/raw_data_sequences.gb`, gravado pelo pipeline.

    Return
    ------
    dict
        ``accession.version -> collection_date`` cru, ou ``None`` quando o
        registro não traz o qualificador.
    """
    with open(gb_path, errors="ignore") as fh:
        texto = fh.read()
    saida: Dict[str, Optional[str]] = {}
    for reg in texto.split("\nLOCUS"):
        acc = re.search(r"^\s*VERSION\s+(\S+)", reg, re.M)
        if not acc:
            continue
        cd = re.search(r'/collection_date="([^"]+)"', reg)
        saida[acc.group(1)] = cd.group(1).strip() if cd else None
    return saida


def converter_data(crua: str, ano_pontual: bool = False) -> Tuple[str, str]:
    """
    Converte `collection_date` para o formato do ``--date`` do IQ-TREE.

    Convenção (E12 §2 — regra 5 do CLAUDE.md, nunca inventar precisão):
    dia → ``YYYY-MM-DD``; mês → ``YYYY-MM`` (o IQ-TREE põe no meio do mês);
    ano-só → intervalo ``YYYY:YYYY+1``, ou ``YYYY`` com `ano_pontual`.

    Return
    ------
    tuple
        ``(valor_para_lsd2, precisao)``, `precisao` em {``dia``, ``mes``, ``ano``}.

    Raises
    ------
    ValueError
        Formato não reconhecido — erro explícito, nunca palpite.
    """
    s = crua.strip()
    if re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", s):
        return s, "dia"
    m = re.fullmatch(r"(\d{1,2})-([A-Za-z]{3})-(\d{4})", s)
    if m and m.group(2).upper() in MESES:
        dia, mes, ano = int(m.group(1)), MESES[m.group(2).upper()], m.group(3)
        return f"{ano}-{mes:02d}-{dia:02d}", "dia"
    m = re.fullmatch(r"([A-Za-z]{3})-(\d{4})", s)
    if m and m.group(1).upper() in MESES:
        return f"{m.group(2)}-{MESES[m.group(1).upper()]:02d}", "mes"
    if re.fullmatch(r"(\d{4})-(\d{2})", s):
        return s, "mes"
    m = re.fullmatch(r"(\d{4})", s)
    if m:
        ano = int(m.group(1))
        return (str(ano) if ano_pontual else f"{ano}:{ano + 1}"), "ano"
    raise ValueError(f"formato de collection_date não reconhecido: {crua!r}")


def _fracao_do_ano(ano: int, mes: int, dia: int) -> float:
    if ano < 1:
        # `datetime` não representa o ano 0 nem negativos, e o LSD2 devolve
        # esses anos quando a raiz é absurda (E12 §5: −5 731 390). O valor é
        # lixo de qualquer jeito; aqui só não pode derrubar a leitura.
        return ano + (mes - 1) / 12 + (dia - 1) / 365.25
    ini = _dt.date(ano, 1, 1)
    fim = _dt.date(ano + 1, 1, 1) if ano < 9999 else _dt.date(ano, 12, 31)
    return ano + (_dt.date(ano, mes, dia) - ini).days / max((fim - ini).days, 1)


def data_decimal(valor_lsd2: str) -> float:
    """
    Valor do arquivo de datas em ano decimal — só para a regressão diagnóstica.

    O intervalo ``YYYY:YYYY+1`` vira o **meio** do intervalo; a datação em si
    recebe o intervalo, não o meio.
    """
    if ":" in valor_lsd2:
        return int(valor_lsd2.split(":")[0]) + 0.5
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", valor_lsd2)
    if m:
        return _fracao_do_ano(*map(int, m.groups()))
    m = re.fullmatch(r"(\d{4})-(\d{2})", valor_lsd2)
    if m:
        ano, mes = map(int, m.groups())
        return _fracao_do_ano(ano, mes, 15)
    return float(valor_lsd2)


def data_anotada_para_decimal(v) -> Optional[float]:
    """
    Anotação de data do LSD2 (``1952-08-21``, ``1963.51``) em ano decimal.

    O LSD2 escreve data de calendário quando a entrada tem precisão de dia e
    decimal quando não tem; os dois casos têm de virar o mesmo tipo, ou o
    consumidor casa só metade dos nós.
    """
    if isinstance(v, (list, tuple)):
        v = v[0] if v else None
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip().strip('"').strip("{}")
    m = re.fullmatch(r"(-?\d{1,7})-(\d{2})-(\d{2})", s)
    if m:
        return _fracao_do_ano(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    try:
        return float(s)
    except ValueError:
        return None


# --------------------------------------------------------------------------- #
# Árvore e alinhamento
# --------------------------------------------------------------------------- #
def ler_rotulos_do_alinhamento(aln_path: str) -> List[str]:
    """Rótulos de um alinhamento FASTA, na ordem do arquivo."""
    rotulos = []
    with open(aln_path) as fh:
        for linha in fh:
            if linha.startswith(">"):
                rotulos.append(linha[1:].strip().split()[0])
    return rotulos


def ler_arvore_de_origem(nexus_path: str):
    """
    Lê uma árvore de `out/Trees` **isoladamente**, como não enraizada.

    Namespace próprio por arquivo: vários Nexus num `TaxonNamespace`
    compartilhado fazem o dendropy abortar (D13).
    """
    import dendropy

    esquema = "nexus" if nexus_path.lower().endswith((".nexus", ".nex")) else "newick"
    arv = dendropy.Tree.get(path=nexus_path, schema=esquema, preserve_underscores=True,
                            terminating_semicolon_required=False)
    arv.is_rooted = False
    for e in arv.postorder_edge_iter():
        if e.length is None:
            e.length = 0.0
    return arv


def _pontas_normalizadas(arv) -> Dict[str, object]:
    return {strip_accession_version(l.taxon.label): l for l in arv.leaf_node_iter()}


def verificar_grupo_externo(arv, declarado: Iterable[str]) -> Dict:
    """
    Confere o grupo externo declarado contra a árvore, **como não enraizada**.

    Monofilia aqui é existir uma aresta que separe exatamente o grupo externo
    do resto — é o que o LSD2 precisa para enraizar nela. Perguntar pela
    monofilia na raiz escrita no arquivo seria perguntar pela convenção de
    escrita da ferramenta (D3).

    Return
    ------
    dict
        ``encontrado``, ``ausente``, ``monofiletico`` e, se não der para
        enraizar, ``motivo``.
    """
    declarado_n = frozenset(strip_accession_version(x) for x in declarado if x and x.strip())
    pontas = _pontas_normalizadas(arv)
    todos = frozenset(pontas)
    encontrado = declarado_n & todos
    ausente = declarado_n - todos
    rel: Dict = {"declarado": len(declarado_n), "encontrado": sorted(encontrado),
                 "ausente": sorted(ausente), "monofiletico": None}
    if not declarado_n:
        rel["motivo"] = ("grupo externo não declarado: o enraizamento da datação tem de "
                         "ser declarado, nunca inferido (E12 §9.1; rooting.py regra 1)")
        return rel
    if not encontrado:
        rel["motivo"] = "nenhum táxon do grupo externo declarado está nesta árvore"
        return rel
    if len(todos - encontrado) < 3:
        rel["motivo"] = "o grupo externo deixa menos de 3 pontas no grupo interno"
        return rel
    if len(encontrado) == 1:
        rel["monofiletico"] = True
        return rel
    lados = set()
    for nd in arv.postorder_node_iter():
        if nd.parent_node is None:
            continue
        lado = frozenset(strip_accession_version(l.taxon.label) for l in nd.leaf_iter())
        lados.add(lado)
        lados.add(todos - lado)
    rel["monofiletico"] = encontrado in lados
    if not rel["monofiletico"]:
        rel["motivo"] = ("grupo externo não é monofilético nesta árvore: nenhuma aresta o "
                         "separa do grupo interno, e o LSD2 abortaria ('The outgroups do "
                         "not form a monophyletic in the tree') — E12 §8.6")
    return rel


def escrever_entrada(arv, rotulos_aln: List[str], destino_nwk: str) -> Dict[str, str]:
    """
    Escreve a árvore em Newick com os rótulos **do alinhamento**.

    IQ-TREE e RAxML-NG truncam ``NC_008030.1`` em ``NC_008030.`` ao gravar a
    árvore (D13); o alinhamento guarda a grafia íntegra, e o ``-te`` do IQ-TREE
    exige que as duas batam. O casamento é pelo acesso normalizado.

    Return
    ------
    dict
        ``acesso normalizado -> rótulo do alinhamento``.

    Raises
    ------
    ValueError
        Ponta da árvore sem sequência no alinhamento.
    """
    por_norm = {strip_accession_version(r): r for r in rotulos_aln}
    faltando = []
    for folha in arv.leaf_node_iter():
        n = strip_accession_version(folha.taxon.label)
        if n not in por_norm:
            faltando.append(folha.taxon.label)
        else:
            folha.taxon.label = por_norm[n]
    if faltando:
        raise ValueError(f"{len(faltando)} ponta(s) da árvore sem sequência no alinhamento "
                         f"(ex.: {faltando[:3]})")
    arv.write(path=destino_nwk, schema="newick", suppress_rooting=True,
              unquoted_underscores=True, suppress_internal_node_labels=True)
    return {strip_accession_version(l.taxon.label): l.taxon.label
            for l in arv.leaf_node_iter()}


def escrever_arquivo_datas(crus: Dict[str, Optional[str]], rotulos: Dict[str, str],
                           destino: str, grupo_interno: FrozenSet[str],
                           ano_pontual: bool = False) -> Dict:
    """
    Escreve ``rótulo<TAB>data`` para toda ponta da árvore que tem data.

    Ponta sem `collection_date` é **omitida** (o LSD2 estima a data dela):
    "sem data" nunca é um número. O arquivo inclui o grupo externo, como na
    corrida validada de E12; a contagem que decide o mínimo é só a do grupo
    interno, porque o grupo externo sai da árvore de tempo.

    Parameters
    ----------
    crus : dict
        Saída de `extrair_collection_date`.
    rotulos : dict
        ``acesso normalizado -> rótulo`` das pontas da árvore.
    grupo_interno : frozenset
        Acessos normalizados do grupo interno.
    """
    crus_n = {strip_accession_version(k): v for k, v in crus.items()}
    linhas: List[str] = []
    prec = {"dia": 0, "mes": 0, "ano": 0}
    sem_data, sem_registro, nao_reconhecidas = [], [], []
    decimais: Dict[str, float] = {}
    for norm, rot in sorted(rotulos.items()):
        if norm not in crus_n:
            sem_registro.append(rot)
            continue
        valor = crus_n[norm]
        if valor is None:
            sem_data.append(rot)
            continue
        try:
            conv, p = converter_data(valor, ano_pontual=ano_pontual)
        except ValueError:
            nao_reconhecidas.append([rot, valor])
            continue
        prec[p] += 1
        linhas.append(f"{rot}\t{conv}")
        decimais[norm] = data_decimal(conv)
    with open(destino, "w") as fh:
        fh.write("\n".join(linhas) + "\n")
    return {
        "pontas_datadas": len(linhas),
        "pontas_datadas_no_grupo_interno": sum(1 for n in decimais if n in grupo_interno),
        "precisao": prec,
        "sem_collection_date": sem_data,
        "sem_registro_genbank": sem_registro,
        "collection_date_nao_reconhecido": nao_reconhecidas,
        "convencao_ano_so": "pontual (1-jan)" if ano_pontual else "intervalo YYYY:YYYY+1",
        "_decimais": decimais,
    }


def minimo_de_pontas_datadas(n_pontas: int) -> int:
    """``max(3, 10 % das pontas)`` — a regra de E12 §8.1, contra o caso VARV-49."""
    return max(3, int(0.1 * n_pontas))


# --------------------------------------------------------------------------- #
# Sinal temporal
# --------------------------------------------------------------------------- #
def _regressao(xs: List[float], ys: List[float]) -> Dict:
    n = len(xs)
    if n < 3:
        return {"n": n, "inclinacao_subs_sitio_ano": None, "r2": None,
                "data_da_raiz_x_intercepto": None,
                "aviso": "menos de 3 pontas datadas — regressão indefinida"}
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx == 0 or syy == 0:
        return {"n": n, "inclinacao_subs_sitio_ano": None, "r2": None,
                "data_da_raiz_x_intercepto": None,
                "aviso": "variância nula — regressão indefinida"}
    s = sxy / sxx
    return {
        "n": n,
        "inclinacao_subs_sitio_ano": s,
        "r2": sxy ** 2 / (sxx * syy),
        "data_da_raiz_x_intercepto": (-(my - s * mx) / s) if s != 0 else None,
    }


def sinal_temporal(arv, decimais: Dict[str, float], grupo_externo: FrozenSet[str]) -> Dict:
    """
    Regressão raiz-ponta no grupo interno, com a raiz na aresta do grupo externo.

    É o teste de sinal temporal que precede qualquer afirmação de relógio
    (E12 §4, no espírito do TempEst). A árvore recebida não é modificada.
    """
    import copy

    a = copy.deepcopy(arv)
    a.is_rooted = False
    todos = frozenset(_pontas_normalizadas(a))
    aresta = None
    for nd in a.postorder_node_iter():
        if nd.parent_node is None:
            continue
        lado = frozenset(strip_accession_version(l.taxon.label) for l in nd.leaf_iter())
        if lado == grupo_externo or todos - lado == grupo_externo:
            aresta = nd.edge
            break
    if aresta is None:
        return {"n": 0, "inclinacao_subs_sitio_ano": None, "r2": None,
                "data_da_raiz_x_intercepto": None,
                "aviso": "aresta do grupo externo não encontrada"}
    a.reroot_at_edge(aresta, update_bipartitions=False)
    # O reenraizamento cria o nó da raiz com aresta sem comprimento.
    for e in a.postorder_edge_iter():
        if e.length is None:
            e.length = 0.0
    a.calc_node_root_distances()
    pares = [
        (decimais[n], l.root_distance)
        for l in a.leaf_node_iter()
        if (n := strip_accession_version(l.taxon.label)) in decimais and n not in grupo_externo
    ]
    return _regressao([p[0] for p in pares], [p[1] for p in pares])


# --------------------------------------------------------------------------- #
# LSD2
# --------------------------------------------------------------------------- #
def interpretar_linha_de_resultado(linha: Optional[str]) -> Dict:
    """
    Decompõe a linha ``rate … , tMRCA …`` do LSD2.

    ``rate 0.0005 [0.0004; 0.0006], tMRCA 1944-12-05 [1932-04-16; 1951-06-21], …``.
    ``:NA`` na linha é o LSD2 dizendo que a solução não é única (E12 §8.2).
    """
    vazio = {"taxa_subs_por_sitio_por_ano": None, "intervalo_da_taxa": None,
             "tmrca": None, "tmrca_lsd2": None, "intervalo_do_tmrca": None,
             "solucao_unica": False}
    if not linha:
        return vazio
    m = re.match(r"\s*rate ([^,\s\[]+)\s*(\[[^\]]*\])?,\s*tMRCA ([^\s,\[]+)\s*(\[[^\]]*\])?", linha)
    if not m:
        return vazio
    taxa_b, ic_taxa_b, tmrca_b, ic_tmrca_b = m.groups()

    def _intervalo(bruto, conv):
        if not bruto:
            return None
        partes = [p.strip() for p in bruto.strip("[]").split(";")]
        if len(partes) != 2:
            return None
        vals = [conv(p) for p in partes]
        return None if any(v is None for v in vals) else vals

    def _float(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None

    unica = ":NA" not in linha and "NA:" not in linha
    return {
        "taxa_subs_por_sitio_por_ano": _float(taxa_b) if unica else None,
        "intervalo_da_taxa": _intervalo(ic_taxa_b, _float) if unica else None,
        "tmrca": data_anotada_para_decimal(tmrca_b) if unica else None,
        "tmrca_lsd2": tmrca_b,
        "intervalo_do_tmrca": _intervalo(ic_tmrca_b, data_anotada_para_decimal) if unica else None,
        "solucao_unica": unica,
    }


def rodar_lsd2(binario: str, alinhamento: str, newick: str, datas: str, prefixo: str,
               grupo_externo_rotulos: List[str], modelo: str = "GTR+G", ci: int = 100,
               threads: int = 1, semente: Optional[int] = None,
               reotimizar_ramos: bool = True, registrar=None) -> Dict:
    """
    Chama ``iqtree --date`` com topologia fixa, grupo externo excluído da árvore de tempo.

    Parameters
    ----------
    grupo_externo_rotulos : list of str
        Rótulos **do alinhamento**. Vão em ``-o`` com ``--date-no-outgroup``:
        o LSD2 enraíza pelo grupo externo e o remove da árvore de tempo — o
        protocolo adotado em E12 §9.1 (os outros dois dão tMRCA 2521 e 5774).
    threads, semente
        ``-T 1 -seed S`` é o que torna a corrida reprodutível (ver o docstring
        do módulo). Mais threads é permitido, e fica registrado.
    registrar : callable, optional
        ``registrar(cmd, saida)`` — o controlador passa `tool_runs.registrar`.
        Chamado **antes** de rodar, como o resto do pipeline (tool_runs, Notes).

    Notes
    -----
    ``--keep-ident`` é obrigatório: sem ele o IQ-TREE descarta sequências
    idênticas e cada ponta descartada perde a calibração em silêncio (18 de 60
    no oráculo de E12). Deduplicação é decisão de método (B-10), não efeito
    colateral da datação.
    """
    cmd = [
        binario, "-s", alinhamento, "-te", newick, "-m", modelo, "--keep-ident",
        "--date", datas, "--prefix", prefixo, "-T", str(threads), "-redo",
        "-o", ",".join(sorted(grupo_externo_rotulos)), "--date-no-outgroup",
    ]
    if semente is not None:
        cmd += ["-seed", str(semente)]
    if not reotimizar_ramos:
        cmd.append("-blfix")
    if ci:
        cmd += ["--date-ci", str(ci)]
    timetree = f"{prefixo}.timetree.nex"
    if registrar is not None:
        registrar(cmd, timetree)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    saida = proc.stdout
    linhas = [l.strip() for l in saida.splitlines() if l.strip().startswith("rate ")]
    bruta = linhas[-1] if linhas else None
    colapso = re.search(r"([\d.]+)% internal branches were collapsed", saida)
    gaps = re.search(r"WARNING: (\d+) sequences contain more than 50% gaps", saida)
    erro = re.findall(r"^(?:ERROR|Error):?\s*(.+)$", saida + "\n" + proc.stderr, re.M)
    res = {
        "comando": cmd,
        "retorno": proc.returncode,
        "linha_de_resultado_lsd2": bruta,
        "aviso_sem_informacao_de_data_suficiente":
            "not enough input date information" in saida,
        "percentual_de_ramos_internos_colapsados_lsd2":
            float(colapso.group(1)) if colapso else None,
        "sequencias_com_mais_de_50pc_gap": int(gaps.group(1)) if gaps else 0,
        "erro_lsd2": erro[-1].strip() if erro else None,
        "timetree": timetree if os.path.exists(timetree) else None,
    }
    res.update(interpretar_linha_de_resultado(bruta))
    if res["aviso_sem_informacao_de_data_suficiente"]:
        res["solucao_unica"] = False
    return res


# --------------------------------------------------------------------------- #
# Datas por clado, pela bipartição canônica
# --------------------------------------------------------------------------- #
def _texto_timetree_legivel(caminho: str) -> str:
    """
    Tira as aspas de ``CI_date="{a,b}"`` antes do dendropy ler.

    O parser de comentário do dendropy corta o valor na vírgula de dentro das
    aspas e devolve ``'"{1951.41'`` (medido; o oráculo de E12 contornava
    reconstruindo de `CI_height`). Sem aspas, ``{a,b}`` vira lista.
    """
    with open(caminho) as fh:
        texto = fh.read()
    return re.sub(r'CI_date="\{([^}]*)\}"', r"CI_date={\1}", texto)


def datas_por_clado(timetree: str, todos: FrozenSet[str]) -> Tuple[Dict[str, Dict], Dict]:
    """
    Data de cada nó interno do timetree, indexada pelo ``clade_id`` canônico.

    Parameters
    ----------
    timetree : str
        ``*.timetree.nex`` do LSD2 (enraizado, sem o grupo externo).
    todos : frozenset of str
        Acessos normalizados de **todas** as pontas da árvore de origem em
        `out/Trees`, grupo externo incluído. É contra esse conjunto que a
        bipartição é calculada — o mesmo que `suporte_de_ramo.py` usa ao ler
        o mesmo arquivo —, e é por isso que o `clade_id` daqui casa com o de
        lá, com o de `metadata.json` e com o do FPMax.

    Return
    ------
    tuple
        ``(datas, contagem)``. `datas` mapeia ``str(clade_id)`` (chave de JSON
        é string; o inteiro vai também dentro do valor) em::

            {"clade_id": int, "digest": str, "n_taxa": int,
             "data_estimada": float, "data_estimada_lsd2": str,
             "intervalo": [lo, hi] | None,
             "n_taxa_do_clado_datado": int,
             "clado_datado_e_o_lado_canonico": bool}

    Notes
    -----
    **O que a data significa.** A árvore de tempo é enraizada; cada aresta
    liga um nó ao pai, e a data que o LSD2 dá é a do nó do lado das pontas —
    o MRCA das pontas desse lado. A bipartição canônica é o lado **menor**
    (`03-metricas §2.2`), que nem sempre é esse: quando não é,
    `clado_datado_e_o_lado_canonico` é ``False`` e a data continua sendo a do
    MRCA do lado **oposto** ao da raiz. `n_taxa` é o tamanho do lado
    canônico, como em `suporte_de_ramo.py`; `n_taxa_do_clado_datado`, o do
    clado cuja origem foi datada.

    **Ramo colapsado não tem data.** O LSD2 colapsa ramo interno com menos de
    meia substituição esperada (47 % em ZIKV, E12 §8.3); a bipartição some da
    árvore de tempo e, portanto, daqui. Ausência não é zero (regra 5).
    """
    import dendropy

    arv = dendropy.Tree.get(data=_texto_timetree_legivel(timetree), schema="nexus",
                            preserve_underscores=True, extract_comment_metadata=True)
    datas: Dict[str, Dict] = {}
    colisoes = 0
    nos_internos = 0
    for no in arv.postorder_node_iter():
        if no.is_leaf():
            continue
        nos_internos += 1
        lado = frozenset(strip_accession_version(l.taxon.label)
                         for l in no.leaf_iter() if l.taxon is not None)
        bip = canonical_bipartition(lado, todos)
        if bip is None:
            continue
        meta = no.annotations.values_as_dict() if no.annotations else {}
        bruto = meta.get("date")
        data = data_anotada_para_decimal(bruto)
        if data is None:
            continue
        intervalo = None
        ci = meta.get("CI_date")
        if isinstance(ci, (list, tuple)) and len(ci) == 2:
            vals = [data_anotada_para_decimal(x) for x in ci]
            if all(v is not None for v in vals):
                intervalo = sorted(vals)
        cid = canonical_item_id(bip)
        chave = str(cid)
        if chave in datas:
            # Com o grupo externo excluído não acontece (os filhos da raiz
            # separam bipartições distintas); fica contado para não sumir.
            colisoes += 1
            continue
        datas[chave] = {
            "clade_id": cid,
            "digest": canonical_digest(bip),
            "n_taxa": len(bip),
            "data_estimada": round(data, 4),
            "data_estimada_lsd2": str(bruto[0] if isinstance(bruto, (list, tuple)) else bruto),
            "intervalo": [round(v, 4) for v in intervalo] if intervalo else None,
            "n_taxa_do_clado_datado": len(lado),
            "clado_datado_e_o_lado_canonico": bip == lado,
        }
    return datas, {"nos_internos_na_arvore_de_tempo": nos_internos,
                   "clados_datados": len(datas), "colisoes_de_biparticao": colisoes}


def biparticoes_do_grupo_interno(arv, grupo_externo: FrozenSet[str]) -> Set[FrozenSet[str]]:
    """
    Bipartições não triviais da árvore de origem que caem **dentro** do grupo interno.

    Mais a aresta grupo interno | grupo externo. É o denominador do colapso:
    toda bipartição daqui que falta em `datas_por_clado` foi colapsada pelo LSD2
    (ou já era ramo de comprimento zero na origem).
    """
    todos = frozenset(_pontas_normalizadas(arv))
    saida: Set[FrozenSet[str]] = set()
    for nd in arv.postorder_node_iter():
        if nd.parent_node is None:
            continue
        lado = frozenset(strip_accession_version(l.taxon.label) for l in nd.leaf_iter())
        for s in (lado, todos - lado):
            if s.isdisjoint(grupo_externo) or s == todos - grupo_externo:
                bip = canonical_bipartition(s, todos)
                if bip is not None:
                    saida.add(bip)
    return saida
