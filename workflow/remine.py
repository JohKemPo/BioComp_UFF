"""
Remineração com proveniência — refaz a mineração de um `out/` já executado,
sem reinferir as árvores.

Nasceu do D29: o item do FPMax gravado era o clado como escrito no arquivo, não
a bipartição canônica, e corrigir o minerador não corrige os artefatos que ele
já gravou. As árvores de `out/Trees` não mudam, então o que se refaz é só o que
depende da mineração: `all_results_fpmax.csv`, `metadata.json`, `metadata.csv`,
o `.cql` e os gráficos de suporte. `workflow.py` não serve: roda sempre o
`TreeBuilderController` antes, e reinferir mudaria o que não está em questão.

O trabalho acontece num diretório de preparo, e os arquivos novos só entram em
`out/` depois que a mineração terminou inteira — uma falha no meio não deixa o
projeto com metade dos artefatos de cada código. O nome do projeto é o do
diretório (`neo4jProcessing.project_name_from_output_path`), por isso o preparo
repete `<projeto>/out` e o `.cql` sai com o nome certo.

Três salvaguardas que o driver improvisado do parecer G9 não tinha:
  * nenhuma sequência é baixada do NCBI (o `.gb` já existente é a fonte e o
    download, que o abre com ``"w"``, levanta erro em vez de truncá-lo);
  * o metadado GenBank de cada terminal é o que o `metadata.json` anterior já
    tinha, não uma releitura do `.gb`. O `.gb` que sobra em `out/outputs/` é o
    da **última** rodada de download — `download_sequences` o reabre com
    ``"w"`` a cada árvore que traz acesso novo —, e em VARV-52 tem 3 dos 51
    registros. Relê-lo trocaria metadado de terminal por ``{"error": ...}`` em
    centenas de nós, e o D29 não tem nada a ver com isso;
  * o manifesto existente recebe o registro da remineração
    (`manifest.anexar_remineracao`), sem ser regravado do zero. Projeto sem
    manifesto é recusado: não se inventa proveniência.

Uso::

    python -m workflow.remine --out projects/<projeto>/out --motivo D29
"""

import argparse
import contextlib
import datetime
import decimal
import glob
import json
import os
import shutil
import sys
import tempfile
from typing import Dict, List, Optional
from unittest import mock

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import ijson
import mlxtend
import pandas
from Bio import Phylo

from workflow.utils import manifest
from workflow.utils.treeUtils import extract_terminais_list


_SEM_METADADO = object()


def _proibido(*args, **kwargs):
    raise RuntimeError("a remineração não baixa sequências: o raw_data_sequences.gb do projeto "
                       "é a fonte, e há um terminal que ele não cobre")


def _digests(arquivos: List[str], raiz: str) -> Dict[str, Optional[str]]:
    return {os.path.relpath(a, raiz): manifest.file_digest(a) for a in sorted(arquivos)}


def _arquivos_de(diretorio: str) -> List[str]:
    return [os.path.join(r, f) for r, _, fs in os.walk(diretorio) for f in fs]


def _sem_decimal(valor):
    """O `ijson` devolve `Decimal` para número não inteiro; o resto do pipeline espera `float`."""
    if isinstance(valor, decimal.Decimal):
        return float(valor)
    if isinstance(valor, dict):
        return {k: _sem_decimal(v) for k, v in valor.items()}
    if isinstance(valor, list):
        return [_sem_decimal(v) for v in valor]
    return valor


def _metadados_anteriores(caminho: str) -> Dict[str, Dict[str, dict]]:
    """
    ``{árvore: {terminal normalizado: metadata}}`` do `metadata.json` em disco.

    Lido em fluxo (os arquivos chegam a 3,2 GB) e reduzido a um registro por
    terminal e árvore: dentro de uma árvore o `.gb` não muda, então o metadado de
    um terminal é o mesmo em todo nó que o contém.

    Nada de ``use_float=True``: o backend C do `ijson` converte com ``strtod``, que
    obedece ao locale do processo, e o matplotlib com Qt (que o pipeline importa)
    troca o locale para o do sistema ao criar a primeira figura — em pt_BR,
    ``1.5`` vira ``1.0`` e o metadado preservado deixaria de ser o original. O
    caminho padrão (``Decimal``) lê o texto e não depende do locale.
    """
    anteriores: Dict[str, Dict[str, dict]] = {}
    if not os.path.isfile(caminho):
        return anteriores
    with open(caminho, "rb") as f:
        inicio = f.read(64).lstrip()
        f.seek(0)
        # modo OFST grava uma lista de grupos, cada um uma lista de árvores
        prefixo = "item.item" if inicio[1:].lstrip()[:1] == b"[" else "item"
        for bloco in ijson.items(f, prefixo):
            for arvore, subarvores in bloco.items():
                por_terminal = anteriores.setdefault(arvore, {})
                for subarvore in subarvores.values():
                    for terminal in subarvore["data_terminals"]:
                        por_terminal.setdefault(terminal["newick"], _sem_decimal(terminal["metadata"]))
    return anteriores


def _arvores(dados):
    """No modo OFST o controlador devolve uma lista de grupos, cada um uma lista de árvores."""
    for d in dados:
        if isinstance(d, list):
            yield from _arvores(d)
        else:
            yield d


def _instalar(origem: str, destino: str) -> None:
    """Copia e troca por `os.replace`: o arquivo antigo nunca fica pela metade."""
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    provisorio = destino + ".remine.tmp"
    shutil.copyfile(origem, provisorio)
    os.replace(provisorio, destino)


def remineirar(out_dir: str, motivo: str) -> Dict:
    """
    Remineira o projeto de `out_dir` e anexa o registro ao manifesto.

    Parameters
    ----------
    out_dir : str
        O diretório ``out/`` do projeto.
    motivo : str
        Por que se remineirou (ex. ``"D29"``); vai para o manifesto.

    Returns
    -------
    dict
        Resumo: padrões gravados, itens por árvore, o que mudou em `Subtrees/`
        e o caminho do registro.

    Raises
    ------
    FileNotFoundError
        Falta manifesto, `config_backup.json`, `out/Trees` ou o `.gb`.
    RuntimeError
        A mineração não gerou `all_results_fpmax.csv` (não aplicável).
    """
    from workflow.controller.subtreeBuilderController import SubtreeBuilderController
    from workflow.subtree_construction import builder as builder_mod

    out_dir = os.path.abspath(out_dir)
    raiz = os.path.dirname(out_dir)
    outputs = os.path.join(out_dir, "outputs")
    trees_dir = os.path.join(out_dir, "Trees")
    gb = os.path.join(outputs, "raw_data_sequences.gb")
    backup = os.path.join(outputs, "config_backup.json")
    for exigido in (os.path.join(outputs, manifest.MANIFEST_FILENAME), backup, gb):
        if not os.path.isfile(exigido):
            raise FileNotFoundError(
                f"{os.path.relpath(exigido, raiz)} não existe: sem ele a remineração não tem "
                "de onde partir nem onde registrar a proveniência.")
    if not os.path.isdir(trees_dir) or not os.listdir(trees_dir):
        raise FileNotFoundError(f"{os.path.relpath(trees_dir, raiz)} não existe ou está vazio.")

    with open(backup, encoding="utf-8") as fh:
        cfg = json.load(fh)["subtree_config"]

    entradas = _digests([os.path.join(trees_dir, f) for f in os.listdir(trees_dir)] + [gb, backup], raiz)

    nome = os.path.basename(raiz)
    preparo = tempfile.mkdtemp(prefix=".remine_", dir=os.path.dirname(raiz))
    preparo_out = os.path.join(preparo, nome, "out")
    os.makedirs(os.path.join(preparo_out, "outputs", "Plots"))
    shutil.copyfile(gb, os.path.join(preparo_out, "outputs", "raw_data_sequences.gb"))

    cfg = dict(cfg, input_path=trees_dir, output_path=preparo_out,
               subtree_miner_configs=dict(cfg["subtree_miner_configs"], output_path=preparo_out))
    log_stdout = os.path.join(preparo_out, "outputs", "output_log.txt")
    anteriores = _metadados_anteriores(os.path.join(outputs, "metadata.json"))
    arvore_atual = []
    preservados, lidos_do_gb = set(), set()
    calcular_original = builder_mod.calculate_tree_hash

    def calcular_preservando(data, is_terminal=False, gbk_file=None):
        if not is_terminal:
            return calcular_original(data, is_terminal, gbk_file=gbk_file)
        # `os.devnull`: o resultado do `.gb` seria descartado, não vale parsear de novo.
        # A função original continua sendo quem monta `newick` e `terminal_hash`.
        resultado = calcular_original(data, True, gbk_file=os.devnull)
        previo = anteriores.get(arvore_atual[-1], {}).get(resultado["newick"], _SEM_METADADO)
        if previo is _SEM_METADADO:
            lidos_do_gb.add((arvore_atual[-1], resultado["newick"]))
            return calcular_original(data, True, gbk_file=gbk_file)
        preservados.add((arvore_atual[-1], resultado["newick"]))
        return dict(resultado, metadata=previo)

    with open(log_stdout, "w") as saida, contextlib.redirect_stdout(saida), \
            mock.patch.object(builder_mod, "download_sequences", _proibido), \
            mock.patch.object(builder_mod, "calculate_tree_hash", calcular_preservando):
        controlador = SubtreeBuilderController(**cfg)
        # `os.listdir` não tem ordem estável entre máquinas, e a ordem das árvores muda o
        # `.cql` (o deduplicador de `Metadata` emite a versão completa na primeira
        # ocorrência completa) e o `metadata.json`. A ordem da execução original é a que
        # está no `metadata.json` anterior; o que ele não conhece vai depois, em ordem alfabética.
        ordem = {nome_arvore: i for i, nome_arvore in enumerate(anteriores)}
        controlador.files = sorted(
            controlador.files, key=lambda f: (ordem.get(str(f.rsplit(".", 1)[0]), len(ordem)), f))
        construir = controlador.builder

        def construir_com_contexto(path, name):
            arvore_atual.append(str(name.rsplit(".", 1)[0]))
            return construir(path, name)
        controlador.builder = construir_com_contexto
        controlador.downloaded_sequences_cache = {
            t.name for f in os.listdir(trees_dir)
            for t in Phylo.read(os.path.join(trees_dir, f), cfg["input_format"]).get_terminals()}
        dados = controlador()

    novo = os.path.join(preparo_out, "outputs")
    if not os.path.isfile(os.path.join(novo, "all_results_fpmax.csv")):
        raise RuntimeError(f"a mineração não gerou all_results_fpmax.csv (não aplicável); "
                           f"nada foi instalado. Preparo mantido em {preparo}")

    # Subtrees/: o conteúdo não depende do item (só `List_terminals_hash` mudou, e ele
    # não está nos arquivos), então o esperado é "igual". O que difere não se sobrescreve.
    subtrees = {"iguais": 0, "novos": 0, "divergentes": []}
    for origem in sorted(glob.glob(os.path.join(preparo_out, "Subtrees", "*"))):
        destino = os.path.join(out_dir, "Subtrees", os.path.basename(origem))
        if not os.path.exists(destino):
            _instalar(origem, destino)
            subtrees["novos"] += 1
        elif manifest.file_digest(origem) == manifest.file_digest(destino):
            subtrees["iguais"] += 1
        else:
            subtrees["divergentes"].append(os.path.basename(origem))

    rotulo = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    reescritos = {os.path.join(novo, "all_results_fpmax.csv"): "all_results_fpmax.csv",
                  os.path.join(novo, "metadata.json"): "metadata.json",
                  os.path.join(novo, "metadata.csv"): "metadata.csv"}
    for cql in glob.glob(os.path.join(novo, "neo4j_commands_*.cql")):
        reescritos[cql] = os.path.basename(cql)
    for png in glob.glob(os.path.join(novo, "Plots", "*")):
        reescritos[png] = os.path.join("Plots", os.path.basename(png))
    registro_dir = os.path.join(outputs, "remine", rotulo)
    registros = {log_stdout: os.path.join(registro_dir, "output_log.txt")}
    for log in glob.glob(os.path.join(novo, "log_*.log")):
        registros[log] = os.path.join(registro_dir, os.path.basename(log))

    anteriores = _digests([os.path.join(outputs, r) for r in reescritos.values()
                           if os.path.isfile(os.path.join(outputs, r))], raiz)
    for origem, relativo in reescritos.items():
        _instalar(origem, os.path.join(outputs, relativo))
    for origem, destino in registros.items():
        _instalar(origem, destino)

    saidas = _digests([os.path.join(outputs, r) for r in reescritos.values()] + list(registros.values()), raiz)
    itens_por_arvore = {nome_arvore: len(extract_terminais_list(arvore))
                        for arvore in _arvores(dados) for nome_arvore in arvore}
    with open(os.path.join(outputs, "all_results_fpmax.csv")) as tabela:
        padroes = sum(1 for _ in tabela) - 1
    detalhes = {
        "comando": f"python -m workflow.remine --out out --motivo {motivo}",
        "arvores": len(itens_por_arvore),
        "itens_por_arvore": itens_por_arvore,
        "padroes_gravados": padroes,
        "subtrees": subtrees,
        "metadado_de_terminal": {"preservado_do_metadata_json_anterior": len(preservados),
                                 "lido_do_gb": len(lidos_do_gb)},
        "versoes": {"python": sys.version.split()[0], "mlxtend": mlxtend.__version__,
                    "pandas": pandas.__version__},
    }
    anexado = manifest.anexar_remineracao(
        out_dir, motivo, entradas, anteriores, saidas,
        repos={"BioComp_UFF": BASE_PATH, "PhyloTreeMiner": os.path.dirname(BASE_PATH)},
        detalhes=detalhes)
    shutil.rmtree(preparo, ignore_errors=True)
    return dict(detalhes, manifesto_atualizado=anexado, registro=os.path.relpath(registro_dir, raiz))


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Refaz a mineração de um out/ já executado, sem reinferir as árvores, "
                    "e anexa o registro ao manifesto.")
    ap.add_argument("--out", required=True, help="diretório out/ do projeto")
    ap.add_argument("--motivo", required=True, help="por que remineirar (ex. D29); vai para o manifesto")
    a = ap.parse_args(argv)
    try:
        resumo = remineirar(a.out, a.motivo)
    except (FileNotFoundError, RuntimeError) as e:
        print(f"remineração recusada: {e}", file=sys.stderr)
        return 2
    print(json.dumps({k: v for k, v in resumo.items() if k != "itens_por_arvore"},
                     indent=2, ensure_ascii=False))
    print("itens por árvore:", json.dumps(resumo["itens_por_arvore"], ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
