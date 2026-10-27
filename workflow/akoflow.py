"""
Workflow granular para o AkôFlow — uma CLI por atividade (meta G4b).

Cada subcomando é uma atividade do DAG de ``docs/akoflow/02-desenho-granular.md``
e **chama as mesmas funções de biblioteca** que o ``workflow.py`` monolítico:
nenhuma fórmula, parâmetro de ferramenta, semente ou política muda. O que muda
é a orquestração — o estado em memória entre etapas vira arquivo, e cada
atividade grava um fragmento de proveniência (``out/outputs/prov/<atividade>.json``)
que ``consolidar`` junta num ``manifest.json`` da mesma forma do monolito.

Os parâmetros científicos vêm do ``config_backup.json``; a linha de comando só
diz a chave da ativação. Os caminhos do arquivo são ignorados: valem os de
``--out``, porque o diretório do *workspace* muda de uma ativação para outra.

Uso::

    python -m workflow.akoflow preparar --config cfg.json --entrada data/X --out P/out
    python -m workflow.akoflow alinhar --config P/out/outputs/config_backup.json \\
        --alinhador mafft --out P/out

D17: RAxML-NG sempre com ``--threads N --workers 1`` (feito pelo ``TreeBuilder``).
"""

import argparse
import contextlib
import datetime
import glob
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import types
import uuid
from pathlib import Path
from typing import Dict, List, Optional
from unittest import mock

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_PATH not in sys.path:
    sys.path.insert(0, BASE_PATH)

from Bio import AlignIO, Phylo

from workflow.utils import external_tools, manifest as manifest_mod, run_logging, tool_runs
from workflow.utils.manifest import ExecutionManifest, file_digest

#: Fragmentos e saídas por atividade (design 02 §1.4: nenhum caminho é escrito por duas ativações).
SUBDIRS_REGISTRO = ("prov", "stdout")
RAIZ_TOKEN = "<RAIZ>"

BASICOS = ("nj_distance", "nj_parsimony", "upgma_distance", "upgma_parsimony")
AVANCADOS = ("iqtree", "fasttree", "raxml", "mrbayes")
COM_MODELO = ("iqtree", "raxml", "mrbayes")


# --------------------------------------------------------------------------- #
# Configuração e caminhos
# --------------------------------------------------------------------------- #

def _agora() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _substituir_projeto(valor, nome: str):
    """`#` → nome do projeto, em todo valor textual (como ``workflow.py``)."""
    if isinstance(valor, str):
        return valor.replace("#", nome)
    if isinstance(valor, dict):
        return {k: _substituir_projeto(v, nome) for k, v in valor.items()}
    if isinstance(valor, list):
        return [_substituir_projeto(v, nome) for v in valor]
    return valor


def _carregar(config: str, out: str):
    """
    Lê o `config_backup.json` e aponta os caminhos para `out`.

    Return
    ------
    tuple
        ``(cfg completo, tree_config, subtree_config, out absoluto)``.
    """
    out = os.path.abspath(out)
    with open(config, encoding="utf-8") as fh:
        cfg = json.load(fh)
    tree = dict(cfg["tree_config"], input_path=os.path.join(out, "entrada"), output_path=out)
    sub = dict(cfg["subtree_config"], input_path=os.path.join(out, "Trees"), output_path=out)
    sub["subtree_miner_configs"] = dict(sub.get("subtree_miner_configs", {}), output_path=out)
    return cfg, tree, sub, out


def _controlador(tree_cfg: dict):
    """
    `TreeBuilderController` sem o efeito colateral do construtor sobre a entrada.

    O construtor chama ``clean_NoPipe(input_path)``; aqui a entrada é
    ``out/entrada``, onde o FASTA deduplicado (``<stem>_NoPipe``) é um dado de
    outra atividade, não lixo.
    """
    from workflow.controller import treeBuilderController as tbc
    with mock.patch.object(tbc, "clean_NoPipe", lambda *_: None):
        return tbc.TreeBuilderController(**tree_cfg)


def _stems(out: str) -> List[str]:
    return sorted(Path(f).stem for f in glob.glob(os.path.join(out, "entrada", "*.fasta")))


def _stem(out: str, pedido: Optional[str]) -> str:
    if pedido:
        return pedido
    stems = _stems(out)
    if len(stems) != 1:
        raise ValueError(f"informe --stem: out/entrada tem {len(stems)} FASTA ({stems})")
    return stems[0]


def _fasta_usado(out: str, stem: str) -> str:
    """O FASTA que o monolito passaria adiante: o deduplicado, quando houve deduplicação."""
    sem_pipe = os.path.join(out, "entrada", f"{stem}_NoPipe")
    return sem_pipe if os.path.isfile(sem_pipe) else os.path.join(out, "entrada", f"{stem}.fasta")


def _nome_atividade(*partes: str) -> str:
    return "-".join(p.replace("_", "-") for p in partes)


def _alinhamento(out: str, stem: str, alinhador: str) -> str:
    return os.path.join(out, "tmp", f"{stem}_{alinhador}.aln")


# --------------------------------------------------------------------------- #
# Proveniência por atividade
# --------------------------------------------------------------------------- #

def _instantaneo(out: str) -> Dict[str, tuple]:
    resultado = {}
    for raiz, subdirs, arquivos in os.walk(out):
        relativa = os.path.relpath(raiz, out)
        if relativa == "outputs":
            subdirs[:] = [d for d in subdirs if d not in SUBDIRS_REGISTRO]
        for nome in arquivos:
            if nome.startswith("log_setup_"):
                continue
            caminho = os.path.join(raiz, nome)
            st = os.stat(caminho)
            resultado[os.path.relpath(caminho, out)] = (st.st_size, st.st_mtime_ns)
    return resultado


def _cota_cpu() -> Optional[float]:
    """Cota de CPU do cgroup v2; ``os.cpu_count()`` num contêiner devolve as CPUs do host."""
    try:
        quota, periodo = open("/sys/fs/cgroup/cpu.max").read().split()
        return None if quota == "max" else round(int(quota) / int(periodo), 2)
    except (OSError, ValueError):
        return None


def _higienizar(valor, raiz: str):
    """Fragmentos relativos à raiz do projeto e sem caminho de fora dela (D15)."""
    m = ExecutionManifest(raiz, {})
    def caminho(p):
        return m._fora_da_raiz(p) if isinstance(p, str) else p
    if isinstance(valor, dict):
        limpo = {}
        for k, v in valor.items():
            if k in ("saida", "alinhamento", "relatorio"):
                limpo[k] = caminho(v)
            elif k == "command" and isinstance(v, list):
                limpo[k] = m._sanitizar_comando(v)
            elif k == "comando" and isinstance(v, list):
                limpo[k] = m._sanitizar_comando(v)
            elif k == "motivo" and isinstance(v, str):
                limpo[k] = m._sanitizar_texto(v)
            else:
                limpo[k] = _higienizar(v, raiz)
        return limpo
    if isinstance(valor, list):
        return [_higienizar(v, raiz) for v in valor]
    return valor


def _reidratar(valor, raiz: str):
    """Inverso de `_higienizar` para a consolidação: o que é relativo volta a ser da raiz."""
    def junta(p):
        return p if (not isinstance(p, str) or os.path.isabs(p) or os.sep not in p) \
            else os.path.join(raiz, p)
    if isinstance(valor, dict):
        novo = {}
        for k, v in valor.items():
            if k in ("saida", "alinhamento", "relatorio"):
                novo[k] = junta(v)
            elif k in ("command", "comando") and isinstance(v, list):
                novo[k] = [junta(t) for t in v]
            else:
                novo[k] = _reidratar(v, raiz)
        return novo
    if isinstance(valor, list):
        return [_reidratar(v, raiz) for v in valor]
    return valor


def _chave_de(caminho: str, out: str) -> str:
    """Relativo a `out`; fora dele, só o nome (D15)."""
    relativo = os.path.relpath(caminho, out)
    return os.path.basename(caminho) if relativo.startswith(os.pardir) else relativo


def _log_atual() -> Optional[str]:
    return run_logging.caminho_atual()


@contextlib.contextmanager
def _atividade(nome: str, out: str, chave: dict, entradas=()):
    """
    Envolve uma atividade: registros zerados, stdout para arquivo, fragmento ao fim.

    O fragmento é gravado também em falha, com o desfecho: o diagnóstico de uma
    execução que morreu é o que o manifesto parcial já fazia no monolito.
    """
    outputs = os.path.join(out, "outputs")
    for d in SUBDIRS_REGISTRO:
        os.makedirs(os.path.join(outputs, d), exist_ok=True)
    raiz = os.path.dirname(out)
    tool_runs.limpar()
    antes = _instantaneo(out)
    inicio = _agora()
    run_logging.configurar(outputs, uuid.uuid4().hex)
    estado = {"desfecho": ("concluido", None), "extra": {}}
    with open(os.path.join(outputs, "stdout", f"{nome}.txt"), "w") as saida, \
            contextlib.redirect_stdout(saida):
        try:
            yield estado
        except BaseException as e:
            estado["desfecho"] = ("falhou", f"{type(e).__name__}: {e}")
            raise
        finally:
            depois = _instantaneo(out)
            novos = sorted(p for p, v in depois.items() if antes.get(p) != v)
            fragmento = {
                "fragmento_versao": 1,
                "atividade": nome,
                "chave": chave,
                "inicio_utc": inicio,
                "fim_utc": _agora(),
                "desfecho": {"status": estado["desfecho"][0], "motivo": estado["desfecho"][1]},
                "entradas_sha256": {_chave_de(e, out): file_digest(e) for e in entradas},
                "saidas_sha256": {p: file_digest(os.path.join(out, p)) for p in novos},
                "tools_invoked": tool_runs.execucoes(),
                "inference_methods": tool_runs.metodos(),
                "model_selection": tool_runs.selecoes_modelo(),
                "log_file": (os.path.relpath(_log_atual(), raiz) if _log_atual() else None),
                "tools_available": dict(manifest_mod.tool_versions()),
                "environment": dict(manifest_mod._environment(), cpu_quota_cgroup=_cota_cpu()),
                **estado["extra"],
            }
            fragmento = _higienizar(fragmento, raiz)
            with open(os.path.join(outputs, "prov", f"{nome}.json"), "w", encoding="utf-8") as fh:
                json.dump(fragmento, fh, indent=2, ensure_ascii=False, sort_keys=True)
                fh.write("\n")


def _ler_fragmento(out: str, nome: str) -> Optional[dict]:
    caminho = os.path.join(out, "outputs", "prov", f"{nome}.json")
    if not os.path.isfile(caminho):
        return None
    with open(caminho, encoding="utf-8") as fh:
        return json.load(fh)


# --------------------------------------------------------------------------- #
# Atividades
# --------------------------------------------------------------------------- #

def preparar(a) -> None:
    """`preparar-entrada`: valida e deduplica o FASTA, grava a configuração e os acessos."""
    from workflow.controller.treeBuilderController import eh_fasta
    from workflow.utils.dataValidation import (deduplicar_por_sequencia, duplicate_names,
                                               duplicate_seq, validate_sequences)
    from Bio import SeqIO

    out = os.path.abspath(a.out)
    with open(a.config, encoding="utf-8") as fh:
        params = json.load(fh)
    nome_projeto = a.projeto_nome or params.get("project_name") or os.path.basename(os.path.dirname(out))
    params["project_name"] = nome_projeto
    params = _substituir_projeto(params, nome_projeto)
    params["output_log"] = "out"
    params["tree_config"].update(input_path="out/entrada", output_path="out")
    params["subtree_config"].update(input_path="out/Trees", output_path="out")
    params["subtree_config"].setdefault("subtree_miner_configs", {})["output_path"] = "out"

    fontes = [os.path.join(a.entrada, n) for n in sorted(os.listdir(a.entrada))
              if eh_fasta(os.path.join(a.entrada, n))]
    if not fontes:
        raise FileNotFoundError(f"nenhum FASTA em {a.entrada}")
    with _atividade("preparar-entrada", out, {"entrada": os.path.basename(a.entrada)}, fontes):
        entrada = os.path.join(out, "entrada")
        os.makedirs(entrada, exist_ok=True)
        acessos = []
        for fonte in fontes:
            stem = Path(fonte).stem
            copia = os.path.join(entrada, f"{stem}.fasta")
            shutil.copyfile(fonte, copia)
            usado = copia
            if not (not duplicate_names(copia) and not duplicate_seq(copia)[0]
                    and validate_sequences(copia)):
                logging.warning(f"Arquivo {stem} contém duplicatas ou erros, deduplicando por sequência.")
                usado, descartados = deduplicar_por_sequencia(stem, copia, entrada)
                if descartados:
                    tool_runs.registrar(
                        "deduplicacao", ["deduplicar_por_sequencia", stem], saida=usado,
                        descartados=[{"descartado": d, "mantido": m} for d, m in descartados],
                        nota=("acessos com sequência idêntica; o mantido é o primeiro do "
                              "arquivo — ver D23"))
            acessos.extend(r.id for r in SeqIO.parse(usado, "fasta"))
        with open(os.path.join(entrada, "acessos.txt"), "w") as fh:
            fh.write("\n".join(acessos) + "\n")
        os.makedirs(os.path.join(out, "outputs"), exist_ok=True)
        with open(os.path.join(out, "outputs", "config_backup.json"), "w", encoding="utf-8") as fh:
            json.dump(params, fh, indent=4)


def baixar_metadados(a) -> None:
    """`baixar-metadados`: o `.gb` completo, baixado uma vez (ou copiado, no modo offline)."""
    out = os.path.abspath(a.out)
    destino = os.path.join(out, "outputs", "raw_data_sequences.gb")
    entradas = [a.acessos] + ([a.gb_congelado] if a.gb_congelado else [])
    with _atividade("baixar-metadados", out, {"modo": "congelado" if a.gb_congelado else "online"},
                    entradas):
        os.makedirs(os.path.dirname(destino), exist_ok=True)
        if a.gb_congelado:
            shutil.copyfile(a.gb_congelado, destino)
            return
        from workflow.utils.treeUtils import download_sequences
        with open(a.acessos) as fh:
            ids = [linha.strip() for linha in fh if linha.strip()]
        resposta = download_sequences(None, [types.SimpleNamespace(name=i) for i in ids],
                                      destino, log=logging)
        if isinstance(resposta, dict) and "error" in resposta:
            raise RuntimeError(f"download do GenBank falhou: {resposta['error']}")
        if not os.path.isfile(destino) or os.path.getsize(destino) == 0:
            raise RuntimeError("o download do GenBank não gravou nenhum registro")


def alinhar(a) -> None:
    """`alinhar-<alg>`: um alinhamento, nomeado pelo alinhador pedido."""
    cfg, tree, _, out = _carregar(a.config, a.out)
    stem = _stem(out, a.stem)
    fasta = _fasta_usado(out, stem)
    nome = _nome_atividade("alinhar", a.alinhador)
    with _atividade(nome, out, {"stem": stem, "alinhador": a.alinhador}, [fasta]):
        ctrl = _controlador(tree)
        caminhos = ctrl._prepare_output_paths(stem)
        saida = _alinhamento(out, stem, a.alinhador)
        efetivo, _ = ctrl._resolver_alinhador(a.alinhador, fasta)
        ctrl._alinhar(efetivo, fasta, saida, caminhos["align_html"])
        os.makedirs(os.path.join(out, "Align"), exist_ok=True)
        shutil.copyfile(saida, os.path.join(out, "Align", os.path.basename(saida)))


def selecionar_modelo(a) -> None:
    """`modelo-<alg>`: ModelFinder uma vez por alinhamento; a tradução vai para `selecao.json`."""
    from workflow.tree_construction.modelo_substituicao import SelecaoDeModeloFalhou
    cfg, tree, _, out = _carregar(a.config, a.out)
    stem = _stem(out, a.stem)
    aln = _alinhamento(out, stem, a.alinhador)
    nome = _nome_atividade("modelo", a.alinhador)
    with _atividade(nome, out, {"stem": stem, "alinhador": a.alinhador}, [aln]) as at:
        ctrl = _controlador(tree)
        alng = AlignIO.read(aln, "fasta")
        selecao = {"alinhamento": os.path.relpath(aln, out), "estado": "concluida",
                   "traducao": None, "motivo": None}
        try:
            selecao["traducao"] = ctrl._modelo_do_alinhamento(aln, alng, a.alinhador)
        except SelecaoDeModeloFalhou as e:
            selecao.update(estado="falhou", motivo=str(e))
        destino = os.path.join(out, "tmp", f"modelfinder_{stem}_{a.alinhador}", "selecao.json")
        os.makedirs(os.path.dirname(destino), exist_ok=True)
        with open(destino, "w", encoding="utf-8") as fh:
            json.dump(selecao, fh, indent=2, ensure_ascii=False)
        at["extra"]["selecao_modelo"] = selecao["estado"]


def _carregar_selecao(ctrl, out: str, stem: str, alinhador: str, aln: str) -> None:
    """Põe no controlador a seleção de modelo já feita, para ele não rodar outra."""
    from workflow.tree_construction.modelo_substituicao import SelecaoDeModeloFalhou
    caminho = os.path.join(out, "tmp", f"modelfinder_{stem}_{alinhador}", "selecao.json")
    if not os.path.isfile(caminho):
        raise FileNotFoundError(f"{os.path.relpath(caminho, out)} não existe: rode `selecionar-modelo` antes")
    with open(caminho, encoding="utf-8") as fh:
        selecao = json.load(fh)
    ctrl._selecoes_modelo[aln] = (SelecaoDeModeloFalhou(selecao["motivo"])
                                  if selecao["estado"] == "falhou" else selecao["traducao"])


def inferir(a) -> None:
    """`arvore-<alg>-<método>`: uma árvore, pelo caminho do monolito (`_process_*_tree_method`)."""
    cfg, tree, _, out = _carregar(a.config, a.out)
    stem = _stem(out, a.stem)
    fasta = _fasta_usado(out, stem)
    aln = _alinhamento(out, stem, a.alinhador)
    nome = _nome_atividade("arvore", a.alinhador, a.metodo)
    with _atividade(nome, out, {"stem": stem, "alinhador": a.alinhador, "metodo": a.metodo},
                    [aln]) as at:
        ctrl = _controlador(tree)
        caminhos = ctrl._prepare_output_paths(stem)
        multi = ctrl._initialize_multi_trees_structure()
        if a.metodo in BASICOS:
            construtor, tipo = a.metodo.split("_")
            ctrl.construct_tree_method = construtor
            if ctrl._should_ignore_method(tipo):
                ctrl._registrar_basico_ignorado(stem, a.alinhador, construtor, tipo)
                return
            executar = lambda: ctrl._process_tree_method(
                stem, fasta, caminhos, a.alinhador, construtor, tipo, multi)
        elif a.metodo in AVANCADOS:
            if ctrl._should_ignore_method(a.metodo):
                ctrl._registrar_ignorado(a.metodo, ctrl._saida_avancado(stem, a.alinhador, a.metodo),
                                         a.alinhador, f"ignore_mode contém {a.metodo!r}")
                return
            if a.metodo in COM_MODELO and ctrl.model_selection != "nenhuma":
                _carregar_selecao(ctrl, out, stem, a.alinhador, aln)
            executar = lambda: ctrl._process_advanced_tree_method(
                stem, fasta, caminhos, a.alinhador, a.metodo, multi)
        else:
            raise ValueError(f"método desconhecido: {a.metodo!r}")
        # `on_method_failure: continue` devolve 0 e registra `tentado_e_falhou`; `fail` levanta.
        at["extra"]["arvores_construidas"] = executar()


def subarvores(a) -> None:
    """`subarvores-<alg>-<método>`: subárvores e metadados de uma árvore, com o `.gb` já baixado."""
    from workflow.subtree_construction.builder import SubtreeBuilder
    cfg, _, sub, out = _carregar(a.config, a.out)
    arvore = os.path.abspath(a.arvore)
    nome_arquivo = os.path.basename(arvore)
    nome_arvore = nome_arquivo.rsplit(".", 1)[0]
    gb = os.path.join(out, "outputs", "raw_data_sequences.gb")
    stem = _stem(out, a.stem)
    nome = "subarvores-" + nome_arvore[len(f"tree_{stem}_"):].replace("_", "-")
    with _atividade(nome, out, {"arvore": nome_arvore}, [arvore, gb]) as at:
        if not os.path.isfile(arvore):
            motivo = _motivo_sem_arvore(out, nome_arquivo)
            at["extra"]["no_op"] = motivo
            logging.warning(f"subarvores no-op para {nome_arvore}: {motivo}")
            return
        os.makedirs(os.path.join(out, "Subtrees"), exist_ok=True)
        cache = {t.name for t in Phylo.read(arvore, sub["input_format"]).get_terminals()}
        dados = SubtreeBuilder(cache, **sub).subtree_constructor(arvore, nome_arquivo)
        destino = os.path.join(out, "outputs", "subarvores", f"{nome_arvore}.json")
        os.makedirs(os.path.dirname(destino), exist_ok=True)
        with open(destino, "w") as fh:
            json.dump(dados, fh)


def _motivo_sem_arvore(out: str, nome_arquivo: str) -> str:
    """Árvore ausente só é *no-op* se alguma atividade declarou por quê; senão é erro."""
    for caminho in sorted(glob.glob(os.path.join(out, "outputs", "prov", "arvore-*.json"))):
        with open(caminho, encoding="utf-8") as fh:
            for m in json.load(fh).get("inference_methods", []):
                if os.path.basename(m["saida"]) == nome_arquivo and \
                        m["estado"] in ("tentado_e_falhou", "ignorado_por_configuracao"):
                    return f"{m['estado']}: {m.get('motivo', '')}".strip()
    raise FileNotFoundError(f"{nome_arquivo} não existe e nenhuma atividade declarou por quê")


def heatmap_rf(a) -> None:
    """`heatmap-rf`: distâncias RF entre as árvores, com o `multi_trees` do laço do monolito."""
    from workflow.utils.metrics import plot_heatmap_distances, process_rf_distance
    import numpy as np
    cfg, tree, _, out = _carregar(a.config, a.out)
    stem = _stem(out, a.stem)
    arvores = sorted(glob.glob(os.path.join(out, "Trees", f"tree_{stem}_*.nexus")))
    with _atividade("heatmap-rf", out, {"stem": stem}, arvores):
        ctrl = _controlador(tree)
        multi = ctrl._initialize_multi_trees_structure()
        # Reproduz as visitas de `_process_advanced_mode`: os avançados são visitados
        # uma vez por método de distância, e a árvore em disco é lida a cada visita.
        for alg in ctrl.aligners:
            for construtor in ("nj", "upgma"):
                for tipo in ("distance", "parsimony"):
                    t = ctrl._load_existing_tree(ctrl._saida_basico(stem, alg, construtor, tipo),
                                                 ctrl.output_format)
                    if t:
                        multi[alg][tipo][construtor].append(t)
                for metodo in ctrl.METODOS_AVANCADOS:
                    t = ctrl._load_existing_tree(ctrl._saida_avancado(stem, alg, metodo),
                                                 ctrl.output_format)
                    if t:
                        multi[alg][metodo].append(t)
        plots = os.path.join(out, "outputs", "Plots")
        os.makedirs(plots, exist_ok=True)
        scores = process_rf_distance(multi)
        matriz = ctrl.somarMatrizes(
            np.zeros((8, 8)),
            plot_heatmap_distances(data_dict=multi, scores=scores, base_name=stem, path=plots))
        plot_heatmap_distances(data_dict=multi, base_name="Acumulate", path=plots,
                               distance_matrix=matriz)


def _carregar_subarvores(out: str) -> List[dict]:
    """Os `result_dict` por árvore, em ordem canônica (nome do arquivo)."""
    arquivos = sorted(glob.glob(os.path.join(out, "outputs", "subarvores", "*.json")))
    if not arquivos:
        raise FileNotFoundError("out/outputs/subarvores/ está vazio: nenhuma árvore chegou à extração")
    dados = []
    for arquivo in arquivos:
        with open(arquivo) as fh:
            dados.append(json.load(fh))
    return dados


def minerar(a) -> None:
    """`minerar-fpmax`: FPMax sobre as bipartições de todas as árvores (`process_group` do monolito)."""
    import pandas as pd
    from workflow.subtree_mining.miner import SubtreeMiner
    cfg, _, sub, out = _carregar(a.config, a.out)
    arquivos = sorted(glob.glob(os.path.join(out, "outputs", "subarvores", "*.json")))
    with _atividade("minerar-fpmax", out, {"modo": sub["subtree_miner_configs"].get("mode")}, arquivos):
        dados = _carregar_subarvores(out)
        minerador = SubtreeMiner(**sub["subtree_miner_configs"])
        if minerador.mode == "OFST":
            resultado = [minerador.process_group(grupo, base)
                         for base, grupo in minerador.group_data_by_tree_base(dados).items()]
        else:
            resultado = minerador.process_group(dados)
        outputs = os.path.join(out, "outputs")
        pd.DataFrame(resultado).to_csv(os.path.join(outputs, "metadata.csv"))
        with open(os.path.join(outputs, "metadata.json"), "w") as fh:
            fh.write(json.dumps(resultado, indent=2))


def _metadata(out: str):
    with open(os.path.join(out, "outputs", "metadata.json")) as fh:
        return json.load(fh)


def gerar_cql(a) -> None:
    """`gerar-cql`: o `.cql` a partir do `metadata.json` minerado (D-a: separado da mineração)."""
    from workflow.utils.neo4jProcessing import parse_tree
    cfg, _, sub, out = _carregar(a.config, a.out)
    entrada = os.path.join(out, "outputs", "metadata.json")
    with _atividade("gerar-cql", out, {}, [entrada]):
        dados = _metadata(out)
        # `parse_tree` deriva o nome do projeto do caminho; o preparo repete `<projeto>/out`
        # para o `.cql` sair com o nome do config, e não com o do diretório do workspace.
        with tempfile.TemporaryDirectory() as preparo:
            destino = os.path.join(preparo, cfg["project_name"], "out", "outputs")
            os.makedirs(destino)
            if sub["subtree_miner_configs"].get("mode") == "OFST":
                for resultado in dados:
                    base = "_".join(next(iter(resultado[0])).split("_")[:2])
                    parse_tree(trees=resultado, path=destino, mode=base)
            else:
                parse_tree(trees=dados, path=destino)
            for cql in glob.glob(os.path.join(destino, "*.cql")):
                shutil.copyfile(cql, os.path.join(out, "outputs", os.path.basename(cql)))


def graficos_suporte(a) -> None:
    """`graficos-suporte`: histogramas de suporte sobre o `metadata.json`."""
    from workflow.utils.metrics import process_histogram_frequence
    cfg, _, _, out = _carregar(a.config, a.out)
    with _atividade("graficos-suporte", out, {}, [os.path.join(out, "outputs", "metadata.json")]):
        plots = os.path.join(out, "outputs", "Plots")
        os.makedirs(plots, exist_ok=True)
        process_histogram_frequence(_metadata(out), plots)


def estabilidade(a) -> None:
    """`estabilidade` (extensão, fora do monolito e da equivalência): CLI que já existe."""
    cfg, _, _, out = _carregar(a.config, a.out)
    raiz = os.path.dirname(out)
    with _atividade("estabilidade", out, {}, [os.path.join(out, "Trees")]):
        cmd = [sys.executable, "-m", "workflow.stability.case_study", "--project", raiz]
        for fasta in _stems(out):
            cmd += ["--fasta", os.path.join(out, "entrada", f"{fasta}.fasta")]
        subprocess.run(cmd, check=True, cwd=BASE_PATH)


# --------------------------------------------------------------------------- #
# Consolidação
# --------------------------------------------------------------------------- #

def _ordem_metodos() -> List[str]:
    """A ordem em que o monolito *visita* os pipelines de um alinhador (a 1ª visita decide)."""
    ordem = []
    for construtor in ("nj", "upgma"):
        for m in (f"{construtor}_distance", f"{construtor}_parsimony") + AVANCADOS:
            if m not in ordem:
                ordem.append(m)
    return ordem


def ordem_de_visita(aligners) -> List[tuple]:
    """
    ``[(atividade, metodo ou None)]`` na ordem em que o monolito executaria os passos
    por alinhador: alinhamento, árvores de distância e parcimônia do NJ, a seleção de
    modelo (preguiçosa: roda antes do primeiro método com modelo), e o resto.
    """
    passos = []
    for alg in aligners:
        passos.append((_nome_atividade("alinhar", alg), None))
        modelo_feito = False
        for metodo in _ordem_metodos():
            if metodo in COM_MODELO and not modelo_feito:
                passos.append((_nome_atividade("modelo", alg), None))
                modelo_feito = True
            passos.append((_nome_atividade("arvore", alg, metodo), (alg, metodo)))
    return passos


def consolidar(a) -> None:
    """`consolidar-manifesto`: recompõe `manifest.json`, o log e o `output_log.txt` do monolito."""
    from workflow.controller.molecularClockController import MolecularClockController
    from workflow.controller.treeBuilderController import MODOS_BASICOS
    from workflow.tree_construction.builder import reproducibility_settings

    cfg, tree, _, out = _carregar(a.config, a.out)
    if tree.get("run_molecular_clock"):
        raise NotImplementedError("a datação granular (datar-*) não faz parte desta entrega")
    raiz = os.path.dirname(out)
    outputs = os.path.join(out, "outputs")
    tool_runs.limpar()
    ctrl = _controlador(tree)
    stems = _stems(out)
    passos = ordem_de_visita(ctrl.aligners)

    # Ordem de visita; o que não é árvore (alinhar, modelo) não tem método.
    ordem = ["preparar-entrada", "baixar-metadados"] + [n for n, _ in passos]
    ordem += ["heatmap-rf"]
    ordem += [Path(p).stem for p in sorted(glob.glob(os.path.join(outputs, "prov", "subarvores-*.json")))]
    ordem += ["minerar-fpmax", "gerar-cql", "graficos-suporte", "estabilidade"]
    fragmentos = {n: _ler_fragmento(out, n) for n in ordem}
    ausentes = []

    m = ExecutionManifest(
        project_root=raiz, params=cfg,
        repos={"BioComp_UFF": BASE_PATH, "PhyloTreeMiner": os.path.dirname(BASE_PATH)})
    presentes = [f for f in fragmentos.values() if f]
    if presentes:
        m.started_at = datetime.datetime.fromisoformat(min(f["inicio_utc"] for f in presentes))
    log = run_logging.configurar(outputs, m.run_id)
    m.register_log(log)
    for stem in stems:
        m.register_input(os.path.join(out, "entrada", f"{stem}.fasta"))
    m.register_reproducibility(reproducibility_settings(tree))
    modo = tree.get("mode", "advanced")
    mapa = {"iqtree": "iqtree", "fasttree": "fasttree", "raxml": "raxml-ng", "mrbayes": "mrbayes"}
    disponiveis = [mt for mt, ferr in mapa.items() if external_tools.resolve_tool(ferr)]
    executados = [] if modo in MODOS_BASICOS else [x for x in disponiveis if x not in ctrl.ignore_methods]
    m.register_execution_mode(modo, disponiveis, executados)

    execucoes: Dict[str, dict] = {}
    metodos: Dict[str, dict] = {}
    selecoes: List[dict] = []
    for nome in ordem:
        f = fragmentos.get(nome)
        if f is None:
            continue
        f = _reidratar(f, raiz)
        for ferramenta, dados in f["tools_invoked"].items():
            alvo = execucoes.setdefault(ferramenta, {"runs": []})
            alvo["runs"].extend(dados.get("runs", []))
            alvo.update({k: v for k, v in dados.items() if k != "runs"})
        # Cada seleção é registrada pelo `modelo-*`; a árvore que a leu só a repete.
        if nome.startswith("modelo-"):
            selecoes.extend(f["model_selection"])
        for reg in f["inference_methods"]:
            metodos.setdefault(reg["saida"], reg)

    # Pipelines sem atividade no DAG: só os ignorados por configuração são legítimos.
    stem0 = stems[0]
    ordenados = []
    for alg in ctrl.aligners:
        for metodo in _ordem_metodos():
            if metodo in BASICOS:
                construtor, tipo = metodo.split("_")
                saida, rotulo, chave_ignore = ctrl._saida_basico(stem0, alg, construtor, tipo), metodo, tipo
            else:
                saida, rotulo, chave_ignore = ctrl._saida_avancado(stem0, alg, metodo), metodo, metodo
            reg = metodos.get(saida)
            if reg is None:
                if ctrl._should_ignore_method(chave_ignore):
                    reg = {"metodo": rotulo, "estado": "ignorado_por_configuracao", "saida": saida,
                           "alinhador": alg, "motivo": f"ignore_mode contém {chave_ignore!r}"}
                else:
                    ausentes.append(_nome_atividade("arvore", alg, metodo))
                    continue
            ordenados.append(reg)
    m.absorver_registros(execucoes, ordenados, selecoes)

    tool_runs.limpar()
    MolecularClockController(**tree)()
    m.absorver_registros(tool_runs.execucoes(), tool_runs.metodos(), [])

    falhas = [n for n, f in fragmentos.items() if f and f["desfecho"]["status"] != "concluido"]
    obrigatorias = ["preparar-entrada", "baixar-metadados", "minerar-fpmax", "gerar-cql", "graficos-suporte"]
    ausentes += [n for n in obrigatorias if fragmentos.get(n) is None]
    desfecho = ("falhou", f"atividades ausentes: {ausentes}; falhas: {falhas}") \
        if (ausentes or falhas) else ("concluido", None)
    m.register_outcome(*desfecho)

    # Logs e stdout por atividade, na ordem de visita, como se fossem de um processo só.
    logging.info("Manifesto consolidado a partir de %d fragmentos", len(presentes))
    logging.shutdown()
    proprio = open(log, encoding="utf-8").read()
    with open(log, "w", encoding="utf-8") as fh:
        for nome in ordem:
            f = fragmentos.get(nome)
            if f and f.get("log_file") and os.path.isfile(os.path.join(raiz, f["log_file"])):
                fh.write(open(os.path.join(raiz, f["log_file"]), encoding="utf-8").read())
        fh.write(proprio)
    with open(os.path.join(outputs, "output_log.txt"), "w") as fh:
        for nome in ordem:
            caminho = os.path.join(outputs, "stdout", f"{nome}.txt")
            if os.path.isfile(caminho):
                fh.write(open(caminho).read())

    m.register_outputs(out)
    m.register_outputs(os.path.join(outputs, "molecular_clock"), suffixes=(".nex", ".json"))
    if presentes:
        m.finished_at = datetime.datetime.fromisoformat(max(f["fim_utc"] for f in presentes))
    caminho = m.write()
    with open(caminho, encoding="utf-8") as fh:
        dados = json.load(fh)
    # Campo aditivo (o mesmo critério dos de M7.6): `manifest_version` continua 2.
    dados["execucao_granular"] = {
        "atividades": [{"atividade": n, "inicio_utc": f["inicio_utc"], "fim_utc": f["fim_utc"],
                        "desfecho": f["desfecho"], "ambiente": f["environment"]}
                       for n in ordem if (f := fragmentos.get(n))],
    }
    with open(caminho, "w", encoding="utf-8") as fh:
        json.dump(dados, fh, indent=2, ensure_ascii=False, sort_keys=True)
    if desfecho[0] == "falhou":
        raise SystemExit(f"consolidação com pendências: {desfecho[1]}")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="python -m workflow.akoflow",
                                 description="Atividades do workflow granular (AkôFlow).")
    sub = ap.add_subparsers(dest="comando", required=True)

    def comum(p, config=True):
        if config:
            p.add_argument("--config", required=True, help="config_backup.json")
        p.add_argument("--out", required=True, help="diretório out/ do projeto")
        p.add_argument("--stem", help="stem do FASTA (padrão: o único em out/entrada)")
        return p

    p = sub.add_parser("preparar", help="valida e deduplica o FASTA; grava config_backup.json")
    p.add_argument("--config", required=True, help="configuração-modelo (zika21.json)")
    p.add_argument("--entrada", required=True, help="diretório com os FASTA")
    p.add_argument("--out", required=True)
    p.add_argument("--projeto-nome", help="substitui project_name (e o '#' do modelo)")
    p.set_defaults(func=preparar)

    p = comum(sub.add_parser("baixar-metadados", help="baixa (ou copia) o .gb completo"))
    p.add_argument("--acessos", required=True, help="acessos.txt")
    p.add_argument("--gb-congelado", help="copia este .gb em vez de baixar (modo offline)")
    p.set_defaults(func=baixar_metadados)

    p = comum(sub.add_parser("alinhar", help="um alinhamento"))
    p.add_argument("--alinhador", required=True)
    p.set_defaults(func=alinhar)

    p = comum(sub.add_parser("selecionar-modelo", help="ModelFinder de um alinhamento"))
    p.add_argument("--alinhador", required=True)
    p.set_defaults(func=selecionar_modelo)

    p = comum(sub.add_parser("inferir", help="uma árvore (alinhador x método)"))
    p.add_argument("--alinhador", required=True)
    p.add_argument("--metodo", required=True, choices=BASICOS + AVANCADOS)
    p.set_defaults(func=inferir)

    p = comum(sub.add_parser("subarvores", help="subárvores e metadados de uma árvore"))
    p.add_argument("--arvore", required=True)
    p.set_defaults(func=subarvores)

    for nome, func, ajuda in (("heatmap-rf", heatmap_rf, "distâncias RF entre as árvores"),
                              ("minerar", minerar, "FPMax"),
                              ("gerar-cql", gerar_cql, "comandos Cypher"),
                              ("graficos-suporte", graficos_suporte, "histogramas de suporte"),
                              ("estabilidade", estabilidade, "extensão opcional"),
                              ("consolidar", consolidar, "manifesto, log e output_log")):
        comum(sub.add_parser(nome, help=ajuda)).set_defaults(func=func)
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    args = _parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
