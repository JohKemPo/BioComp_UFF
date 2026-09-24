"""
Manifesto de execução — o que é preciso para reexecutar e obter o mesmo resultado.

Uma figura só é reproduzível se for possível dizer, sem ambiguidade, **qual
código**, **qual entrada** e **qual ambiente** a produziram. Hoje nada disso é
registrado: os `config_backup.json` guardam parâmetros, mas não a versão das
ferramentas, não o commit, não as sementes efetivas e não o hash das entradas.
É o defeito D11, e é o que torna inatingível o item "cada figura reproduzível
por script + commit + hash" do checklist de submissão.

Este módulo coleta esses fatos e os grava em `out/outputs/manifest.json`.

**Regra de privacidade.** O manifesto vai para o repositório e pode ir para o
material suplementar de um artigo. Ele **não** registra nome de usuário,
*hostname* nem caminho absoluto — foi assim que D15 vazou `/home/<usuário>` de
um terceiro pela API. Todo caminho é relativo à raiz do projeto.
"""

from __future__ import annotations

import datetime
import functools
import hashlib
import json
import logging
import os
import platform
import re
import subprocess
import sys
import uuid
from typing import Dict, List, Optional

from workflow.utils.external_tools import resolve_tool
from workflow.utils import tool_runs

__all__ = [
    "ExecutionManifest",
    "tool_versions",
    "file_digest",
    "MANIFEST_FILENAME",
    "anexar_execucao_isolada",
]

MANIFEST_FILENAME = "manifest.json"

#: Ferramentas externas cuja versão muda o resultado, e como perguntá-la.
#: A saída de cada uma é filtrada por uma regex, porque quase nenhuma respeita
#: `--version` da mesma forma — o FastTree, por exemplo, imprime a versão numa
#: mensagem de uso e sai com código diferente de zero.
#: `chave em external_tools.CANDIDATOS -> (argumentos de versão, regex)`.
#: O binário não é fixado aqui: quem o encontra é `resolve_tool`, porque o nome
#: muda com a versão do pacote — `iqtree` 3.x não instala `iqtree2`.
_TOOLS = {
    "mafft": (["--version"], r"v[\d.]+"),
    "clustalo": (["--version"], r"[\d.]+"),
    "muscle": (["-version"], r"v?[\d.]+"),
    "fasttree": ([], r"Version\s+([\d.]+)"),
    "iqtree": (["--version"], r"version\s+([\d.]+)"),
    "raxml-ng": (["--version"], r"v\.\s*([\d.]+)"),
    # `mb -h` imprime só o uso, sem versão; ela sai no banner de abertura, que
    # aparece quando o binário roda com stdin fechado (D20).
    "mrbayes": ([], r"MrBayes\s+v?([\d.]+)"),
}


def _executar(cmd: List[str], timeout: int = 10) -> str:
    """
    Roda um comando e devolve stdout+stderr, ou '' se a ferramenta não existe.

    `stdin=DEVNULL` não é detalhe: o FastTree invocado sem argumentos **lê a
    entrada padrão** e fica bloqueado até o timeout. Sem isso, coletar as
    versões custaria dezenas de segundos no início de toda execução do
    pipeline — medido: 20 s por chamada.
    """
    try:
        resultado = subprocess.run(cmd, capture_output=True, text=True,
                                   stdin=subprocess.DEVNULL,
                                   timeout=timeout, check=False)
        return (resultado.stdout or "") + (resultado.stderr or "")
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return ""


@functools.lru_cache(maxsize=1)
def tool_versions() -> Dict[str, Optional[str]]:
    """
    Versão de cada ferramenta externa do pipeline.

    Return
    ------
    dict
        ``nome -> versão`` ou ``nome -> None`` quando a ferramenta não está no
        PATH. **Ausente é `None`, nunca string vazia nem "desconhecida"**: uma
        ferramenta que não existe é um fato, e o manifesto tem de dizê-lo.
    """
    versoes: Dict[str, Optional[str]] = {}
    for nome, (args, padrao) in _TOOLS.items():
        caminho = resolve_tool(nome)
        if caminho is None:
            versoes[nome] = None
            continue
        saida = _executar([caminho, *args])
        achado = re.search(padrao, saida) if saida else None
        versoes[nome] = (achado.group(achado.lastindex or 0) if achado else None)
    return versoes


def file_digest(path: str, chunk: int = 1024 * 1024) -> Optional[str]:
    """
    SHA-256 de um arquivo, lido em blocos.

    Os `metadata.json` chegam a 3,2 GB: ler tudo em memória para hashear seria
    trocar um problema por outro.

    Return
    ------
    str or None
        Digest hexadecimal, ou ``None`` se o arquivo não existe.
    """
    if not os.path.isfile(path):
        return None
    digestor = hashlib.sha256()
    with open(path, "rb") as handle:
        for bloco in iter(lambda: handle.read(chunk), b""):
            digestor.update(bloco)
    return digestor.hexdigest()


def _relativo_a(path: str, raiz: str) -> str:
    """Versão livre de `ExecutionManifest._relativo` — mesma regra, sem instância."""
    try:
        return os.path.relpath(os.path.abspath(path), raiz)
    except ValueError:
        return os.path.basename(path)


def _fora_da_raiz_de(path: str, raiz: str) -> str:
    """Versão livre de `ExecutionManifest._fora_da_raiz` — mesma regra, sem instância."""
    relativo = _relativo_a(path, raiz)
    return os.path.basename(path) if relativo.startswith(os.pardir) else relativo


def anexar_execucao_isolada(out_dir: str, diretorio_saida: str) -> bool:
    """
    Acrescenta ao `manifest.json` já existente o desfecho de uma etapa
    disparada **fora** de `workflow.py` (E12.11: relógio molecular via CLI/API
    isolada, sem reprocessar o resto do projeto).

    Drena `tool_runs` exatamente como `ExecutionManifest.drain_tool_runs` drena
    numa execução normal — mesma fonte, mesmos campos —, mas anexa a um
    manifesto que a execução original já fechou, em vez de construir um novo.
    Um manifesto inteiro (``git``, ``environment``, ``params``, ``inputs_sha256``
    da árvore original) descreve a execução que **construiu** o projeto; gerar
    outro do zero para uma etapa isolada sobrescreveria essa proveniência com
    uma que nunca existiu — pior que a lacuna que isto fecha.

    Se o projeto não tem manifesto (anterior a M2.5, ou nunca gerado), não
    inventa um: devolve `False` e deixa `tool_runs` como estava, para o
    chamador decidir se isso é motivo de aviso.

    Parameters
    ----------
    out_dir : str
        O diretório ``out/`` do projeto (mesmo `--out` do CLI da etapa).
    diretorio_saida : str
        Onde a etapa isolada gravou seus próprios arquivos (ex.
        ``out/outputs/molecular_clock/``) — hasheados e somados a
        `outputs_sha256`, sem tocar o que já estava lá.

    Returns
    -------
    bool
        `True` se um manifesto existia e foi atualizado; `False` se não havia
        manifesto para anexar.
    """
    project_root = os.path.dirname(os.path.abspath(out_dir))
    caminho = os.path.join(out_dir, "outputs", MANIFEST_FILENAME)
    if not os.path.isfile(caminho):
        return False

    with open(caminho, "r", encoding="utf-8") as f:
        dados = json.load(f)

    agora = datetime.datetime.now(datetime.timezone.utc).isoformat()

    metodos_existentes = dados.get("inference_methods") or []
    for registro in tool_runs.metodos():
        entrada = {
            "metodo": registro["metodo"],
            "estado": registro["estado"],
            "saida": _fora_da_raiz_de(registro["saida"], project_root),
            "via": "isolado",
            "registrado_em": agora,
        }
        if registro.get("alinhador") is not None:
            entrada["alinhador"] = registro["alinhador"]
        if registro.get("motivo"):
            entrada["motivo"] = ExecutionManifest._CAMINHO_EM_TEXTO.sub(
                lambda m: _fora_da_raiz_de(m.group(1), project_root), registro["motivo"])
        metodos_existentes.append(entrada)
    dados["inference_methods"] = metodos_existentes

    tools_existentes = dados.get("tools_invoked") or {}
    for ferramenta, execucao in tool_runs.execucoes().items():
        entrada = tools_existentes.setdefault(ferramenta, {"runs": []})
        for chamada in execucao.get("runs", []):
            comando = [_fora_da_raiz_de(t, project_root) if os.sep in t else t
                       for t in chamada.get("command", [])]
            registro_chamada = {"command": comando, "via": "isolado", "registrado_em": agora}
            if chamada.get("saida"):
                registro_chamada["saida"] = _fora_da_raiz_de(chamada["saida"], project_root)
            entrada["runs"].append(registro_chamada)
    dados["tools_invoked"] = tools_existentes

    outputs_existentes = dados.get("outputs_sha256") or {}
    if os.path.isdir(diretorio_saida):
        for raiz, subdirs, arquivos in os.walk(diretorio_saida):
            subdirs[:] = sorted(subdirs)
            for arquivo in sorted(arquivos):
                caminho_arquivo = os.path.join(raiz, arquivo)
                outputs_existentes[_relativo_a(caminho_arquivo, project_root)] = (
                    file_digest(caminho_arquivo))
    dados["outputs_sha256"] = outputs_existentes

    with open(caminho, "w", encoding="utf-8") as f:
        json.dump(dados, f, indent=2, ensure_ascii=False)
        f.write("\n")
    return True


def _git(repo: str, *args: str) -> Optional[str]:
    saida = _executar(["git", "-C", repo, *args]).strip()
    return saida or None


def _git_state(repo: str) -> Dict[str, Optional[object]]:
    """Commit, ramo e se havia mudança não commitada no momento da execução."""
    if not os.path.isdir(repo):
        return {"commit": None, "branch": None, "dirty": None}
    sujo = _git(repo, "status", "--porcelain")
    return {
        "commit": _git(repo, "rev-parse", "HEAD"),
        "branch": _git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(sujo) if sujo is not None else None,
    }


def _environment() -> Dict[str, object]:
    """
    Ambiente de execução, **sem identificar a máquina nem o usuário**.

    Registra o que muda resultado ou desempenho — sistema, arquitetura, número
    de núcleos, memória, versão do Python — e omite `hostname`, usuário e
    caminhos absolutos, que são dado de terceiro num artefato publicado (D15).
    """
    try:
        memoria_gb = round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9, 1)
    except (ValueError, OSError, AttributeError):
        memoria_gb = None

    return {
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "processor": platform.processor() or None,
        "cpu_count_logical": os.cpu_count(),
        "memory_gb": memoria_gb,
        "python": sys.version.split()[0],
    }


class ExecutionManifest:
    """
    Manifesto de uma execução do workflow.

    Uso::

        manifesto = ExecutionManifest(project_root, params)
        manifesto.register_input(caminho_do_fasta)
        manifesto.write()                       # parcial, já no início
        ...                                     # o pipeline roda
        manifesto.register_outputs(diretorio_de_saida)
        manifesto.finish()                      # grava a versão final

    Gravar cedo é deliberado: se a execução morrer no meio — e
    [D17](../../../docs/science/02-defeitos-que-alteram-resultado.md#d17) mostra
    que morre —, o manifesto parcial diz em que ambiente ela morreu.

    Parameters
    ----------
    project_root : str
        Raiz do projeto (`out/` fica dentro dela). Todo caminho gravado é
        relativo a ela.
    params : dict
        Configuração do workflow, como recebida por `workflow.py`.
    repos : dict of str -> str, optional
        ``nome -> caminho`` dos repositórios cujo commit deve ser registrado.
    """

    def __init__(self, project_root: str, params: Dict,
                 repos: Optional[Dict[str, str]] = None) -> None:
        self.project_root = os.path.abspath(project_root)
        self.params = params
        self.repos = repos or {}
        self._inputs: Dict[str, Optional[str]] = {}
        self._outputs: Dict[str, Optional[str]] = {}
        self._tools_effective: Dict[str, Dict] = {}
        self._reproducibility: Dict[str, int] = {}
        self._execution_mode: Dict[str, object] = {}
        self._inference_methods: List[Dict] = []
        self._model_selection: List[Dict] = []
        self._outcome: Optional[Dict[str, Optional[str]]] = None
        self._log_file: Optional[str] = None

        self.run_id = uuid.uuid4().hex
        self.started_at = datetime.datetime.now(datetime.timezone.utc)
        self.finished_at: Optional[datetime.datetime] = None

    # ------------------------------------------------------------------ #
    # Coleta
    # ------------------------------------------------------------------ #

    def _relativo(self, path: str) -> str:
        """Caminho relativo à raiz do projeto; nunca absoluto (D15)."""
        return _relativo_a(path, self.project_root)

    def _fora_da_raiz(self, path: str) -> str:
        """
        Como `_relativo`, mas reduz a **nome do arquivo** o que estiver fora do
        projeto.

        `os.path.relpath` de um caminho fora da raiz devolve `../../..` até
        alcançá-lo, e o resto do caminho absoluto vai junto: o binário do
        RAxML-NG mora em `<home do usuário>/miniconda3/envs/...`, de modo que
        relativizá-lo grava o nome do usuário no manifesto. É exatamente
        [D15](../../../docs/science/02-defeitos-que-alteram-resultado.md#d15),
        que já vazou `/home/<usuário>` de um terceiro pela API.

        Dentro do projeto, caminho relativo — é informação do experimento.
        Fora, só o nome — a versão de cada ferramenta já está em
        `tools_available`, e *onde* ela estava instalada não é reproduzível
        noutra máquina de qualquer forma.
        """
        return _fora_da_raiz_de(path, self.project_root)

    def _sanitizar_comando(self, comando: List[str]) -> List[str]:
        """
        Linha de comando sem caminho absoluto, preservando o que ela informa.

        Um token é tratado como caminho quando tem separador de diretório: os
        parâmetros que decidem o resultado (`--threads`, `--seed`, `GTR+G`)
        nunca têm, e passam intactos.
        """
        return [self._fora_da_raiz(token) if os.sep in token else token
                for token in comando]

    #: Caminho absoluto dentro de texto livre: começa em `/` depois de início
    #: de linha, espaço ou delimitador, e vai até o próximo espaço/aspa.
    _CAMINHO_EM_TEXTO = re.compile(r"(?:^|(?<=[\s'\"(\[=,:]))(/[^\s'\"\]),]+)", re.MULTILINE)

    def _sanitizar_texto(self, texto: str) -> str:
        """
        Texto livre sem caminho absoluto.

        O motivo de uma falha vem de exceção e de `stderr` de ferramenta, e os
        dois citam caminho absoluto — `FileNotFoundError` do `require_tool`
        cita o PATH, e `CalledProcessError` cita a linha de comando inteira,
        com o binário dentro do env conda do usuário. Gravar cru é
        [D15](../../../docs/science/02-defeitos-que-alteram-resultado.md#d15).
        """
        return self._CAMINHO_EM_TEXTO.sub(lambda m: self._fora_da_raiz(m.group(1)), texto)

    def _sanitizar_params(self, valor):
        """
        Configuração sem caminho absoluto, em qualquer profundidade.

        O módulo promete, desde a primeira linha, que "todo caminho é relativo à
        raiz do projeto" — e `params` era gravado cru, com `input_path` e
        `output_path` absolutos. O manifesto declarava não vazar D15 e vazava,
        e a conferência não pegava porque só varre as chaves de
        `inputs_sha256`/`outputs_sha256`.

        Um valor é tratado como caminho quando é absoluto: um parâmetro comum
        (`"advanced"`, `"GTR+G"`, `"nexus"`) nunca é.
        """
        if isinstance(valor, dict):
            return {k: self._sanitizar_params(v) for k, v in valor.items()}
        if isinstance(valor, list):
            return [self._sanitizar_params(v) for v in valor]
        if isinstance(valor, str) and os.path.isabs(valor):
            return self._fora_da_raiz(valor)
        return valor

    def register_input(self, path: str, suffixes=(".fasta", ".fa", ".fna", ".csv", ".gb")) -> None:
        """
        Registra uma entrada e seu SHA-256.

        Aceita arquivo **ou diretório**: o `input_path` do `tree_config` é um
        diretório de FASTA, e registrar só o caminho dele não diria nada sobre o
        conteúdo — que é o que precisa ser idêntico para a execução ser a mesma.
        """
        if os.path.isdir(path):
            for raiz, subdirs, arquivos in os.walk(path):
                subdirs[:] = sorted(subdirs)
                for arquivo in sorted(arquivos):
                    if arquivo.endswith(suffixes):
                        caminho = os.path.join(raiz, arquivo)
                        self._inputs[self._relativo(caminho)] = file_digest(caminho)
        else:
            self._inputs[self._relativo(path)] = file_digest(path)

    def register_outputs(self, directory: str, suffixes=(".nexus", ".nwk", ".aln", ".csv")) -> None:
        """
        Registra o SHA-256 de toda saída relevante encontrada em `directory`.

        `metadata.json` e `raw_data_sequences.gb` ficam de fora por padrão: têm
        gigabytes e seu conteúdo é derivado das árvores e do GenBank, que já
        estão cobertos.
        """
        for raiz, subdirs, arquivos in os.walk(directory):
            # Poda `tmp/` na travessia, comparando o caminho RELATIVO ao
            # diretório varrido: testar `"/tmp" in caminho_absoluto` casaria com
            # qualquer projeto guardado sob um diretório chamado `tmp` — o que
            # inclui todo diretório temporário de teste.
            subdirs[:] = [d for d in sorted(subdirs) if d != "tmp"]
            for arquivo in sorted(arquivos):
                if arquivo.endswith(suffixes):
                    caminho = os.path.join(raiz, arquivo)
                    self._outputs[self._relativo(caminho)] = file_digest(caminho)

    def register_log(self, path: str) -> None:
        """
        Registra o arquivo de log desta execução.

        É o que faz manifesto e log apontarem um para o outro: a pergunta "que
        log produziu esta árvore" passa a ter resposta exata, em vez de "o mais
        recente por data de modificação" — que escolhia entre execuções
        diferentes sem dizer qual (D22).
        """
        self._log_file = self._relativo(path) if path else None

    def register_reproducibility(self, settings: Dict[str, int]) -> None:
        """
        Registra semente e paralelização efetivas da execução.

        Deve receber o resultado de `builder.reproducibility_settings`, e não
        valores recalculados: se o manifesto declarar um número e o pipeline
        usar outro, ele deixa de ser manifesto e passa a ser ficção.
        """
        self._reproducibility = dict(settings)

    def register_execution_mode(self, requested: str, methods_advanced_available: List[str],
                                 methods_advanced_executed: List[str]) -> None:
        """
        Registra o modo pedido contra os métodos avançados de fato executados.

        [D18](../../../docs/science/02-defeitos-que-alteram-resultado.md#d18):
        o modo `auto`/`basic` roda só distância e parcimônia, e nada no
        artefato dizia isso além do `mode` cru em `config_backup.json` — a
        mesma distinção "disponível não é executado" de `tools_invoked`
        (`tool_runs.py`), agora também para a escolha de método, não só de
        ferramenta. `M`, o denominador de todo suporte metodológico, deixa de
        ser um número sem proveniência.

        Parameters
        ----------
        requested : str
            O `mode` cru como pedido (`"auto"`, `"basic"` ou `"advanced"`).
        methods_advanced_available : list of str
            Nomes de **método do controlador** (`iqtree`, `fasttree`, `raxml`,
            `mrbayes` — não a chave de ferramenta `external_tools.CANDIDATOS`,
            que para o RAxML é `"raxml-ng"`). Vem de `external_tools.resolve_tool`
            por método, não de suposição; quem procurar `"raxml-ng"` aqui não
            acha — `workflow.py` já faz essa tradução antes de chamar isto.
        methods_advanced_executed : list of str
            Os que de fato rodaram. Vazio em modo básico; igual a
            `methods_advanced_available` menos `ignore_methods` em `advanced`.
        """
        self._execution_mode = {
            "mode_solicitado": requested,
            "metodos_avancados_disponiveis": list(methods_advanced_available),
            "metodos_avancados_executados": list(methods_advanced_executed),
            "metodos_avancados_pulados": sorted(
                set(methods_advanced_available) - set(methods_advanced_executed)
            ),
        }

    def register_tool_run(self, tool: str, command: List[str],
                          saida: Optional[str] = None, **extra) -> None:
        """
        Registra a linha de comando efetiva de uma chamada de ferramenta.

        É o que responde "com que semente e com que paralelização esta árvore
        foi feita" — e, depois de D17, sabe-se que a paralelização **muda a
        topologia** mesmo com a semente fixa.

        Chamar mais de uma vez para a mesma ferramenta **acumula**: um
        delineamento com dois alinhadores invoca o RAxML-NG duas vezes, e
        guardar só a última chamada seria declarar como único o comando que
        produziu metade das árvores.

        Os caminhos são higienizados por `_sanitizar_comando` antes de entrar:
        o manifesto vai para o repositório e pode ir para material suplementar.
        """
        entrada = self._tools_effective.setdefault(tool, {"runs": []})
        entrada.update({k: v for k, v in extra.items() if v is not None})
        chamada: Dict = {"command": self._sanitizar_comando(command)}
        if saida is not None:
            chamada["saida"] = self._fora_da_raiz(saida)
        entrada["runs"].append(chamada)

    def drain_tool_runs(self) -> None:
        """
        Recolhe em `tool_runs` o que o pipeline registrou durante a execução.

        Chamado no encerramento, inclusive quando a execução falhou: saber qual
        comando estava rodando é metade do diagnóstico.
        """
        for ferramenta, dados in tool_runs.execucoes().items():
            chamadas = dados.pop("runs", [])
            for chamada in chamadas:
                extras_da_chamada = {k: v for k, v in chamada.items()
                                     if k not in ("command", "saida")}
                self.register_tool_run(ferramenta, chamada["command"],
                                       saida=chamada.get("saida"), **dados)
                # Fatos anotados depois da chamada (`tool_runs.anotar`) são da
                # chamada, não da ferramenta: o ASDSF de uma cadeia do MrBayes
                # não é o da outra.
                if extras_da_chamada:
                    self._tools_effective[ferramenta]["runs"][-1].update(extras_da_chamada)

        for registro in tool_runs.metodos():
            self.register_method_state(**registro)

        for selecao in tool_runs.selecoes_modelo():
            self.register_model_selection(selecao)

    def register_model_selection(self, selecao: Dict) -> None:
        """
        Registra a seleção de modelo de um alinhamento — M7.3.

        `alinhamento` e `relatorio` são caminhos (relativizados, D15); o
        `comando` é higienizado como os de `tools_invoked`; o `motivo` de uma
        falha, como o de `inference_methods`.
        """
        registro = dict(selecao)
        for chave in ("alinhamento", "relatorio"):
            if registro.get(chave):
                registro[chave] = self._fora_da_raiz(registro[chave])
        if registro.get("comando"):
            registro["comando"] = self._sanitizar_comando(registro["comando"])
        if registro.get("motivo"):
            registro["motivo"] = self._sanitizar_texto(registro["motivo"])
        self._model_selection.append(registro)

    def register_method_state(self, metodo: str, estado: str, saida: str,
                              alinhador: Optional[str] = None,
                              motivo: Optional[str] = None) -> None:
        """
        Registra o desfecho de um pipeline — M7.6.

        `execution_mode` é o **plano**, calculado antes de rodar; este é o
        **fato**. Uma execução que pediu RAxML-NG e viu o RAxML-NG morrer tem
        o método em `metodos_avancados_executados` do plano e aqui como
        ``tentado_e_falhou``, com o motivo.
        """
        if estado not in tool_runs.ESTADOS_METODO:
            raise ValueError(f"Estado de método desconhecido: {estado!r}")
        registro: Dict[str, str] = {"metodo": metodo, "estado": estado,
                                    "saida": self._fora_da_raiz(saida)}
        if alinhador is not None:
            registro["alinhador"] = alinhador
        if motivo:
            registro["motivo"] = self._sanitizar_texto(motivo)
        self._inference_methods.append(registro)

    def register_outcome(self, status: str, motivo: Optional[str] = None) -> None:
        """
        Desfecho da execução inteira: ``concluido`` ou ``falhou``.

        O manifesto é fechado no `finally` de `workflow.py`, então
        `finished_at_utc` existe também quando a execução morreu — e até aqui
        quem separava as duas coisas era o log (`execution_state.py` do
        backend diz isso num comentário). Um manifesto que conclui sem dizer
        se deu certo é meia declaração.
        """
        if status not in ("concluido", "falhou"):
            raise ValueError(f"Desfecho desconhecido: {status!r}")
        self._outcome = {"status": status,
                         "motivo": self._sanitizar_texto(motivo) if motivo else None}

    # ------------------------------------------------------------------ #
    # Serialização
    # ------------------------------------------------------------------ #

    def to_dict(self) -> Dict:
        # Versão 2: `tools_invoked` deixou de ser `ferramenta -> {command, ...}`
        # e passou a ser `ferramenta -> {parâmetros, "runs": [...]}`. A forma
        # antiga guardava uma chamada por ferramenta e perdia as demais — e,
        # como o campo nunca chegou a ser populado (DEC-045), nenhum artefato
        # em disco tem a forma antiga.
        return {
            "manifest_version": 2,
            "run_id": self.run_id,
            "started_at_utc": self.started_at.isoformat(),
            "finished_at_utc": self.finished_at.isoformat() if self.finished_at else None,
            "git": {nome: _git_state(caminho) for nome, caminho in self.repos.items()},
            "environment": _environment(),
            "tools_available": dict(tool_versions()),
            "tools_invoked": self._tools_effective,
            "log_file": self._log_file,
            "reproducibility": dict(self._reproducibility),
            "execution_mode": dict(self._execution_mode) or None,
            # M7.6 — campos acrescentados sem mudar a forma de nenhum outro,
            # por isso `manifest_version` continua 2: um leitor da v2 que não
            # os conhece segue lendo o resto igual. `None` = a execução ainda
            # não chegou lá (manifesto parcial), não "nenhum método".
            "inference_methods": ([dict(m) for m in self._inference_methods]
                                  if self._inference_methods else None),
            "outcome": dict(self._outcome) if self._outcome else None,
            # M7.3 — uma entrada por alinhamento: modelo escolhido pelo BIC,
            # candidatos e a tradução para cada ferramenta. Acrescentado, como
            # os de M7.6; `None` = nenhum método baseado em modelo chegou a
            # precisar de seleção (ou manifesto parcial).
            "model_selection": ([dict(m) for m in self._model_selection]
                                if self._model_selection else None),
            "params": self._sanitizar_params(self.params),
            "inputs_sha256": self._inputs,
            "outputs_sha256": self._outputs,
        }

    def path(self) -> str:
        return os.path.join(self.project_root, "out", "outputs", MANIFEST_FILENAME)

    def write(self) -> str:
        """Grava o manifesto no estado atual. Idempotente."""
        destino = self.path()
        os.makedirs(os.path.dirname(destino), exist_ok=True)
        with open(destino, "w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2, ensure_ascii=False, sort_keys=True)
        return destino

    def finish(self) -> str:
        """Fecha o manifesto com o horário de término e grava."""
        self.finished_at = datetime.datetime.now(datetime.timezone.utc)
        destino = self.write()
        logging.info(f"Manifesto de execução gravado em {self._relativo(destino)} "
                     f"(run_id {self.run_id})")
        return destino
