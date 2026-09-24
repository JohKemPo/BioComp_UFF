"""
Registro das chamadas de ferramenta externa feitas durante uma execução.

O manifesto sempre soube dizer **o que estava instalado** (`tools_available`) e
nunca **o que rodou**. `ExecutionManifest.register_tool_run` existe desde M2.5,
com teste de unidade, e nenhum ponto do pipeline a chamava: o campo
`tools_invoked` saía vazio de toda execução (DEC-045). É a mesma distinção que
[D18](../../../docs/science/02-defeitos-que-alteram-resultado.md#d18) cobra do
modo `auto` — disponível não é executado.

Este módulo é o coletor. Quem invoca uma ferramenta registra aqui; o manifesto
drena no fim. A separação existe porque as duas pontas têm donos diferentes:

- **Quem chama registra o fato bruto.** O `TreeBuilder` é construído seis vezes
  por execução, do fundo da pilha do controlador, e não tem — nem deve ter —
  uma referência ao manifesto.
- **O manifesto decide o que pode ser gravado.** A linha de comando efetiva
  contém caminhos absolutos, e o binário resolvido mora dentro do ambiente
  conda do usuário: gravá-la crua reintroduz
  [D15](../../../docs/science/02-defeitos-que-alteram-resultado.md#d15). A
  higienização é de lá, não daqui.

O registro é **de processo**: uma execução do pipeline é um processo, e é o
escopo natural do manifesto que ela grava. Em teste, chame `limpar()` no
`setUp` — é o preço explícito dessa escolha.
"""

from __future__ import annotations

import threading
from typing import Dict, List, Optional

__all__ = ["registrar", "execucoes", "anotar", "registrar_metodo", "metodos",
           "registrar_selecao_modelo", "selecoes_modelo", "ESTADOS_SELECAO_MODELO",
           "limpar", "ESTADOS_METODO"]

#: `ferramenta -> {parâmetros comuns, "runs": [chamada, ...]}`.
_EXECUCOES: Dict[str, Dict] = {}

#: M7.6 — desfecho de cada método de inferência em cada braço do delineamento.
#: `tools_invoked` diz o que foi **chamado**; isto diz o que cada pipeline
#: **deu**. Os dois não coincidem: uma chamada registrada pode ter morrido, e
#: um método ignorado nunca chega a chamar ferramenta nenhuma — e era aí que
#: "excluído de propósito" e "quebrou e foi excluído depois" se confundiam
#: (D18, DM-11).
_METODOS: List[Dict] = []

#: M7.3 — seleção de modelo de substituição, uma por alinhamento. É um fato
#: **do alinhamento**, não de um método: IQ-TREE, RAxML-NG e MrBayes do mesmo
#: braço leem a mesma escolha, e é isso que permite dizer que rodaram o mesmo
#: modelo.
_SELECOES_MODELO: List[Dict] = []

#: ``concluida`` (modelo escolhido e traduzido), ``falhou`` (motivo
#: obrigatório; os métodos que dependiam dela ficam ``tentado_e_falhou``) ou
#: ``desligada_por_configuracao`` (`model_selection: "nenhuma"` — GTR+G fixo).
ESTADOS_SELECAO_MODELO = ("concluida", "falhou", "desligada_por_configuracao")

#: Os quatro desfechos possíveis de um pipeline numa execução.
#:
#: - ``executado``: a árvore foi produzida **nesta** execução.
#: - ``reaproveitado``: a árvore já estava em disco e foi lida, não produzida.
#:   Não é o mesmo que executado — a chamada registrada em `tools_invoked` de
#:   outra execução é que a produziu, com a semente e a versão de então.
#: - ``ignorado_por_configuracao``: o experimento pediu para não rodar
#:   (`ignore_mode`, ou `mode` básico para os métodos avançados).
#: - ``tentado_e_falhou``: rodou e não produziu árvore utilizável.
ESTADOS_METODO = ("executado", "reaproveitado", "ignorado_por_configuracao",
                  "tentado_e_falhou")

#: Estados que só fazem sentido acompanhados de um porquê.
_EXIGEM_MOTIVO = ("ignorado_por_configuracao", "tentado_e_falhou", "reaproveitado")

#: Teto do motivo. Um `stderr` de ferramenta pode ter megabytes; o manifesto
#: precisa do diagnóstico, não do log inteiro — o log da execução está ao lado.
_MOTIVO_MAX = 600

#: O pipeline é sequencial hoje, mas `parsl` está nas dependências e a
#: paralelização por pipeline é um caminho previsto. Um dicionário perdendo
#: escrita por corrida seria uma chamada some do manifesto sem nenhum sinal.
_TRAVA = threading.Lock()


def registrar(ferramenta: str, comando: List[str],
              saida: Optional[str] = None, **parametros) -> None:
    """
    Registra uma chamada de ferramenta externa.

    Parameters
    ----------
    ferramenta : str
        Chave de `external_tools.CANDIDATOS` — ``"raxml-ng"``, ``"iqtree"``,
        ``"mafft"``… É a mesma chave de `tools_available`, para que as duas
        metades do manifesto se cruzem sem tradução.
    comando : list of str
        Linha de comando efetiva, como entregue ao `subprocess`. Crua: a
        higienização de caminho é do manifesto.
    saida : str or None, optional
        Arquivo que esta chamada produz. É o que liga a chamada ao resultado —
        sem ele, "com que semente esta árvore foi feita" continua sem resposta
        quando a mesma ferramenta roda em vários braços do delineamento.
    **parametros
        Parâmetros que decidem o resultado e não se leem da linha de comando
        sem interpretá-la: ``seed``, ``threads``, ``workers``, ``model``.
        Ficam no nível da ferramenta porque são invariantes da execução.

    Notes
    -----
    Registra-se **no momento da invocação**, não depois do sucesso. O manifesto
    de uma execução que morreu no meio é o que permite diagnosticá-la, e
    [D17](../../../docs/science/02-defeitos-que-alteram-resultado.md#d17)
    mostrou que elas morrem — foi um `SIGSEGV` do RAxML-NG.
    """
    with _TRAVA:
        entrada = _EXECUCOES.setdefault(ferramenta, {"runs": []})
        entrada.update({k: v for k, v in parametros.items() if v is not None})
        chamada: Dict = {"command": list(comando)}
        if saida is not None:
            chamada["saida"] = saida
        entrada["runs"].append(chamada)


def execucoes() -> Dict[str, Dict]:
    """Cópia do que foi registrado até agora, para o manifesto drenar."""
    with _TRAVA:
        return {
            ferramenta: {**dados, "runs": [dict(r) for r in dados["runs"]]}
            for ferramenta, dados in _EXECUCOES.items()
        }


def anotar(ferramenta: str, saida: str, **campos) -> bool:
    """
    Acrescenta fatos a uma chamada **já registrada**, identificada pela saída.

    Existe porque alguns fatos só se sabem depois que a ferramenta termina: o
    ASDSF do MrBayes (M7.4), ou que o RAxML-NG caiu no `.raxml.bestTree` sem
    suporte porque o bootstrap não gravou `.raxml.support`. `registrar` é
    chamado **antes** de rodar, de propósito (ver as Notes dele) — e por isso
    não pode carregar esses números.

    Return
    ------
    bool
        ``True`` se achou a chamada. ``False`` não é erro silencioso: quem
        chama decide se isso importa, e o valor de retorno diz que não achou.
    """
    with _TRAVA:
        entrada = _EXECUCOES.get(ferramenta)
        if not entrada:
            return False
        for chamada in reversed(entrada["runs"]):
            if chamada.get("saida") == saida:
                chamada.update({k: v for k, v in campos.items() if v is not None})
                return True
        return False


def registrar_metodo(metodo: str, estado: str, saida: str,
                     alinhador: Optional[str] = None,
                     motivo: Optional[str] = None) -> bool:
    """
    Registra o desfecho de um pipeline (método x alinhador x arquivo).

    Parameters
    ----------
    metodo : str
        Rótulo do método como aparece no nome do arquivo de árvore:
        ``iqtree``, ``raxml``, ``nj_distance``, ``upgma_parsimony``...
    estado : str
        Um de `ESTADOS_METODO`. Outro valor é ``ValueError`` — um estado fora
        da enumeração é exatamente a ambiguidade que este registro existe para
        eliminar.
    saida : str
        Arquivo de árvore que o pipeline produz (ou produziria). É a chave:
        identifica o pipeline sem ambiguidade entre alinhadores e arquivos de
        entrada.
    alinhador : str, optional
    motivo : str, optional
        **Obrigatório** para ``ignorado_por_configuracao``,
        ``tentado_e_falhou`` e ``reaproveitado``.

    Return
    ------
    bool
        ``False`` quando o pipeline já tinha desfecho nesta execução — o
        primeiro prevalece. O modo `advanced` visita cada método avançado uma
        vez por método de distância (NJ, UPGMA): sem isto, toda árvore
        apareceria como "executado" e, logo depois, "reaproveitado".
    """
    if estado not in ESTADOS_METODO:
        raise ValueError(f"Estado de método desconhecido: {estado!r}. "
                         f"Válidos: {ESTADOS_METODO}")
    if estado in _EXIGEM_MOTIVO and not (motivo and motivo.strip()):
        raise ValueError(f"O estado {estado!r} exige motivo — é o que o distingue "
                         f"dos outros (M7.6).")
    with _TRAVA:
        if any(m["saida"] == saida for m in _METODOS):
            return False
        registro: Dict = {"metodo": metodo, "estado": estado, "saida": saida}
        if alinhador is not None:
            registro["alinhador"] = alinhador
        if motivo:
            motivo = motivo.strip()
            if len(motivo) > _MOTIVO_MAX:
                motivo = motivo[:_MOTIVO_MAX] + " [...]"
            registro["motivo"] = motivo
        _METODOS.append(registro)
        return True


def metodos() -> List[Dict]:
    """Cópia dos desfechos registrados, na ordem em que aconteceram."""
    with _TRAVA:
        return [dict(m) for m in _METODOS]


def registrar_selecao_modelo(alinhamento: str, estado: str, **campos) -> bool:
    """
    Registra a seleção de modelo de um alinhamento — M7.3.

    Parameters
    ----------
    alinhamento : str
        Arquivo de alinhamento (a chave; o primeiro registro prevalece, como em
        `registrar_metodo`).
    estado : str
        Um de `ESTADOS_SELECAO_MODELO`. ``falhou`` e
        ``desligada_por_configuracao`` exigem ``motivo``.
    **campos
        ``modelo_escolhido``, ``traducao``, ``candidatos``, ``relatorio``,
        ``motivo``… — o manifesto higieniza caminhos.
    """
    if estado not in ESTADOS_SELECAO_MODELO:
        raise ValueError(f"Estado de seleção de modelo desconhecido: {estado!r}. "
                         f"Válidos: {ESTADOS_SELECAO_MODELO}")
    motivo = campos.get("motivo")
    if estado != "concluida" and not (motivo and str(motivo).strip()):
        raise ValueError(f"O estado {estado!r} exige motivo (M7.3).")
    with _TRAVA:
        if any(s["alinhamento"] == alinhamento for s in _SELECOES_MODELO):
            return False
        registro = {"alinhamento": alinhamento, "estado": estado}
        registro.update({k: v for k, v in campos.items() if v is not None})
        if isinstance(registro.get("motivo"), str) and len(registro["motivo"]) > _MOTIVO_MAX:
            registro["motivo"] = registro["motivo"][:_MOTIVO_MAX] + " [...]"
        _SELECOES_MODELO.append(registro)
        return True


def selecoes_modelo() -> List[Dict]:
    """Cópia das seleções de modelo registradas, na ordem em que aconteceram."""
    with _TRAVA:
        return [dict(s) for s in _SELECOES_MODELO]


def limpar() -> None:
    """Esvazia o registro. Para uso em teste."""
    with _TRAVA:
        _EXECUCOES.clear()
        _METODOS.clear()
        _SELECOES_MODELO.clear()
