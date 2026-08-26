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

__all__ = ["registrar", "execucoes", "limpar"]

#: `ferramenta -> {parâmetros comuns, "runs": [chamada, ...]}`.
_EXECUCOES: Dict[str, Dict] = {}

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


def limpar() -> None:
    """Esvazia o registro. Para uso em teste."""
    with _TRAVA:
        _EXECUCOES.clear()
