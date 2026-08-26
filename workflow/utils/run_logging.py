"""
Um arquivo de log por execução — item 4 de
[D22](../../../docs/science/02-defeitos-que-alteram-resultado.md#d22).

O log chamava-se `log_setup_{ano}_{mês}_{dia}.log` e era aberto em modo
*append*. Duas execuções no mesmo dia caíam no **mesmo arquivo**, e quem
tentasse medir a duração lendo do primeiro ao último carimbo somava as duas
mais o intervalo ocioso entre elas: 1 960 s onde a última execução levara
396 s. Nos artefatos em disco há dois arquivos com **dois**
`Completed successfully!` dentro.

O consumidor foi ensinado a recortar o log nas fronteiras de execução
([DEC-048](../../../docs/automation/07-log-de-execucao.md)), mas isso é
heurística: separa duas execuções quando a primeira **concluiu** e não separa
quando ela morreu no meio. A separação de verdade é na origem, e é esta —
**o nome do arquivo carrega o `run_id`**, de modo que duas execuções nunca
disputam o mesmo arquivo.

O `run_id` é o mesmo do `manifest.json`: um arquivo de log e um manifesto
apontam um para o outro, e a pergunta "que log produziu esta árvore" passa a
ter resposta exata em vez de "o mais recente por data de modificação".

**Por que centralizar.** As mesmas três linhas de `logging.basicConfig` estavam
copiadas em **sete** módulos, cada uma recalculando o nome do arquivo. Sete
lugares que precisam concordar sobre um nome é a forma de
[D5](../../../docs/science/02-defeitos-que-alteram-resultado.md#d5) noutro
assunto — e, pior, `basicConfig` só tem efeito na **primeira** chamada, então
qual das sete vencia dependia da ordem de importação.
"""

from __future__ import annotations

import datetime
import logging
import os
from typing import Optional

__all__ = ["configurar", "garantir", "caminho_atual", "FORMATO"]

FORMATO = "%(asctime)s - %(levelname)s - %(message)s"

#: Caminho do arquivo em uso, para o manifesto registrar. `None` até configurar.
_ATUAL: Optional[str] = None


def caminho_atual() -> Optional[str]:
    """Arquivo de log em uso nesta execução, ou `None`."""
    return _ATUAL


def _nome(run_id: Optional[str]) -> str:
    """
    `log_setup_{AAAA-MM-DD}_{run_id}.log`.

    A data continua no nome porque é o que torna um diretório de logs legível
    por quem o abre; o `run_id` é o que torna o nome **único**. Sem `run_id` —
    módulo usado fora do workflow — cai no nome antigo, que é o comportamento
    que já existia.
    """
    hoje = datetime.date.today().isoformat()
    if not run_id:
        agora = datetime.datetime.now()
        return f"log_setup_{agora.year}_{agora.month}_{agora.day}.log"
    return f"log_setup_{hoje}_{run_id[:12]}.log"


def _instalar(caminho: str) -> None:
    """
    Põe um `FileHandler` na raiz, substituindo o que houver.

    Não usa `logging.basicConfig`: ele é no-op quando a raiz já tem handler, e
    era justamente isso que fazia a configuração depender de quem importasse
    primeiro. Aqui a substituição é explícita.
    """
    raiz = logging.getLogger()
    for handler in list(raiz.handlers):
        raiz.removeHandler(handler)
        try:
            handler.close()
        except Exception:      # handler de terceiro que não fecha limpo
            pass
    manipulador = logging.FileHandler(caminho, encoding="utf-8")
    manipulador.setFormatter(logging.Formatter(FORMATO))
    raiz.addHandler(manipulador)
    raiz.setLevel(logging.INFO)


def configurar(outputs_dir: str, run_id: Optional[str] = None) -> str:
    """
    Abre o log **desta** execução e o torna o destino de todo `logging`.

    Parameters
    ----------
    outputs_dir : str
        Diretório `out/outputs` do projeto.
    run_id : str or None
        O mesmo do manifesto. Sem ele o nome cai no formato antigo, por data.

    Return
    ------
    str
        Caminho do arquivo aberto.
    """
    global _ATUAL
    os.makedirs(outputs_dir, exist_ok=True)
    caminho = os.path.join(outputs_dir, _nome(run_id))
    _instalar(caminho)
    _ATUAL = caminho
    return caminho


def garantir(outputs_dir: str) -> str:
    """
    Configura o log **só se ninguém o tiver configurado**.

    É o que os módulos chamam no lugar das sete cópias de `basicConfig`. Dentro
    de uma execução do workflow, `configurar` já rodou e esta função não toca em
    nada — o log da execução não é trocado no meio por um módulo importado
    depois. Fora dela, cada módulo continua tendo para onde escrever.
    """
    if _ATUAL is not None or logging.getLogger().handlers:
        return _ATUAL or ""
    return configurar(outputs_dir)
