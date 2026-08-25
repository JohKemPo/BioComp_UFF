"""
Resolução dos binários das ferramentas externas.

**O nome do binário não acompanha o nome do pacote de forma estável.** O pacote
`iqtree` do bioconda instalava `iqtree2` na série 2.x e passou a instalar
`iqtree` e `iqtree3` na 3.x — sem `iqtree2`. O pipeline chamava `iqtree2` fixo,
então quem seguisse o `environment.yml` do projeto instalava o IQ-TREE com
sucesso e mesmo assim não conseguia rodar: nada encontrava o binário.

O mesmo vale para o FastTree, que aparece como `FastTree` ou `fasttree`
conforme a distribuição, e para o MrBayes, que é `mb` na maioria das e
`mrbayes` em algumas.

Este módulo existe para que **haja um lugar só** que saiba os nomes possíveis
de cada ferramenta, consultado tanto pelo pipeline quanto pelo manifesto.
Duas listas de nomes divergindo seria o defeito D5 em outro assunto.
"""

from __future__ import annotations

import os
import shutil
from typing import Dict, Optional, Sequence, Tuple

__all__ = ["CANDIDATOS", "resolve_tool", "require_tool", "resolved_tools"]


#: Nomes possíveis de cada ferramenta, em ordem de preferência.
#:
#: A ordem importa: prefere-se a série mais nova quando ambas existem, porque é
#: a que o `environment.yml` instala. Um binário antigo no PATH do sistema não
#: deve ganhar do binário do ambiente do projeto.
CANDIDATOS: Dict[str, Sequence[str]] = {
    "iqtree": ("iqtree3", "iqtree2", "iqtree"),
    "fasttree": ("FastTree", "fasttree", "FastTreeMP"),
    "raxml-ng": ("raxml-ng",),
    "mrbayes": ("mb", "mrbayes"),
    "mafft": ("mafft",),
    "clustalo": ("clustalo",),
    "muscle": ("muscle",),
}


def resolve_tool(tool: str, prefer_dir: Optional[str] = None) -> Optional[str]:
    """
    Encontra o binário de uma ferramenta entre os nomes possíveis.

    Parameters
    ----------
    tool : str
        Chave de `CANDIDATOS` — ``"iqtree"``, ``"fasttree"``, ``"mrbayes"``…
    prefer_dir : str or None, optional
        Diretório consultado **antes** do PATH. Serve para privilegiar o `bin`
        do ambiente do projeto: nesta base, o env tem FastTree 2.2.0 enquanto o
        PATH resolvia 2.1.11 de `/usr/bin`, e medir o segundo chamando-o de
        "versão do projeto" produziu um registro errado que durou dias.

    Return
    ------
    str or None
        Caminho do executável, ou `None` quando nenhum candidato existe.
    """
    nomes = CANDIDATOS.get(tool, (tool,))

    if prefer_dir and os.path.isdir(prefer_dir):
        for nome in nomes:
            caminho = os.path.join(prefer_dir, nome)
            if os.path.isfile(caminho) and os.access(caminho, os.X_OK):
                return caminho

    for nome in nomes:
        caminho = shutil.which(nome)
        if caminho:
            return caminho

    return None


def require_tool(tool: str, prefer_dir: Optional[str] = None) -> str:
    """
    Como `resolve_tool`, mas levanta com uma mensagem acionável.

    Raises
    ------
    FileNotFoundError
        Nomeando os candidatos procurados e como instalar. Um "command not
        found" cru não diz a quem o lê que o pacote pode estar instalado sob
        outro nome de binário.
    """
    caminho = resolve_tool(tool, prefer_dir)
    if caminho:
        return caminho

    nomes = ", ".join(CANDIDATOS.get(tool, (tool,)))
    raise FileNotFoundError(
        f"{tool} não encontrado. Procurei por: {nomes}. "
        f"Instale o ambiente do projeto com `bash scripts/setup_env.sh` e ative-o, "
        f"ou confira com `bash scripts/check_dependencies.sh`."
    )


def resolved_tools(prefer_dir: Optional[str] = None) -> Dict[str, Optional[str]]:
    """``ferramenta -> caminho`` de todas, para registro no manifesto."""
    return {tool: resolve_tool(tool, prefer_dir) for tool in CANDIDATOS}
