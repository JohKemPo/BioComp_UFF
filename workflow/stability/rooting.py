"""
Enraizamento explícito por grupo externo declarado — lote M2.3.

FastTree, IQ-TREE, RAxML e NJ emitem topologia **não enraizada**, escrita em
Newick com raiz trifurcante. A raiz que se lê do arquivo é convenção de escrita
da ferramenta, não hipótese biológica; o UPGMA é o único método do conjunto que
produz raiz, e a produz impondo um relógio molecular que os dados violam.

[D3](../../../docs/science/02-defeitos-que-alteram-resultado.md#d3) resolveu a
metade defensiva do problema: a unidade de comparação passou a ser a
**bipartição**, que ignora onde a raiz caiu. Este módulo resolve a outra metade
— quando a análise enraizada é o que se quer, ela exige **enraizamento
explícito e comum a todos os métodos, pelo grupo externo declarado**, nunca a
raiz arbitrária do arquivo.

Três regras governam tudo aqui, e todas existem porque enraizar errado é pior
que não enraizar:

1. **O grupo externo é declarado, nunca inferido.** Ele vem da configuração do
   experimento. Adivinhar qual táxon é externo a partir do próprio dado é
   circular.
2. **Grupo externo não monofilético não enraíza.** Se os táxons do grupo externo
   não formam um clado na árvore, existe mais de uma aresta candidata a raiz e a
   escolha seria arbitrária. O resultado é `None` com o motivo declarado — não
   uma árvore enraizada em algum lugar plausível.
3. **Ou todos, ou nenhum.** Se um método não pode ser enraizado, a análise
   enraizada daquele conjunto não é comparável. `root_tree_set` devolve o
   relatório completo para que a decisão seja tomada com o quadro à vista.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Callable, Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

from Bio.Phylo.BaseTree import Tree

from workflow.stability.clade_identity import strip_accession_version

__all__ = [
    "RootingReport",
    "root_at_outgroup",
    "root_tree_set",
    "outgroup_from_classifier",
]


@dataclass
class RootingReport:
    """
    O que aconteceu ao tentar enraizar uma árvore.

    Attributes
    ----------
    rooted : bool
        Se a árvore foi efetivamente enraizada.
    reason : str or None
        Motivo da recusa, quando `rooted` é False. `None` em caso de sucesso.
    outgroup_found : frozenset of str
        Táxons do grupo externo presentes na árvore, já normalizados.
    outgroup_missing : frozenset of str
        Táxons do grupo externo declarados e **ausentes** da árvore.
    monophyletic : bool or None
        Se os táxons do grupo externo formam um clado. `None` quando não havia
        táxons suficientes para a pergunta fazer sentido.
    """

    rooted: bool
    reason: Optional[str] = None
    outgroup_found: FrozenSet[str] = field(default_factory=frozenset)
    outgroup_missing: FrozenSet[str] = field(default_factory=frozenset)
    monophyletic: Optional[bool] = None

    def summary(self) -> Dict[str, object]:
        """Resumo achatado, pronto para o manifesto ou uma tabela."""
        return {
            "rooted": self.rooted,
            "reason": self.reason,
            "outgroup_found": sorted(self.outgroup_found),
            "outgroup_missing": sorted(self.outgroup_missing),
            "monophyletic": self.monophyletic,
        }


def _normalizar(nomes: Iterable[str]) -> FrozenSet[str]:
    return frozenset(strip_accession_version(n) for n in nomes if n)


def outgroup_from_classifier(tree: Tree,
                             classifier: Callable[[str], str],
                             ingroup_codes: Sequence[str]) -> FrozenSet[str]:
    """
    Deriva o grupo externo a partir de um classificador de espécie.

    Serve ao caso dos experimentos de *Variola*, em que o grupo externo é
    "tudo que não é VARV" — camelpox, cowpox, taterapox. O que se **declara** é
    o grupo **interno**; o externo é o complemento, e assim a declaração
    continua sendo do pesquisador.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        Árvore de onde sair os terminais.
    classifier : callable
        Mapeia nome de terminal em sigla de espécie — ver
        `workflow.stability.case_study.build_classifier`.
    ingroup_codes : sequence of str
        Siglas que compõem o grupo **interno**.

    Return
    ------
    frozenset of str
        Nomes normalizados dos terminais que **não** pertencem ao grupo interno.
    """
    internos = set(ingroup_codes)
    return frozenset(
        strip_accession_version(t.name)
        for t in tree.get_terminals()
        if t.name and classifier(t.name) not in internos
    )


def root_at_outgroup(tree: Tree,
                     outgroup: Iterable[str],
                     require_monophyletic: bool = True) -> Tuple[Optional[Tree], RootingReport]:
    """
    Enraíza uma árvore no grupo externo declarado.

    A árvore original **não é modificada**: devolve-se uma cópia enraizada, para
    que a mesma árvore possa ser analisada nas duas formas sem que uma
    contamine a outra.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        Árvore a enraizar.
    outgroup : iterable of str
        Nomes dos terminais do grupo externo. Normalizados internamente, de modo
        que a grafia truncada de IQ-TREE e RAxML case com a íntegra (D13).
    require_monophyletic : bool, optional
        Se True (padrão), recusa enraizar quando o grupo externo não é um clado.
        Desligar isto significa escolher arbitrariamente uma entre várias arestas
        candidatas — só faz sentido para exploração, nunca para um resultado.

    Return
    ------
    tuple
        ``(árvore enraizada ou None, RootingReport)``.
    """
    declarado = _normalizar(outgroup)
    presentes = _normalizar(t.name for t in tree.get_terminals())

    if not declarado:
        return None, RootingReport(
            rooted=False,
            reason="grupo externo vazio: o enraizamento tem de ser declarado, não inferido")

    encontrados = declarado & presentes
    ausentes = declarado - presentes

    if not encontrados:
        return None, RootingReport(
            rooted=False,
            reason="nenhum táxon do grupo externo está nesta árvore",
            outgroup_missing=ausentes)

    if encontrados == presentes:
        return None, RootingReport(
            rooted=False,
            reason="o grupo externo é a árvore inteira; não sobra grupo interno",
            outgroup_found=encontrados, outgroup_missing=ausentes)

    copia = copy.deepcopy(tree)
    terminais = {strip_accession_version(t.name): t
                 for t in copia.get_terminals() if t.name}
    alvos = [terminais[n] for n in sorted(encontrados)]

    monofiletico: Optional[bool] = None
    if len(alvos) > 1:
        monofiletico = bool(copia.is_monophyletic(alvos))
        if require_monophyletic and not monofiletico:
            return None, RootingReport(
                rooted=False,
                reason=("grupo externo não é monofilético nesta árvore: há mais de uma "
                        "aresta candidata a raiz e a escolha seria arbitrária"),
                outgroup_found=encontrados, outgroup_missing=ausentes,
                monophyletic=False)

    copia.root_with_outgroup(*alvos)
    copia.rooted = True

    return copia, RootingReport(
        rooted=True,
        outgroup_found=encontrados,
        outgroup_missing=ausentes,
        monophyletic=monofiletico)


def root_tree_set(trees: Dict[str, Tree],
                  outgroup: Iterable[str],
                  require_monophyletic: bool = True
                  ) -> Tuple[Dict[str, Tree], Dict[str, RootingReport]]:
    """
    Enraíza todas as árvores de um conjunto no **mesmo** grupo externo.

    Enraizamento comum é o que torna a comparação por clado legítima: árvores
    enraizadas em pontos diferentes não são comparáveis por clado, e foi isso
    que D3 descreveu.

    Devolve as duas coisas separadamente — as árvores que enraizaram e o
    relatório de **todas** —, porque o que decide se a análise enraizada é
    publicável não é quantas enraizaram, e sim se **alguma** ficou de fora.

    Parameters
    ----------
    trees : dict of str -> Tree
        Árvores indexadas pelo rótulo do pipeline.
    outgroup : iterable of str
        Grupo externo declarado, comum a todas.
    require_monophyletic : bool, optional
        Repassado a `root_at_outgroup`.

    Return
    ------
    tuple
        ``(árvores enraizadas, relatórios por pipeline)``.
    """
    enraizadas: Dict[str, Tree] = {}
    relatorios: Dict[str, RootingReport] = {}

    for nome, arvore in trees.items():
        enraizada, relatorio = root_at_outgroup(arvore, outgroup, require_monophyletic)
        relatorios[nome] = relatorio
        if enraizada is not None:
            enraizadas[nome] = enraizada

    return enraizadas, relatorios
