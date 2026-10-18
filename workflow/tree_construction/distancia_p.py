"""
D28 — distância e parcimônia sem tratar lacuna como caráter.

Achado (`docs/science/02-defeitos-que-alteram-resultado.md#d28`, medido em
2026-09-23 durante M7.5): `Bio.Phylo.TreeConstruction.DistanceCalculator`,
quando `scoring_matrix` é `None` (o caso do modelo `'identity'` que
`builder.distance_matrix` usa), calcula o denominador como
``max_score = len(seq1)`` **incondicionalmente** — mesmo quando
``skip_letters`` está preenchido. `skip_letters` só afeta o numerador
(quantas posições contam como identidade), nunca o denominador. Confirmado
lendo o fonte instalado (`Bio.Phylo.TreeConstruction.DistanceCalculator._pairwise`),
não só documentação. O resultado: pares com lacuna medem "quanto do
alinhamento inteiro é igual", não "quanto das posições comparáveis é igual" —
e a `ParsimonyScorer` (Fitch) do mesmo módulo trata `'-'` como um estado de
caráter comum, a mesma classe de erro.

Este módulo dá duas coisas ao `builder.py`, sem tocar o Biopython instalado:

1. `matriz_distancia_p` — p-distância nos sítios em que **ambas** as
   sequências têm base determinada, denominador ajustado ao par (não ao
   alinhamento inteiro). Lacuna e caracteres totalmente ambíguos (`n`/`N`)
   não contam nem no numerador nem no denominador.
2. `FitchAmbiguidadeIupac` — mesmo algoritmo de Fitch da `ParsimonyScorer`,
   mas o estado de uma folha para um código IUPAC (`-`, `n`, `r`, `y`, ...)
   é o conjunto de bases que o código representa, não o caractere literal.
   Lacuna e `n` viram "qualquer base" (interseção com o estado do irmão nunca
   fica vazia por causa da lacuna, current: 0 passos de parcimônia cobrados);
   códigos parciais (`r` = A ou G) intersectam só com as bases que
   representam. Colunas onde os únicos caracteres distintos são lacuna/`n`
   continuam contribuindo 0 (a checagem "coluna não-informativa" da
   `ParsimonyScorer`, feita por igualdade literal de caractere, já cobre esse
   caso: uma coluna só de `-`/`n` é igual ao seu primeiro caractere em toda
   posição).

Ambas seguem a regra 5 do projeto: um par sem nenhum sítio comparável (as
duas sequências não têm nenhuma posição em que ambas tenham base) não vira
`1.0` como se fosse a distância máxima medida — levanta `SemSitioComparavel`,
porque "sem dado" não é o mesmo número que "totalmente diferente".

Isso acontece com dado real, não só em teste: um conjunto que mistura
fragmentos de genes diferentes (em `Zika479_Test_large`, um fragmento de NS3
contra quatro de NS5) não tem sítio homólogo entre esses pares. A p-distância
é indefinida ali, e NJ/UPGMA sobre o alinhamento inteiro também (DEC-152).
"""

from __future__ import annotations

from typing import Dict, FrozenSet, Iterable, List, Sequence, Tuple

import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, ParsimonyScorer

#: Tabela IUPAC de ambiguidade de nucleotídeo, minúscula (as sequências do
#: pipeline chegam em minúsculo — conferido em dados reais de
#: `Zika_21seq_validacao`). `-` (lacuna) e `n` são tratados como "qualquer
#: base": a ambiguidade mais ampla possível, e a mesma semântica de "dado
#: ausente" que os métodos de ML (IQ-TREE/RAxML-NG/FastTree) já usam.
BASES = frozenset("acgt")
IUPAC_AMBIGUIDADE: Dict[str, FrozenSet[str]] = {
    "a": frozenset("a"),
    "c": frozenset("c"),
    "g": frozenset("g"),
    "t": frozenset("t"),
    "u": frozenset("t"),  # RNA lido como DNA, mesma convenção do resto do pipeline
    "r": frozenset("ag"),
    "y": frozenset("ct"),
    "s": frozenset("gc"),
    "w": frozenset("at"),
    "k": frozenset("gt"),
    "m": frozenset("ac"),
    "b": frozenset("cgt"),
    "d": frozenset("agt"),
    "h": frozenset("act"),
    "v": frozenset("acg"),
    "n": BASES,
    "-": BASES,
}

#: Letras que não carregam informação de identidade — lacuna e ambiguidade
#: total. Usadas só pela distância p (a parcimônia usa a tabela inteira acima).
LETRAS_SEM_INFORMACAO = frozenset({"-", "n"})


#: Quantos pares a mensagem da exceção nomeia. O motivo vai para o manifesto,
#: que corta em 600 caracteres (`tool_runs._MOTIVO_MAX`); a lista inteira fica
#: em `SemSitioComparavel.pares`.
_PARES_NA_MENSAGEM = 5


class SemSitioComparavel(ValueError):
    """
    Nenhum sítio em que as duas sequências de um par têm base determinada.

    Attributes
    ----------
    pares : list of tuple(str, str)
        Todos os pares sem sítio comparável. Vazio quando a exceção vem de
        `p_distancia_par`, que só conhece as duas sequências, não os nomes.
    """

    def __init__(self, mensagem: str, pares: Sequence[Tuple[str, str]] = ()) -> None:
        super().__init__(mensagem)
        self.pares = list(pares)

    @classmethod
    def para_pares(cls, pares: Sequence[Tuple[str, str]],
                   n_sequencias: int) -> "SemSitioComparavel":
        nomeados = "; ".join(f"{a} × {b}" for a, b in pares[:_PARES_NA_MENSAGEM])
        resto = len(pares) - _PARES_NA_MENSAGEM
        if resto > 0:
            nomeados += f"; e mais {resto}"
        envolvidos = len({t for par in pares for t in par})
        return cls(
            f"{len(pares)} par(es) sem nenhum sítio em que as duas sequências têm "
            f"base determinada ({envolvidos} de {n_sequencias} sequências envolvidas) "
            f"— p-distância indefinida, NJ/UPGMA não se calculam sobre este "
            f"alinhamento (D28): {nomeados}. Com on_method_failure='continue' o "
            f"método fica tentado_e_falhou e os demais seguem.",
            pares,
        )


def _normalizar(seq) -> str:
    return str(seq).lower()


def pares_sem_sitio_comparavel(alignment) -> List[Tuple[str, str]]:
    """
    Pares do alinhamento sem nenhum sítio em que ambas têm base determinada.

    Mesmo critério de `p_distancia_par` (fora de `LETRAS_SEM_INFORMACAO`),
    calculado de uma vez para todos os pares: sobreposição = M·Mᵀ, onde
    M[i, c] = 1 se a sequência i tem base na coluna c.

    Return
    ------
    list of tuple(str, str)
        `(id_i, id_j)` com `j < i`, na ordem em que `matriz_distancia_p`
        visita os pares.
    """
    nomes = [r.id for r in alignment]
    if len(nomes) < 2:
        return []
    sem_info = np.frombuffer("".join(sorted(LETRAS_SEM_INFORMACAO)).encode("ascii"),
                             dtype=np.uint8)
    # float32 conta exato até 2^24 colunas — muito acima de qualquer genoma
    # viral alinhado (o maior aqui, VARV, tem ~214 mil colunas).
    tem_base = np.stack([
        ~np.isin(np.frombuffer(_normalizar(r.seq).encode("ascii", "replace"),
                               dtype=np.uint8), sem_info)
        for r in alignment
    ]).astype(np.float32)
    sobreposicao = tem_base @ tem_base.T
    return [(nomes[i], nomes[j])
            for i in range(len(nomes)) for j in range(i)
            if sobreposicao[i, j] == 0]


def p_distancia_par(seq1, seq2) -> Tuple[float, int]:
    """
    p-distância entre duas sequências alinhadas, nos sítios em que ambas
    têm base determinada (fora de `LETRAS_SEM_INFORMACAO`).

    Return
    ------
    (distancia, n_sitios_comparaveis) : Tuple[float, int]
        `distancia` é a fração de sítios comparáveis que diferem. Levanta
        `SemSitioComparavel` se `n_sitios_comparaveis == 0` — a regra 5 do
        projeto proíbe devolver um número (`0.0` ou `1.0`) onde a métrica é
        indefinida.
    """
    s1, s2 = _normalizar(seq1), _normalizar(seq2)
    if len(s1) != len(s2):
        raise ValueError(f"Sequências de comprimento diferente: {len(s1)} != {len(s2)}")
    comparaveis = 0
    diferentes = 0
    for c1, c2 in zip(s1, s2):
        if c1 in LETRAS_SEM_INFORMACAO or c2 in LETRAS_SEM_INFORMACAO:
            continue
        comparaveis += 1
        if c1 != c2:
            diferentes += 1
    if comparaveis == 0:
        raise SemSitioComparavel(
            "nenhum sítio em que as duas sequências têm base determinada — "
            "p-distância indefinida para este par (D28)"
        )
    return diferentes / comparaveis, comparaveis


def matriz_distancia_p(alignment) -> DistanceMatrix:
    """
    Matriz de distâncias compatível com `Bio.Phylo.TreeConstruction`
    (mesma interface que `DistanceCalculator.get_distance` devolve), usando
    `p_distancia_par` em vez do `'identity'` do Biopython — corrige o
    denominador fixo em `len(seq1)` (D28).

    Parameters
    ----------
    alignment : Bio.Align.MultipleSeqAlignment

    Raises
    ------
    SemSitioComparavel
        Se algum par do alinhamento não tiver nenhum sítio comparável —
        não silenciado como distância máxima. Levantada antes de calcular
        qualquer distância, com **todos** os pares em `.pares`.
    """
    sem_sitio = pares_sem_sitio_comparavel(alignment)
    if sem_sitio:
        raise SemSitioComparavel.para_pares(sem_sitio, len(alignment))
    nomes = [r.id for r in alignment]
    dm = DistanceMatrix(nomes)
    for i in range(len(alignment)):
        for j in range(i):
            dist, _n = p_distancia_par(alignment[i].seq, alignment[j].seq)
            dm[nomes[i], nomes[j]] = dist
    return dm


class FitchAmbiguidadeIupac(ParsimonyScorer):
    """
    `ParsimonyScorer` (Fitch) tratando lacuna e códigos IUPAC de ambiguidade
    como o conjunto de bases que representam, não como um estado de caráter
    literal — corrige a leitura de `'-'` como um 5º estado (D28).

    Reimplementa `get_score` a partir do fonte do Biopython instalado
    (`Bio.Phylo.TreeConstruction.ParsimonyScorer.get_score`), mudando só a
    linha de inicialização dos estados de folha. O restante do algoritmo
    (interseção/união bottom-up, contagem de passos) é idêntico — conferido
    linha a linha contra o fonte em 2026-09-23.
    """

    def get_score(self, tree, alignment):
        if not tree.is_bifurcating():
            raise ValueError("The tree provided should be bifurcating.")
        if not tree.rooted:
            tree.root_at_midpoint()
        terms = tree.get_terminals()
        terms.sort(key=lambda term: term.name)
        alignment.sort()
        if not all(t.name == a.id for t, a in zip(terms, alignment)):
            raise ValueError(
                "Taxon names of the input tree should be the same with the alignment."
            )
        if self.matrix:
            # Sankoff (matriz de penalidade) não é usado pelo pipeline hoje
            # (parsimony_constructor sempre chama ParsimonyScorer() sem
            # matriz) — fora do escopo de D28, delega ao Biopython original.
            return super().get_score(tree, alignment)

        score = 0
        for i in range(len(alignment[0])):
            column_i = _normalizar(alignment[:, i])
            if column_i == len(column_i) * column_i[0]:
                continue  # coluna não-informativa (todas as letras iguais)
            clade_states = dict(
                zip(terms, [IUPAC_AMBIGUIDADE.get(c, frozenset(c)) for c in column_i])
            )
            score_i = 0
            for clade in tree.get_nonterminals(order="postorder"):
                esquerdo, direito = clade.clades
                estado_esquerdo = clade_states[esquerdo]
                estado_direito = clade_states[direito]
                estado = estado_esquerdo & estado_direito
                if not estado:
                    estado = estado_esquerdo | estado_direito
                    score_i += 1
                clade_states[clade] = estado
            score += score_i
        return score
