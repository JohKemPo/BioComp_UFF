"""
Identidade de clados para mineração de padrões filogenéticos.

Um clado só pode ser tratado como *item* em mineração de padrões frequentes se
duas árvores que recuperam a mesma partição de táxons produzirem exatamente o
mesmo identificador. A identidade utilizada por `workflow.utils.treeUtils`
(`encode_list_to_int`) não satisfaz essa condição por dois motivos:

1. **Dependência de ordem** -- o hash é calculado sobre a *representação textual
   da lista* de terminais na ordem de travessia. Duas árvores com a mesma
   topologia, porém com filhos em ordem distinta, geram identificadores
   diferentes para o mesmo clado.
2. **Espaço de 16 bits** -- apenas os 4 primeiros dígitos hexadecimais do MD5
   são usados (65 536 valores). Com alguns milhares de clados, colisões deixam
   de ser raras (aniversário: ~2 colisões esperadas para 1 000 clados).

Ambos os efeitos atuam na mesma direção perigosa: o primeiro *subestima* o
suporte de clados verdadeiros (fragmenta um clado em vários itens) e o segundo
*fabrica* suporte (funde clados distintos em um item).

Este módulo define a identidade canônica -- `frozenset` de nomes de terminais,
opcionalmente resumido em um digest de 128 bits -- e mantém a implementação
legada para fins de auditoria e reprodutibilidade dos resultados anteriores.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from Bio.Phylo.BaseTree import Clade, Tree

__all__ = [
    "CladeIdentity",
    "canonical_clade_id",
    "canonical_digest",
    "canonical_item_id",
    "canonical_bipartition",
    "strip_accession_version",
    "ITEM_ID_BITS",
    "legacy_terminal_hash",
    "legacy_clade_hash",
    "clade_identities",
    "audit_legacy_identity",
    "LegacyAudit",
]


#: Largura do identificador inteiro de item (`canonical_item_id`).
#:
#: 52 bits, e o limite é o consumidor mais estreito da cadeia: o valor viaja no
#: JSON da API até o navegador, e `Number` do JavaScript só representa inteiros
#: exatamente até 2^53 - 1. Um identificador maior seria **arredondado no
#: cliente** — trocar a colisão de 16 bits por um arredondamento silencioso não
#: seria correção. Cabe também no inteiro de 64 bits com sinal do Neo4j.
#:
#: A margem continua enorme: 4,5 * 10^15 valores contra os 65 536 do esquema
#: legado. Com 10^5 clados distintos — três ordens de grandeza acima do maior
#: experimento atual, VARV-121 com 270 — a probabilidade de qualquer colisão
#: fica na casa de 10^-6.
ITEM_ID_BITS = 52

_VERSION_SUFFIX = re.compile(r"\.\d*$")


def strip_accession_version(name: str) -> str:
    """
    Normaliza o rótulo de um terminal removendo o sufixo de versão do acesso.

    IQ-TREE e RAxML-NG reescrevem rótulos ao emitir Newick/Nexus e, para acessos
    RefSeq, truncam o dígito de versão (``NC_008030.1`` torna-se ``NC_008030.``).
    Como o dígito é o único caractere perdido, remover o sufixo inteiro reconcilia
    as duas grafias sem ambiguidade -- dois acessos do mesmo registro em versões
    diferentes não coexistem em um alinhamento.

    Vive neste módulo, e não em `stability`, porque normalizar o rótulo **é** parte
    de decidir a identidade: dois clados só são o mesmo item se seus terminais
    forem os mesmos, e isso depende de como o nome é lido (D5, D13).

    Parameters
    ----------
    name : str
        Rótulo do terminal como lido da árvore.

    Return
    ------
    str
        Rótulo sem o sufixo de versão.
    """
    return _VERSION_SUFFIX.sub("", name.strip().strip("'\""))


# --------------------------------------------------------------------------- #
# Identidade canônica
# --------------------------------------------------------------------------- #

def canonical_clade_id(clade: Clade) -> FrozenSet[str]:
    """
    Retorna a identidade canônica de um clado.

    A identidade é o conjunto (não ordenado) dos nomes de seus terminais, o que
    a torna invariante à ordem de travessia, à rotação de nós e ao enraizamento
    arbitrário de subárvores.

    Parameters
    ----------
    clade : Bio.Phylo.BaseTree.Clade
        Clado do qual se deseja a identidade.

    Return
    ------
    frozenset of str
        Conjunto dos nomes dos terminais descendentes do clado.
    """
    return frozenset(terminal.name for terminal in clade.get_terminals())


def canonical_digest(clade_id: Iterable[str]) -> str:
    """
    Resume uma identidade canônica em um digest hexadecimal de 128 bits.

    Os nomes são ordenados antes do hashing, de modo que o digest herda a
    invariância à ordem. O espaço de 128 bits torna colisões desprezíveis na
    escala de qualquer análise filogenômica praticável.

    Parameters
    ----------
    clade_id : iterable of str
        Nomes dos terminais que compõem o clado.

    Return
    ------
    str
        Digest MD5 completo (32 dígitos hexadecimais).
    """
    payload = "\n".join(sorted(clade_id))
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def canonical_bipartition(taxa: Iterable[str],
                          all_taxa: FrozenSet[str]) -> Optional[FrozenSet[str]]:
    """
    Reduz um clado à bipartição não enraizada canônica (D3).

    Remover uma aresta interna de uma árvore não enraizada parte ``X`` em
    ``A | X∖A``; a bipartição é o par **não ordenado** ``{A, X∖A}``. FastTree,
    IQ-TREE, RAxML e NJ escrevem topologia não enraizada com raiz trifurcante:
    o clado que se lê do Newick depende de onde o arquivo pôs a raiz, que é
    convenção de escrita e não hipótese biológica. Comparar clados entre essas
    árvores mede a convenção; comparar bipartições mede a topologia.

    Forma canônica adotada (`03-metricas §2.2`): o lado de menor cardinalidade;
    em empate, o menor na ordem lexicográfica dos nomes ordenados. Assim, duas
    escritas da mesma topologia produzem o mesmo objeto.

    Parameters
    ----------
    taxa : iterable of str
        Terminais descendentes do clado.
    all_taxa : frozenset of str
        Conjunto completo de terminais da árvore.

    Return
    ------
    frozenset of str or None
        Lado menor da bipartição, ou ``None`` se ela for trivial — isto é, se
        algum lado tiver menos de 2 táxons, caso em que a aresta é externa e
        não carrega informação topológica.
    """
    side = frozenset(taxa)
    other = frozenset(all_taxa) - side
    if len(side) < 2 or len(other) < 2:
        return None
    return min((side, other), key=lambda s: (len(s), sorted(s)))


def canonical_item_id(clade_id: Iterable[str], bits: int = ITEM_ID_BITS) -> int:
    """
    Identidade canônica de um clado como inteiro, para uso como *item* de mineração.

    É o `canonical_digest` truncado em `bits`. Serve ao pipeline de produção, onde
    o identificador precisa ser um inteiro: é chave em `metadata.json`, item do
    FPMax e propriedade de nó no Neo4j, que só guarda inteiro de 64 bits com sinal.

    Herda do digest as duas propriedades que o esquema legado não tinha: é
    invariante à ordem de travessia e tem espaço suficiente para que colisões
    sejam desprezíveis (ver `ITEM_ID_BITS`).

    Os nomes devem chegar já normalizados por `strip_accession_version` — é o que
    faz o mesmo clado ter a mesma identidade numa árvore de IQ-TREE (rótulo
    truncado) e numa de FastTree (rótulo íntegro).

    Parameters
    ----------
    clade_id : iterable of str
        Nomes dos terminais que compõem o clado.
    bits : int, optional
        Largura do identificador. O padrão, `ITEM_ID_BITS`, cabe em 64 bits com sinal.

    Return
    ------
    int
        Identificador inteiro, invariante à ordem dos nomes.
    """
    if bits <= 0 or bits > 128:
        raise ValueError("bits deve estar entre 1 e 128")
    # `frozenset` antes do digest: a identidade é o CONJUNTO de terminais, como em
    # `CladeIdentity.from_clade`. Sem isso, um nome repetido — possível quando dois
    # rótulos truncados normalizam para o mesmo acesso (D13) — mudaria o digest.
    digest = canonical_digest(frozenset(clade_id))
    return int(digest, 16) >> (128 - bits)


@dataclass(frozen=True)
class CladeIdentity:
    """
    Identidade completa de um clado observado em uma árvore.

    Attributes
    ----------
    taxa : frozenset of str
        Identidade canônica (conjunto de terminais).
    size : int
        Número de terminais.
    digest : str
        Digest canônico de 128 bits.
    legacy : int
        Identificador legado de 16 bits, dependente da ordem de travessia.
    """

    taxa: FrozenSet[str]
    size: int
    digest: str
    legacy: int

    @classmethod
    def from_clade(cls, clade: Clade) -> "CladeIdentity":
        """
        Constrói a identidade a partir de um clado do Biopython.

        Parameters
        ----------
        clade : Bio.Phylo.BaseTree.Clade
            Clado a ser identificado.

        Return
        ------
        CladeIdentity
            Identidade canônica e legada do clado.
        """
        ordered = [terminal.name for terminal in clade.get_terminals()]
        taxa = frozenset(ordered)
        return cls(
            taxa=taxa,
            size=len(taxa),
            digest=canonical_digest(taxa),
            legacy=legacy_clade_hash([legacy_terminal_hash(name) for name in ordered]),
        )


# --------------------------------------------------------------------------- #
# Identidade legada (reimplementação fiel, para auditoria)
# --------------------------------------------------------------------------- #

def legacy_terminal_hash(name: str) -> int:
    """
    Reproduz `workflow.utils.treeUtils.calculate_tree_hash` para um terminal.

    Parameters
    ----------
    name : str
        Nome do terminal.

    Return
    ------
    int
        Inteiro de 16 bits derivado dos 4 primeiros dígitos do MD5 do nome.
    """
    return int(hashlib.md5(name.encode()).hexdigest()[:4], 16)


def legacy_clade_hash(terminal_hashes: Sequence[int]) -> int:
    """
    Reproduz `workflow.utils.treeUtils.encode_list_to_int`.

    O argumento é consumido na ordem em que é fornecido -- é justamente essa
    dependência de ordem que a identidade canônica elimina.

    Parameters
    ----------
    terminal_hashes : sequence of int
        Hashes dos terminais, na ordem de travessia da árvore.

    Return
    ------
    int
        Inteiro de 16 bits identificando o clado sob o esquema legado.
    """
    return int(hashlib.md5(str(list(terminal_hashes)).encode()).hexdigest()[:4], 16)


# --------------------------------------------------------------------------- #
# Extração e auditoria
# --------------------------------------------------------------------------- #

def clade_identities(tree: Tree, include_root: bool = False) -> List[CladeIdentity]:
    """
    Extrai a identidade de todos os clados internos informativos de uma árvore.

    Clados triviais são descartados: folhas (um terminal) e, por padrão, o clado
    que contém todos os táxons, que não carrega informação topológica.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        Árvore a ser percorrida.
    include_root : bool, optional
        Se True, mantém o clado universal (todos os terminais).

    Return
    ------
    list of CladeIdentity
        Identidades dos clados internos, sem repetição.
    """
    total = len(tree.get_terminals())
    seen: Set[FrozenSet[str]] = set()
    identities: List[CladeIdentity] = []

    for clade in tree.find_clades():
        identity = CladeIdentity.from_clade(clade)
        if identity.size <= 1:
            continue
        if not include_root and identity.size == total:
            continue
        if identity.taxa in seen:
            continue
        seen.add(identity.taxa)
        identities.append(identity)

    return identities


@dataclass
class LegacyAudit:
    """
    Resultado da comparação entre a identidade canônica e a legada.

    Attributes
    ----------
    n_canonical : int
        Número de clados distintos sob a identidade canônica.
    n_legacy : int
        Número de itens distintos sob a identidade legada.
    fragmented : dict
        Clados canônicos mapeados para mais de um identificador legado
        (subestimação de suporte por dependência de ordem).
    collided : dict
        Identificadores legados compartilhados por clados canônicos distintos
        (superestimação de suporte por colisão de 16 bits).
    """

    n_canonical: int = 0
    n_legacy: int = 0
    fragmented: Dict[FrozenSet[str], Set[int]] = field(default_factory=dict)
    collided: Dict[int, Set[FrozenSet[str]]] = field(default_factory=dict)

    @property
    def fragmentation_rate(self) -> float:
        """Fração dos clados canônicos fragmentados em múltiplos itens legados."""
        return len(self.fragmented) / self.n_canonical if self.n_canonical else 0.0

    @property
    def collision_rate(self) -> float:
        """Fração dos itens legados que fundem clados canônicos distintos."""
        return len(self.collided) / self.n_legacy if self.n_legacy else 0.0

    def summary(self) -> Dict[str, float]:
        """Resume a auditoria em um dicionário achatado, pronto para tabulação."""
        return {
            "canonical_clades": self.n_canonical,
            "legacy_items": self.n_legacy,
            "fragmented_clades": len(self.fragmented),
            "fragmentation_rate": round(self.fragmentation_rate, 4),
            "colliding_items": len(self.collided),
            "collision_rate": round(self.collision_rate, 4),
        }


def audit_legacy_identity(trees: Dict[str, Tree]) -> LegacyAudit:
    """
    Quantifica a divergência entre a identidade canônica e a legada.

    Parameters
    ----------
    trees : dict of str -> Bio.Phylo.BaseTree.Tree
        Conjunto de árvores construídas sobre o mesmo conjunto de táxons.

    Return
    ------
    LegacyAudit
        Contagens de fragmentação (ordem) e de colisão (16 bits).
    """
    canonical_to_legacy: Dict[FrozenSet[str], Set[int]] = {}
    legacy_to_canonical: Dict[int, Set[FrozenSet[str]]] = {}

    for tree in trees.values():
        for identity in clade_identities(tree):
            canonical_to_legacy.setdefault(identity.taxa, set()).add(identity.legacy)
            legacy_to_canonical.setdefault(identity.legacy, set()).add(identity.taxa)

    return LegacyAudit(
        n_canonical=len(canonical_to_legacy),
        n_legacy=len(legacy_to_canonical),
        fragmented={k: v for k, v in canonical_to_legacy.items() if len(v) > 1},
        collided={k: v for k, v in legacy_to_canonical.items() if len(v) > 1},
    )
