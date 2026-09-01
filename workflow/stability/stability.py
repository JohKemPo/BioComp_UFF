"""
Estabilidade de clados entre pipelines filogenéticos.

A ideia central é tratar cada *pipeline* (combinação alinhador x método de
inferência) como uma transação e cada clado como um item. O suporte de um clado
passa a ser a fração de pipelines que o recuperam -- uma medida de robustez
**metodológica**, complementar (e ortogonal) ao bootstrap, que mede robustez
**amostral** dentro de um único pipeline.

Como o número de pipelines M é pequeno (tipicamente 8 a 10), o reticulado de
conjuntos de clados fechados pode ser enumerado exatamente em O(2**M * |C|),
dispensando heurísticas de mineração. `StabilityAnalyzer.maximal_patterns`
devolve, portanto, os conjuntos maximais exatos -- não uma aproximação.
"""

from __future__ import annotations

import itertools
import os
import re
from dataclasses import dataclass
from typing import Callable, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from workflow.alignment.aligners import ALIGNERS as _ALINHADORES_REGISTRADOS
from workflow.stability.clade_identity import (
    CladeIdentity,
    LegacyAudit,
    audit_legacy_identity,
    canonical_bipartition,
    clade_identities,
    strip_accession_version,
)

__all__ = ["PipelineLabel", "TreeSet", "StabilityAnalyzer", "CladeRecord", "Pattern"]

#: Métodos de inferência reconhecidos nos nomes de arquivo produzidos pelo workflow.
INFERENCE_METHODS = ("iqtree", "raxml", "fasttree", "mrbayes",
                     "nj_distance", "upgma_distance",
                     "nj_parsimony", "upgma_parsimony", "parsimony")

#: Alinhadores reconhecidos nos nomes de arquivo produzidos pelo workflow —
#: a mesma fonte que o controlador usa para gerar esses nomes
#: (`workflow.alignment.aligners.ALIGNERS`), não uma cópia paralela. Cópia
#: paralela foi a causa do bug encontrado em 2026-09-01: quando DEC-050
#: acrescentou "mafft_iterative" ao registro real, esta tupla — que só
#: existia aqui, desatualizada — ficou parada em `("mafft", "clustalo")`, e
#: todo pipeline do braço iterativo colidia com o do braço progressivo.
ALIGNERS = tuple(_ALINHADORES_REGISTRADOS.keys())


#: Sufixo de versão de acesso GenBank, possivelmente truncado pelo inferidor.
# `strip_accession_version` vive em `clade_identity`: normalização e identidade
# são a mesma decisão e não podem divergir (D5). Reexportado aqui porque é a
# porta de entrada histórica do módulo.


@dataclass(frozen=True)
class PipelineLabel:
    """
    Decomposição do nome de uma árvore nos fatores do delineamento experimental.

    Attributes
    ----------
    name : str
        Rótulo completo do pipeline (ex.: ``mafft_iqtree``).
    aligner : str
        Alinhador múltiplo utilizado.
    inference : str
        Método de construção da árvore.
    """

    name: str
    aligner: str
    inference: str

    @classmethod
    def parse(cls, filename: str, prefix: str = "tree_") -> "PipelineLabel":
        """
        Extrai alinhador e método de inferência do nome de arquivo da árvore.

        Parameters
        ----------
        filename : str
            Nome do arquivo, com ou sem extensão (ex.:
            ``tree_dataset_final_mafft_iqtree.nexus``).
        prefix : str, optional
            Prefixo a ser ignorado no início do nome.

        Return
        ------
        PipelineLabel
            Rótulo decomposto; campos não reconhecidos recebem ``"unknown"``.
        """
        stem = os.path.splitext(os.path.basename(filename))[0]
        if stem.startswith(prefix):
            stem = stem[len(prefix):]

        # Delimitado por "_" dos dois lados (a própria posição no stem é livre
        # — o antigo `tokens = set(stem.split("_"))` também não exigia
        # posição), e o MAIS LONGO vence, não o primeiro: "mafft" é substring
        # de "mafft_iterative", e seria sempre encontrado primeiro. O split
        # em tokens soltos nunca reconhecia o composto, porque o próprio
        # split já separava os dois em tokens distintos antes da comparação.
        stem_delimitado = f"_{stem}_"
        candidatos_alinhador = [
            a for a in ALIGNERS if f"_{a}_" in stem_delimitado
        ]
        aligner = max(candidatos_alinhador, key=len) if candidatos_alinhador else "unknown"
        # Sufixo MAIS LONGO, não o primeiro que casar: `clustalo_nj_parsimony`
        # termina tanto em `parsimony` quanto em `nj_parsimony`, e escolher o
        # curto fundia NJ com UPGMA num único pipeline `clustalo_parsimony` —
        # com uma árvore sobrescrevendo a outra em silêncio (D19). Ordenar por
        # tamanho torna a regra independente da ordem da tupla.
        candidatos = [m for m in INFERENCE_METHODS if stem.endswith(m)]
        inference = max(candidatos, key=len) if candidatos else "unknown"
        name = f"{aligner}_{inference}"
        return cls(name=name, aligner=aligner, inference=inference)


@dataclass(frozen=True)
class CladeRecord:
    """
    Um clado observado, com o conjunto de pipelines que o recuperam.

    Attributes
    ----------
    taxa : frozenset of str
        Identidade canônica do clado.
    size : int
        Número de terminais.
    pipelines : frozenset of str
        Pipelines que recuperaram o clado.
    support : float
        Fração de pipelines que recuperaram o clado.
    """

    taxa: FrozenSet[str]
    size: int
    pipelines: FrozenSet[str]
    support: float


@dataclass(frozen=True)
class Pattern:
    """
    Conjunto maximal de clados co-recuperados por um mesmo grupo de pipelines.

    Attributes
    ----------
    clades : Tuple[FrozenSet[str], ...]
        Clados que compõem o padrão.
    pipelines : frozenset of str
        Pipelines que contêm todos os clados do padrão.
    support : float
        Fração de pipelines de suporte.
    """

    clades: Tuple[FrozenSet[str], ...]
    pipelines: FrozenSet[str]
    support: float

    @property
    def size(self) -> int:
        """Número de clados no padrão."""
        return len(self.clades)


class TreeSet:
    """
    Conjunto de árvores inferidas a partir do mesmo alinhamento-fonte.

    Attributes
    ----------
    trees : dict of str -> Bio.Phylo.BaseTree.Tree
        Árvores indexadas pelo rótulo do pipeline.
    labels : dict of str -> PipelineLabel
        Decomposição de cada rótulo nos fatores do delineamento.
    taxa : list of str
        Nomes dos terminais, na ordem da primeira árvore.
    mangled_labels : dict of str -> list of str
        Por pipeline, os táxons cuja grafia bruta divergia da majoritária --
        isto é, os rótulos adulterados pelo programa de inferência.
    """

    def __init__(
        self,
        trees: Dict[str, Tree],
        labels: Dict[str, PipelineLabel],
        normalizer: Optional[Callable[[str], str]] = strip_accession_version,
    ) -> None:
        if not trees:
            raise ValueError("TreeSet exige ao menos uma árvore.")

        self.trees = trees
        self.labels = labels
        self.mangled_labels: Dict[str, List[str]] = {name: [] for name in trees}

        if normalizer is not None:
            raw: Dict[str, Dict[str, str]] = {}
            for name, tree in trees.items():
                raw[name] = {}
                for terminal in tree.get_terminals():
                    normalized = normalizer(terminal.name)
                    raw[name][normalized] = terminal.name
                    terminal.name = normalized

            # Um rótulo é considerado adulterado quando sua grafia bruta diverge
            # da grafia majoritária entre os pipelines para o mesmo táxon. Em caso
            # de empate vence a grafia mais longa: adulteração aqui é truncamento,
            # portanto a forma íntegra nunca é a mais curta.
            modal: Dict[str, str] = {}
            for key in {k for mapping in raw.values() for k in mapping}:
                spellings = [mapping[key] for mapping in raw.values() if key in mapping]
                modal[key] = max(set(spellings), key=lambda s: (spellings.count(s), len(s)))
            for name, mapping in raw.items():
                self.mangled_labels[name] = sorted(
                    key for key, spelling in mapping.items() if spelling != modal[key]
                )

        first = next(iter(trees.values()))
        self.taxa = [terminal.name for terminal in first.get_terminals()]

        expected = set(self.taxa)
        for name, tree in trees.items():
            observed = {terminal.name for terminal in tree.get_terminals()}
            if observed != expected:
                missing = expected ^ observed
                raise ValueError(
                    f"Árvore '{name}' não compartilha o mesmo conjunto de táxons "
                    f"(diferença: {sorted(missing)[:5]}...)."
                )

    @classmethod
    def from_directory(
        cls,
        path: str,
        pattern: str = ".nexus",
        tree_format: str = "nexus",
        prefix: str = "tree_dataset_final_",
        normalizer: Optional[Callable[[str], str]] = strip_accession_version,
    ) -> "TreeSet":
        """
        Carrega todas as árvores de um diretório `out/Trees`.

        Parameters
        ----------
        path : str
            Diretório contendo as árvores.
        pattern : str, optional
            Sufixo dos arquivos a carregar.
        tree_format : str, optional
            Formato aceito por `Bio.Phylo.parse`.
        prefix : str, optional
            Prefixo removido do nome do arquivo antes da análise do rótulo.
        normalizer : callable or None, optional
            Normalizador de rótulos de terminais; None desativa a reconciliação
            e faz o construtor rejeitar árvores com rótulos divergentes.

        Return
        ------
        TreeSet
            Conjunto de árvores pronto para análise.
        """
        trees: Dict[str, Tree] = {}
        labels: Dict[str, PipelineLabel] = {}

        origem: Dict[str, str] = {}

        for filename in sorted(os.listdir(path)):
            if not filename.endswith(pattern):
                continue
            label = PipelineLabel.parse(filename, prefix=prefix)
            if label.name in trees:
                # Dois arquivos mapeados ao mesmo pipeline: antes, o segundo
                # sobrescrevia o primeiro sem uma linha de log, e `M` — o
                # denominador de todo suporte metodológico — vinha menor que o
                # número de árvores em disco. Perder pipeline em silêncio é
                # inaceitável; recusar é (D19).
                raise ValueError(
                    f"Dois arquivos designam o mesmo pipeline '{label.name}': "
                    f"'{origem[label.name]}' e '{filename}'. O rótulo não distingue "
                    f"os dois métodos — acrescente o sufixo a INFERENCE_METHODS."
                )
            origem[label.name] = filename
            trees[label.name] = next(Phylo.parse(os.path.join(path, filename), tree_format))
            labels[label.name] = label

        return cls(trees, labels, normalizer=normalizer)

    def __len__(self) -> int:
        return len(self.trees)

    @property
    def n_taxa(self) -> int:
        """Número de táxons compartilhado por todas as árvores."""
        return len(self.taxa)


class StabilityAnalyzer:
    """
    Analisador de estabilidade de clados entre pipelines.

    Parameters
    ----------
    tree_set : TreeSet
        Conjunto de árvores a comparar.
    include_root : bool, optional
        Se True, o clado universal é mantido entre os itens. Só tem efeito na
        análise enraizada; como bipartição, o clado universal é trivial.
    rooted : bool, optional
        Se False (padrão), a unidade de comparação é a **bipartição canônica**
        não enraizada — ver D3. Se True, são os clados enraizados, o que só é
        legítimo quando todas as árvores têm raiz biologicamente significativa
        e comum. No pipeline atual, apenas o UPGMA tem: FastTree, IQ-TREE,
        RAxML e NJ emitem topologia não enraizada escrita com raiz trifurcante,
        e comparar seus clados mede a convenção de escrita do arquivo, não a
        topologia. Em VARV-6 isso reportava 75% de discordância entre três
        métodos que produzem a **mesma** topologia.

    Attributes
    ----------
    clade_sets : dict of str -> set of frozenset
        Unidade de comparação por pipeline: bipartições canônicas (padrão) ou
        clados enraizados (``rooted=True``).
    """

    def __init__(self, tree_set: TreeSet, include_root: bool = False,
                 rooted: bool = False) -> None:
        self.tree_set = tree_set
        self.include_root = include_root
        self.rooted = rooted

        if rooted:
            self.clade_sets: Dict[str, Set[FrozenSet[str]]] = {
                name: {identity.taxa for identity in clade_identities(tree, include_root)}
                for name, tree in tree_set.trees.items()
            }
        else:
            all_taxa = frozenset(tree_set.taxa)
            self.clade_sets = {
                name: {
                    bipartition
                    for identity in clade_identities(tree, include_root=True)
                    for bipartition in (canonical_bipartition(identity.taxa, all_taxa),)
                    if bipartition is not None
                }
                for name, tree in tree_set.trees.items()
            }

    def bipartition_counts(self) -> Dict[str, int]:
        """
        Número de bipartições não triviais recuperadas por cada pipeline.

        `03-metricas §3` exige reportar ``|B(T)|`` ao lado de toda RF: uma
        politomia reduz ``|B(T)|``, e a normalização por ``n - 3`` — que supõe
        árvore binária — passa a **subestimar** a distância. Sem esse número ao
        lado, um valor baixo é ambíguo entre "topologias parecidas" e "árvore
        malresolvida".

        Return
        ------
        dict
            ``pipeline -> |B(T)|``.
        """
        return {name: len(clades) for name, clades in self.clade_sets.items()}

    # ------------------------------------------------------------------ #
    # Suporte de clados
    # ------------------------------------------------------------------ #

    def clade_records(self) -> List[CladeRecord]:
        """
        Tabula cada clado distinto com o conjunto de pipelines que o recuperam.

        Return
        ------
        list of CladeRecord
            Registros ordenados por suporte decrescente e, em empate, por tamanho.
        """
        occurrences: Dict[FrozenSet[str], Set[str]] = {}
        for pipeline, clades in self.clade_sets.items():
            for taxa in clades:
                occurrences.setdefault(taxa, set()).add(pipeline)

        total = len(self.tree_set)
        records = [
            CladeRecord(
                taxa=taxa,
                size=len(taxa),
                pipelines=frozenset(pipelines),
                support=len(pipelines) / total,
            )
            for taxa, pipelines in occurrences.items()
        ]
        records.sort(key=lambda r: (-r.support, -r.size))
        return records

    def support_profile(self) -> Dict[int, Dict[str, int]]:
        """
        Conta, para cada nível k, quantos clados são recuperados por exatamente
        k pipelines e por pelo menos k pipelines.

        Return
        ------
        dict
            Mapeamento ``k -> {"exact": n, "cumulative": n}``.
        """
        records = self.clade_records()
        total = len(self.tree_set)
        profile: Dict[int, Dict[str, int]] = {}

        for k in range(total, 0, -1):
            counts = [len(r.pipelines) for r in records]
            profile[k] = {
                "exact": sum(1 for c in counts if c == k),
                "cumulative": sum(1 for c in counts if c >= k),
            }
        return profile

    def consensus_clades(self, min_support: float = 1.0, min_size: int = 2) -> List[CladeRecord]:
        """
        Retorna os clados cujo suporte metodológico atinge o limiar informado.

        Parameters
        ----------
        min_support : float, optional
            Suporte mínimo, em fração de pipelines.
        min_size : int, optional
            Número mínimo de terminais.

        Return
        ------
        list of CladeRecord
            Clados estáveis, do maior para o menor.
        """
        return [
            record
            for record in self.clade_records()
            if record.support >= min_support - 1e-9 and record.size >= min_size
        ]

    # ------------------------------------------------------------------ #
    # Mineração exata de padrões
    # ------------------------------------------------------------------ #

    def closed_patterns(self) -> List[Pattern]:
        """
        Enumera exatamente todos os conjuntos fechados de clados.

        Para cada subconjunto não vazio S de pipelines, o conjunto de clados
        comuns I(S) é fechado; seu suporte é o conjunto de todos os pipelines
        que contêm I(S). O custo é O(2**M * |C|), praticável porque M é o número
        de pipelines, não de táxons.

        Return
        ------
        list of Pattern
            Padrões fechados, sem repetição, ordenados por suporte decrescente.
        """
        pipelines = list(self.clade_sets)
        total = len(pipelines)
        seen: Dict[FrozenSet[FrozenSet[str]], Pattern] = {}

        for size in range(1, total + 1):
            for subset in itertools.combinations(pipelines, size):
                common: Set[FrozenSet[str]] = set.intersection(
                    *(self.clade_sets[name] for name in subset)
                )
                if not common:
                    continue
                key = frozenset(common)
                if key in seen:
                    continue
                carriers = frozenset(
                    name for name in pipelines if common <= self.clade_sets[name]
                )
                seen[key] = Pattern(
                    clades=tuple(sorted(common, key=lambda t: (-len(t), sorted(t)))),
                    pipelines=carriers,
                    support=len(carriers) / total,
                )

        patterns = list(seen.values())
        patterns.sort(key=lambda p: (-p.support, -p.size))
        return patterns

    def maximal_patterns(self, min_support: float = 0.5) -> List[Pattern]:
        """
        Retorna os conjuntos de clados maximais acima de um suporte mínimo.

        Um padrão é maximal quando nenhum outro padrão frequente o contém
        propriamente -- é o análogo exato do resultado do FPMax, obtido aqui por
        enumeração completa do reticulado.

        Parameters
        ----------
        min_support : float, optional
            Suporte mínimo, em fração de pipelines.

        Return
        ------
        list of Pattern
            Padrões maximais frequentes.
        """
        frequent = [p for p in self.closed_patterns() if p.support >= min_support - 1e-9]
        maximal: List[Pattern] = []

        for candidate in frequent:
            clades = set(candidate.clades)
            if any(clades < set(other.clades) for other in frequent):
                continue
            maximal.append(candidate)

        return maximal

    # ------------------------------------------------------------------ #
    # Distâncias topológicas e decomposição da variância
    # ------------------------------------------------------------------ #

    def rf_matrix(self, normalized: bool = True) -> Dict[str, Dict[str, Optional[float]]]:
        """
        Calcula a matriz de distâncias Robinson-Foulds entre os pipelines.

        A unidade de comparação é a que `clade_sets` guarda: bipartição canônica
        (padrão) ou clado enraizado (``rooted=True``).

        Parameters
        ----------
        normalized : bool, optional
            Se True, divide pelo máximo possível, que **depende do enraizamento**
            (`03-metricas §3`):

            ==================================  ==========================
            Caso                                Denominador
            ==================================  ==========================
            Não enraizada, binária              ``2 * (n - 3)``
            Enraizada, binária (raiz excluída)  ``2 * (n - 2)``
            ==================================  ==========================

            Usar ``2 * (n - 2)`` sobre bipartições desloca todo valor
            normalizado — era o segundo erro de D3, somado ao primeiro.

        Return
        ------
        dict
            Matriz aninhada ``pipeline -> pipeline -> distância``. O valor é
            ``None`` quando a distância **não está definida**: a RF não
            enraizada exige ``n >= 4`` e a enraizada, ``n >= 3``. A diagonal é
            sempre ``0``. "Não aplicável" nunca é um número
            (`04-rigor-cientifico §3`) — devolver ``0`` faria topologias
            incomparáveis passarem por idênticas.

        See Also
        --------
        bipartition_counts : ``|B(T)|`` por pipeline, exigido ao lado da RF
            quando há politomia.
        """
        names = list(self.clade_sets)
        n = self.tree_set.n_taxa
        minimo = 3 if self.rooted else 4
        definida = n >= minimo

        if normalized:
            denominator = 2 * (n - 2) if self.rooted else 2 * (n - 3)
        else:
            denominator = 1

        matrix: Dict[str, Dict[str, Optional[float]]] = {}
        for a in names:
            matrix[a] = {}
            for b in names:
                if a == b:
                    matrix[a][b] = 0.0
                    continue
                if not definida or denominator <= 0:
                    matrix[a][b] = None
                    continue
                raw = len(self.clade_sets[a] - self.clade_sets[b]) + len(
                    self.clade_sets[b] - self.clade_sets[a]
                )
                matrix[a][b] = raw / denominator
        return matrix

    def factor_effects(self, normalized: bool = True) -> Dict[str, Dict[str, float]]:
        """
        Decompõe a discordância topológica nos dois fatores do delineamento.

        Compara a distância RF média entre pares que diferem apenas no alinhador
        (mesmo método de inferência) com a média entre pares que diferem apenas
        no método de inferência (mesmo alinhador).

        Parameters
        ----------
        normalized : bool, optional
            Repassado a `rf_matrix`.

        Return
        ------
        dict
            Estatísticas por fator: ``n``, ``mean``, ``min`` e ``max``.
        """
        matrix = self.rf_matrix(normalized=normalized)
        labels = self.tree_set.labels
        buckets: Dict[str, List[float]] = {"aligner": [], "inference": [], "both": []}

        for a, b in itertools.combinations(matrix, 2):
            same_aligner = labels[a].aligner == labels[b].aligner
            same_inference = labels[a].inference == labels[b].inference
            distancia = matrix[a][b]
            if distancia is None:
                # RF indefinida (n pequeno demais) não entra em média nenhuma.
                continue
            if same_inference and not same_aligner:
                buckets["aligner"].append(distancia)
            elif same_aligner and not same_inference:
                buckets["inference"].append(distancia)
            elif not same_aligner and not same_inference:
                buckets["both"].append(distancia)

        return {
            factor: {
                "n": len(values),
                "mean": sum(values) / len(values) if values else None,
                "min": min(values) if values else None,
                "max": max(values) if values else None,
            }
            for factor, values in buckets.items()
        }

    # ------------------------------------------------------------------ #
    # Coerência taxonômica
    # ------------------------------------------------------------------ #

    def taxonomic_purity(
        self,
        classifier: Callable[[str], str],
        min_size: int = 3,
    ) -> List[Dict[str, object]]:
        """
        Avalia se clados estáveis correspondem a grupos taxonômicos coerentes.

        Parameters
        ----------
        classifier : callable
            Função que mapeia o nome de um terminal em um rótulo taxonômico.
        min_size : int, optional
            Tamanho mínimo do clado a ser avaliado.

        Return
        ------
        list of dict
            Um registro por clado, com composição, rótulo majoritário e pureza.
        """
        rows: List[Dict[str, object]] = []
        for record in self.clade_records():
            if record.size < min_size:
                continue
            composition: Dict[str, int] = {}
            for taxon in record.taxa:
                label = classifier(taxon)
                composition[label] = composition.get(label, 0) + 1
            majority, majority_count = max(composition.items(), key=lambda kv: kv[1])
            rows.append(
                {
                    "size": record.size,
                    "support": record.support,
                    "n_pipelines": len(record.pipelines),
                    "composition": composition,
                    "majority_label": majority,
                    "purity": majority_count / record.size,
                    "monophyletic_group": len(composition) == 1,
                }
            )
        return rows

    # ------------------------------------------------------------------ #
    # Auditoria da identidade legada
    # ------------------------------------------------------------------ #

    def legacy_audit(self) -> LegacyAudit:
        """
        Compara a identidade canônica com a de 16 bits dependente de ordem.

        Return
        ------
        LegacyAudit
            Contagens de fragmentação e colisão.
        """
        return audit_legacy_identity(self.tree_set.trees)

    def legacy_support_profile(self) -> Dict[int, int]:
        """
        Perfil de suporte reconstruído sob a identidade legada.

        Serve de contrafactual: mostra quantos clados o esquema anterior
        reconheceria em cada nível de suporte.

        Return
        ------
        dict
            Mapeamento ``k -> número de itens legados com suporte >= k``.
        """
        occurrences: Dict[int, Set[str]] = {}
        for pipeline, tree in self.tree_set.trees.items():
            for identity in clade_identities(tree, self.include_root):
                occurrences.setdefault(identity.legacy, set()).add(pipeline)

        total = len(self.tree_set)
        return {
            k: sum(1 for pipelines in occurrences.values() if len(pipelines) >= k)
            for k in range(total, 0, -1)
        }
