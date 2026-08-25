"""
Módulo de análise de estabilidade de clados entre pipelines (cross-pipeline clade stability).

Este pacote complementa o `subtree_mining`, oferecendo:

* `clade_identity`  -- identidade canônica de clados, invariante à ordem de travessia,
                       e auditoria da identidade legada de 16 bits;
* `stability`       -- perfil de suporte, mineração exata de conjuntos de clados
                       maximais/fechados, matriz Robinson-Foulds e decomposição da
                       variância topológica entre alinhador e método de inferência;
* `report`          -- tabelas e figuras derivadas da análise.

O ponto de entrada de alto nível é `StabilityAnalyzer`.
"""

from workflow.stability.clade_identity import (
    CladeIdentity,
    canonical_clade_id,
    legacy_clade_hash,
    legacy_terminal_hash,
    audit_legacy_identity,
)
from workflow.stability.stability import (
    StabilityAnalyzer,
    TreeSet,
    PipelineLabel,
)

__all__ = [
    "CladeIdentity",
    "canonical_clade_id",
    "legacy_clade_hash",
    "legacy_terminal_hash",
    "audit_legacy_identity",
    "StabilityAnalyzer",
    "TreeSet",
    "PipelineLabel",
]
