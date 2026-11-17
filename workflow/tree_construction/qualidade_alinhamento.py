"""
P2 / DM-10 — perfil de lacuna e filtro de colunas do alinhamento.

`docs/science/19-protocolo-reexecucao-baselines.md` §3.2: coluna com **mais de**
`LACUNA_MAX_COLUNA` de lacuna sai antes da construção de árvore; táxon que,
depois disso, fica com mais de `LACUNA_MAX_TAXON` de lacuna é **reportado, nunca
removido** — a composição do conjunto é a do baseline e não muda por qualidade
de alinhamento. Se o filtro removeria mais de `LIMITE_APLICACAO` das colunas, ele
não se aplica: tanta remoção redefine a região analisada em vez de limpar ruído
(Zika-21: 93,0%; VARV-52: 13,0% a 13,3%). Só `-` conta como lacuna; `n` é ambíguo, não ausente (diferente
de U4, que não conta `n` como base).

O filtro é opt-in, mas o perfil é medido sempre: o usuário é informado antes de
qualquer aplicação.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

LACUNA_MAX_COLUNA = 0.50
LACUNA_MAX_TAXON = 0.90
LIMITE_APLICACAO = 0.25

MOTIVO_DESLIGADO = "desligado_por_configuracao"
MOTIVO_ACIMA_DO_TETO = "acima_do_teto"

_LACUNA = ord("-")


def _matriz_lacuna(alignment) -> np.ndarray:
    n, comprimento = len(alignment), alignment.get_alignment_length() if len(alignment) else 0
    if n == 0 or comprimento == 0:
        return np.zeros((n, comprimento), dtype=bool)
    bruto = "".join(str(r.seq) for r in alignment).encode("ascii", "replace")
    return np.frombuffer(bruto, dtype=np.uint8).reshape(n, comprimento) == _LACUNA


def medir_perfil(alignment, lacuna_max_coluna: float = LACUNA_MAX_COLUNA,
                 lacuna_max_taxon: float = LACUNA_MAX_TAXON) -> Tuple[np.ndarray, Dict]:
    """
    Mede o perfil de lacuna sem alterar o alinhamento.

    Parameters
    ----------
    alignment : Bio.Align.MultipleSeqAlignment
    lacuna_max_coluna : float
        Fração de lacuna acima da qual a coluna é removível (estritamente maior).
    lacuna_max_taxon : float
        Fração de lacuna, nas colunas que sobram, acima da qual o táxon é
        reportado (estritamente maior).

    Return
    ------
    (manter, perfil) : Tuple[np.ndarray, dict]
        `manter` é a máscara booleana das colunas que o filtro conservaria.
        `perfil` traz `colunas_antes`, `colunas_removiveis`,
        `fracao_removivel` (`None` com 0 colunas), `colunas_depois` (se o
        filtro for aplicado) e `taxons_acima_limiar`, calculada sobre as colunas
        conservadas; a fração de um táxon sem coluna restante é `None`.
    """
    for nome, valor in (("lacuna_max_coluna", lacuna_max_coluna),
                        ("lacuna_max_taxon", lacuna_max_taxon)):
        if not 0.0 <= valor <= 1.0:
            raise ValueError(f"{nome} deve estar em [0, 1]: {valor!r}")
    lacuna = _matriz_lacuna(alignment)
    n, colunas = lacuna.shape
    # Divisão, não `limiar * n`: a comparação com o literal do limiar fica exata.
    fracao_coluna = lacuna.sum(axis=0) / n if n else np.zeros(colunas)
    manter = ~(fracao_coluna > lacuna_max_coluna)
    depois = int(manter.sum())

    acima: List[Dict] = []
    if depois:
        fracao_taxon = lacuna[:, manter].sum(axis=1) / depois
        for registro, fracao in zip(alignment, fracao_taxon):
            if fracao > lacuna_max_taxon:
                acima.append({"id": registro.id, "fracao_lacuna": float(fracao)})
    return manter, {
        "colunas_antes": colunas,
        "colunas_removiveis": colunas - depois,
        "fracao_removivel": (colunas - depois) / colunas if colunas else None,
        "colunas_depois": depois,
        "taxons_acima_limiar": acima,
    }


def filtrar_colunas(alignment, manter: np.ndarray) -> MultipleSeqAlignment:
    """Alinhamento só com as colunas de `manter`; todos os táxons, na mesma ordem."""
    indices = np.flatnonzero(manter)
    registros = []
    for r in alignment:
        texto = str(r.seq)
        registros.append(SeqRecord(Seq("".join(texto[i] for i in indices)),
                                   id=r.id, name=r.name, description=r.description))
    return MultipleSeqAlignment(registros)


def avaliar_alinhamento(alignment, ligado: bool = False,
                        lacuna_max_coluna: float = LACUNA_MAX_COLUNA,
                        lacuna_max_taxon: float = LACUNA_MAX_TAXON,
                        teto_aplicacao: float = LIMITE_APLICACAO
                        ) -> Tuple[MultipleSeqAlignment, Dict]:
    """
    Mede o perfil e, só se `ligado` e dentro do teto, devolve o alinhamento filtrado.

    Sem aplicação, o alinhamento sai idêntico e o relatório diz `aplicado: False`,
    com `motivo` (`desligado_por_configuracao` ou `acima_do_teto`) e
    `colunas_depois` igual a `colunas_antes`; `colunas_removiveis` informa o que
    o filtro removeria. O teto vale quando a fração removível é **maior** que
    `teto_aplicacao`.
    """
    if not 0.0 <= teto_aplicacao <= 1.0:
        raise ValueError(f"teto_aplicacao deve estar em [0, 1]: {teto_aplicacao!r}")
    manter, perfil = medir_perfil(alignment, lacuna_max_coluna, lacuna_max_taxon)
    fracao = perfil["fracao_removivel"]
    acima = fracao is not None and fracao > teto_aplicacao
    aplicar = bool(ligado) and not acima
    motivo = None if aplicar else (MOTIVO_ACIMA_DO_TETO if ligado else MOTIVO_DESLIGADO)
    relatorio = {"aplicado": aplicar, "motivo": motivo,
                 "lacuna_max_coluna": lacuna_max_coluna,
                 "lacuna_max_taxon": lacuna_max_taxon,
                 "teto_aplicacao": teto_aplicacao, **perfil}
    if not aplicar:
        relatorio["colunas_depois"] = perfil["colunas_antes"]
        relatorio["fracao_removida"] = 0.0 if perfil["colunas_antes"] else None
        return alignment, relatorio
    relatorio["fracao_removida"] = fracao
    return filtrar_colunas(alignment, manter), relatorio
