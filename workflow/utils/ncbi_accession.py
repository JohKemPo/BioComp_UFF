"""
Preferência RefSeq/GenBank na aquisição — [D23](../../../docs/science/02-defeitos-que-alteram-resultado.md#d23).

Os conjuntos de *Variola* baixam, para o mesmo genoma, o registro original do
GenBank e a cópia curada do RefSeq (`NC_008291.1` == `DQ437594.1`,
`NC_003391.1` == `AF438165.1`) — dois acessos, uma sequência. Até
[DEC-082](../../../docs/automation/07-log-de-execucao.md#dec-082--2026-09-04--d23--decisão-do-usuário-entre-as-opções-a-d-a-preferir-refseq-relabelando-na-posição-de-hoje)
a deduplicação por conteúdo guardava a **primeira ocorrência no arquivo**, e a
ordem de download decidia qual acesso representava o táxon — o mesmo genoma
virava `DQ437594.1` num experimento e `NC_008291.1` noutro.

DEC-082 decidiu a opção **(A)**: quando um grupo de sequências idênticas
contém um acesso RefSeq, ele é o rótulo do sobrevivente — na MESMA posição em
que o sobrevivente por ordem de chegada já estava, sem reordenar o conjunto (a
opção (B), que reordenava, foi descartada por mudar a entrada do MAFFT sem
medição do efeito). `efetch` nos três pares conhecidos confirmou `source`
idêntico (mesmo `geo_loc_name`/`collection_date`, quando presentes): preferir
RefSeq não perde metadado filogeográfico nestes casos.
"""

from __future__ import annotations

import re

__all__ = ["eh_refseq"]

# Nomenclatura de acesso do NCBI: RefSeq é 1-3 letras + "_" + dígitos
# (com ou sem sufixo de versão, ex. "NC_008291" / "NC_008291.1"). Acesso
# GenBank/INSDC nunca tem "_" — é letras seguidas diretamente por dígitos
# (ex. "DQ437594", "AF438165").
_PADRAO_REFSEQ = re.compile(r"^[A-Za-z]{1,3}_\d")


def eh_refseq(accession: str) -> bool:
    """
    Devolve se `accession` segue o formato de acesso RefSeq.

    Parameters
    ----------
    accession : str
        Acesso a conferir, com ou sem sufixo de versão.

    Return
    ------
    bool
    """
    if not accession:
        return False
    return bool(_PADRAO_REFSEQ.match(accession))
