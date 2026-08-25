"""
Filtro taxonômico declarado e verificação pós-download — lote M2.2.

A consulta de aquisição nunca restringiu o táxon, e o conjunto de dados foi
contaminado ([D6](../../../docs/science/02-defeitos-que-alteram-resultado.md#d6)):
`Nile crocodilepox virus` está em VARV-6, VARV-52 e VARV-121; `Yoka poxvirus` e
dois `Saltwater crocodilepox` estão em VARV-121. Nenhum deles é *Orthopoxvirus*
— são gêneros diferentes, e sua presença distorce o alinhamento, a topologia e
todo painel derivado.

Duas defesas, e a segunda existe porque a primeira não basta:

1. **Filtro na consulta** — `txid10242[Organism:exp]` restringe a busca ao
   clado na origem.
2. **Verificação pós-download** — confere a linhagem de cada registro **baixado**
   contra o clado declarado. Necessária porque o filtro da consulta não cobre
   todos os caminhos de entrada: o `download_from_csv` e um FASTA fornecido à
   mão nunca passam por uma consulta. Uma verificação que só olha a consulta
   confia no que deveria conferir.

A verificação é **offline**: a linhagem vem de `annotations['taxonomy']`, que
todo registro GenBank carrega. Isso a torna aplicável aos conjuntos que já
existem em disco, sem nova consulta ao NCBI.

**O filtro é declarado, nunca presumido.** O padrão é `None` — sem filtro —, e
nesse caso o manifesto registra que não houve filtro. Um projeto sobre outro
gênero declara outro táxon; o que não pode existir é a ausência silenciosa.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "TaxonFilter",
    "ORTHOPOXVIRUS",
    "ORTHOFLAVIVIRUS",
    "entrez_term",
    "record_lineage",
    "within_taxon",
    "TaxonomyAudit",
    "audit_records",
    "audit_genbank",
]


@dataclass(frozen=True)
class TaxonFilter:
    """
    Um clado declarado, pelo identificador do NCBI e pelo nome na linhagem.

    Os dois são necessários e servem a coisas diferentes: o `taxid` restringe a
    **consulta** ao NCBI; o `name` é o que aparece em `annotations['taxonomy']`
    e permite conferir um registro **já baixado**, sem rede.

    Attributes
    ----------
    taxid : str
        Identificador taxonômico, na forma ``"txid10242"``.
    name : str
        Nome do clado como o GenBank o escreve na linhagem (``"Orthopoxvirus"``).
    """

    taxid: str
    name: str

    def entrez_clause(self) -> str:
        """Cláusula de consulta que restringe a busca a este clado e seus descendentes."""
        return f"{self.taxid}[Organism:exp]"


#: *Orthopoxvirus* — o gênero de VARV, CMLV, CPXV, TATV e MPXV. É o clado do
#: baseline de Li *et al.* (2007) e o filtro que faltava em toda a aquisição.
ORTHOPOXVIRUS = TaxonFilter(taxid="txid10242", name="Orthopoxvirus")

#: *Orthoflavivirus* — o gênero do vírus Zika, na nomenclatura atual do ICTV
#: (o antigo *Flavivirus*). Existe para deixar explícito que **o filtro é do
#: experimento, não do projeto**: rodar os conjuntos de Zika contra
#: `ORTHOPOXVIRUS` reprova todos os 20 táxons, corretamente.
ORTHOFLAVIVIRUS = TaxonFilter(taxid="txid11051", name="Orthoflavivirus")


def entrez_term(query: str, taxon: Optional[TaxonFilter] = None) -> str:
    """
    Compõe o termo de busca do Entrez com o filtro taxonômico declarado.

    Parameters
    ----------
    query : str
        Consulta do experimento, como o pesquisador a escreveu.
    taxon : TaxonFilter or None, optional
        Clado a que restringir. ``None`` devolve a consulta intacta — o que é
        legítimo, desde que registrado como decisão.

    Return
    ------
    str
        Termo composto, ou a consulta original quando não há filtro.
    """
    consulta = (query or "").strip()
    if taxon is None or not consulta:
        return consulta
    if taxon.taxid in consulta:
        # Já declarado à mão na consulta; não duplicar.
        return consulta
    return f"({consulta}) AND {taxon.entrez_clause()}"


def record_lineage(record) -> Tuple[str, ...]:
    """
    Linhagem taxonômica de um registro GenBank, do reino ao gênero.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord
        Registro lido de um arquivo GenBank.

    Return
    ------
    tuple of str
        Nomes da linhagem. Vazia quando o registro não traz a anotação — o que
        é um fato a registrar, não um motivo para presumir pertencimento.
    """
    anotacoes = getattr(record, "annotations", None) or {}
    return tuple(anotacoes.get("taxonomy") or ())


def within_taxon(record, taxon: TaxonFilter) -> Optional[bool]:
    """
    Diz se um registro pertence ao clado declarado.

    Return
    ------
    bool or None
        ``True``/``False`` conforme a linhagem; ``None`` quando o registro **não
        tem** anotação de taxonomia. Ausência de informação não é o mesmo que
        pertencer, e nem que não pertencer — devolver ``False`` aqui descartaria
        registros por falta de metadado, e devolver ``True`` os aceitaria pelo
        mesmo motivo.
    """
    linhagem = record_lineage(record)
    if not linhagem:
        return None
    return taxon.name in linhagem


@dataclass
class TaxonomyAudit:
    """
    Resultado da conferência taxonômica de um conjunto de registros.

    Attributes
    ----------
    taxon : TaxonFilter
        Clado contra o qual se conferiu.
    within : list of str
        Acessos dentro do clado.
    outside : dict
        ``acesso -> {"organism": ..., "lineage": (...)}`` dos que estão fora.
    unknown : list of str
        Acessos sem anotação de taxonomia — indecidíveis sem consultar o NCBI.
    """

    taxon: TaxonFilter
    within: List[str] = field(default_factory=list)
    outside: Dict[str, Dict[str, object]] = field(default_factory=dict)
    unknown: List[str] = field(default_factory=list)

    @property
    def clean(self) -> bool:
        """Verdadeiro só quando nada está fora **e** nada é indecidível."""
        return not self.outside and not self.unknown

    @property
    def total(self) -> int:
        return len(self.within) + len(self.outside) + len(self.unknown)

    def summary(self) -> Dict[str, object]:
        """Resumo achatado, para o manifesto ou uma tabela."""
        return {
            "taxon": {"taxid": self.taxon.taxid, "name": self.taxon.name},
            "total": self.total,
            "within": len(self.within),
            "outside": {a: v["organism"] for a, v in sorted(self.outside.items())},
            "unknown": sorted(self.unknown),
            "clean": self.clean,
        }

    def raise_if_contaminated(self) -> None:
        """
        Levanta `ValueError` se houver registro fora do clado.

        O modo estrito existe para o dataset de referência, onde um táxon fora
        do clado invalida o alinhamento inteiro e não deve poder passar por
        distração.
        """
        if self.outside:
            itens = ", ".join(f"{a} ({v['organism']})" for a, v in sorted(self.outside.items()))
            raise ValueError(
                f"{len(self.outside)} registro(s) fora de {self.taxon.name}: {itens}")


def _acesso(record) -> str:
    bruto = getattr(record, "id", "") or getattr(record, "name", "") or "?"
    return bruto.split(".")[0]


def audit_records(records: Iterable, taxon: TaxonFilter = ORTHOPOXVIRUS) -> TaxonomyAudit:
    """
    Confere uma coleção de registros contra o clado declarado.

    Acessos repetidos — o `raw_data_sequences.gb` guarda o mesmo registro uma
    vez por árvore — contam uma vez só.

    Parameters
    ----------
    records : iterable
        Registros do Biopython.
    taxon : TaxonFilter, optional
        Clado exigido. O padrão é *Orthopoxvirus*.

    Return
    ------
    TaxonomyAudit
    """
    auditoria = TaxonomyAudit(taxon=taxon)
    vistos = set()

    for record in records:
        acesso = _acesso(record)
        if acesso in vistos:
            continue
        vistos.add(acesso)

        pertence = within_taxon(record, taxon)
        if pertence is None:
            auditoria.unknown.append(acesso)
        elif pertence:
            auditoria.within.append(acesso)
        else:
            anotacoes = getattr(record, "annotations", None) or {}
            auditoria.outside[acesso] = {
                "organism": anotacoes.get("organism", "desconhecido"),
                "lineage": record_lineage(record),
            }

    auditoria.within.sort()
    auditoria.unknown.sort()
    return auditoria


def audit_genbank(path: str, taxon: TaxonFilter = ORTHOPOXVIRUS) -> TaxonomyAudit:
    """
    Confere um arquivo GenBank em disco, sem consultar o NCBI.

    É o que permite auditar os conjuntos que **já existem** — a linhagem está
    dentro do próprio arquivo.

    Parameters
    ----------
    path : str
        Caminho do `.gb`.
    taxon : TaxonFilter, optional
        Clado exigido.

    Return
    ------
    TaxonomyAudit
    """
    from Bio import SeqIO

    with open(path, encoding="utf-8", errors="replace") as handle:
        return audit_records(SeqIO.parse(handle, "genbank"), taxon)
