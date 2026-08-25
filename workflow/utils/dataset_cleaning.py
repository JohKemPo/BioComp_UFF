"""
Limpeza taxonômica de conjuntos de dados — D6, insumo de M2.6.

Os conjuntos de *Variola* estão contaminados: `Nile crocodilepox` em VARV-6,
VARV-52 e VARV-121; mais `Yokapox` e dois `Saltwater crocodilepox` em VARV-121.
Nenhum é *Orthopoxvirus*.

**Nada é apagado.** O conjunto original permanece intacto — ele é o histórico de
como o workflow evoluiu, e serve como subamostra do conjunto completo para
teste. A limpeza **cria um conjunto novo, ao lado**, com:

* o FASTA sem os táxons fora do clado declarado;
* um `PROVENIENCIA.md` que diz de onde veio, o que saiu, por quê, e com que
  comando reproduzir.

Isso é deliberado e não é conservadorismo: um conjunto limpo sem proveniência é
tão indefensável quanto um conjunto contaminado. O que torna o novo conjunto
publicável não é ele estar limpo — é **ser possível provar** o que foi retirado.
"""

from __future__ import annotations

import datetime
import hashlib
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from workflow.utils.taxonomy import ORTHOPOXVIRUS, TaxonFilter, TaxonomyAudit

__all__ = ["CleaningReport", "clean_dataset"]


@dataclass
class CleaningReport:
    """
    O que a limpeza fez, em termos verificáveis.

    Attributes
    ----------
    source : str
        Conjunto de origem.
    destination : str
        Conjunto criado.
    taxon : TaxonFilter
        Clado exigido.
    kept : list of str
        Acessos mantidos.
    removed : dict
        ``acesso -> motivo``.
    unresolved : list of str
        Acessos presentes no FASTA para os quais **não havia** registro GenBank
        que permitisse decidir. Ficam no conjunto e são declarados: retirá-los
        seria descartar dado por falta de metadado.
    """

    source: str
    destination: str
    taxon: TaxonFilter
    kept: List[str] = field(default_factory=list)
    removed: Dict[str, str] = field(default_factory=dict)
    unresolved: List[str] = field(default_factory=list)

    @property
    def total_before(self) -> int:
        return len(self.kept) + len(self.removed) + len(self.unresolved)

    def as_markdown(self, comando: str) -> str:
        """O `PROVENIENCIA.md` que acompanha o conjunto limpo."""
        agora = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
        linhas = [
            f"# Proveniência — {os.path.basename(self.destination)}",
            "",
            "Conjunto derivado por remoção taxonômica. **O conjunto de origem não foi alterado.**",
            "",
            "| Campo | Valor |",
            "|---|---|",
            f"| Origem | `{self.source}` |",
            f"| Clado exigido | {self.taxon.name} (`{self.taxon.taxid}`) |",
            f"| Sequências antes | {self.total_before} |",
            f"| Sequências depois | {len(self.kept) + len(self.unresolved)} |",
            f"| Removidas | {len(self.removed)} |",
            f"| Sem linhagem, mantidas | {len(self.unresolved)} |",
            f"| Gerado em (UTC) | {agora} |",
            "",
        ]

        if self.removed:
            linhas += [
                "## Removidas",
                "",
                "| Acesso | Motivo |",
                "|---|---|",
            ]
            linhas += [f"| `{a}` | {m} |" for a, m in sorted(self.removed.items())]
            linhas.append("")

        if self.unresolved:
            linhas += [
                "## Mantidas sem decisão",
                "",
                "Sem registro GenBank que permitisse decidir o clado. **Permanecem no conjunto**: "
                "retirá-las seria descartar dado por falta de metadado, que é uma decisão diferente "
                "e precisa ser tomada explicitamente.",
                "",
                "".join(f"- `{a}`\n" for a in sorted(self.unresolved)),
            ]

        linhas += [
            "## Reproduzir",
            "",
            "```bash",
            comando,
            "```",
            "",
            "A conferência independente é `docs/science/scripts/auditar_taxonomia.py`.",
        ]
        return "\n".join(linhas)


def _acesso(cabecalho: str) -> str:
    return cabecalho.split()[0].lstrip(">").split(".")[0]


def clean_dataset(source_fasta: str,
                  destination_fasta: str,
                  genbank: Optional[str] = None,
                  taxon: TaxonFilter = ORTHOPOXVIRUS) -> CleaningReport:
    """
    Escreve uma cópia do FASTA sem os táxons fora do clado declarado.

    A decisão vem do arquivo GenBank correspondente, que traz a linhagem de cada
    registro. Sem GenBank, nada é removido — não se decide taxonomia por nome de
    arquivo.

    Parameters
    ----------
    source_fasta : str
        FASTA de origem. **Não é modificado.**
    destination_fasta : str
        Caminho do FASTA limpo a criar.
    genbank : str or None
        `.gb` com os registros correspondentes. Sem ele, todos os acessos ficam
        como `unresolved` e nada é removido.
    taxon : TaxonFilter
        Clado exigido.

    Return
    ------
    CleaningReport
    """
    from Bio import SeqIO

    from workflow.utils.taxonomy import audit_genbank

    auditoria: Optional[TaxonomyAudit] = None
    if genbank and os.path.exists(genbank):
        auditoria = audit_genbank(genbank, taxon)

    dentro = set(auditoria.within) if auditoria else set()
    fora = dict(auditoria.outside) if auditoria else {}

    relatorio = CleaningReport(source=source_fasta, destination=destination_fasta, taxon=taxon)

    os.makedirs(os.path.dirname(os.path.abspath(destination_fasta)), exist_ok=True)
    mantidos = []

    for record in SeqIO.parse(source_fasta, "fasta"):
        acesso = _acesso(record.id)
        if acesso in fora:
            info = fora[acesso]
            linhagem = " > ".join(info["lineage"][-2:]) if info["lineage"] else "linhagem ausente"
            relatorio.removed[acesso] = f"{info['organism']} — {linhagem}, fora de {taxon.name}"
            continue
        if acesso in dentro:
            relatorio.kept.append(acesso)
        else:
            relatorio.unresolved.append(acesso)
        mantidos.append(record)

    SeqIO.write(mantidos, destination_fasta, "fasta")
    relatorio.kept.sort()
    relatorio.unresolved.sort()
    return relatorio
