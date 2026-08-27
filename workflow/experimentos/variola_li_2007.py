"""
Aquisição do baseline de Li *et al.* (2007) — M2.1.

O experimento existia **comentado** no rodapé de `workflow_dataAcquisition.py`,
com as 48 accessions escritas à mão, o e-mail do Entrez em `"seu_email@dominio.com"`
e todos os parâmetros embutidos. Um bloco comentado não é reprodutível: não roda,
não é testado, e ninguém sabe se ainda corresponde ao que gerou os artefatos.

Aqui ele vira código executável e parametrizado. Três regras que o bloco antigo
violava:

1. **O e-mail do Entrez vem do ambiente**, nunca do código. O NCBI o exige para
   identificar quem consulta, e é dado pessoal: não entra em repositório
   ([`04-rigor`](../../../docs/automation/04-rigor-cientifico.md), regra 7).
2. **As accessions ficam num arquivo versionado**, não numa string de 48 linhas
   dentro de um `if __name__`. Elas *são* a definição do experimento.
3. **Os parâmetros são declarados num lugar só** e vão para o manifesto.

⚠️ **Este módulo baixa os dados; ele não corrige [D23](../../../docs/science/02-defeitos-que-alteram-resultado.md#d23).**
A consulta de origem traz, para o mesmo genoma, o registro do GenBank **e** a
cópia curada do RefSeq (`NC_008291` = `DQ437594`, `NC_003391` = `AF438165`), e o
pipeline descarta um dos dois adiante. A decisão do usuário em 2026-08-26 foi
**declarar agora e corrigir a aquisição depois** — corrigir muda a composição do
conjunto e portanto toda árvore publicada. `ACESSOS_DUPLICADOS_CONHECIDOS`
existe para que a checagem abaixo acuse o problema em vez de deixá-lo passar.
"""

from __future__ import annotations

import logging
import os
import pathlib
from typing import Dict, List, Optional, Tuple

__all__ = [
    "PARAMETROS",
    "OUTGROUP_QUERY",
    "ACESSOS_DUPLICADOS_CONHECIDOS",
    "email_do_ambiente",
    "carregar_acessos",
    "conferir_duplicatas",
    "montar_workflow",
]

#: Arquivo com um acesso por linha. As accessions **são** a definição do
#: experimento: versioná-las é o que permite dizer que duas execuções partiram
#: do mesmo conjunto.
ARQUIVO_ACESSOS = pathlib.Path(__file__).with_name("variola_li_2007_acessos.txt")

#: Parâmetros de Li *et al.* (2007), como reconstruídos em
#: [`08-ficha-de-fatos §4`](../../../docs/automation/08-ficha-de-fatos.md).
#: O genoma do VARV tem ~186 kb; os cortes de comprimento existem para descartar
#: submissões parciais. Vírus não têm UTR no sentido eucarioto, daí os `None`.
PARAMETROS: Dict[str, object] = {
    "initial_min_length": 180_000,
    "refined_min_length": 183_000,
    "utr5_end": None,
    "utr3_start": None,
    "similarity_threshold": 0.999,
    "retmax": 200,
}

#: Grupo externo declarado: *Taterapox* e *Camelpox*, as duas espécies irmãs de
#: *Variola* dentro de *Orthopoxvirus*. A janela de datas reproduz a do estudo.
OUTGROUP_QUERY = (
    '("Taterapox virus"[Organism] OR "Camelpox virus"[Organism]) '
    'AND (complete[All Fields] AND genome[All Fields]) '
    'AND "1900"[PDAT] : "2007"[PDAT]'
)

#: Pares em que a consulta devolve o mesmo genoma sob dois acessos — o do
#: GenBank e a cópia curada do RefSeq. Ver D23. Declarados aqui para que a
#: conferência os reconheça e os **nomeie**, em vez de o pipeline descartá-los
#: em silêncio adiante.
ACESSOS_DUPLICADOS_CONHECIDOS: Tuple[Tuple[str, str], ...] = (
    ("NC_008291", "DQ437594"),   # Taterapox virus
    ("NC_003391", "AF438165"),   # Camelpox virus
)


def email_do_ambiente(variavel: str = "NCBI_EMAIL") -> str:
    """
    E-mail do Entrez, lido do ambiente.

    Raises
    ------
    RuntimeError
        Com a instrução de como definir. O NCBI **exige** o e-mail e bloqueia
        quem consulta sem ele; falhar aqui, com a razão, é melhor do que falhar
        no meio de um download de horas.
    """
    email = os.environ.get(variavel, "").strip()
    if not email:
        raise RuntimeError(
            f"{variavel} não está definida. O NCBI exige um e-mail de contato "
            f"para consultas ao Entrez, e ele é dado pessoal: não pode ser "
            f"embutido no código.\n"
            f"  export {variavel}='voce@instituicao.br'"
        )
    return email


def carregar_acessos(caminho: Optional[os.PathLike] = None) -> List[str]:
    """
    Lê os acessos do arquivo versionado, ignorando comentários e linhas vazias.
    """
    caminho = pathlib.Path(caminho or ARQUIVO_ACESSOS)
    acessos = []
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        linha = linha.split("#", 1)[0].strip().rstrip(",")
        if linha:
            acessos.append(linha)
    return acessos


def conferir_duplicatas(acessos: List[str]) -> List[Tuple[str, str]]:
    """
    Aponta os pares RefSeq/GenBank presentes na lista.

    Não remove nada: a decisão de 2026-08-26 foi declarar, não corrigir. O que
    esta função impede é o caso que [D23](../../../docs/science/02-defeitos-que-alteram-resultado.md#d23)
    documenta — o par entrar na consulta sem que ninguém perceba, e o pipeline
    descartar um dos dois adiante, mudando `n` sem declarar.
    """
    presentes = {a.split(".")[0] for a in acessos}
    encontrados = [(a, b) for a, b in ACESSOS_DUPLICADOS_CONHECIDOS
                   if a in presentes and b in presentes]
    for refseq, genbank in encontrados:
        logging.warning(
            f"D23 — {refseq} e {genbank} são o mesmo genoma (RefSeq e GenBank). "
            f"Ambos estão na lista de acessos; o pipeline vai descartar um deles "
            f"adiante e `n` será menor que o número de acessos."
        )
    return encontrados


def montar_workflow(work_dir: str, email: Optional[str] = None, **sobrescritas):
    """
    Instancia o workflow de aquisição com os parâmetros do experimento.

    Parameters
    ----------
    work_dir : str
        Diretório de trabalho do download.
    email : str or None
        Sendo `None`, lido de `NCBI_EMAIL`.
    **sobrescritas
        Sobrepõem `PARAMETROS`. Toda sobrescrita muda o experimento e precisa
        ir para o manifesto — não use para "só testar".
    """
    from workflow.workflow_dataAcquisition import workflowAquisitionDatasetNCBI

    parametros = {**PARAMETROS, **sobrescritas}
    return workflowAquisitionDatasetNCBI(
        email=email or email_do_ambiente(),
        work_dir=work_dir,
        taxon_filter="Orthopoxvirus",   # D6 — filtro taxonômico declarado (M2.2)
        **parametros,
    )


def main(work_dir: str = "replication-RetMax200-ITRs") -> None:
    """Executa a aquisição completa. Precisa de rede e de `NCBI_EMAIL`."""
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")

    acessos = carregar_acessos()
    conferir_duplicatas(acessos)
    logging.info(f"{len(acessos)} acessos declarados para o baseline de Li et al. (2007)")

    workflow = montar_workflow(work_dir)
    workflow.run_workflow(
        query=",".join(acessos),
        outgroup_query_or_file=OUTGROUP_QUERY,
        download_method="query",
    )
    workflow.generate_fasta(
        os.path.join(work_dir, "dataset_with_outgroup.gb"),
        os.path.join(work_dir, "dataset_final.fasta"),
    )


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "replication-RetMax200-ITRs")
