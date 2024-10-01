from Bio import AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MafftCommandline

import os

class AlignmentSeqs():
    """
    Classe responsável por alinhar sequências de DNA ou proteínas utilizando ferramentas de alinhamento como ClustalW e MAFFT.

    Esta classe fornece métodos para alinhar sequências utilizando ClustalW e MAFFT, retornando os alinhamentos em 
    formatos apropriados para posterior processamento filogenético.
    """
    def __init__(self) -> None:
        """
        Inicializa a instância da classe AlignmentSeqs.

        Configura os atributos a partir dos argumentos fornecidos via kwargs, cria os diretórios
        de saída necessários e executa as limpezas de diretórios temporários.

        Return
        ------
        None
        """

    def align_sequences_clustalw(self,
                                 fasta_path: str,
                                 path_dnd: str,
                                 output_path_dnd: str,
                                 output_path_align: str):
        """
        Alinha sequências utilizando o ClustalW.

        Executa o comando ClustalW para alinhar as sequências presentes no arquivo FASTA fornecido.
        Gera dois arquivos de saída: um contendo o alinhamento no formato Clustal e outro contendo 
        informações sobre o agrupamento hierárquico.

        Return
        ------
        Bio.Align.MultipleSeqAlignment
            Objeto de alinhamento resultante do ClustalW.
        """
        clustalw_cline = ClustalwCommandline("clustalw", infile=fasta_path, outfile=output_path_align)
        clustalw_cline()

        os.rename(path_dnd, output_path_dnd)
        return AlignIO.read(output_path_align, "clustal")

    def align_sequences_mafft(self,
                              fasta_path: str):
        """
        Alinha sequências utilizando o MAFFT.

        Executa o comando MAFFT para alinhar as sequências presentes no arquivo FASTA fornecido.
        Retorna o alinhamento como uma string.

        Return
        ------
        str
            Alinhamento gerado pelo MAFFT como uma string no formato padrão de saída.
        """
        mafft_cline = MafftCommandline(input=fasta_path)
        stdout, stderr = mafft_cline()
        return stdout