from Bio import AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MafftCommandline
import html
import os, re
from pathlib import Path
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
                                 output_path_align: str,
                                 output_path_html: str):
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
        alignment = AlignIO.read(output_path_align, "clustal")

        self.save_alignment_as_html_from_file(output_path_align, output_path_html)

        return alignment

    def align_sequences_mafft(self,
                              fasta_path: str,
                              output_path_html: str):
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
    
    # def save_alignment_as_html(self, alignment, output_html_path):
    #     with open(output_html_path, "w") as html_file:
    #         html_file.write("<html><head><style>")
    #         html_file.write("pre {font-family: monospace; font-size: 14px;}")
    #         html_file.write("</style></head><body>")
    #         html_file.write("<h2>Alinhamento de Sequências</h2>")
    #         html_file.write("<pre>")
    #         for record in alignment:
    #             html_file.write(f"{html.escape(record.id):<15} {html.escape(str(record.seq))}\n")
    #         html_file.write("</pre>")
    #         html_file.write("</body></html>")
    def parse_clustal_blocks(self, aln_path):
        """
        Lê um arquivo Clustal (.aln) com blocos e linha de consenso.
        Retorna um dicionário de sequências e uma lista de linhas de consenso por bloco.
        """
        with open(aln_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        sequence_data = {}
        consensus_lines = []
        current_block_consensus = []

        for line in lines:
            if line.startswith("CLUSTAL") or line.strip() == "":
                continue

            if re.match(r"^\s+[.*: ]", line):
                # linha de consenso, preserva os espaços à esquerda
                current_block_consensus.append(line.rstrip("\n"))
                continue

            parts = line.rstrip("\n").split(maxsplit=1)
            if len(parts) != 2:
                continue

            id_, seq = parts
            sequence_data[id_] = sequence_data.get(id_, "") + seq

            # detecta fim do bloco ao chegar ao final das sequências
            if len(current_block_consensus) > 0:
                consensus_lines.append(current_block_consensus.pop(0))

        return sequence_data, consensus_lines

    def build_alignment_blocks_html(self, sequences, consensus_lines, block_size=60):
        """
        Gera os blocos HTML para substituição no template.
        """
        ids = list(sequences.keys())
        total_length = len(next(iter(sequences.values())))
        num_blocks = (total_length + block_size - 1) // block_size
        blocks_html = []

        for block_idx in range(num_blocks):
            start = block_idx * block_size
            end = min(start + block_size, total_length)
            block_html = ['<div class="block">']

            # Sequências
            for id_ in ids:
                segment = sequences[id_][start:end]
                seq_html = ''.join(f'<span class="res">{html.escape(res)}</span>' for res in segment)
                block_html.append(
                    f'<div class="line">'
                    f'<span class="seq-id">{html.escape(id_)}</span>'
                    f'<span class="residues">{seq_html}</span>'
                    f'<span class="pos" style="margin-left: 20px;">{start + len(segment)}</span>'
                    f'</div>'
                )

            # Linha de consenso (se houver)
            if block_idx < len(consensus_lines):
                consensus_segment = consensus_lines[block_idx][start:end]
                consensus_html = ''.join(f'<span class="res">{html.escape(c)}</span>' for c in consensus_segment)
                block_html.append(
                    f'<div class="line consensus">'
                    f'<span class="seq-id"></span>'
                    f'<span class="residues">{consensus_html}</span>'
                    f'</div>'
                )

            block_html.append('</div>')
            blocks_html.append('\n'.join(block_html))

        return '\n'.join(blocks_html)


    def save_alignment_as_html_from_file(self, aln_path, output_html_path):
        """
        Lê um arquivo Clustal (.aln), aplica ao template HTML e salva como arquivo final.
        """
        sequences, consensus_lines = self.parse_clustal_blocks(aln_path)
        html_blocks = self.build_alignment_blocks_html(sequences, consensus_lines)

        with open("workflow/alignment/utils/template/alignmentViz.html", "r", encoding="utf-8") as f:
            template = f.read()

        
        title = Path(aln_path).stem

        result_html = template.replace("{{alignment_blocks}}", html_blocks)
        result_html = result_html.replace("{{alignment_title}}", title)

        Path(output_html_path).write_text(result_html, encoding="utf-8")
        return output_html_path
