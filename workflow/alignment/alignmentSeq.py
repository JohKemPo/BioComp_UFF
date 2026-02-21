from Bio import AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline, MafftCommandline, ClustalOmegaCommandline
import html
import os, re
import tempfile
import subprocess
import random
import datetime
import logging
import shutil
from pathlib import Path
import psutil  

class AlignmentSeqs():
    """
    Classe responsável por alinhar sequências de DNA ou proteínas utilizando ferramentas de alinhamento como ClustalW e MAFFT.

    Esta classe fornece métodos para alinhar sequências utilizando ClustalW e MAFFT, retornando os alinhamentos em 
    formatos apropriados para posterior processamento filogenético.
    """
    def __init__(self, config: dict) -> None:
        """
        Inicializa a instância da classe AlignmentSeqs.

        Configura os atributos a partir dos argumentos fornecidos via kwargs, cria os diretórios
        de saída necessários e executa as limpezas de diretórios temporários.

        Return
        ------
        None
        """
        
        date = datetime.datetime.now()
        logging.basicConfig(level=logging.INFO, 
                    filename=config.get('logfile_path', 1),
                    format='%(asctime)s - %(levelname)s - %(message)s')
        
        
        self.num_threads = config.get('num_threads', 1)
        self.max_memory_gb = config.get('max_memory_gb', 4)  # Limite de memória em GB
        self.max_sequences = config.get('max_sequences', 500)  # Número máximo de sequências para alinhamento completo
        logging.info(f" Alinhamento configurado para usar até {self.num_threads} thread(s).")
        logging.info(f" Limite de memória: {self.max_memory_gb} GB")
        logging.info(f" Número máximo de sequências: {self.max_sequences}")
    
    def _check_memory_available(self, required_gb=2):
        """
        Verifica se há memória disponível suficiente.
        """
        try:
            memory = psutil.virtual_memory()
            available_gb = memory.available / (1024**3)
            logging.info(f" Memória disponível: {available_gb:.2f} GB")
            return available_gb >= required_gb
        except:
            # Se psutil não estiver disponível, assume que há memória
            return True
    
    def _count_sequences(self, fasta_path):
        """
        Conta o número de sequências em um arquivo FASTA.
        """
        count = 0
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    
    def _get_sequence_stats(self, fasta_path, sample_size=100):
        """
        Obtém estatísticas das sequências para estimar requisitos de memória.
        """
        lengths = []
        count = 0
        current_len = 0
        
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_len > 0:
                        lengths.append(current_len)
                        count += 1
                    current_len = 0
                else:
                    current_len += len(line.strip())
            if current_len > 0:
                lengths.append(current_len)
                count += 1
        
        if not lengths:
            return 0, 0, 0
        
        avg_length = sum(lengths) / len(lengths)
        max_length = max(lengths)
        
        # Estimar memória necessária (aproximadamente: n² * comprimento médio * 4 bytes)
        estimated_memory_gb = (count ** 2 * avg_length * 4) / (1024**3)
        
        return count, avg_length, estimated_memory_gb
    
    def _subsample_sequences(self, fasta_path, output_path, max_sequences):
        """
        Subamostra sequências aleatoriamente se houver muitas.
        """
        sequences = []
        current_header = None
        current_seq = []
        
        # Ler todas as sequências
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_header and current_seq:
                        sequences.append((current_header, ''.join(current_seq)))
                    current_header = line.strip()
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_header and current_seq:
                sequences.append((current_header, ''.join(current_seq)))
        
        total_sequences = len(sequences)
        
        if total_sequences <= max_sequences:
            # Se já está dentro do limite, copia o arquivo original
            shutil.copy2(fasta_path, output_path)
            return output_path, total_sequences
        
        # Selecionar sequências aleatoriamente
        selected = random.sample(sequences, max_sequences)
        
        # Escrever sequências selecionadas
        with open(output_path, 'w') as f:
            for header, seq in selected:
                f.write(f"{header}\n")
                # Quebrar sequência em linhas de 60 caracteres
                for i in range(0, len(seq), 60):
                    f.write(f"{seq[i:i+60]}\n")
        
        logging.info(f" Subamostragem: {total_sequences} -> {max_sequences} sequências")
        return output_path, max_sequences
    
    def _run_clustalo_with_timeout(self, cmd, timeout_minutes=600):
        """
        Executa ClustalO com timeout para evitar processos travados.
        """
        try:
            result = subprocess.run(
                cmd,
                timeout=timeout_minutes * 60,
                capture_output=True,
                text=True
            )
            return result
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"ClustalO excedeu o tempo limite de {timeout_minutes} minutos")
    
    def align_sequences_clustalo(self,
                                 fasta_path: str,
                                 output_path_align: str,
                                 output_path_html: str,
                                 force_subsample: bool = False):
        """
        Alinha sequências utilizando o Clustal Omega com suporte a paralelização e 
        gerenciamento de memória.

        Executa o comando Clustal Omega utilizando o número de threads definido
        na inicialização da classe.

        Args:
            fasta_path: Caminho para o arquivo FASTA de entrada
            output_path_align: Caminho para o arquivo de alinhamento de saída
            output_path_html: Caminho para o arquivo HTML de visualização
            force_subsample: Forçar subamostragem mesmo se memória for suficiente

        Return:
            Bio.Align.MultipleSeqAlignment: Objeto de alinhamento resultante
        """
        # Verificar memória disponível
        if not self._check_memory_available(required_gb=self.max_memory_gb):
            logging.warning("A  VISO: Memória baixa detectada. Usando configurações conservadoras.")
            force_subsample = True
        
        # Contar sequências
        num_sequences = self._count_sequences(fasta_path)
        logging.info(f" Número de sequências no arquivo: {num_sequences}")
        
        # Estimar requisitos de memória
        count, avg_len, est_memory = self._get_sequence_stats(fasta_path)
        logging.info(f" Comprimento médio das sequências: {avg_len:.0f} pb")
        logging.info(f" Memória estimada necessária: {est_memory:.2f} GB")
        
        # Criar arquivo temporário para subamostragem se necessário
        temp_dir = None
        input_file = fasta_path
        
        try:
            if force_subsample or num_sequences > self.max_sequences or est_memory > self.max_memory_gb:
                temp_dir = tempfile.mkdtemp()
                subsampled_path = os.path.join(temp_dir, "subsampled.fasta")
                input_file, actual_sequences = self._subsample_sequences(
                    fasta_path, 
                    subsampled_path, 
                    self.max_sequences
                )
                logging.info(f" Usando subamostragem com {actual_sequences} sequências")
            
            # Configurar parâmetros otimizados para memória
            clustalo_cmd = [
                "clustalo",
                "-i", input_file,
                "-o", output_path_align,
                "--outfmt", "fasta",
                "--auto",
                "--threads", str(self.num_threads),
                "--force",
                "--max-guidetree-iterations", "2",  # Reduzir iterações para economizar memória
                "--max-hmm-iterations", "2"
            ]
            
            # Executar ClustalO
            logging.info("  Executando Clustal Omega...")
            logging.info(f" Comando: {' '.join(clustalo_cmd)}")
            
            # Usar subprocess para melhor controle
            result = self._run_clustalo_with_timeout(clustalo_cmd)
            
            if result.returncode != 0:
                if result.returncode == 137:  # Killed (OOM)
                    error_msg = (
                        f"ClustalO foi morto por falta de memória (código 137).\n"
                        f"Tente reduzir ainda mais o número de sequências ou aumentar a memória disponível.\n"
                        f"Saída de erro: {result.stderr}"
                    )
                    raise RuntimeError(error_msg)
                else:
                    raise RuntimeError(f"Erro no Clustal Omega (código {result.returncode}): {result.stderr}")
            
            if result.stderr:
                logging.warning(f"  Status/Avisos do Clustal Omega:\n{result.stderr}")
            
            # Ler o alinhamento
            alignment = AlignIO.read(output_path_align, "fasta")
            
            
            
            return alignment
            
        except Exception as e:
            logging.error(f"    Erro durante o alinhamento: {e}")
            raise
            
        finally:
            # Limpar diretório temporário
            if temp_dir and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
    
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
        # Verificar número de sequências
        num_sequences = self._count_sequences(fasta_path)
        if num_sequences > 200:  # ClustalW é mais pesado
            logging.warning(f"  AVISO: ClustalW com {num_sequences} sequências pode consumir muita memória")
        
        clustalw_cline = ClustalwCommandline(
            "clustalw", 
            infile=fasta_path, 
            outfile=output_path_align
        )
        clustalw_cline()

        os.rename(path_dnd, output_path_dnd)
        alignment = AlignIO.read(output_path_align, "clustal")

        return alignment

    def align_sequences_mafft(self,
                              fasta_path: str,
                              output_path_align: str,
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
        # Verificar número de sequências
        num_sequences = self._count_sequences(fasta_path)
        
        # Configurar parâmetros MAFFT baseado no tamanho do dataset
        if num_sequences > 500:
            logging.warning(f"  Dataset grande ({num_sequences} sequências). Usando modo rápido do MAFFT.")
            # Usar estratégia FFT-NS-2 para datasets grandes
            mafft_cline = MafftCommandline(
                input=fasta_path, 
                thread=self.num_threads,
                auto=True,
                maxiterate=2
            )
        else:
            mafft_cline = MafftCommandline(
                input=fasta_path, 
                thread=self.num_threads, 
                auto=True
            )
        
        stdout, stderr = mafft_cline()
        
        if stderr:
            logging.error("     Erro durante a execução do MAFFT:")
            logging.error(f'    {stderr}')   
        
        with open(output_path_align, "w") as f:
            f.write(stdout)
        
        alignment = AlignIO.read(output_path_align, "fasta")
        
        return alignment
    
   