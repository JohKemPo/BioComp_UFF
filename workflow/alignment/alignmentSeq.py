from Bio import AlignIO, SeqIO
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

from workflow.utils import tool_runs

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
        
        
        self.num_threads = config.get('num_threads', psutil.cpu_count(logical=True))
        self.max_memory_gb = config.get('max_memory_gb', 4)  # Limite de memória em GB
        self.max_sequences = config.get('max_sequences', 500)  # Número máximo de sequências para alinhamento completo
        logging.info(f" Alinhamento configurado para usar até {self.num_threads} thread(s).")
        logging.info(f" Limite de memória: {self.max_memory_gb} GB")
        logging.info(f" Número máximo de sequências: {self.max_sequences}")
        
        tmp_dir = os.path.join(config.get('output_path'),'tmp','tmpFilesAlignment')
        os.makedirs(tmp_dir, exist_ok=True)
        self.temp_base_dir = tmp_dir if os.path.exists('/dev/shm') else tempfile.gettempdir()
    
    
    def _get_sequence_stats_stream(self, fasta_path):
        """Usa geradores para extrair estatísticas sem carregar o arquivo na RAM."""
        count = 0
        total_length = 0
        for record in SeqIO.parse(fasta_path, "fasta"):
            count += 1
            total_length += len(record.seq)
            
        avg_length = total_length / count if count > 0 else 0
        est_memory_gb = (count ** 2 * avg_length * 4) / (1024**3)
        return count, avg_length, est_memory_gb
    
    
    def _subsample_sequences_stream(self, fasta_path, output_path, max_sequences):
        """
        Subamostragem usando duas passagens contínuas.
        Impede carregamento na memória de datasets grandes.
        """
        total_sequences = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
        
        if total_sequences <= max_sequences:
            shutil.copy2(fasta_path, output_path)
            return output_path, total_sequences

        # F1: Define quais índices manter (Reservoir/Random Sampling)
        indices_to_keep = set(random.sample(range(total_sequences), max_sequences))
        
        # F2: Escreve apenas os índices selecionados sob demanda
        with open(output_path, "w") as out_f:
            records_to_write = (
                rec for i, rec in enumerate(SeqIO.parse(fasta_path, "fasta")) 
                if i in indices_to_keep
            )
            SeqIO.write(records_to_write, out_f, "fasta")
            
        logging.info(f"Subamostragem concluída: {total_sequences} -> {max_sequences} sequências")
        return output_path, max_sequences

        
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
        logging.info("STEP: Aligning with CLUSTALO")
        
        num_seqs, _, est_memory = self._get_sequence_stats_stream(fasta_path)
        temp_dir = None
        input_file = fasta_path

        try:
            
            if force_subsample or num_seqs > self.max_sequences or est_memory > self.max_memory_gb:
                temp_dir = tempfile.mkdtemp(dir=self.temp_base_dir)
                subsampled_path = os.path.join(temp_dir, "subsampled.fasta")
                input_file, _ = self._subsample_sequences_stream(fasta_path, subsampled_path, self.max_sequences)

            clustalo_cmd = [
                "clustalo",
                "-i", input_file,
                "-o", output_path_align,
                "--outfmt", "fasta",
                "--threads", str(self.num_threads),
                "--force"
            ]

            # Datasets muito grandes requerem menos iterações no HMM para não estourar tempo/RAM
            if num_seqs > 5000:
                clustalo_cmd.extend(["--max-guidetree-iterations", "1", "--max-hmm-iterations", "1"])

            tool_runs.registrar('clustalo', clustalo_cmd, saida=output_path_align,
                                threads=self.num_threads, n_sequencias=num_seqs)
            result = subprocess.run(clustalo_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"Falha no ClustalO: {result.stderr}")
                
            return AlignIO.read(output_path_align, "fasta")

        finally:
            if temp_dir and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
    
    def align_sequences_muscle(self,
                               fasta_path: str,
                               output_path_align: str,
                               output_path_html: str):
        """
        Alinha sequências utilizando o MUSCLE.

        Suporta as **duas gerações de linha de comando**, porque elas são
        incompatíveis e as duas circulam:

        ==========  ==========================================
        Versão      Invocação
        ==========  ==========================================
        3.8.x       ``muscle -in entrada -out saida``
        5.x         ``muscle -align entrada -output saida``
        ==========  ==========================================

        A versão instalada é detectada uma vez e registrada no log — não se
        adivinha por tentativa e erro, porque um erro de sintaxe e uma falha de
        alinhamento produzem o mesmo código de saída, e confundi-los esconderia
        a segunda.

        Parameters
        ----------
        fasta_path : str
            FASTA de entrada.
        output_path_align : str
            Caminho do alinhamento de saída, em FASTA.
        output_path_html : str
            Mantido por simetria com os outros alinhadores; o MUSCLE não emite HTML.

        Return
        ------
        Bio.Align.MultipleSeqAlignment
        """
        logging.info("STEP: Aligning with MUSCLE")

        num_seqs, avg_len, _ = self._get_sequence_stats_stream(fasta_path)
        versao = self._muscle_major_version()

        if versao >= 5:
            cmd = ["muscle", "-align", fasta_path, "-output", output_path_align]
        else:
            cmd = ["muscle", "-in", fasta_path, "-out", output_path_align]
            # O refinamento iterativo do MUSCLE 3.8 domina o custo em conjuntos
            # grandes. `-maxiters 2` é a recomendação do próprio manual acima de
            # algumas centenas de sequências; acima disso o padrão (16) não
            # termina em tempo útil.
            if num_seqs > 500:
                cmd.extend(["-maxiters", "2"])
                logging.info("MUSCLE 3.8: -maxiters 2 (conjunto grande)")

        logging.info(f"Executando comando: {' '.join(cmd)}")
        # A versão maior escolhe a sintaxe (`-align/-output` × `-in/-out`) e não
        # se lê da linha de comando: vai explícita. É a diferença que fez a
        # medição de custo do MUSCLE 3.8 não valer para o 5.3 (DEC-044).
        tool_runs.registrar('muscle', cmd, saida=output_path_align,
                            versao_maior=versao, n_sequencias=num_seqs)
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"Falha no MUSCLE: {result.stderr[-2000:]}")

        return AlignIO.read(output_path_align, "fasta")

    def _muscle_major_version(self) -> int:
        """
        Versão maior do MUSCLE instalado, para escolher a sintaxe.

        Devolve 3 quando não consegue determinar: a sintaxe antiga é a que está
        instalada nesta base, e errar para o lado conhecido é melhor que errar
        para o lado suposto.
        """
        if getattr(self, "_muscle_version_cache", None) is not None:
            return self._muscle_version_cache

        maior = 3
        for args in (["-version"], ["--version"]):
            try:
                r = subprocess.run(["muscle", *args], capture_output=True, text=True,
                                   stdin=subprocess.DEVNULL, timeout=10, check=False)
            except (OSError, subprocess.TimeoutExpired):
                continue
            achado = re.search(r"(\d+)\.\d+", (r.stdout or "") + (r.stderr or ""))
            if achado:
                maior = int(achado.group(1))
                break

        logging.info(f"MUSCLE detectado: versão maior {maior}")
        self._muscle_version_cache = maior
        return maior

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

        Parametrização dinâmica baseada em heurísticas de tamanho do dataset.
        
        Return
        ------
        str
            Alinhamento gerado pelo MAFFT como uma string no formato padrão de saída.
        """
        num_seqs, avg_len, _ = self._get_sequence_stats_stream(fasta_path)
        
        # Estratégia de Alinhamento Dinâmica
        if num_seqs < 100 and avg_len < 1000:
            # Alta precisão para datasets pequenos/médios
            strategy = ["--localpair", "--maxiterate", "1000"]
            logging.info("Estratégia MAFFT: L-INS-i (Alta Precisão)")
            logging.info("STEP: MAFFT Strategy: L-INS-i")
        elif num_seqs < 10000:
            # Balanceado
            strategy = ["--auto"]
            logging.info("Estratégia MAFFT: FFT-NS-1/2 (Auto)")
            logging.info("STEP: MAFFT Strategy: FFT-NS-1/2")
        else:
            # Escalonamento massivo usando construção de árvores particionadas
            strategy = ["--parttree", "--retree", "1", "--partsize", "1000"]
            logging.info("Estratégia MAFFT: PartTree (Escalonamento para Datasets Massivos)")
            logging.info("STEP: MAFFT Strategy: PartTree")
            

        cmd = ["mafft", "--thread", str(self.num_threads)] + strategy + [fasta_path]
        
        logging.info(f"Executando comando: {' '.join(cmd)}")
        
        # A estratégia (`--auto` × `--parttree`) é escolhida pelo tamanho do
        # conjunto e muda o alinhamento — logo, muda a árvore. Até aqui ela só
        # existia no log da execução, que não acompanha o artefato.
        tool_runs.registrar('mafft', cmd, saida=output_path_align,
                            threads=self.num_threads,
                            estrategia=' '.join(strategy))
        # Redirecionamento direto de stdout para evitar overhead de strings no Python
        with open(output_path_align, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"Falha no MAFFT: {result.stderr}")
            
        return AlignIO.read(output_path_align, "fasta")
    
   