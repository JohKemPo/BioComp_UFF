
# from Bio import SeqIO
# from Bio import Entrez

# Entrez.email = 'joaovitormoraesjp@gmail.com'
# search_handler = Entrez.esearch(db="nucleotide", term="SPOV")

# search_records = Entrez.read(search_handler)

# for idx, record_id in enumerate(search_records['IdList']):

#   fetch_handler = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
#   fetch_records = SeqIO.parse(fetch_handler, 'genbank')

#   for record in fetch_records:

#     print(f'\nId: {idx} -------------------------------')
#     print(f'Record Id: {record_id}')
#     print('Record accession: ', record.id)
#     print('Record description: ', record.description)
#     print('Record sequence length: ', len(record.seq))
#     print('Record features count: ', len(record.features))

import os
import subprocess
import logging
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

class ZikaWorkflow:
    def __init__(self, email, work_dir="workflow_dataAcquisition",
                 initial_min_length=700, refined_min_length=9000,
                 utr5_end=None, utr3_start=None, similarity_threshold=0.99, retmax=1000):
        """
        Inicializa os parâmetros do workflow.
        
        Parameters
        ----------

        email : str
            E-mail para consulta ao GenBank (obrigatório para Entrez).
        work_dir : str
            Diretório base para salvar os arquivos intermediários e finais.
        initial_min_length : int
            Comprimento mínimo (em pb) para filtrar sequências inicialmente.
        refined_min_length : int
            Comprimento mínimo para o dataset refinado.
        utr5_end : int or str
            Posição final do UTR 5' (se conhecido); se None, não será removido.
        utr3_start : int or str
            Posição inicial do UTR 3' (se conhecido); se None, não será removido.
        similarity_threshold: float
            Limite para considerar duas sequências como idênticas (para remoção de duplicatas/overrepresentation).
        retmax : int
            Limite de sequências que serão baixadas.
        """
        Entrez.email = email
        self.work_dir = work_dir
        os.makedirs(self.work_dir, exist_ok=True)
        
        self.logger = logging.getLogger("ZikaWorkflow")
        self.logger.setLevel(logging.INFO)
        log_path = os.path.join(self.work_dir, "workflow.log")
        file_handler = logging.FileHandler(log_path)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        if not self.logger.handlers:
            self.logger.addHandler(file_handler)
        
        self.logger.info("Workflow iniciado com os seguintes parâmetros:")
        self.logger.info(f"email: {email}")
        self.logger.info(f"work_dir: {work_dir}")
        self.logger.info(f"initial_min_length: {initial_min_length}")
        self.logger.info(f"refined_min_length: {refined_min_length}")
        self.logger.info(f"utr5_end: {utr5_end}")
        self.logger.info(f"utr3_start: {utr3_start}")
        self.logger.info(f"similarity_threshold: {similarity_threshold}")
        self.logger.info(f"retmax: {retmax}")

        self.initial_min_length = initial_min_length
        self.refined_min_length = refined_min_length
        self.utr5_end = utr5_end
        self.utr3_start = utr3_start
        self.similarity_threshold = similarity_threshold
        self.retmax = retmax

    def download_sequences(self, query, output_file):
        """
        Passo 1: Baixa sequências do GenBank usando o Biopython.
        
        Parameters
        ----------

        query: str
            String de consulta para o GenBank.
        output_file: str
            Caminho para salvar as sequências baixadas (formato GenBank).
        """
        self.logger.info(f"Baixando sequências com query: {query}")
        print(f"Baixando sequências com query: {query}")
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=self.retmax)
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
            self.logger.info(f"Número de IDs encontrados: {len(id_list)}")
            handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
            with open(output_file, "w") as f:
                f.write(handle.read())
            handle.close()
            self.logger.info(f"Sequências salvas em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro no download de sequências: {e}")

    def filter_sequences(self, input_file, output_file):
        """
        Passo 2: Filtra as sequências para manter apenas aquelas com metadados e
        comprimento mínimo (e remove duplicatas).
        
        Parameters
        ----------

        input_file: str
            Arquivo de entrada com sequências (formato GenBank).
        output_file: str
            Arquivo de saída com sequências filtradas.
        """
        print(f"Iniciando filtragem de sequências: {input_file}")
        self.logger.info(f"Iniciando filtragem de sequências: {input_file}")
        try:
            records = list(SeqIO.parse(input_file, "genbank"))
            filtered = []
            seen_seqs = set()
            for rec in records:
                # Verifica se a sequência tem metadados (ex.: data, local) e comprimento mínimo
                if len(rec.seq) < self.initial_min_length:
                    continue
                seq_str = str(rec.seq).upper()
                if seq_str in seen_seqs:
                    continue
                seen_seqs.add(seq_str)
                filtered.append(rec)
            SeqIO.write(filtered, output_file, "genbank")
            self.logger.info(f"Total de sequências filtradas: {len(filtered)}")
            self.logger.info(f"Arquivo filtrado salvo em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro na filtragem de sequências: {e}")

    def remove_utrs(self, input_file, output_file):
        """
        Passo 3: Remove regiões UTR se as posições forem fornecidas.
        Caso contrário, apenas copia o arquivo de entrada.
        Parameters
        ----------

        input_file: str
            Arquivo de entrada com as sequências (formato GenBank).
        output_file: str
            Arquivo de saída com as sequências sem UTRs.
        """
        print("Removendo UTRs (se as posições estiverem definidas)...")
        self.logger.info("Removendo UTRs (se as posições estiverem definidas)...")
        try:
            records = []
            for rec in SeqIO.parse(input_file, "genbank"):
                if self.utr5_end is not None and self.utr3_start is not None:
                    # Considerando que a CDS esteja entre utr5_end+1 e utr3_start-1
                    cds_seq = rec.seq[rec.features[1].location.start:rec.features[1].location.end]
                    new_rec = SeqRecord(cds_seq, id=rec.id, name=rec.name,
                                        description=rec.description, annotations=rec.annotations)
                    records.append(new_rec)
                else:
                    records.append(rec)
            SeqIO.write(records, output_file, "genbank")
            self.logger.info(f"Sequências sem UTRs salvas em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro ao remover UTRs: {e}")

    def refine_dataset(self, input_file, output_file):
        """
        Passo 4: Refina o conjunto de dados:
          - Filtra por um comprimento mínimo mais alto (ex.: 9000 nt).
          - (Opcional) Remove sequências muito similares para reduzir super-representação.
        
        Parameters
        ----------

        input_file: str
            Arquivo de entrada com as sequências (após remoção de UTRs).
        output_file: str
            Arquivo de saída com o dataset refinado.
        """
        print("Refinando o dataset...")
        self.logger.info("Refinando o dataset...")
        try:
            records = list(SeqIO.parse(input_file, "genbank"))
            refined = []
            for rec in records:
                if len(rec.seq) < self.refined_min_length:
                    continue
                duplicate = False
                for existing in refined:
                    # Similaridade simples
                    matches = sum(1 for a, b in zip(str(rec.seq).upper(), str(existing.seq).upper()) if a == b)
                    similarity = matches / min(len(rec.seq), len(existing.seq))
                    if similarity >= self.similarity_threshold:
                        duplicate = True
                        break
                if not duplicate:
                    refined.append(rec)
            SeqIO.write(refined, output_file, "genbank")
            self.logger.info(f"Total de sequências no dataset refinado: {len(refined)}")
            self.logger.info(f"Dataset refinado salvo em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro ao refinar o dataset: {e}")

    def add_outgroup(self, input_file, outgroup_file, output_file):
        """
        Passo 5: Adiciona a(s) sequência(s) de outgroup ao dataset refinado.
        
        Parameters
        ----------

        input_file: str
            Arquivo do dataset refinado (formato GenBank).
        outgroup_file: str
            Arquivo com a(s) sequência(s) do outgroup (formato GenBank).
        output_file: str
            Arquivo final combinando dataset e outgroup.
        """
        print("Adicionando outgroup...")
        self.logger.info("Adicionando outgroup...")
        try:
            records_dataset = list(SeqIO.parse(input_file, "genbank"))
            records_outgroup = list(SeqIO.parse(outgroup_file, "genbank"))
            combined = records_dataset + records_outgroup
            SeqIO.write(combined, output_file, "genbank")
            self.logger.info(f"Dataset com outgroup salvo em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro ao adicionar outgroup: {e}")

    def align_sequences(self, input_file, output_file, mafft_path="mafft"):
        """
        Passo 6: Alinha as sequências utilizando MAFFT.
        
        Parameters
        ----------

        input_file: str
            Arquivo de entrada (formato Genbank ou FASTA).
        output_file: str
            Arquivo de saída com o alinhamento (formato FASTA).
        mafft_path: str
            Caminho para o executável do MAFFT (assumindo que esteja no PATH, por padrão "mafft").
        """
        print("Alinhando sequências com MAFFT...")
        self.logger.info("Alinhando sequências com MAFFT...")
        try:
            fasta_temp = os.path.join(self.work_dir, "temp_sequences.fasta")
            records = list(SeqIO.parse(input_file, "genbank"))
            SeqIO.write(records, fasta_temp, "fasta")
            
            cmd = [mafft_path, "--auto", fasta_temp]
            with open(output_file, "w") as outf:
                subprocess.run(cmd, stdout=outf, check=True)
            self.logger.info(f"Alinhamento salvo em: {output_file}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Erro ao executar MAFFT: {e}")
        except Exception as e:
            self.logger.error(f"Erro inesperado no alinhamento: {e}")

    def run_workflow(self, query, outgroup_query_or_file):
        """
        Executa o workflow completo:
          1. Baixar sequências.
          2. Filtrar sequências.
          3. Remover UTRs.
          4. Refinar o dataset.
          5. Adicionar outgroup.
          6. Alinhar sequências.
        
        Parameters
        ----------

        query: str
            Query para baixar as sequências do GenBank.
        outgroup_query_or_file: str
            Pode ser uma query para o outgroup ou um arquivo já existente.
        """
        self.logger.info("Iniciando o workflow completo...")
        raw_file = os.path.join(self.work_dir, "raw_sequences.gb")
        filtered_file = os.path.join(self.work_dir, "filtered_sequences.gb")
        no_utrs_file = os.path.join(self.work_dir, "no_utrs_sequences.gb")
        refined_file = os.path.join(self.work_dir, "refined_dataset.gb")
        dataset_outgroup_file = os.path.join(self.work_dir, "dataset_with_outgroup.gb")
        alignment_file = os.path.join(self.work_dir, "final_alignment.fasta")
        
        # Passo 1: Baixar sequências 
        self.download_sequences(query, raw_file)
        
        # Passo 2: Filtrar sequências
        self.filter_sequences(raw_file, filtered_file)
        
        # Passo 3: Remover UTRs (se parâmetros definidos)
        self.remove_utrs(filtered_file, no_utrs_file)
        
        # Passo 4: Refinar o dataset
        self.refine_dataset(no_utrs_file, refined_file)
        
        # Passo 5: Adicionar outgroup
        if os.path.isfile(outgroup_query_or_file):
            outgroup_file = outgroup_query_or_file
            self.logger.info("Outgroup informado por arquivo.")
        else:
            outgroup_file = os.path.join(self.work_dir, "outgroup.gb")
            self.logger.info(f"Baixando outgroup com query: {outgroup_query_or_file}")
            self.download_sequences(outgroup_query_or_file, outgroup_file)
        
        self.add_outgroup(refined_file, outgroup_file, dataset_outgroup_file)
        
        # Passo 6: Alinhar sequências
        self.align_sequences(dataset_outgroup_file, alignment_file)
        
        self.logger.info("Workflow concluído com sucesso.")

    def generate_fasta(self, input_file, output_file):
        """
        Converte um arquivo de sequências (ex.: GenBank) em um arquivo FASTA.
        
        input_file: str
            Caminho do arquivo de entrada (formato GenBank).
        output_file: str
            Caminho do arquivo de saída (formato FASTA) que conterá todas as sequências.
        """
        self.logger.info(f"Gerando arquivo FASTA a partir de: {input_file}")
        try:
            records = list(SeqIO.parse(input_file, "genbank"))
            SeqIO.write(records, output_file, "fasta")
            self.logger.info(f"Arquivo FASTA gerado em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro ao gerar FASTA: {e}")

    def slice_file(self, input_file, output_prefix, slice_size=50):
        """
        Divide um arquivo de sequências (FASTA) em múltiplos arquivos, cada um contendo 
        no máximo 'slice_size' sequências.
        
        input_file: str
            Caminho para o arquivo FASTA de entrada.
        output_prefix: str
            Prefixo para os arquivos de saída; serão gerados arquivos com nome
                              <output_prefix>_1.fasta, <output_prefix>_2.fasta, etc.
        slice_size: str
            Número máximo de sequências por arquivo (default: 50).
        """
        self.logger.info(f"Iniciando a divisão do arquivo: {input_file}")
        try:
            records = list(SeqIO.parse(input_file, "fasta"))
            total_records = len(records)
            self.logger.info(f"Número total de sequências encontradas: {total_records}")
            
            output_dir = os.path.join(os.path.dirname(input_file), "slices")
            os.makedirs(output_dir, exist_ok=True)
            
            slice_number = 1
            for i in range(0, total_records, slice_size):
                slice_records = records[i:i+slice_size]
                output_file = os.path.join(output_dir, f"{output_prefix}_{slice_number}.fasta")
                SeqIO.write(slice_records, output_file, "fasta")
                self.logger.info(f"Arquivo gerado: {output_file} (Sequências: {len(slice_records)})")
                slice_number += 1
        except Exception as e:
            self.logger.error(f"Erro ao dividir o arquivo: {e}")
    
if __name__ == "__main__":
    path = "workflow_dataAcquisition_1"
    workflow = ZikaWorkflow(work_dir=path,
                            email="joaovitormoraesjp@gmail.com",
                            utr5_end=True,  
                            utr3_start=True,  
                            similarity_threshold=0.99)
    
    # QUERY: Primeira tentativa 270hits
    zika_query = "Zika virus[Organism] AND complete genome"
    outgroup_query = "Spondweni virus[Organism] AND complete genome"
    
    # QUERY: Segunda tentativa 128hits
    # zika_query = '"Zika virus"[Organism] OR Zika virus[All Fields]'
    # outgroup_query = 'spondweni[All Fields] AND ("Viruses"[Organism] OR viruses[All Fields])'
    
    workflow.run_workflow(query=zika_query, outgroup_query_or_file=outgroup_query)
    
    
    input_genbank = f"{path}/dataset_with_outgroup.gb"
    output_fasta = f"{path}/dataset_final.fasta"
    workflow.generate_fasta(input_genbank, output_fasta)
    
    input_fasta = f"{path}/dataset_final.fasta"
    workflow.slice_file(input_fasta, output_prefix="dataset_slice", slice_size=50)