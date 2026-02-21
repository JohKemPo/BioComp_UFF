import os
import subprocess
import logging
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio.Align import PairwiseAligner
import hashlib
import gc

class workflowAquisitionDatasetNCBI:
    def __init__(self, email, work_dir="workflow_dataAcquisition",
                 initial_min_length=700, refined_min_length=700,
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
        if initial_min_length and refined_min_length and initial_min_length > refined_min_length:
            raise ValueError("initial_min_length não pode ser maior que refined_min_length")
        
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
        
        self.custom_filters = []

    def download_sequences(self, query, output_file):
        """
        Baixa sequências do GenBank usando o Biopython.
        
        Parameters
        ----------
        query : str
            String de consulta para o GenBank. Se vazia, cria arquivo vazio.
        output_file : str
            Caminho para salvar as sequências baixadas (formato GenBank).
            
        Raises
        ------
        ValueError
            Se os parâmetros forem inválidos.
        IOError
            Se não for possível escrever o arquivo de saída.
        """
        self.logger.info(f"Baixando sequências com query: {query}")
        print(f"Baixando sequências com query: {query}")
        try:
            if not query or query.strip() == "":
                self.logger.warning("Query vazia, criando arquivo vazio")
                with open(output_file, "w") as f:
                    f.write("")
                return
            
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=self.retmax)
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                self.logger.warning("Nenhum resultado encontrado para a query")
                with open(output_file, "w") as f:
                    f.write("")
                return
            
            id_list = record["IdList"]
            self.logger.info(f"Número de IDs encontrados: {len(id_list)}")
            handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
            with open(output_file, "w") as f:
                f.write(handle.read())
            handle.close()
            self.logger.info(f"Sequências salvas em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro no download: {e}")
            # Criar arquivo vazio para não quebrar o pipeline
            with open(output_file, "w") as f:
                f.write("")

    
    def download_from_csv(self, csv_path, output_file):
        """
        Baixa sequências com base em accession numbers listados em um CSV.
        
        Parameters
        ----------
        csv_path : str
            Caminho para o arquivo CSV contendo coluna 'Accession'.
        output_file : str
            Caminho do arquivo de saída (formato GenBank).
        """
        self.logger.info(f"Lendo CSV: {csv_path}")
        try:
            df = pd.read_csv(csv_path)
            if 'Accession' not in df.columns:
                raise ValueError("O CSV precisa conter uma coluna 'Accession'")
            
            accession_list = df['Accession'].dropna().astype(str).tolist()
            self.logger.info(f"Total de accessions encontrados: {len(accession_list)}")

            handle = Entrez.efetch(db="nucleotide", id=accession_list, rettype="gb", retmode="text")
            with open(output_file, "w") as f:
                f.write(handle.read())
            handle.close()
            self.logger.info(f"Sequências baixadas e salvas em: {output_file}")
        except Exception as e:
            self.logger.error(f"Erro ao processar CSV: {e}")

            
    def filter_sequences(self, input_file, output_file):
        """
        Versão otimizada com processamento em lote e limpeza de memória.
        """
        print(f"Iniciando filtragem otimizada: {input_file}")
        self.logger.info(f"Iniciando filtragem otimizada: {input_file}")
        
        if not os.path.exists(input_file) or os.path.getsize(input_file) == 0:
            with open(output_file, "w") as f:
                f.write("")
            return []
        
        try:
            filtered = []
            seen_seqs = set()
            batch_size = 5  # Processar em lotes pequenos
            
            # Usar gerador, NÃO list()
            for i, rec in enumerate(SeqIO.parse(input_file, "genbank")):
                try:
                    if rec.seq is None or len(rec.seq) == 0:
                        continue
                        
                    if self.initial_min_length is not None and len(rec.seq) < self.initial_min_length:
                        continue
                        
                    seq_str = str(rec.seq).upper()
                    if seq_str in seen_seqs:
                        continue
                        
                    seen_seqs.add(seq_str)
                    filtered.append(rec)
                    
                    # A cada batch_size registros, forçar limpeza de memória
                    if len(filtered) % batch_size == 0:
                        gc.collect()  # Forçar garbage collector
                        
                except Exception as e:
                    self.logger.warning(f"Erro ao processar sequência {rec.id}: {e}")
                    continue
            
            # Salvar resultado
            SeqIO.write(filtered, output_file, "genbank")
            self.logger.info(f"Total de sequências filtradas: {len(filtered)}")
            
            # Limpar referências grandes
            del seen_seqs
            del filtered
            gc.collect()
            
            return filtered
            
        except Exception as e:
            self.logger.error(f"Erro na filtragem: {e}")
            with open(output_file, "w") as f:
                f.write("")
            return []

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
                    # Encontrar feature CDS explicitamente
                    cds_features = [feat for feat in rec.features if feat.type == "CDS"]
                    if cds_features:
                        cds = cds_features[0]
                        cds_seq = cds.extract(rec.seq)
                        new_rec = SeqRecord(cds_seq, id=rec.id, name=rec.name,
                                        description=rec.description + " | CDS only", 
                                        annotations=rec.annotations)
                        records.append(new_rec)
                    else:
                        self.logger.warning(f"Sequência {rec.id} não tem feature CDS, mantendo original")
                        records.append(rec)
                else:
                    records.append(rec)
            SeqIO.write(records, output_file, "genbank")
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
            seq_hashes = set()
            
            for rec in records:
                if self.refined_min_length is not None and len(rec.seq) < self.refined_min_length:
                    continue
                
                # Usar hash para detecção rápida de duplicatas exatas
                seq_hash = hashlib.md5(str(rec.seq).encode()).hexdigest()
                if seq_hash in seq_hashes:
                    continue
                    
                # Para similaridade, usar amostragem ou métodos mais eficientes
                # if not self.is_similar_to_existing(rec, refined):
                refined.append(rec)
                seq_hashes.add(seq_hash)
                    
            SeqIO.write(refined, output_file, "genbank")
            self.logger.info(f"Dataset refinado com {len(refined)} sequências")
        except Exception as e:
            self.logger.error(f"Erro ao Refinar dataset: {e}")


    def is_similar_to_existing(self, record, existing_records, sample_size=5):
        """
        Verifica se uma sequência é similar às sequências existentes.
        Usa amostragem para melhor performance com grandes datasets.
        
        Parameters
        ----------
        record : SeqRecord
            Sequência a ser verificada.
        existing_records : list
            Lista de SeqRecords já incluídos.
        sample_size : int
            Número máximo de sequências para comparar (amostragem).
            
        Returns
        -------
        bool
            True se similar acima do threshold, False caso contrário.
        """
        if not existing_records:
            return False
            
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        
        # Amostra aleatória para melhor performance
        import random
        sample_records = random.sample(existing_records, min(sample_size, len(existing_records)))
        
        for existing_rec in sample_records:
            try:
                alignment = aligner.align(record.seq, existing_rec.seq)
                best_score = alignment[0].score
                max_possible = max(len(record.seq), len(existing_rec.seq))
                similarity = best_score / max_possible
                
                if similarity >= self.similarity_threshold:
                    return True
            except Exception as e:
                self.logger.warning(f"Erro no alinhamento entre {record.id} e {existing_rec.id}: {e}")
                continue
                
        return False

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

    def run_workflow(self, outgroup_query_or_file="", query=None ,csv_path = None, 
                     download_method = "query", expected_accessions_file=None):
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
        expected_accessions_file : str, optional
            Arquivo com lista de accession numbers esperados do estudo original
        """
        self.logger.info("Iniciando o workflow completo...")
        raw_file = os.path.join(self.work_dir, "raw_sequences.gb")
        filtered_file = os.path.join(self.work_dir, "filtered_sequences.gb")
        no_utrs_file = os.path.join(self.work_dir, "no_utrs_sequences.gb")
        refined_file = os.path.join(self.work_dir, "refined_dataset.gb")
        dataset_outgroup_file = os.path.join(self.work_dir, "dataset_with_outgroup.gb")
        alignment_file = os.path.join(self.work_dir, "final_alignment.fasta")
        
        # Passo 1: Baixar sequências 
        if download_method == "query":
            self.download_sequences(query, raw_file)
        elif download_method == "csv":
            self.download_from_csv(csv_path, raw_file)
        
        # if expected_accessions_file:
        raw_stats = self.verify_downloaded_accessions(
            raw_file, 
            expected_file=expected_accessions_file,
            stage="raw"
        )
        
        if len(raw_stats['missing']) > 0:
            missing_pct = (len(raw_stats['missing']) / 
                        (len(raw_stats['missing']) + len(raw_stats['accessions']))) * 100
            self.logger.warning(f"Dataset incompleto: {missing_pct:.1f}% dos accessions esperados faltando")
        
        # Passo 2: Filtrar sequências
        filtered = self.filter_sequences(raw_file, filtered_file)
        if len(filtered) == 0:
            filtered_file = raw_file
            
        if expected_accessions_file and len(filtered) > 0:
            self.verify_downloaded_accessions(
                filtered_file, 
                expected_file=expected_accessions_file,
                stage="filtered"
            )
            
        # Passo 3: Remover UTRs (se parâmetros definidos)
        if self.utr5_end is not None and self.utr3_start is not None:
            self.remove_utrs(filtered_file, no_utrs_file)
        else:
            self.logger.info("Parâmetros UTR não definidos, pulando remoção")
            no_utrs_file = filtered_file  # Usa o mesmo arquivo
        
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
        # self.align_sequences(dataset_outgroup_file, alignment_file)
        
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
            
    def verify_downloaded_accessions(self, genbank_file, expected_file=None, stage="raw"):
        """
        Verifica os accession numbers baixados e gera relatório.
        
        Parameters
        ----------
        genbank_file : str
            Arquivo GenBank baixado
        expected_file : str, optional
            Arquivo com lista de accessions esperados (um por linha)
        stage : str
            Estágio do workflow (raw, filtered, refined)
        
        Returns
        -------
        dict
            Dicionário com estatísticas da verificação
        """
        self.logger.info(f"Verificando accessions em: {genbank_file} (estágio: {stage})")
        
        stats = {
            'total': 0,
            'accessions': [],
            'missing': [],
            'unexpected': [],
            'file_exists': os.path.exists(genbank_file),
            'file_size': os.path.getsize(genbank_file) if os.path.exists(genbank_file) else 0
        }
        
        try:
            if not stats['file_exists'] or stats['file_size'] == 0:
                self.logger.warning(f"Arquivo vazio ou não existe: {genbank_file}")
                return stats
            
            records = list(SeqIO.parse(genbank_file, "genbank"))
            stats['total'] = len(records)
            stats['accessions'] = [rec.id for rec in records]
            
            self.logger.info(f"Accessions encontrados: {stats['total']}")
            
            # Se tiver um arquivo de referência com os accessions esperados
            if expected_file and os.path.exists(expected_file):
                with open(expected_file, 'r') as f:
                    expected = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                
                stats['missing'] = sorted(set(expected) - set(stats['accessions']))
                stats['unexpected'] = sorted(set(stats['accessions']) - set(expected))
                
                if stats['missing']:
                    self.logger.warning(f"Accessions esperados NÃO encontrados ({len(stats['missing'])}):")
                    for acc in stats['missing'][:10]:  # Mostra só os 10 primeiros
                        self.logger.warning(f"  - {acc}")
                    if len(stats['missing']) > 10:
                        self.logger.warning(f"  ... e mais {len(stats['missing'])-10}")
                
                if stats['unexpected']:
                    self.logger.info(f"Accessions adicionais encontrados ({len(stats['unexpected'])}):")
                    for acc in stats['unexpected'][:10]:
                        self.logger.info(f"  + {acc}")
            
            # Salva relatório detalhado
            report_file = os.path.join(self.work_dir, f"accession_report_{stage}.txt")
            with open(report_file, 'w') as f:
                f.write(f"=== Relatório de Accessions - Estágio: {stage} ===\n")
                f.write(f"Arquivo: {genbank_file}\n")
                f.write(f"Total de sequências: {stats['total']}\n")
                f.write(f"Data da verificação: {pd.Timestamp.now()}\n\n")
                
                f.write("=== Accessions Encontrados ===\n")
                for acc in sorted(stats['accessions']):
                    f.write(f"{acc}\n")
                
                if stats['missing']:
                    f.write("\n=== Accessions Esperados NÃO Encontrados ===\n")
                    for acc in stats['missing']:
                        f.write(f"{acc}\n")
                
                if stats['unexpected']:
                    f.write("\n=== Accessions Adicionais Encontrados ===\n")
                    for acc in stats['unexpected']:
                        f.write(f"{acc}\n")
            
            self.logger.info(f"Relatório salvo em: {report_file}")
            return stats
            
        except Exception as e:
            self.logger.error(f"Erro na verificação: {e}")
            return stats
    
if __name__ == "__main__":
    # path = "workflow_dataAcquisition_SupplementaryTable_filtered_1"
    # workflow = workflowAquisitionDatasetNCBI(work_dir=path,
    #                         email="email@gmail.com",
    #                         utr5_end=True,  
    #                         utr3_start=True,  
    #                         initial_min_length=700,
    #                         refined_min_length=9000,
    #                         similarity_threshold=0.99)
    
    # QUERY: Primeira tentativa 270hits
    # zika_query = "Zika virus[Organism] AND complete genome"
    # outgroup_query = "Spondweni virus[Organism] AND complete genome"
    
    # QUERY: Segunda tentativa 128hits
    # zika_query = '"Zika virus"[Organism] OR Zika virus[All Fields]'
    # outgroup_query = 'spondweni[All Fields] AND ("Viruses"[Organism] OR viruses[All Fields])'
    
    #QUERY: Bactéria da gastrite
    # zika_query = 'Helicobacter pylori[Organism] AND complete genome'
    # outgroup_query = ""
    
    
    # workflow.run_workflow(csv_path="data/dataset_zikaVirus.csv", outgroup_query_or_file="", download_method="csv")
    
    
    # input_genbank = f"{path}/dataset_with_outgroup.gb"

    # input_genbank = f"{path}/raw_sequences.gb" # processamento do csv.
    # output_fasta = f"{path}/dataset_final.fasta"
    # workflow.generate_fasta(input_genbank, output_fasta)
    
    # input_fasta = f"{path}/dataset_final.fasta"
    # workflow.slice_file(input_fasta, output_prefix="dataset_slice", slice_size=50)
    
    # EXPERIMENTO - CORONAVIRUS ----------------------------------------------------------------------------
    # path = "workflow_dataAcquisition_coronavirus"
    # covid_workflow = workflowAquisitionDatasetNCBI(
    #     email="email@dominio.com",
    #     work_dir=path,
    #     initial_min_length=29000,    # Genomas completos
    #     refined_min_length=29500,    # Filtro mais rigoroso
    #     utr5_end=None,              # Manter UTRs para estudos de regulação
    #     utr3_start=None,
    #     similarity_threshold=0.999,  # Alta similaridade devido à conservação
    #     retmax=100                 # Muitas sequências disponíveis
    # )
    
    # covid_workflow.run_workflow(
    #     query='("Severe acute respiratory syndrome coronavirus 2"[Organism] AND complete genome) AND 2023[PDAT]',
    #     outgroup_query_or_file='"SARS coronavirus"[Organism]',
    #     download_method="query"
    # )
    
    # input_genbank = f"{path}/dataset_with_outgroup.gb"
    # output_fasta = f"{path}/dataset_final.fasta"
    # covid_workflow.generate_fasta(input_genbank, output_fasta)
    
    # # EXPERIMENTO - TUBERCULOSE ----------------------------------------------------------------------------
    # path = "workflow_dataAcquisition_tuberculosis"
    # tb_workflow = workflowAquisitionDatasetNCBI(
    #     email="email@dominio.com",
    #     work_dir=path,
    #     # initial_min_length=40000,   # Genoma bacteriano ~4.4Mb
    #     # refined_min_length=42000,   # Filtro para genomas mais completos
    #     utr5_end=None,               # Geralmente não se remove UTRs em bactérias
    #     utr3_start=None,
    #     similarity_threshold=0.98,    # Mais tolerante devido à diversidade
    #     retmax=100                  # Menos sequências completas disponíveis
    # )
    
    # tb_workflow.run_workflow(
    #     query='"Mycobacterium tuberculosis"[Organism] AND complete genome',
    #     outgroup_query_or_file='"Mycobacterium bovis"[Organism]', 
    #     download_method="query"
    # )
    
    # input_genbank = f"{path}/dataset_with_outgroup.gb"
    # output_fasta = f"{path}/dataset_final.fasta"
    # tb_workflow.generate_fasta(input_genbank, output_fasta)
    
    # EXPERIMENTO - TUBERCULOSE ----------------------------------------------------------------------------
    path = "workflow_dataAcquisition_li_et_al_2007_replication-RetMax100"
    workflow = workflowAquisitionDatasetNCBI(
        email="seu_email@dominio.com",  # Substitua pelo seu email
        work_dir=path,
        initial_min_length=180000,      # Genoma do VARV tem ~186kb, filtro inicial
        refined_min_length=183000,      # Para garantir genomas praticamente completos
        utr5_end=None,                   # Vírus não têm UTRs como eucariotos
        utr3_start=None,                  # Manter genoma completo
        similarity_threshold=0.999,      # Vírus têm alta similaridade (>99.6% entre isolados)
        retmax=100                        # Para garantir que pegue todos os disponíveis
    )
    
    # Query para as 47 amostras do estudo (ou o máximo disponível atualmente)
    variola_query = '''
        ("Variola virus"[Organism] OR "Variola virus"[All Fields]) 
        AND complete genome 
        AND 1000:200000[SLEN]
    '''

    # Query para os outgroups 
    outgroup_query = '''
        ("Taterapox virus"[Organism] OR "Taterapox virus"[All Fields] OR 
        "Camelpox virus"[Organism] OR "Camelpox virus"[All Fields])
        AND complete genome
    '''
    
    workflow.run_workflow(
        query=variola_query,
        outgroup_query_or_file=outgroup_query,  
        download_method="query"
    )
    
    input_genbank = f"{path}/dataset_with_outgroup.gb"
    output_fasta = f"{path}/dataset_final.fasta"
    workflow.generate_fasta(input_genbank, output_fasta)
    
    
