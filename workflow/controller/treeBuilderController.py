import os, time, sys
from pathlib import Path

from tqdm import tqdm
from Bio import Phylo
import matplotlib.pyplot as plt
from Bio import AlignIO
import numpy as np
import logging, datetime


current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../..'))

from workflow.tree_construction.builder import TreeBuilder
from workflow.utils.dataValidation import duplicate_names, duplicate_seq, validate_sequences, remove_pipe
from workflow.utils.dataCleaning import clean_NoPipe, clean_tmp, copiar_arquivos
from workflow.utils.messages import Messages
from workflow.utils.metrics import process_rf_distance, plot_heatmap_distances
from workflow.alignment.alignmentSeq import AlignmentSeqs

class TreeBuilderController:
    """
    Controlador responsável pela construção e manipulação de árvores filogenéticas.

    Esta classe coordena o processo de construção de árvores usando métodos baseados em distâncias
    e parcimônia, integrando funcionalidades de validação, alinhamento e visualização de árvores.

    Attributes
    ----------
    msg : Messages
        Instância da classe Messages para exibir informações no console.
    start : float
        Tempo de início da execução do processo de construção de árvores.
    files : list
        Lista de arquivos de entrada que serão processados.
    count_trees : int
        Contador de árvores construídas.
    count_nodes : list
        Lista de contagens de nós terminais em cada árvore construída.
    list_times : list
        Lista de tempos de execução para cada ciclo de construção de árvore.

    """

    def __init__(self, **kwargs) -> None:
        """
        Inicializa a instância da classe TreeBuilderController.

        Configura os atributos a partir dos argumentos fornecidos via kwargs, cria os diretórios
        de saída necessários e executa as limpezas de diretórios temporários.

        Parameters
        ----------
        **kwargs : dict
            Dicionário contendo pares chave-valor que serão configurados como atributos da instância.

        Return
        ------
        None
        """
        for k, v in kwargs.items():
            setattr(self, k, v)
                        
        self.msg = Messages(logPath=os.path.join(self.output_path,'outputs'))
        self.start = time.time()
        
        print(self.msg.init_message())
        print('               - CONSTRUÇÃO DE ÁRVORES -               \n')
        print('------------------------------------------------------')
        
        os.makedirs(self.output_path, exist_ok=True)
        dirs = ['outputs', 'tmp','Align', 'Trees', 'outputs/Plots']

        for dir in dirs:
            path = os.path.join(self.output_path, dir)
            os.makedirs(path, exist_ok=True)
        
        self.files = sorted(os.listdir(self.input_path))
        self.count_trees = 0
        self.count_nodes = list()
        self.list_times = list()
        self.aligner = AlignmentSeqs({'num_threads': self.num_threads})

        # clean_tmp(self.output_path)
        clean_NoPipe(self.input_path)

        date = datetime.datetime.now()

        logging.basicConfig(level=logging.INFO, 
                    filename=os.path.join(self.output_path,'outputs',f"log_setup_{date.year}_{date.month}_{date.day}.log"),
                    format='%(asctime)s - %(levelname)s - %(message)s')

    def __call__(self):
        """
        Executa o processo de construção das árvores baseado no modo especificado.

        Esta função coordena o fluxo de trabalho para construir árvores filogenéticas utilizando
        os métodos de construção por matriz de distâncias (NJ/UPGMA) ou parcimônia, dependendo do modo selecionado.

        Return
        ------
        None
        """
        if self.mode == "distance":
            print("Iniciando construção das árvores utilizando o método: DISTANCE TREE CONSTRUCTOR + CLUSTALW\n")
            logging.info("Iniciando construção das árvores utilizando o método: DISTANCE TREE CONSTRUCTOR + CLUSTALW")
            logging.info("STEP: construction of trees using the method: DISTANCE TREE CONSTRUCTOR + CLUSTALW")
        elif self.mode == "parsimony":
            print("Iniciando construção das árvores utilizando o método: PARSIMONY + CLUSTALW\n")
            logging.info("Iniciando construção das árvores utilizando o método: PARSIMONY + CLUSTALW") 
            logging.info("STEP: construction of trees using the method: PARSIMONY + CLUSTALW")                           
        elif self.mode == "auto":
            print("Iniciando construção das árvores utilizando ambos os métodos: DISTANCE TREE CONSTRUCTOR e PARSIMONY\n")
            logging.info("Iniciando construção das árvores utilizando ambos os métodos: DISTANCE TREE CONSTRUCTOR e PARSIMONY")
            logging.info("STEP: construction of trees using both methods: DISTANCE TREE CONSTRUCTOR and PARSIMONY")
        elif self.mode == "advanced":
            print("Iniciando construção das árvores utilizando métodos avançados\n")
            logging.info("Iniciando construção das árvores utilizando métodos avançados")
        else:
            logging.error(f"Modo desconhecido: {self.mode}")
            raise ValueError(f"Modo desconhecido: {self.mode}")
        
        
        print(f"Ignorando o método: {self.ignore_mode.upper()}\n")
        logging.info(f"Ignorando o método: {self.ignore_mode.upper()}")
        heatmap_matrix = np.zeros((8, 8))

        multi_trees = {
            "clustalo": {
                "distance": {"nj": [], "upgma": []}, 
                "parsimony": {"nj": [], "upgma": []},
                "iqtree": [],
                "fasttree": [],
                "raxml": [],
                "mrbayes": []
            },
            "mafft": {
                "distance": {"nj": [], "upgma": []}, 
                "parsimony": {"nj": [], "upgma": []},
                "iqtree": [],
                "fasttree": [],
                "raxml": [],
                "mrbayes": []
            }
        }
        
        base_folder = self.input_path.split('/')[-1]

        for file in tqdm(self.files, desc="Construindo árvores...", ascii="░▒█"):
            if '.dnd' in file or '.json' in file:
                logging.debug(f"Arquivo {file} ignorado por ser .dnd/.json")
                continue
            
            logging.info(f"Iniciando processamento do arquivo: {os.path.join(base_folder,file)}")
            
            start_cycle = time.time()
            fasta_path = os.path.join(self.input_path, file)
            output_path_align_base = os.path.join(self.output_path, 'tmp', f'{Path(file).stem}.aln')
            output_path_align_html = os.path.join(self.output_path, 'Align', f'{Path(file).stem}.html')
            output_path_dnd = os.path.join(self.output_path, 'tmp', f'{Path(file).stem}.dnd')
            path_dnd = os.path.join(self.input_path, f'{Path(file).stem}.dnd')
            
            if not(duplicate_names(fasta_path)) and not(duplicate_seq(fasta_path)[0]) and validate_sequences(fasta_path):
                logging.info(f"Arquivo {file} passou nas validações de sequência.")
            else:
                logging.warning(f"Arquivo {file} contém duplicatas ou erros, iniciando remoção de pipes.")
                fasta_path = remove_pipe(Path(file).stem, fasta_path, self.input_path)
                path_dnd = os.path.join(f'{fasta_path}.dnd')
            
            output_path_tree = os.path.join(self.output_path, 'Trees')
            output_path_tree_image = os.path.join(self.output_path, 'Trees')
            
            try:
                if self.mode.lower() == "distance":
                    self.count_trees += 1
                    name = f'tree_{Path(file).stem}_distance.{self.output_format}'
                    logging.info(f"Construindo árvore de distância para o arquivo {file}")
                    tree = self.build_tree_distance_matrix(fasta_path, output_path_align_base, output_path_dnd, path_dnd, os.path.join(output_path_tree, name), self.align_method, output_path_align_html)
                    self.save_tree_image(title=name, tree=[tree], path=os.path.join(output_path_tree_image, name).replace('Trees', 'outputs/Plots'))
                    end_time = time.time()
                elif self.mode.lower() == "parsimony":
                    self.count_trees += 1
                    name = f'tree_{Path(file).stem}_parsimony.{self.output_format}'
                    logging.info(f"Construindo árvore por parcimônia para o arquivo {file}")
                    tree = self.build_tree_parsimony(fasta_path, output_path_align_base, output_path_dnd, path_dnd, os.path.join(output_path_tree, name), self.align_method, output_path_align_html)
                    self.save_tree_image(title=name, tree=[tree], path=os.path.join(output_path_tree_image, name).replace('Trees', 'outputs/Plots'))
                    end_time = time.time()
                elif self.mode.lower() == "auto":
                    
                    logging.info(f"Construindo árvores com múltiplos métodos para o arquivo {file}")
                    for method in ['nj', 'upgma']:
                        for alg in ['clustalo', 'mafft']:
                            self.count_trees += 2
                            name_distance = f'tree_{Path(file).stem}_{alg}_{method}_distance.{self.output_format}'
                            name_parsimony = f'tree_{Path(file).stem}_{alg}_{method}_parsimony.{self.output_format}'
                            output_path_align = output_path_align_base.replace('.aln',f'_{alg}.aln')

                            output_path_tree_distance = os.path.join(self.output_path, 'Trees', name_distance)
                            output_path_tree_parsimony = os.path.join(self.output_path, 'Trees', name_parsimony)

                            if os.path.exists(output_path_tree_distance) and os.path.exists(output_path_tree_parsimony):
                                
                                tree_distance = Phylo.read(output_path_tree_distance, self.output_format)
                                tree_parsimony = Phylo.read(output_path_tree_parsimony, self.output_format)
                                multi_trees[alg]['distance'][method].append(tree_distance)
                                multi_trees[alg]['parsimony'][method].append(tree_parsimony)
                                continue

                            self.construct_tree_method = method

                            if self.ignore_mode.lower() and multi_trees.get(self.ignore_mode.lower()):
                                multi_trees[alg].pop(self.ignore_mode.lower())
                            
                            if not self.ignore_mode.lower() == "distance" and \
                                not os.path.exists(output_path_tree_distance) and \
                                not self.ignore_mode.lower() == "distance":
                                    
                                logging.debug(f"Construindo árvore de distância ({alg} - {method}) para o arquivo {file}")
                                tree_distance = self.build_tree_distance_matrix(fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree_distance, alg, output_path_align_html)
                                self.save_tree_image(title=name_distance, tree=[tree_distance], path=os.path.join(output_path_tree_image, name_distance).replace('Trees', 'outputs/Plots'))
                                multi_trees[alg]['distance'][method].append(tree_distance)
                            
                            if not self.ignore_mode.lower() == "parsimony" and\
                                not os.path.exists(output_path_tree_parsimony) and\
                                not self.ignore_mode.lower() == "parsimony":
                                        
                                logging.debug(f"Construindo árvore por parcimônia ({alg} - {method}) para o arquivo {file}")
                                tree_parsimony = self.build_tree_parsimony(fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree_parsimony, alg, output_path_align_html)
                                self.save_tree_image(title=name_parsimony, tree=[tree_parsimony], path=os.path.join(output_path_tree_image, name_parsimony).replace('Trees', 'outputs/Plots'))
                                multi_trees[alg]['parsimony'][method].append(tree_parsimony)
                        
                    # self.save_tree_image(title=f'tree_{Path(file).stem}', tree=multi_trees, path=os.path.join(output_path_tree_image, f'tree_{Path(file).stem}').replace('Trees', 'outputs/Plots'))
                    end_time = time.time()
                    rf_scores = process_rf_distance(multi_trees)
                    heatmap_matrix = self.somarMatrizes(
                        heatmap_matrix,
                        plot_heatmap_distances(data_dict=multi_trees, scores=rf_scores, base_name=Path(file).stem, path=os.path.join(self.output_path, 'outputs/Plots'))
                    )
                elif self.mode == "advanced":
                    for alg in ['clustalo', 'mafft']:
                        for method in ['nj', 'upgma']:
                            name_distance = f'tree_{Path(file).stem}_{alg}_{method}_distance.{self.output_format}'
                            name_parsimony = f'tree_{Path(file).stem}_{alg}_{method}_parsimony.{self.output_format}'
                            name_iqtree = f'tree_{Path(file).stem}_{alg}_iqtree.{self.output_format}'
                            name_fasttree = f'tree_{Path(file).stem}_{alg}_fasttree.{self.output_format}'
                            name_raxml = f'tree_{Path(file).stem}_{alg}_raxml.{self.output_format}'
                            name_mrbayes = f'tree_{Path(file).stem}_{alg}_mrbayes.{self.output_format}'
                           
                            output_path_align = output_path_align_base.replace('.aln',f'_{alg}.aln')

                            output_path_tree_distance = os.path.join(self.output_path, 'Trees', name_distance)
                            output_path_tree_parsimony = os.path.join(self.output_path, 'Trees', name_parsimony)
                            output_path_tree_iqtree = os.path.join(self.output_path, 'Trees', name_iqtree)
                            output_path_tree_fasttree = os.path.join(self.output_path, 'Trees', name_fasttree)
                            output_path_tree_raxml = os.path.join(self.output_path, 'Trees', name_raxml)
                            output_path_tree_mrbayes = os.path.join(self.output_path, 'Trees', name_mrbayes)

                            if os.path.exists(output_path_tree_distance):

                                tree_distance = Phylo.read(output_path_tree_distance, self.output_format)
                                tree_parsimony = Phylo.read(output_path_tree_parsimony, self.output_format)
                                tree_iqtree = Phylo.read(output_path_tree_iqtree, self.output_format)
                                tree_fasttree = Phylo.read(output_path_tree_fasttree, self.output_format)
                                tree_raxml = Phylo.read(output_path_tree_raxml, self.output_format)
                                tree_mrbayes = Phylo.read(output_path_tree_mrbayes, self.output_format)
                                multi_trees[alg]['distance'][method].append(tree_distance)
                                multi_trees[alg]['parsimony'][method].append(tree_parsimony)
                                multi_trees[alg]['iqtree'].append(tree_iqtree)
                                multi_trees[alg]['fasttree'].append(tree_fasttree)
                                multi_trees[alg]['raxml'].append(tree_raxml)
                                multi_trees[alg]['mrbayes'].append(tree_mrbayes)
                                self.count_trees = 12
                                continue

                            self.construct_tree_method = method

                            if self.ignore_mode.lower() and multi_trees.get(self.ignore_mode.lower()):
                                multi_trees[alg].pop(self.ignore_mode.lower())
                            
                            if not self.ignore_mode.lower() == "distance" and \
                                not os.path.exists(output_path_tree_distance) and \
                                not self.ignore_mode.lower() == "distance":
                                    
                                self.count_trees += 1
                                logging.debug(f"Construindo árvore de distância ({alg} - {method}) para o arquivo {file}")
                                tree_distance = self.build_tree_distance_matrix(fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree_distance, alg, output_path_align_html)
                                self.save_tree_image(title=name_distance, tree=[tree_distance], path=os.path.join(output_path_tree_image, name_distance).replace('Trees', 'outputs/Plots'))
                                multi_trees[alg]['distance'][method].append(tree_distance)
                            
                            if not self.ignore_mode.lower() == "parsimony" and\
                                not os.path.exists(output_path_tree_parsimony) and\
                                not self.ignore_mode.lower() == "parsimony":
                                        
                                self.count_trees += 1
                                logging.debug(f"Construindo árvore por parcimônia ({alg} - {method}) para o arquivo {file}")
                                tree_parsimony = self.build_tree_parsimony(fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree_parsimony, alg, output_path_align_html)
                                self.save_tree_image(title=name_parsimony, tree=[tree_parsimony], path=os.path.join(output_path_tree_image, name_parsimony).replace('Trees', 'outputs/Plots'))
                                multi_trees[alg]['parsimony'][method].append(tree_parsimony)
                            
                            if not self.ignore_mode.lower() == "iqtree" and \
                                not os.path.exists(output_path_tree_iqtree) and \
                                not self.ignore_mode.lower() == "iqtree":

                                self.count_trees += 1
                                logging.debug(f"Construindo árvore por iqtree  ({alg} - {method}) para o arquivo {file}")
                                tree_iqtree = self.build_tree_iqtree(fasta_path, output_path_align, output_path_tree_iqtree, alg, output_path_align_html)
                                self.save_tree_image(title=name_iqtree, tree=[tree_iqtree], path=os.path.join(output_path_tree_image, name_iqtree).replace('Trees', 'outputs/Plots'))

                                multi_trees[alg]['iqtree'].append(tree_iqtree)
                            
                            if not self.ignore_mode.lower() == "fasttree" and \
                                not os.path.exists(output_path_tree_fasttree) and \
                                not self.ignore_mode.lower() == "fasttree":

                                self.count_trees += 1
                                logging.debug(f"Construindo árvore por fasttree  ({alg} - {method}) para o arquivo {file}")
                                tree_fasttree = self.build_tree_fasttree(fasta_path, output_path_align, output_path_tree_fasttree, alg, output_path_align_html)
                                self.save_tree_image(title=name_fasttree, tree=[tree_fasttree], path=os.path.join(output_path_tree_image, name_fasttree).replace('Trees', 'outputs/Plots'))

                                multi_trees[alg]['fasttree'].append(tree_fasttree)
                            
                            if not self.ignore_mode.lower() == "raxml" and \
                                not os.path.exists(output_path_tree_raxml) and \
                                not self.ignore_mode.lower() == "raxml":

                                self.count_trees += 1
                                logging.debug(f"Construindo árvore por raxml  ({alg} - {method}) para o arquivo {file}")
                                tree_raxml = self.build_tree_raxml(fasta_path, output_path_align, output_path_tree_raxml, alg, output_path_align_html)
                                self.save_tree_image(title=name_raxml, tree=[tree_raxml], path=os.path.join(output_path_tree_image, name_raxml).replace('Trees', 'outputs/Plots'))

                                multi_trees[alg]['raxml'].append(tree_raxml)
                            
                            if not self.ignore_mode.lower() == "mrbayes" and \
                                not os.path.exists(output_path_tree_mrbayes) and \
                                not self.ignore_mode.lower() == "mrbayes":

                                self.count_trees += 1
                                logging.debug(f"Construindo árvore por mrbayes  ({alg} - {method}) para o arquivo {file}")
                                tree_mrbayes = self.build_tree_mrbayes(fasta_path, output_path_align, output_path_tree_mrbayes, alg, output_path_align_html)
                                self.save_tree_image(title=name_mrbayes, tree=[tree_mrbayes], path=os.path.join(output_path_tree_image, name_mrbayes).replace('Trees', 'outputs/Plots'))

                                multi_trees[alg]['mrbayes'].append(tree_mrbayes)
                                
                    end_time = time.time()
                    rf_scores = process_rf_distance(multi_trees)
                    heatmap_matrix = self.somarMatrizes(
                        heatmap_matrix,
                        plot_heatmap_distances(data_dict=multi_trees, scores=rf_scores, base_name=Path(file).stem, path=os.path.join(self.output_path, 'outputs/Plots'))
                    )
                else:
                    logging.error(f"Modo não suportado: {self.mode}")
                    continue
                cycle_time = end_time - start_cycle
                self.list_times.append(cycle_time)
                logging.info(f"Arquivo {file} processado com sucesso em {cycle_time:.2f} segundos.")
            except Exception as e:
                logging.error(f"Erro ao processar o arquivo {file}: {e}", exc_info=True)
                exit(1)
                
        if self.mode == "auto" :
            plot_heatmap_distances(data_dict=multi_trees, base_name='Acumulate', 
                               path=os.path.join(self.output_path, 'outputs/Plots'), distance_matrix=heatmap_matrix)
            logging.info("Heatmap acumulado gerado com sucesso.")
        
        clean_NoPipe(self.input_path)
        copiar_arquivos(os.path.join(self.output_path, 'tmp'),os.path.join(self.output_path, 'Align'))
        # clean_tmp(self.output_path)
        logging.info("Diretórios temporários limpos.")

        self.msg.resume_tree(start=self.start, sum_time=self.list_times, num_trees=self.count_trees, 
                             output_format=self.output_format, n_nodes=self.count_nodes, method=self.construct_tree_method)
        logging.info("Processo de construção de árvores finalizado.")
        
    def somarMatrizes(self, matriz1, matriz2):
        """
        Soma duas matrizes elemento por elemento.

        Esta função realiza a soma de duas matrizes de distância, elemento por elemento,
        verificando se as dimensões das matrizes são compatíveis.

        Parameters
        ----------
        matriz1 : np.ndarray
            A primeira matriz de distância.
        matriz2 : np.ndarray
            A segunda matriz de distância.

        Return
        ------
        np.ndarray or None
            A matriz resultante da soma ou None se as dimensões forem incompatíveis.
        """
        if len(matriz1) != len(matriz2) or len(matriz1[0]) != len(matriz2[0]):
            return None
        result = []
        for i in range(len(matriz1)):   
            result.append([])
            for j in range(len(matriz1[0])):
                result[i].append(round(matriz1[i][j] + matriz2[i][j]))
        return result

    def build_tree_distance_matrix(self, fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree, align_method, output_path_align_html):
        """
        Constrói uma árvore filogenética usando matriz de distâncias.

        Alinha as sequências, calcula a matriz de distâncias e constrói uma árvore filogenética
        usando o método de Neighbor-Joining ou UPGMA, dependendo do método especificado.

        Parameters
        ----------
        fasta_path : str
            Caminho para o arquivo FASTA contendo as sequências.
        output_path_align : str
            Caminho para salvar o arquivo de alinhamento gerado.
        output_path_dnd : str
            Caminho para salvar o arquivo .dnd gerado.
        path_dnd : str
            Caminho do arquivo .dnd original.
        output_path_tree : str
            Caminho para salvar a árvore gerada.
        align_method : str
            Método de alinhamento a ser usado: 'clustalw' ou 'mafft'.

        Return
        ------
        Bio.Phylo.BaseTree.Tree
            A árvore filogenética construída.
        """
        logging.info(f"Iniciando construção de árvore por matriz de distância para {fasta_path}")
        logging.info(f"STEP: Construction of distance matrix.")
        builder = TreeBuilder(fasta_path=fasta_path, 
                              output_path_align=output_path_align, 
                              output_path_dnd=output_path_dnd, 
                              path_dnd=path_dnd, 
                              output_path_tree=output_path_tree)
        
        try:
            if os.path.exists(output_path_align):
                logging.info(f"Arquivo de alinhamento já existe: {output_path_align}. Reutilizando.")
                # with open(output_path_align, "rb") as f:
                #     content = f.read().decode("utf-8-sig")
                # tmp_path = os.path.dirname(os.path.abspath(output_path_align))
                # tmp_alg_path = os.path.join(tmp_path,"temp.aln")
                # with open(tmp_alg_path, "w") as f:
                #     f.write(content)
                alng = AlignIO.read(output_path_align, "fasta")
                
                
                
                # alng = AlignIO.read(output_path_align, "clustal" if align_method == "clustalw" else "fasta")
            else:
                if align_method == "clustalo":
                    logging.debug(f"Alinhando sequências com Clustalo para {fasta_path}.")
                    alng = self.aligner.align_sequences_clustalo(
                        fasta_path=fasta_path,
                        # path_dnd=path_dnd,
                        # output_path_dnd=output_path_dnd,
                        output_path_align=output_path_align,
                        output_path_html=output_path_align_html
                    )

                elif align_method == "mafft":
                    logging.debug(f"Alinhando sequências com MAFFT para {fasta_path}.")
                    alng = self.aligner.align_sequences_mafft(
                        fasta_path=fasta_path,
                        output_path_align=output_path_align,
                        output_path_html=output_path_align_html
                    )
                    
                    logging.debug("Arquivo de alinhamento gerado e lido com sucesso (MAFFT).")
                else:
                    logging.error(f"Método de alinhamento desconhecido: {align_method}")
                    raise ValueError("Método de alinhamento não suportado.")
        except Exception as e:
            logging.error(f"Erro no alinhamento das sequências para {fasta_path}: {e}", exc_info=True)
            raise

        try:
            distance_matrix = builder.distance_matrix(alignment=alng)
            logging.debug("Matriz de distâncias calculada com sucesso.")
            tree = builder.distance_constructor(distance_matrix=distance_matrix, 
                                                construct_tree_method=self.construct_tree_method)
            nodes = tree.count_terminals()
            self.count_nodes.append(nodes)
            logging.info(f"Árvore construída com {nodes} nós terminais.")
            builder.save_tree(tree=tree, path=output_path_tree, format=self.output_format)
            logging.info(f"Árvore salva em {output_path_tree}.")
        except Exception as e:
            logging.error(f"Erro na construção da árvore por matriz de distância para {fasta_path}: {e}", exc_info=True)
            raise
        return tree
    
    def build_tree_parsimony(self, fasta_path, output_path_align, output_path_dnd,
                         path_dnd, output_path_tree, align_method, output_path_align_html):
        """
        Constrói uma árvore filogenética usando o método de parcimônia.
        Se o alinhamento já existir, ele será reutilizado.

        Parameters
        ----------
        fasta_path : str
            Caminho para o arquivo FASTA contendo as sequências.
        output_path_align : str
            Caminho para salvar o arquivo de alinhamento gerado.
        output_path_dnd : str
            Caminho para salvar o arquivo .dnd gerado.
        path_dnd : str
            Caminho do arquivo .dnd original.
        output_path_tree : str
            Caminho para salvar a árvore gerada.
        align_method : str
            Método de alinhamento a ser usado: 'clustalw' ou 'mafft'.
        output_path_align_html : str
            Caminho para salvar o alinhamento em HTML.

        Return
        ------
        Bio.Phylo.BaseTree.Tree
            A árvore filogenética construída.
        """
        logging.info(f"Iniciando construção de árvore por parcimônia para {fasta_path}")
        logging.info(f"STEP: Tree Construction with parsimony method.")
        builder = TreeBuilder(fasta_path=fasta_path, 
                            output_path_align=output_path_align, 
                            output_path_dnd=output_path_dnd, 
                            path_dnd=path_dnd, 
                            output_path_tree=output_path_tree)
        
        try:
            if os.path.exists(output_path_align):
                logging.info(f"Arquivo de alinhamento já existe: {output_path_align}. Reutilizando.")
                # with open(output_path_align, "rb") as f:
                #     content = f.read().decode("utf-8-sig")
                # tmp_path = os.path.dirname(os.path.abspath(output_path_align))
                # tmp_alg_path = os.path.join(tmp_path,"temp.aln")
                # with open(tmp_alg_path, "w") as f:
                #     f.write(content)
                alng = AlignIO.read(output_path_align, "fasta")
                
                
                
                # alng = AlignIO.read(output_path_align, "clustal" if align_method == "clustalw" else "fasta")
            else:
                if align_method == "clustalo":
                    logging.debug(f"Alinhando sequências com Clustalo para {fasta_path}.")
                    alng = self.aligner.align_sequences_clustalo(
                        fasta_path=fasta_path,
                        # path_dnd=path_dnd,
                        # output_path_dnd=output_path_dnd,
                        output_path_align=output_path_align,
                        output_path_html=output_path_align_html
                    )
                    # alng = AlignIO.read(output_path_align, "clustal")

                elif align_method == "mafft":
                    logging.debug(f"Alinhando sequências com MAFFT para {fasta_path}.")
                    alng = self.aligner.align_sequences_mafft(
                        fasta_path=fasta_path,
                        output_path_align=output_path_align,
                        output_path_html=output_path_align_html
                    )
                    
                    logging.debug("Arquivo de alinhamento gerado e lido com sucesso (MAFFT).")
                else:
                    logging.error(f"Método de alinhamento desconhecido: {align_method}")
                    raise ValueError("Método de alinhamento não suportado.")
        except Exception as e:
            logging.error(f"Erro no alinhamento das sequências para {fasta_path}: {e}", exc_info=True)
            raise

        try:
            distance_matrix = builder.distance_matrix(alignment=alng)
            logging.debug("Matriz de distâncias calculada para parcimônia com sucesso.")
            tree = builder.parsimony_constructor(
                distance_matrix=distance_matrix,
                alignment=alng,
                construct_tree_method=self.construct_tree_method
            )
            builder.save_tree(tree=tree, path=output_path_tree, format=self.output_format)
            logging.info(f"Árvore de parcimônia salva em {output_path_tree}.")
        except Exception as e:
            logging.error(f"Erro na construção da árvore por parcimônia para {fasta_path}: {e}", exc_info=True)
            raise
        return tree
    
    def build_tree_iqtree(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando IQ-TREE."""
        logging.info(f"Iniciando construção de árvore com IQ-TREE para {fasta_path}")
        logging.info(f"STEP: Tree Construction with IQ-TREE method.")

        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree)
        
        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        tree = builder.iqtree_constructor(alng, output_path_tree)
        self.count_nodes.append(tree.count_terminals())
        return tree
    
    def build_tree_fasttree(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando FastTree."""
        logging.info(f"Iniciando construção de árvore com FastTree para {fasta_path}")
        logging.info(f"STEP: Tree Construction with FastTree method.")

        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree)
        
        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        tree = builder.fasttree_constructor(alng, output_path_tree)
        self.count_nodes.append(tree.count_terminals())
        return tree
    
    def build_tree_raxml(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando RAxML-NG."""
        logging.info(f"Iniciando construção de árvore com RAxML-NG para {fasta_path}")
        logging.info(f"STEP: Tree Construction with RAxML-NG method.")

        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree)
        
        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        tree = builder.raxml_ng_constructor(alng, output_path_tree)
        self.count_nodes.append(tree.count_terminals())
        return tree
    
    def build_tree_mrbayes(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando MrBayes."""
        logging.info(f"Iniciando construção de árvore com MrBayes para {fasta_path}")
        logging.info(f"STEP: Tree Construction with MrBayes method.")

        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree)
        
        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        tree = builder.mrbayes_constructor(alng, output_path_tree)
        self.count_nodes.append(tree.count_terminals())
        return tree
    
    def _get_alignment(self, fasta_path, output_path_align, align_method, output_path_align_html):
        """Método auxiliar para obter alinhamento."""
        try:
            if os.path.exists(output_path_align):
                logging.info(f"Arquivo de alinhamento já existe: {output_path_align}. Reutilizando.")
                return AlignIO.read(output_path_align, "fasta")
            else:
                if align_method == "clustalo":
                    return self.aligner.align_sequences_clustalo(
                        fasta_path=fasta_path,
                        output_path_align=output_path_align,
                        output_path_html=output_path_align_html
                    )
                elif align_method == "mafft":
                    return self.aligner.align_sequences_mafft(
                        fasta_path=fasta_path,
                        output_path_align=output_path_align,
                        output_path_html=output_path_align_html
                    )
        except Exception as e:
            logging.error(f"Erro ao obter alinhamento: {e}")
            raise

    def save_tree_image(self, title, tree, path):
        """
        Salva uma imagem da árvore filogenética em formato gráfico.

        Salva a imagem de uma árvore filogenética (ou múltiplas árvores) em um arquivo PNG,
        ajustando o layout e o título de acordo com os parâmetros fornecidos.

        Parameters
        ----------
        title : str
            O título da imagem.
        tree : Bio.Phylo.BaseTree.Tree or dict
            A árvore filogenética ou um dicionário contendo múltiplas árvores.
        path : str
            Caminho do arquivo onde a imagem será salva.

        Return
        ------
        None
        """
        logging.info(f"Starting to save the tree image: {title}")
        logging.info(f"STEP: saving the tree image.")
        try:
            if isinstance(tree, list):
                fig = plt.figure(figsize=(21, 10))
                ax = fig.add_subplot(1, 1, 1)
                Phylo.draw(tree[0], do_show=False, axes=ax, 
                           label_func=lambda x: None if x.name is None or ('Inner' in x.name if x.name else False) else x.name)
                plt.title(title)
                plt.tight_layout()
                img_path = path.replace(self.output_format, 'png')
                fig.savefig(img_path)
                plt.close(fig)
                logging.info(f"Imagem da árvore {title} salva em {img_path}.")
            else:
                fig, axs = plt.subplots(4, 4, figsize=(42, 36))  
                axs = axs.flatten()
                ax_index = 0

                for alg, methods in tree.items():
                    for method_name, trees in methods.items():
                        if isinstance(trees, dict): 
                            for sub_method, tree_list in trees.items():
                                if ax_index < len(axs) and len(tree_list) > 0:
                                    Phylo.draw(tree_list[0], do_show=False, axes=axs[ax_index])
                                    axs[ax_index].set_title(f"{alg} - {method_name} - {sub_method}", fontsize=8)
                                    ax_index += 1
                        else: 
                            if ax_index < len(axs) and len(trees) > 0:
                                Phylo.draw(trees[0], do_show=False, axes=axs[ax_index])
                                axs[ax_index].set_title(f"{alg} - {method_name}", fontsize=8)
                                ax_index += 1

                plt.tight_layout()
                img_path = path.replace(self.output_format, 'png')
                fig.savefig(img_path)
                plt.close(fig)
                logging.info(f"Imagem composta das árvores {title} salva em {img_path}.")
        except Exception as e:
            logging.error(f"Erro ao salvar a imagem da árvore {title}: {e}", exc_info=True)