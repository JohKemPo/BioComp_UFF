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
from workflow.utils.dataCleaning import clean_NoPipe, clean_tmp
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
        dirs = ['outputs', 'tmp', 'Trees', 'outputs/Plots']

        for dir in dirs:
            path = os.path.join(self.output_path, dir)
            os.makedirs(path, exist_ok=True)
        
        self.files = sorted(os.listdir(self.input_path))
        self.count_trees = 0
        self.count_nodes = list()
        self.list_times = list()
        self.aligner = AlignmentSeqs()

        clean_tmp(self.output_path)
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
        elif self.mode == "parsimony":
            print("Iniciando construção das árvores utilizando o método: PARSIMONY + CLUSTALW\n")
            logging.info("Iniciando construção das árvores utilizando o método: PARSIMONY + CLUSTALW")
        elif self.mode == "auto":
            print("Iniciando construção das árvores utilizando ambos os métodos: DISTANCE TREE CONSTRUCTOR e PARSIMONY\n")
            logging.info("Iniciando construção das árvores utilizando ambos os métodos: DISTANCE TREE CONSTRUCTOR e PARSIMONY")
        else:
            logging.error(f"Modo desconhecido: {self.mode}")
            raise ValueError(f"Modo desconhecido: {self.mode}")

        heatmap_matrix = np.zeros((8, 8))
        for file in tqdm(self.files, desc="Construindo árvores...", ascii="░▒█"):
            logging.info(f"Iniciando processamento do arquivo: {file}")
            if '.dnd' in file:
                logging.debug(f"Arquivo {file} ignorado por ser .dnd")
                continue
            
            start_cycle = time.time()
            fasta_path = os.path.join(self.input_path, file)
            output_path_align = os.path.join(self.output_path, 'tmp', f'{Path(file).stem}.aln')
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
                if self.mode == "distance":
                    self.count_trees += 1
                    name = f'tree_{Path(file).stem}_distance.{self.output_format}'
                    logging.info(f"Construindo árvore de distância para o arquivo {file}")
                    tree = self.build_tree_distance_matrix(fasta_path, output_path_align, output_path_dnd, path_dnd, os.path.join(output_path_tree, name), self.align_method)
                    self.save_tree_image(title=name, tree=[tree], path=os.path.join(output_path_tree_image, name).replace('Trees', 'outputs/Plots'))
                    end_time = time.time()
                elif self.mode == "parsimony":
                    self.count_trees += 1
                    name = f'tree_{Path(file).stem}_parsimony.{self.output_format}'
                    logging.info(f"Construindo árvore por parcimônia para o arquivo {file}")
                    tree = self.build_tree_parsimony(fasta_path, output_path_align, output_path_dnd, path_dnd, os.path.join(output_path_tree, name), self.align_method)
                    self.save_tree_image(title=name, tree=[tree], path=os.path.join(output_path_tree_image, name).replace('Trees', 'outputs/Plots'))
                    end_time = time.time()
                elif self.mode == "auto":
                    multi_trees = {
                        "clustalw": {"distance": {"nj": [], "upgma": []}, "parsimony": {"nj": [], "upgma": []}},
                        "mafft": {"distance": {"nj": [], "upgma": []}, "parsimony": {"nj": [], "upgma": []}}
                    }
                    logging.info(f"Construindo árvores com múltiplos métodos para o arquivo {file}")
                    for method in ['nj', 'upgma']:
                        for alg in ['clustalw', 'mafft']:
                            self.count_trees += 2
                            name_distance = f'tree_{Path(file).stem}_{alg}_{method}_distance.{self.output_format}'
                            name_parsimony = f'tree_{Path(file).stem}_{alg}_{method}_parsimony.{self.output_format}'
                            output_path_tree_distance = os.path.join(self.output_path, 'Trees', name_distance)
                            output_path_tree_parsimony = os.path.join(self.output_path, 'Trees', name_parsimony)
                            self.construct_tree_method = method
                            logging.debug(f"Construindo árvore de distância ({alg} - {method}) para o arquivo {file}")
                            tree_distance = self.build_tree_distance_matrix(fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree_distance, alg)
                            self.save_tree_image(title=name_distance.split('.')[0], tree=[tree_distance], path=os.path.join(output_path_tree_image, name).replace('Trees', 'outputs/Plots'))
                            logging.debug(f"Construindo árvore por parcimônia ({alg} - {method}) para o arquivo {file}")
                            tree_parsimony = self.build_tree_parsimony(fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree_parsimony, alg)
                            self.save_tree_image(title=name_parsimony.split('.')[0], tree=[tree_parsimony], path=os.path.join(output_path_tree_image, name).replace('Trees', 'outputs/Plots'))
                            multi_trees[alg]['distance'][method].append(tree_distance)
                            multi_trees[alg]['parsimony'][method].append(tree_parsimony)
                    
                    self.save_tree_image(title=f'tree_{Path(file).stem}', tree=multi_trees, path=os.path.join(output_path_tree_image, f'tree_{Path(file).stem}').replace('Trees', 'outputs/Plots'))
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
                
        if self.mode == "auto":
            plot_heatmap_distances(data_dict=multi_trees, base_name='Acumulate', 
                                path=os.path.join(self.output_path, 'outputs/Plots'), distance_matrix=heatmap_matrix)
            logging.info("Heatmap acumulado gerado com sucesso.")
        
        clean_tmp(self.output_path)
        clean_NoPipe(self.input_path)
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

    def build_tree_distance_matrix(self, fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree, align_method):
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
        builder = TreeBuilder(fasta_path=fasta_path, 
                              output_path_align=output_path_align, 
                              output_path_dnd=output_path_dnd, 
                              path_dnd=path_dnd, 
                              output_path_tree=output_path_tree)
        
        try:
            if align_method == "clustalw":
                logging.debug(f"Alinhando sequências com ClustalW para {fasta_path}.")
                alng = self.aligner.align_sequences_clustalw(fasta_path=fasta_path,
                                                             path_dnd=path_dnd,
                                                             output_path_dnd=output_path_dnd,
                                                             output_path_align=output_path_align)
            elif align_method == "mafft":
                logging.debug(f"Alinhando sequências com MAFFT para {fasta_path}.")
                alng = self.aligner.align_sequences_mafft(fasta_path=fasta_path)
            else:
                logging.error(f"Método de alinhamento desconhecido: {align_method}")
                raise ValueError("Método de alinhamento não suportado.")
            
            if align_method == "mafft":
                tmp_align = os.path.join(self.output_path, 'tmp', 'align.aln')
                with open(tmp_align, "w") as f:
                    f.write(alng)
                alng = AlignIO.read(tmp_align, "fasta")
                logging.debug("Arquivo de alinhamento gerado e lido com sucesso (MAFFT).")
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
    
    def build_tree_parsimony(self, fasta_path, output_path_align, output_path_dnd, path_dnd, output_path_tree, align_method):
        """
        Constrói uma árvore filogenética usando o método de parcimônia.

        Alinha as sequências, calcula a matriz de distâncias e constrói uma árvore filogenética
        usando o método de parcimônia com uma árvore inicial gerada por NJ ou UPGMA.

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
        logging.info(f"Iniciando construção de árvore por parcimônia para {fasta_path}")
        builder = TreeBuilder(fasta_path=fasta_path, 
                              output_path_align=output_path_align, 
                              output_path_dnd=output_path_dnd, 
                              path_dnd=path_dnd, 
                              output_path_tree=output_path_tree)
        
        try:
            if align_method == "clustalw":
                logging.debug(f"Alinhando sequências com ClustalW para {fasta_path}.")
                alng = self.aligner.align_sequences_clustalw(fasta_path=fasta_path,
                                                             path_dnd=path_dnd,
                                                             output_path_dnd=output_path_dnd,
                                                             output_path_align=output_path_align)
            elif align_method == "mafft":
                logging.debug(f"Alinhando sequências com MAFFT para {fasta_path}.")
                alng = self.aligner.align_sequences_mafft(fasta_path=fasta_path)
            else:
                logging.error(f"Método de alinhamento desconhecido: {align_method}")
                raise ValueError("Método de alinhamento não suportado.")
            
            if align_method == "mafft":
                tmp_align = os.path.join(self.output_path, 'tmp', 'align.aln')
                with open(tmp_align, "w") as f:
                    f.write(alng)
                alng = AlignIO.read(tmp_align, "fasta")
                logging.debug("Arquivo de alinhamento gerado e lido com sucesso (MAFFT).")
        except Exception as e:
            logging.error(f"Erro no alinhamento das sequências para {fasta_path}: {e}", exc_info=True)
            raise

        try:
            distance_matrix = builder.distance_matrix(alignment=alng)
            logging.debug("Matriz de distâncias calculada para parcimônia com sucesso.")
            tree = builder.parsimony_constructor(distance_matrix=distance_matrix, 
                                                  alignment=alng, 
                                                  construct_tree_method=self.construct_tree_method)
            builder.save_tree(tree=tree, path=output_path_tree, format=self.output_format)
            logging.info(f"Árvore de parcimônia salva em {output_path_tree}.")
        except Exception as e:
            logging.error(f"Erro na construção da árvore por parcimônia para {fasta_path}: {e}", exc_info=True)
            raise
        return tree

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
        logging.info(f"Iniciando salvamento da imagem da árvore: {title}")
        try:
            if isinstance(tree, list):
                fig = plt.figure(figsize=(21, 10))
                ax = fig.add_subplot(1, 1, 1)
                Phylo.draw(tree[0], do_show=False, axes=ax, 
                           label_func=lambda x: None if 'Inner' in x.name else x.name)
                plt.title(title)
                plt.tight_layout()
                img_path = path.replace(self.output_format, 'png')
                fig.savefig(img_path)
                plt.close(fig)
                logging.info(f"Imagem da árvore {title} salva em {img_path}.")
            else:
                ax_index = 0
                fig, axs = plt.subplots(2, 4, figsize=(42, 18))
                axs = axs.flatten()

                width, height = fig.get_size_inches()
                fontsize = min(width, height) * 6
                plt.xticks(fontsize=fontsize * 0.2)
                plt.subplots_adjust(wspace=0.2, hspace=0.2)

                for alg, trs in tree.items():
                    for method, trees in trs.items():
                        for mtd, tr in trees.items():
                            if ax_index < len(axs):
                                Phylo.draw(tr[0], do_show=False, axes=axs[ax_index], 
                                           label_func=lambda x: None if x.name is None or 'Inner' in x.name else x.name)
                                axs[ax_index].set_title(f"{title} - {alg.lower()} - {method.lower()} - {mtd}", 
                                                         fontsize=fontsize * 0.15)
                            ax_index += 1

                plt.tight_layout()
                img_path = path.replace(self.output_format, 'png')
                fig.savefig(img_path)
                plt.close(fig)
                logging.info(f"Imagem composta das árvores {title} salva em {img_path}.")
        except Exception as e:
            logging.error(f"Erro ao salvar a imagem da árvore {title}: {e}", exc_info=True)