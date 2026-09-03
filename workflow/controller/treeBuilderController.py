from workflow.utils import run_logging, tool_runs
import os
import time
import sys
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
from Bio import Phylo
import matplotlib.pyplot as plt
from Bio import AlignIO
import numpy as np
import logging
import datetime


current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../..'))

from workflow.tree_construction.builder import TreeBuilder
from workflow.utils.dataValidation import (duplicate_names, duplicate_seq,
                                           validate_sequences, deduplicar_por_sequencia)
from workflow.utils.dataCleaning import clean_NoPipe, clean_tmp, copiar_arquivos
from workflow.utils.messages import Messages
from workflow.utils.metrics import process_rf_distance, plot_heatmap_distances
from workflow.alignment.alignmentSeq import AlignmentSeqs

from workflow.alignment.aligners import ALIGNERS, AlignerPolicy, resolve_aligner

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
    ignore_methods : list
        Lista de métodos a serem ignorados durante a construção de árvores.
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
            # `aligners` é propriedade (valida contra a biblioteca), então o
            # valor do experimento entra pelo campo privado.
            setattr(self, '_aligners' if k == 'aligners' else k, v)

        # D22 — o log é aberto uma vez por execução, com o `run_id` no nome.
        # `garantir` não troca o que o workflow já configurou.
        logfile_path = run_logging.garantir(os.path.join(self.output_path, 'outputs'))
        
        

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
        self.aligner = AlignmentSeqs({
            'num_threads': self.num_threads, 
            # 'num_threads': 16, 
            # 'max_memory_gb':,
            # 'max_sequences':,
            'output_path': self.output_path, 
            'logfile_path': logfile_path 
            })
        
        # Definir métodos a serem ignorados
        self.ignore_methods = self._parse_ignore_methods()

        # clean_tmp(self.output_path)
        clean_NoPipe(self.input_path)

        

    def _parse_ignore_methods(self):
        """
        Parseia os métodos a serem ignorados.
        Agora já recebe como lista do frontend.
        
        Return
        ------
        list
            Lista de métodos a serem ignorados.
        """
        if not hasattr(self, 'ignore_mode') or not self.ignore_mode:
            return []
        
        if isinstance(self.ignore_mode, list):
            return [method.lower() for method in self.ignore_mode if method and method != 'none']
        
        if isinstance(self.ignore_mode, str):
            return [method.strip().lower() for method in self.ignore_mode.split(',') if method.strip() and method.strip() != 'none']
        
        return []

    def _should_ignore_method(self, method):
        """
        Verifica se um método deve ser ignorado.
        
        Parameters
        ----------
        method : str
            Nome do método a ser verificado.
            
        Return
        ------
        bool
            True se o método deve ser ignorado, False caso contrário.
        """
        return method.lower() in self.ignore_methods

    def _initialize_multi_trees_structure(self):
        """
        Inicializa a estrutura para armazenar múltiplas árvores.

        As chaves vêm de `self.aligners`. Eram `"clustalo"` e `"mafft"` escritas
        à mão — um **terceiro** lugar que precisava concordar com a lista de
        alinhadores, depois dos dois laços do controlador. Com um braço novo, a
        estrutura não tinha a chave e as árvores daquele braço se perdiam: com
        `('mafft', 'mafft_iterative')` saíram **8 árvores em vez de 14**, e o
        segundo braço rendeu uma só. Duas listas divergindo é
        [D5](../../../docs/science/02-defeitos-que-alteram-resultado.md#d5) em
        outro assunto; três, é o mesmo defeito com mais chances de acontecer.

        Return
        ------
        dict
            Estrutura inicializada para armazenar árvores, uma entrada por
            alinhador declarado.
        """
        return {
            alinhador: {
                "distance": {"nj": [], "upgma": []},
                "parsimony": {"nj": [], "upgma": []},
                "iqtree": [],
                "fasttree": [],
                "raxml": [],
                "mrbayes": []
            }
            for alinhador in self.aligners
        }

    def _load_existing_tree(self, tree_path, tree_format):
        """
        Carrega uma árvore existente do arquivo.
        
        Parameters
        ----------
        tree_path : str
            Caminho para o arquivo da árvore.
        tree_format : str
            Formato da árvore.
            
        Return
        ------
        Bio.Phylo.BaseTree.Tree or None
            Árvore carregada ou None se não existir.
        """
        if os.path.exists(tree_path):
            try:
                return Phylo.read(tree_path, tree_format)
            except Exception as e:
                logging.warning(f"Erro ao carregar árvore existente {tree_path}: {e}")
        return None

    def __call__(self):
        """
        Executa o processo de construção das árvores baseado no modo especificado.

        Esta função coordena o fluxo de trabalho para construir árvores filogenéticas utilizando
        os métodos de construção por matriz de distâncias (NJ/UPGMA) ou parcimônia, dependendo do modo selecionado.

        Return
        ------
        None
        """
        # D18 — "auto" sugeria "escolhe sozinho o que faz sentido"; o que o
        # modo faz é rodar só distância e parcimônia, nunca os métodos
        # avançados. "basic" é o nome novo, honesto com o comportamento;
        # "auto" continua aceito como alias — os ~20 projetos já em disco com
        # mode: "auto" em config_backup.json não podem quebrar.
        MODOS_BASICOS = ("auto", "basic")
        descricao_basico = (
            "BÁSICO (apenas distância NJ/UPGMA e parcimônia — nenhum método "
            "avançado roda; use mode='advanced' para IQ-TREE/FastTree/RAxML-NG/MrBayes)"
        )
        mode_descriptions = {
            "distance": "DISTANCE TREE CONSTRUCTOR + CLUSTALW",
            "parsimony": "PARSIMONY + CLUSTALW",
            "auto": descricao_basico,
            "basic": descricao_basico,
            "advanced": "Advanced Methods"
        }

        if self.mode in mode_descriptions:
            description = mode_descriptions[self.mode]
            print(f"Iniciando construção das árvores utilizando o método: {description}\n")
            logging.info(f"Iniciando construção das árvores utilizando o método: {description}")
            logging.info(f"STEP: construction of trees using the method: {description}")
        else:
            logging.error(f"Modo desconhecido: {self.mode}")
            raise ValueError(f"Modo desconhecido: {self.mode}")

        if self.mode in MODOS_BASICOS:
            # D18, opção 2 do defeito: aviso explícito, não só descrição —
            # "Completed successfully!" no fim da execução (emitido pelo
            # subtree builder, sempre) não distingue básico de advanced.
            aviso = (
                "Modo básico: métodos avançados (IQ-TREE, FastTree, RAxML-NG, "
                "MrBayes) NÃO serão executados nesta execução. Use "
                "mode='advanced' para incluí-los."
            )
            print(f"AVISO: {aviso}\n")
            logging.warning(aviso)

        print(f"Métodos ignorados: {', '.join([m.upper() for m in self.ignore_methods])}\n")
        logging.info(f"Métodos ignorados: {', '.join([m.upper() for m in self.ignore_methods])}")
        
        heatmap_matrix = np.zeros((8, 8))
        multi_trees = self._initialize_multi_trees_structure()
        base_folder = self.input_path.split('/')[-1]

        for file in tqdm(self.files, desc="Building trees...", ascii="░▒█"):
            if any(ext in file for ext in ['.dnd', '.json']):
                logging.debug(f"Arquivo {file} ignorado por ser .dnd/.json")
                continue
            
            logging.info(f"Iniciando processamento do arquivo: {os.path.join(base_folder,file)}")
            
            start_cycle = time.time()
            fasta_path = os.path.join(self.input_path, file)
            file_stem = Path(file).stem
            
            # Preparar caminhos de saída
            output_paths = self._prepare_output_paths(file_stem)
            
            # Validar e preparar arquivo FASTA
            fasta_path = self._validate_and_prepare_fasta(file, fasta_path, file_stem)
            
            try:
                if self.mode == "distance":
                    trees_built = self._process_single_mode(file_stem, fasta_path, output_paths, "distance")
                    self.count_trees += trees_built
                    
                elif self.mode == "parsimony":
                    trees_built = self._process_single_mode(file_stem, fasta_path, output_paths, "parsimony")
                    self.count_trees += trees_built
                    
                elif self.mode in MODOS_BASICOS:
                    trees_built, file_multi_trees = self._process_auto_mode(file_stem, fasta_path, output_paths)
                    self.count_trees += trees_built
                    
                    # Atualizar multi_trees e processar RF distance
                    for alg in file_multi_trees:
                        for method_type in file_multi_trees[alg]:
                            if isinstance(file_multi_trees[alg][method_type], dict):
                                for sub_method in file_multi_trees[alg][method_type]:
                                    multi_trees[alg][method_type][sub_method].extend(file_multi_trees[alg][method_type][sub_method])
                            else:
                                multi_trees[alg][method_type].extend(file_multi_trees[alg][method_type])
                    
                    rf_scores = process_rf_distance(file_multi_trees)
                    heatmap_matrix = self.somarMatrizes(
                        heatmap_matrix,
                        plot_heatmap_distances(data_dict=file_multi_trees, scores=rf_scores, 
                                             base_name=file_stem, path=os.path.join(self.output_path, 'outputs/Plots'))
                    )
                    
                elif self.mode == "advanced":
                    trees_built, file_multi_trees = self._process_advanced_mode(file_stem, fasta_path, output_paths)
                    self.count_trees += trees_built
                    
                    # Atualizar multi_trees e processar RF distance
                    for alg in file_multi_trees:
                        for method_type in file_multi_trees[alg]:
                            if isinstance(file_multi_trees[alg][method_type], dict):
                                for sub_method in file_multi_trees[alg][method_type]:
                                    multi_trees[alg][method_type][sub_method].extend(file_multi_trees[alg][method_type][sub_method])
                            else:
                                multi_trees[alg][method_type].extend(file_multi_trees[alg][method_type])
                    
                    rf_scores = process_rf_distance(file_multi_trees)
                    heatmap_matrix = self.somarMatrizes(
                        heatmap_matrix,
                        plot_heatmap_distances(data_dict=file_multi_trees, scores=rf_scores, 
                                             base_name=file_stem, path=os.path.join(self.output_path, 'outputs/Plots'))
                    )
                else:
                    logging.error(f"Modo não suportado: {self.mode}")
                    continue
                    
                cycle_time = time.time() - start_cycle
                self.list_times.append(cycle_time)
                logging.info(f"Arquivo {file} processado com sucesso em {cycle_time:.2f} segundos.")
                
            except Exception as e:
                logging.error(f"Erro ao processar o arquivo {file}: {e}", exc_info=True)
                exit(1)
                
        # Gerar heatmap acumulado para modos auto e advanced
        if self.mode in MODOS_BASICOS or self.mode == "advanced":
            plot_heatmap_distances(data_dict=multi_trees, base_name='Acumulate', 
                               path=os.path.join(self.output_path, 'outputs/Plots'), distance_matrix=heatmap_matrix)
            logging.info("Heatmap acumulado gerado com sucesso.")
        
        # Limpeza final
        clean_NoPipe(self.input_path)
        copiar_arquivos(os.path.join(self.output_path, 'tmp'), os.path.join(self.output_path, 'Align'))
        # clean_tmp(self.output_path)
        logging.info("Diretórios temporários limpos.")

        # self.msg.resume_tree(start=self.start, sum_time=self.list_times, num_trees=self.count_trees, 
        #                      output_format=self.output_format, n_nodes=self.count_nodes, 
        #                      method=getattr(self, 'construct_tree_method', 'unknown'))
        logging.info("Processo de construção de árvores finalizado.")

    def _prepare_output_paths(self, file_stem):
        """
        Prepara os caminhos de saída para um arquivo.
        
        Parameters
        ----------
        file_stem : str
            Nome do arquivo sem extensão.
            
        Return
        ------
        dict
            Dicionário com os caminhos de saída.
        """
        return {
            'align_base': os.path.join(self.output_path, 'tmp', f'{file_stem}.aln'),
            'align_html': os.path.join(self.output_path, 'Align', f'{file_stem}.html'),
            'dnd': os.path.join(self.output_path, 'tmp', f'{file_stem}.dnd'),
            'dnd_original': os.path.join(self.input_path, f'{file_stem}.dnd'),
            'tree': os.path.join(self.output_path, 'Trees'),
            'tree_image': os.path.join(self.output_path, 'Trees')
        }

    def _validate_and_prepare_fasta(self, file, fasta_path, file_stem):
        """
        Valida e prepara o arquivo FASTA.
        
        Parameters
        ----------
        file : str
            Nome do arquivo.
        fasta_path : str
            Caminho para o arquivo FASTA.
        file_stem : str
            Nome do arquivo sem extensão.
            
        Return
        ------
        str
            Caminho para o arquivo FASTA válido.
        """
        if not duplicate_names(fasta_path) and not duplicate_seq(fasta_path)[0] and validate_sequences(fasta_path):
            logging.info(f"Arquivo {file} passou nas validações de sequência.")
        else:
            # D23 — o descarte deixa de ser silencioso. `n` é número de Methods:
            # o projeto chamado VARV-49 recebe 52 registros e produz 49 folhas,
            # porque a aquisição baixa RefSeq e GenBank do mesmo genoma. A
            # composição não muda aqui — a decisão do usuário foi declarar agora
            # e corrigir a aquisição depois —, mas quem descartou o quê passa a
            # ficar registrado no log e no manifesto.
            logging.warning(f"Arquivo {file} contém duplicatas ou erros, deduplicando por sequência.")
            fasta_path, descartados = deduplicar_por_sequencia(
                file_stem, fasta_path, self.input_path)
            if descartados:
                tool_runs.registrar(
                    'deduplicacao', ['deduplicar_por_sequencia', file_stem],
                    saida=fasta_path,
                    descartados=[{"descartado": d, "mantido": m} for d, m in descartados],
                    nota=('acessos com sequência idêntica; o mantido é o primeiro do '
                          'arquivo — ver D23'))

        return fasta_path

    def _process_single_mode(self, file_stem, fasta_path, output_paths, mode_type):
        """
        Processa um arquivo no modo single (distance ou parsimony).
        
        Parameters
        ----------
        file_stem : str
            Nome do arquivo sem extensão.
        fasta_path : str
            Caminho para o arquivo FASTA.
        output_paths : dict
            Dicionário com caminhos de saída.
        mode_type : str
            Tipo de modo ('distance' ou 'parsimony').
            
        Return
        ------
        int
            Número de árvores construídas.
        """
        name = f'tree_{file_stem}_{mode_type}.{self.output_format}'
        logging.info(f"Construindo árvore por {mode_type} para o arquivo {file_stem}")
        
        if mode_type == "distance":
            tree = self.build_tree_distance_matrix(
                fasta_path, output_paths['align_base'], output_paths['dnd'], 
                output_paths['dnd_original'], os.path.join(output_paths['tree'], name), 
                self.align_method, output_paths['align_html']
            )
        else:  # parsimony
            tree = self.build_tree_parsimony(
                fasta_path, output_paths['align_base'], output_paths['dnd'],
                output_paths['dnd_original'], os.path.join(output_paths['tree'], name),
                self.align_method, output_paths['align_html']
            )
            
        self.save_tree_image(
            title=name, tree=[tree], 
            path=os.path.join(output_paths['tree_image'], name).replace('Trees', 'outputs/Plots')
        )
        
        return 1

    def _process_auto_mode(self, file_stem, fasta_path, output_paths):
        """
        Processa um arquivo no modo auto.
        
        Parameters
        ----------
        file_stem : str
            Nome do arquivo sem extensão.
        fasta_path : str
            Caminho para o arquivo FASTA.
        output_paths : dict
            Dicionário com caminhos de saída.
            
        Return
        ------
        tuple
            (número de árvores construídas, dicionário com árvores)
        """
        multi_trees = self._initialize_multi_trees_structure()
        trees_built = 0
        
        for method in ['nj', 'upgma']:
            for alg in self.aligners:
                self.construct_tree_method = method
                
                # Processar árvore de distância
                if not self._should_ignore_method("distance"):
                    trees_built += self._process_tree_method(
                        file_stem, fasta_path, output_paths, alg, method, "distance", multi_trees
                    )
                
                # Processar árvore de parcimônia
                if not self._should_ignore_method("parsimony"):
                    trees_built += self._process_tree_method(
                        file_stem, fasta_path, output_paths, alg, method, "parsimony", multi_trees
                    )
        
        return trees_built, multi_trees

    def _process_advanced_mode(self, file_stem, fasta_path, output_paths):
        """
        Processa um arquivo no modo advanced.
        
        Parameters
        ----------
        file_stem : str
            Nome do arquivo sem extensão.
        fasta_path : str
            Caminho para o arquivo FASTA.
        output_paths : dict
            Dicionário com caminhos de saída.
            
        Return
        ------
        tuple
            (número de árvores construídas, dicionário com árvores)
        """
        multi_trees = self._initialize_multi_trees_structure()
        trees_built = 0
        advanced_methods = ['iqtree', 'fasttree', 'raxml', 'mrbayes']
        
        for alg in self.aligners:
            for method in ['nj', 'upgma']:
                self.construct_tree_method = method
                
                # Processar métodos básicos
                if not self._should_ignore_method("distance"):
                    trees_built += self._process_tree_method(
                        file_stem, fasta_path, output_paths, alg, method, "distance", multi_trees
                    )
                
                if not self._should_ignore_method("parsimony"):
                    trees_built += self._process_tree_method(
                        file_stem, fasta_path, output_paths, alg, method, "parsimony", multi_trees
                    )
                
                # Processar métodos avançados
                for adv_method in advanced_methods:
                    if not self._should_ignore_method(adv_method):
                        trees_built += self._process_advanced_tree_method(
                            file_stem, fasta_path, output_paths, alg, adv_method, multi_trees
                        )
        
        return trees_built, multi_trees

    def _process_tree_method(self, file_stem, fasta_path, output_paths, alg, method, method_type, multi_trees):
        """
        Processa um método de construção de árvore específico.
        
        Parameters
        ----------
        file_stem : str
            Nome do arquivo sem extensão.
        fasta_path : str
            Caminho para o arquivo FASTA.
        output_paths : dict
            Dicionário com caminhos de saída.
        alg : str
            Algoritmo de alinhamento ('clustalo' ou 'mafft').
        method : str
            Método de construção ('nj' ou 'upgma').
        method_type : str
            Tipo de método ('distance' ou 'parsimony').
        multi_trees : dict
            Dicionário para armazenar as árvores.
            
        Return
        ------
        int
            1 se uma árvore foi construída, 0 caso contrário.
        """
        name = f'tree_{file_stem}_{alg}_{method}_{method_type}.{self.output_format}'
        output_path_align = output_paths['align_base'].replace('.aln', f'_{alg}.aln')
        output_path_tree = os.path.join(self.output_path, 'Trees', name)
        
        # Verificar se a árvore já existe
        existing_tree = self._load_existing_tree(output_path_tree, self.output_format)
        if existing_tree:
            multi_trees[alg][method_type][method].append(existing_tree)
            return 0
        
        # Construir nova árvore
        logging.debug(f"Construindo árvore de {method_type} ({alg} - {method}) para o arquivo {file_stem}")
        
        if method_type == "distance":
            tree = self.build_tree_distance_matrix(
                fasta_path, output_path_align, output_paths['dnd'], 
                output_paths['dnd_original'], output_path_tree, alg, output_paths['align_html']
            )
        else:  # parsimony
            tree = self.build_tree_parsimony(
                fasta_path, output_path_align, output_paths['dnd'],
                output_paths['dnd_original'], output_path_tree, alg, output_paths['align_html']
            )
            
        self.save_tree_image(
            title=name, tree=[tree], 
            path=os.path.join(output_paths['tree_image'], name).replace('Trees', 'outputs/Plots')
        )
        
        multi_trees[alg][method_type][method].append(tree)
        return 1

    def _process_advanced_tree_method(self, file_stem, fasta_path, output_paths, alg, method, multi_trees):
        """
        Processa um método avançado de construção de árvore.
        
        Parameters
        ----------
        file_stem : str
            Nome do arquivo sem extensão.
        fasta_path : str
            Caminho para o arquivo FASTA.
        output_paths : dict
            Dicionário com caminhos de saída.
        alg : str
            Algoritmo de alinhamento ('clustalo' ou 'mafft').
        method : str
            Método avançado ('iqtree', 'fasttree', 'raxml', 'mrbayes').
        multi_trees : dict
            Dicionário para armazenar as árvores.
            
        Return
        ------
        int
            1 se uma árvore foi construída, 0 caso contrário.
        """
        name = f'tree_{file_stem}_{alg}_{method}.{self.output_format}'
        output_path_align = output_paths['align_base'].replace('.aln', f'_{alg}.aln')
        output_path_tree = os.path.join(self.output_path, 'Trees', name)
        
        # Verificar se a árvore já existe
        existing_tree = self._load_existing_tree(output_path_tree, self.output_format)
        if existing_tree:
            multi_trees[alg][method].append(existing_tree)
            return 0
        
        # Construir nova árvore
        logging.debug(f"Construindo árvore por {method} ({alg}) para o arquivo {file_stem}")
        
        builder_methods = {
            'iqtree': self.build_tree_iqtree,
            'fasttree': self.build_tree_fasttree, 
            'raxml': self.build_tree_raxml,
            'mrbayes': self.build_tree_mrbayes
        }
        
        if method in builder_methods:
            tree = builder_methods[method](
                fasta_path, output_path_align, output_path_tree, alg, output_paths['align_html']
            )
            
            self.save_tree_image(
                title=name, tree=[tree],
                path=os.path.join(output_paths['tree_image'], name).replace('Trees', 'outputs/Plots')
            )
            
            multi_trees[alg][method].append(tree)
            return 1
        
        return 0

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
                logging.info(f"STEP: Reusing Aligning...")
                
                logging.info(f"Arquivo de alinhamento já existe: {output_path_align}. Reutilizando.")
                alng = AlignIO.read(output_path_align, "fasta")
            else:
                logging.info(f"STEP: Aligning seqs...")
                align_method, _ = self._resolver_alinhador(align_method, fasta_path)
                logging.debug(f"Alinhando sequências com {align_method} para {fasta_path}.")
                alng = self._alinhar(align_method, fasta_path,
                                     output_path_align, output_path_align_html)
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
                alng = AlignIO.read(output_path_align, "fasta")
            else:
                logging.info(f"STEP: Aligning seqs...")
                align_method, _ = self._resolver_alinhador(align_method, fasta_path)
                logging.debug(f"Alinhando sequências com {align_method} para {fasta_path}.")
                alng = self._alinhar(align_method, fasta_path,
                                     output_path_align, output_path_align_html)
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
    
    def _dimensoes_do_conjunto(self, fasta_path):
        """Número de sequências e comprimento da MAIOR delas, em pares de base.

        É o comprimento máximo, e não a média: uma sequência só é o bastante
        para estourar a memória do alinhador."""
        n = 0
        maior = 0
        for record in SeqIO.parse(fasta_path, "fasta"):
            n += 1
            maior = max(maior, len(record.seq))
        return n, maior

    def _resolver_alinhador(self, align_method, fasta_path):
        """
        Decide qual alinhador vai rodar, e **devolve o nome do que rodou**.

        Substitui `_isExecutableByClustalO`, que trocava Clustal Omega por MAFFT
        e devolvia a troca sem que o chamador nomeasse o arquivo pelo alinhador
        efetivo. É [D1](../../../docs/science/02-defeitos-que-alteram-resultado.md#d1):
        nos experimentos de Variola, metade dos "pipelines" são cópias byte a
        byte de MAFFT com nome de `clustalo`, e o fator alinhador não existe.

        O padrão agora é **falhar** com o motivo. A substituição só acontece se
        o experimento a autorizar (`aligner_on_unavailable="fallback"`), e nesse
        caso quem chama **tem de usar o nome devolvido** para nomear a saída.

        Return
        ------
        tuple
            ``(alinhador efetivo, motivo da substituição ou None)``.
        """
        n_seqs, maior_bp = self._dimensoes_do_conjunto(fasta_path)
        politica = AlignerPolicy(
            on_unavailable=getattr(self, "aligner_on_unavailable", "fail"))
        efetivo, motivo = resolve_aligner(align_method, n_seqs, maior_bp, politica)

        if motivo:
            logging.warning(f"      SUBSTITUIÇÃO DE ALINHADOR: {motivo}")
            logging.warning(f"      A saída será nomeada como '{efetivo}', não '{align_method}'.")
        return efetivo, motivo

    #: Braços do fator alinhador quando o experimento não os declara.
    #:
    #: Era `['clustalo', 'mafft']` **fixo no código**, em dois lugares. Em
    #: genoma de poxvírus o Clustal Omega não termina (medido: 1 h sem concluir,
    #: com pico de 220 MB — é limite de tempo, não de memória), então o braço
    #: `clustalo` acabava sendo MAFFT com outro nome: é a
    #: [D1](../../../docs/science/02-defeitos-que-alteram-resultado.md#d1).
    #:
    #: Decisão do usuário em 2026-08-26: o fator passa a ser **duas estratégias
    #: do MAFFT**. Mesma ferramenta, mesma versão, mesmo binário; o que muda é o
    #: algoritmo — progressivo contra iterativo. É o único par que existe tanto
    #: em *Variola* quanto em Zika, porque o MUSCLE 5.3 recusa genoma longo por
    #: projeto e o Clustal Omega não termina.
    ALINHADORES_PADRAO = ('mafft', 'mafft_iterative')

    @property
    def aligners(self):
        """
        Braços do fator alinhador, declarados pelo experimento.

        Vem de `tree_config.aligners`; sem ele, `ALINHADORES_PADRAO`. Um
        alinhador desconhecido é erro, não aviso — um braço que não existe
        produziria árvore com nome de um método que nunca rodou, que é a forma
        de D1.
        """
        escolhidos = tuple(getattr(self, '_aligners', None) or self.ALINHADORES_PADRAO)
        desconhecidos = [a for a in escolhidos if a not in ALIGNERS]
        if desconhecidos:
            raise ValueError(
                f"Alinhador(es) desconhecido(s) em `aligners`: {desconhecidos}. "
                f"Disponíveis: {sorted(ALIGNERS)}")
        return escolhidos

    def _alinhar(self, align_method, fasta_path, output_path_align, output_path_align_html):
        """Executa o alinhador pedido. Um método desconhecido é erro, não aviso."""
        if align_method == "clustalo":
            return self.aligner.align_sequences_clustalo(
                fasta_path=fasta_path, output_path_align=output_path_align,
                output_path_html=output_path_align_html)
        if align_method in ("mafft", "mafft_iterative"):
            # A estratégia vem da biblioteca, não de uma tabela local: dois
            # braços do fator alinhador usam o mesmo binário e diferem só nela.
            return self.aligner.align_sequences_mafft(
                fasta_path=fasta_path, output_path_align=output_path_align,
                output_path_html=output_path_align_html,
                estrategia=ALIGNERS[align_method].estrategia)
        if align_method == "muscle":
            return self.aligner.align_sequences_muscle(
                fasta_path=fasta_path, output_path_align=output_path_align,
                output_path_html=output_path_align_html)
        raise ValueError(
            f"Método de alinhamento não suportado: '{align_method}'. "
            f"Disponíveis: {sorted(ALIGNERS)}")

    def _get_alignment(self, fasta_path, output_path_align, align_method, output_path_align_html):
        """
        Alinhamento do arquivo, reaproveitando o que já estiver em disco.

        **Existir não é servir.** A checagem era `os.path.exists` e nada mais.
        O alinhador escreve por redirecionamento de stdout, então uma execução
        interrompida no meio do alinhamento deixa um arquivo de **0 byte** — e a
        execução seguinte o encontrava, anunciava `Reutilizando` e devolvia um
        alinhamento vazio. Medido ao interromper o MAFFT iterativo no VARV-49.

        É a terceira vez que este projeto tropeça na mesma forma: foi assim com
        o binário resolvido por `command -v` que não executava, e com o `pnpm`
        que existia no PATH e abortava a cada chamada.
        """
        try:
            if os.path.exists(output_path_align):
                if os.path.getsize(output_path_align) == 0:
                    # Um arquivo vazio é o rastro de um alinhamento interrompido.
                    # Reaproveitá-lo silenciosamente é pior do que refazê-lo.
                    logging.warning(
                        f"Alinhamento existente está vazio (0 byte) e será refeito: "
                        f"{os.path.basename(output_path_align)}. É o rastro de uma "
                        f"execução interrompida.")
                    os.remove(output_path_align)
                else:
                    logging.info(f"Arquivo de alinhamento já existe: {output_path_align}. Reutilizando.")
                    return AlignIO.read(output_path_align, "fasta")

            logging.info(f"STEP: Aligning seqs...")
            align_method, _ = self._resolver_alinhador(align_method, fasta_path)
            return self._alinhar(align_method, fasta_path,
                                 output_path_align, output_path_align_html)
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