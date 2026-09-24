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

from workflow.tree_construction.builder import (TreeBuilder, MRBAYES_DEFAULTS,
                                                reproducibility_settings)
from workflow.tree_construction.modelo_substituicao import (
    MODELO_LEGADO, POLITICA_SELECAO_PADRAO, POLITICAS_SELECAO_MODELO,
    TAXAS_CANDIDATAS, SelecaoDeModeloFalhou, selecionar_modelo)
from workflow.utils.dataValidation import (duplicate_names, duplicate_seq,
                                           validate_sequences, deduplicar_por_sequencia)
from workflow.utils.dataCleaning import clean_NoPipe, clean_tmp, copiar_arquivos
from workflow.utils.messages import Messages
from workflow.utils.metrics import process_rf_distance, plot_heatmap_distances
from workflow.alignment.alignmentSeq import AlignmentSeqs

from workflow.alignment.aligners import ALIGNERS, AlignerPolicy, resolve_aligner

#: D18 — "basic" é o nome honesto; "auto" continua aceito como alias legado
#: (projetos já em disco). Nível de módulo, importável por `workflow.py` — que
#: precisa da mesma lista para computar `execution_mode` do manifesto **antes**
#: do controlador rodar. Duas cópias da mesma tupla já divergiram uma vez neste
#: projeto (D5, D19); aqui é a fonte única.
MODOS_BASICOS = ("auto", "basic")

#: M7.6 — o que fazer quando um método de inferência não produz árvore.
#:
#: ``fail`` (padrão) é o comportamento de sempre: a exceção sobe e a execução
#: para. ``continue`` segue com os demais pipelines. Nos dois casos o manifesto
#: registra o pipeline como ``tentado_e_falhou``, com o motivo — o que muda é
#: só se a execução continua. É a mesma forma de `aligner_on_unavailable`:
#: o padrão é falhar, e seguir adiante é decisão declarada do experimento.
#:
#: Sem ``continue``, a saída para um método que quebra era tirá-lo à mão via
#: `ignore_mode` e rodar de novo — e aí "quebrou" e "excluído de propósito"
#: ficavam indistinguíveis no artefato (D18, DM-11).
POLITICAS_FALHA_METODO = ("fail", "continue")


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

        # M7.6 — validado aqui, antes de qualquer alinhamento: um valor
        # desconhecido descoberto só na primeira falha seria uma segunda falha.
        self.on_method_failure = getattr(self, 'on_method_failure', None) or 'fail'
        if self.on_method_failure not in POLITICAS_FALHA_METODO:
            raise ValueError(
                f"on_method_failure desconhecido: {self.on_method_failure!r}. "
                f"Válidos: {POLITICAS_FALHA_METODO}")
        # Pipelines que já falharam nesta execução. O modo `advanced` visita
        # cada método avançado uma vez por método de distância (NJ e UPGMA);
        # com `continue`, sem isto o mesmo método quebrado rodaria de novo.
        self._falharam = set()

        # M7.3 — política de modelo, validada na largada pelo mesmo motivo.
        self.model_selection = (getattr(self, 'model_selection', None)
                                or POLITICA_SELECAO_PADRAO)
        if self.model_selection not in POLITICAS_SELECAO_MODELO:
            raise ValueError(
                f"model_selection desconhecido: {self.model_selection!r}. "
                f"Válidos: {POLITICAS_SELECAO_MODELO}")
        # Uma seleção por alinhamento: `caminho do alinhamento -> tradução`,
        # ou a exceção, se falhou (não se tenta de novo para o método seguinte).
        self._selecoes_modelo = {}

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

    # ------------------------------------------------------------------ #
    # M7.6 — desfecho de cada pipeline, no manifesto
    # ------------------------------------------------------------------ #

    def _registrar_ignorado(self, rotulo, saida, alinhador, motivo):
        tool_runs.registrar_metodo(rotulo, "ignorado_por_configuracao", saida,
                                   alinhador=alinhador, motivo=motivo)

    def _construir_com_desfecho(self, rotulo, saida, alinhador, construir):
        """
        Roda `construir()` e registra o desfecho do pipeline.

        Toda exceção vira ``tentado_e_falhou`` **antes** de qualquer decisão
        sobre continuar: mesmo com a política ``fail``, o `finally` de
        `workflow.py` drena o registro, e o manifesto da execução que morreu
        diz qual pipeline a matou e por quê.

        Return
        ------
        tree or None
            ``None`` só com a política ``continue`` e depois de uma falha.
        """
        if saida in self._falharam:
            return None
        try:
            tree = construir()
        except Exception as e:
            self._falharam.add(saida)
            tool_runs.registrar_metodo(rotulo, "tentado_e_falhou", saida,
                                       alinhador=alinhador, motivo=self._motivo(e))
            if self.on_method_failure == 'fail':
                raise
            logging.error(
                f"Pipeline {alinhador}/{rotulo} FALHOU e foi deixado de fora desta "
                f"execução (on_method_failure='continue'): {self._motivo(e)}")
            return None
        tool_runs.registrar_metodo(rotulo, "executado", saida, alinhador=alinhador)
        return tree

    @staticmethod
    def _motivo(e):
        """Tipo e mensagem da exceção, mais a cauda do `stderr` quando houver."""
        motivo = f"{type(e).__name__}: {e}"
        stderr = getattr(e, 'stderr', None)
        if isinstance(stderr, bytes):
            stderr = stderr.decode('utf-8', errors='ignore')
        if stderr and stderr.strip():
            motivo += f" | stderr: {stderr.strip()[-300:]}"
        return motivo

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
        # mode: "auto" em config_backup.json não podem quebrar. MODOS_BASICOS
        # é constante de módulo (topo do arquivo) — `workflow.py` usa a mesma.
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
        
        saida = os.path.join(output_paths['tree'], name)
        if mode_type == "distance":
            construir = lambda: self.build_tree_distance_matrix(
                fasta_path, output_paths['align_base'], output_paths['dnd'], 
                output_paths['dnd_original'], saida,
                self.align_method, output_paths['align_html']
            )
        else:  # parsimony
            construir = lambda: self.build_tree_parsimony(
                fasta_path, output_paths['align_base'], output_paths['dnd'],
                output_paths['dnd_original'], saida,
                self.align_method, output_paths['align_html']
            )
        tree = self._construir_com_desfecho(mode_type, saida, self.align_method, construir)
        if tree is None:
            return 0
            
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
                else:
                    self._registrar_basico_ignorado(file_stem, alg, method, "distance")
                
                # Processar árvore de parcimônia
                if not self._should_ignore_method("parsimony"):
                    trees_built += self._process_tree_method(
                        file_stem, fasta_path, output_paths, alg, method, "parsimony", multi_trees
                    )
                else:
                    self._registrar_basico_ignorado(file_stem, alg, method, "parsimony")

        # D18/M7.6 — o modo básico não roda método avançado nenhum. Isso é
        # configuração, e o manifesto diz: cada pipeline avançado aparece como
        # ignorado, com o `mode` como motivo — não simplesmente ausente.
        for alg in self.aligners:
            for adv_method in self.METODOS_AVANCADOS:
                self._registrar_ignorado(
                    adv_method, self._saida_avancado(file_stem, alg, adv_method), alg,
                    f"mode={self.mode!r} é básico e não executa métodos avançados (D18)")
        
        return trees_built, multi_trees

    #: Os quatro métodos avançados, na ordem em que o modo `advanced` os roda.
    METODOS_AVANCADOS = ('iqtree', 'fasttree', 'raxml', 'mrbayes')

    def _saida_basico(self, file_stem, alg, method, method_type):
        name = f'tree_{file_stem}_{alg}_{method}_{method_type}.{self.output_format}'
        return os.path.join(self.output_path, 'Trees', name)

    def _saida_avancado(self, file_stem, alg, method):
        name = f'tree_{file_stem}_{alg}_{method}.{self.output_format}'
        return os.path.join(self.output_path, 'Trees', name)

    def _registrar_basico_ignorado(self, file_stem, alg, method, method_type):
        self._registrar_ignorado(
            f'{method}_{method_type}', self._saida_basico(file_stem, alg, method, method_type),
            alg, f"ignore_mode contém {method_type!r}")

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
        advanced_methods = list(self.METODOS_AVANCADOS)
        
        for alg in self.aligners:
            for method in ['nj', 'upgma']:
                self.construct_tree_method = method
                
                # Processar métodos básicos
                if not self._should_ignore_method("distance"):
                    trees_built += self._process_tree_method(
                        file_stem, fasta_path, output_paths, alg, method, "distance", multi_trees
                    )
                else:
                    self._registrar_basico_ignorado(file_stem, alg, method, "distance")
                
                if not self._should_ignore_method("parsimony"):
                    trees_built += self._process_tree_method(
                        file_stem, fasta_path, output_paths, alg, method, "parsimony", multi_trees
                    )
                else:
                    self._registrar_basico_ignorado(file_stem, alg, method, "parsimony")
                
                # Processar métodos avançados
                for adv_method in advanced_methods:
                    if not self._should_ignore_method(adv_method):
                        trees_built += self._process_advanced_tree_method(
                            file_stem, fasta_path, output_paths, alg, adv_method, multi_trees
                        )
                    else:
                        self._registrar_ignorado(
                            adv_method, self._saida_avancado(file_stem, alg, adv_method), alg,
                            f"ignore_mode contém {adv_method!r}")
        
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
        rotulo = f'{method}_{method_type}'
        
        # Verificar se a árvore já existe
        existing_tree = self._load_existing_tree(output_path_tree, self.output_format)
        if existing_tree:
            tool_runs.registrar_metodo(
                rotulo, "reaproveitado", output_path_tree, alinhador=alg,
                motivo="árvore já estava em disco; lida, não produzida nesta execução")
            multi_trees[alg][method_type][method].append(existing_tree)
            return 0
        
        # Construir nova árvore
        logging.debug(f"Construindo árvore de {method_type} ({alg} - {method}) para o arquivo {file_stem}")
        
        if method_type == "distance":
            construir = lambda: self.build_tree_distance_matrix(
                fasta_path, output_path_align, output_paths['dnd'], 
                output_paths['dnd_original'], output_path_tree, alg, output_paths['align_html']
            )
        else:  # parsimony
            construir = lambda: self.build_tree_parsimony(
                fasta_path, output_path_align, output_paths['dnd'],
                output_paths['dnd_original'], output_path_tree, alg, output_paths['align_html']
            )
        tree = self._construir_com_desfecho(rotulo, output_path_tree, alg, construir)
        if tree is None:
            return 0
            
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
            tool_runs.registrar_metodo(
                method, "reaproveitado", output_path_tree, alinhador=alg,
                motivo="árvore já estava em disco; lida, não produzida nesta execução")
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
            tree = self._construir_com_desfecho(
                method, output_path_tree, alg,
                lambda: builder_methods[method](
                fasta_path, output_path_align, output_path_tree, alg, output_paths['align_html']
                ))
            if tree is None:
                return 0
            
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
    
    #: D26 — os quatro métodos avançados instanciavam `TreeBuilder` sem
    #: repassar `random_seed`/`raxml_threads`/`iqtree_threads`, então
    #: `reproducibility_settings` (`builder.py`) sempre caía nos defaults do
    #: módulo (12345/4/4) — não importa o que o `tree_config` do experimento
    #: pedisse. O manifesto (`workflow.py:100`) chama a mesma função com o
    #: `tree_config` de verdade e por isso *declarava* o valor pedido enquanto
    #: a chamada real usava outro. `getattr` com fallback ausente (não `None`)
    #: porque `reproducibility_settings` faz `config.get(chave, padrao)`: uma
    #: chave presente com valor `None` quebraria em `int(None)`.
    def _reproducibility_kwargs(self):
        """`tree_config` do experimento, restrito ao que reprodutibilidade lê."""
        chaves = ('random_seed', 'raxml_threads', 'iqtree_threads')
        return {chave: getattr(self, chave) for chave in chaves if hasattr(self, chave)}

    def _mrbayes_kwargs(self):
        """
        `mrbayes_*` do `tree_config` — M7.4. Mesma lição de D26: um parâmetro
        que o builder sabe ler e o controlador não repassa é constante com
        nome de configuração.

        Repassa a chave **presente**, inclusive com valor `None`: para
        `mrbayes_asdsf_max`, `null` explícito significa "gate desligado" e é
        diferente de ausente (padrão 0,01) — `mrbayes_settings` distingue os
        dois por `in`, então o `None` não pode ser filtrado aqui.
        """
        return {chave: getattr(self, chave) for chave in MRBAYES_DEFAULTS
                if hasattr(self, chave)}

    def _modelo_do_alinhamento(self, output_path_align, alng, alinhador):
        """
        Modelo de substituição do alinhamento — M7.3, opção C′.

        A seleção (`iqtree -m MF --mset mrbayes -mrate E,I,G,I+G`, só
        ModelFinder, sem árvore) roda **uma vez por alinhamento**, na primeira
        vez em que um método baseado em modelo (IQ-TREE, RAxML-NG, MrBayes)
        precisa dela — e antes dele. Preguiçosa de propósito: se as três
        árvores já estão em disco (`reaproveitado`), nada as usa, e uma seleção
        rodada à toa entraria no manifesto como se tivesse decidido alguma
        coisa. O FastTree não passa por aqui: não segue modelo escolhido
        (GTR+CAT20, declarado como aproximação).

        Falha na seleção (binário ausente, tempo, relatório sem modelo, modelo
        sem tradução) é guardada e **relevantada para cada método que
        dependia dela**: cada um vira ``tentado_e_falhou`` com o motivo (M7.6),
        e a política ``on_method_failure`` decide se a execução continua — o
        FastTree e os métodos de distância seguem com ``continue``. Não há
        recuo para GTR+G: seria substituição silenciosa (D1). Para rodar o
        literal antigo, `model_selection: "nenhuma"`, declarado.

        Return
        ------
        dict
            Tradução do modelo para cada ferramenta
            (`modelo_substituicao.traduzir`), ou `MODELO_LEGADO`.

        Raises
        ------
        SelecaoDeModeloFalhou
        """
        selecoes = self.__dict__.setdefault('_selecoes_modelo', {})
        if output_path_align in selecoes:
            anterior = selecoes[output_path_align]
            if isinstance(anterior, SelecaoDeModeloFalhou):
                raise SelecaoDeModeloFalhou(
                    f"seleção de modelo deste alinhamento já falhou: {anterior}")
            return anterior

        politica = getattr(self, 'model_selection', None) or POLITICA_SELECAO_PADRAO
        if politica == 'nenhuma':
            tool_runs.registrar_selecao_modelo(
                output_path_align, 'desligada_por_configuracao', alinhador=alinhador,
                politica=politica, traducao=dict(MODELO_LEGADO),
                motivo=("model_selection='nenhuma': literais anteriores a M7.3 "
                        "(IQ-TREE/RAxML-NG -m GTR+G, MrBayes lset nst=6 rates=gamma)"))
            selecoes[output_path_align] = MODELO_LEGADO
            return MODELO_LEGADO

        diretorio = os.path.join(self.output_path, 'tmp',
                                 f'modelfinder_{Path(output_path_align).stem}')
        semente = reproducibility_settings(self._reproducibility_kwargs())['random_seed']
        logging.info(f"STEP: seleção de modelo (ModelFinder, BIC) para "
                     f"{os.path.basename(output_path_align)}")
        try:
            resultado = selecionar_modelo(
                alng, diretorio, semente,
                timeout_s=getattr(self, 'model_selection_timeout_s', None))
        except SelecaoDeModeloFalhou as e:
            selecoes[output_path_align] = e
            tool_runs.registrar_selecao_modelo(
                output_path_align, 'falhou', alinhador=alinhador, politica=politica,
                motivo=f"SelecaoDeModeloFalhou: {e}")
            logging.error(f"Seleção de modelo FALHOU para "
                          f"{os.path.basename(output_path_align)}: {e}")
            raise

        traducao = resultado['traducao']
        tool_runs.registrar_selecao_modelo(
            output_path_align, 'concluida', alinhador=alinhador, politica=politica,
            criterio='BIC', conjunto_candidatos=f'--mset mrbayes -mrate {TAXAS_CANDIDATAS}',
            modelo_escolhido=resultado['modelo_escolhido'], traducao=traducao,
            candidatos=resultado['candidatos'], relatorio=resultado['relatorio'],
            comando=resultado['comando'], seed=semente)
        logging.info(f"Modelo escolhido por BIC: {resultado['modelo_escolhido']} -> "
                     f"IQ-TREE {traducao['iqtree']}, RAxML-NG {traducao['raxml-ng']}, "
                     f"MrBayes '{traducao['mrbayes_lset']}'")
        selecoes[output_path_align] = traducao
        return traducao

    def build_tree_iqtree(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando IQ-TREE."""
        logging.info(f"Iniciando construção de árvore com IQ-TREE para {fasta_path}")
        logging.info(f"STEP: Tree Construction with IQ-TREE method.")

        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree,
                              modelo_substituicao=self._modelo_do_alinhamento(
                                  output_path_align, alng, align_method),
                              **self._reproducibility_kwargs())
        tree = builder.iqtree_constructor(alng, output_path_tree)
        self.count_nodes.append(tree.count_terminals())
        return tree

    def build_tree_fasttree(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando FastTree."""
        logging.info(f"Iniciando construção de árvore com FastTree para {fasta_path}")
        logging.info(f"STEP: Tree Construction with FastTree method.")

        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree,
                              **self._reproducibility_kwargs())

        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        tree = builder.fasttree_constructor(alng, output_path_tree)
        self.count_nodes.append(tree.count_terminals())
        return tree

    def build_tree_raxml(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando RAxML-NG."""
        logging.info(f"Iniciando construção de árvore com RAxML-NG para {fasta_path}")
        logging.info(f"STEP: Tree Construction with RAxML-NG method.")

        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree,
                              modelo_substituicao=self._modelo_do_alinhamento(
                                  output_path_align, alng, align_method),
                              **self._reproducibility_kwargs())
        tree = builder.raxml_ng_constructor(alng, output_path_tree)
        self.count_nodes.append(tree.count_terminals())
        return tree

    def build_tree_mrbayes(self, fasta_path, output_path_align, output_path_tree, align_method, output_path_align_html):
        """Constrói árvore usando MrBayes."""
        logging.info(f"Iniciando construção de árvore com MrBayes para {fasta_path}")
        logging.info(f"STEP: Tree Construction with MrBayes method.")

        alng = self._get_alignment(fasta_path, output_path_align, align_method, output_path_align_html)
        builder = TreeBuilder(fasta_path=fasta_path, output_path_tree=output_path_tree,
                              modelo_substituicao=self._modelo_do_alinhamento(
                                  output_path_align, alng, align_method),
                              **self._reproducibility_kwargs(), **self._mrbayes_kwargs())
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