import os, sys, time, json, logging, datetime
from tqdm import tqdm
import pandas as pd
from typing import Any, Dict, List, Union

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../..'))

from workflow.utils.messages import Messages
from workflow.subtree_construction.builder import SubtreeBuilder
from workflow.controller.subtreeMinerController import SubtreeMiner
from workflow.utils.metrics import process_histogram_frequence


class SubtreeBuilderController:
    """
    Controlador responsável pela construção e mineração de subárvores em árvores filogenéticas.

    A classe coordena a construção das subárvores, o armazenamento dos dados gerados,
    e a mineração das subárvores frequentes, se configurada para tal.

    Attributes
    ----------
    output_path : str
        Caminho para o diretório de saída.
    input_path : str
        Caminho para o diretório de entrada, onde as árvores estão localizadas.
    save_metadata : bool
        Indica se os metadados das subárvores devem ser salvos em um arquivo JSON.
    subtree_miner : bool
        Indica se a mineração de subárvores deve ser realizada após a construção.
    subtree_miner_configs : Dict[str, Any]
        Configurações para o minerador de subárvores.
    count_trees : int
        Número de árvores processadas.
    count_subtrees : int
        Número de subárvores construídas.
    raw_data : List[Dict[str, Any]]
        Lista de dados brutos das árvores e subárvores.
    matrix_subtree : List[List[Any]]
        Matriz contendo as subárvores extraídas.
    list_times : List[float]
        Lista dos tempos de execução por árvore.
    count_nodes : List[int]
        Lista com o número de nós em cada árvore processada.
    """

    def __init__(self, **kwargs) -> None:
        """
        Inicializa a classe SubtreeBuilderController, configurando os atributos com base nas opções fornecidas via kwargs.

        Parameters
        ----------
        **kwargs : dict
            Argumentos passados para configurar os atributos da classe, como caminhos de entrada/saída e opções de mineração.
        """
        
        
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        date = datetime.datetime.now()
        logging.basicConfig(level=logging.INFO, 
                    filename=os.path.join(self.output_path,'outputs',f"log_setup_{date.year}_{date.month}_{date.day}.log"),
                    format='%(asctime)s - %(levelname)s - %(message)s')
            
        self.start = time.time()
        self.msg = Messages(logPath=os.path.join(self.output_path,'outputs'))
        print(self.msg.init_message())
        print('             - CONTRUÇÃO DE SUBÁRVORES -             \n')
        print('------------------------------------------------------')
        logging.info("Inicializando SubtreeBuilderController.")
        
        os.makedirs(self.output_path, exist_ok=True)
        dirs = ['Subtrees', 'outputs']

        for dir in dirs:
            path = os.path.join(self.output_path, dir)
            os.makedirs(path, exist_ok=True)
        
        self.files = os.listdir(self.input_path)
        self.count_trees = len(self.files)
        self.count_nodes = []
        self.list_times = []
        self.raw_data = []
        self.matrix_subtree = []
        self.subtree_kwargs = kwargs
        self.count_subtrees = 0
        logging.info(f"{self.count_trees} arquivo(s) encontrado(s) em {self.input_path}.")

    def __call__(self) -> List[Dict[str, Any]]:
        """
        Executa o processo de construção das subárvores e, opcionalmente, mineração das subárvores frequentes.

        Este método é responsável por iterar sobre os arquivos de entrada, construir subárvores usando o SubtreeBuilder,
        armazenar os resultados em um arquivo JSON, e executar a mineração de subárvores frequentes, se configurado para tal.

        Returns
        -------
        List[Dict[str, Any]]
            Lista contendo os dados das subárvores processadas e mineradas.
        """
        logging.info("Início da construção das subárvores.")
        logging.info("STEP: Construction of Subtrees.")
        for name in tqdm(self.files, desc="Gerando subárvores....", ascii="░▒█"):
            path = os.path.join(self.input_path, name)
            logging.debug(f"Iniciando construção da subárvore para o arquivo: {name}")
            self.raw_data.append(self.builder(path, name))
            logging.debug(f"Subárvore construída para: {name}")

        json_result = json.dumps(self.raw_data, indent=2)
        logging.info("Construção das subárvores concluída.")
        
        if self.save_metadata:
            metadata_path = os.path.join(self.output_path, 'outputs', 'metadata.json')
            try:
                with open(metadata_path, 'w') as output_file:
                    output_file.write(json_result)
                logging.info(f"Metadados salvos em JSON: {metadata_path}")
            except Exception as e:
                logging.error(f"Erro ao salvar metadados em JSON: {e}", exc_info=True)


        self.msg.resume_subtree(start=self.start,
                                num_trees=self.count_trees,
                                output_format=self.output_format,
                                num_subtrees=self.count_subtrees)
        data = self.raw_data
        
        if self.subtree_miner:
            logging.info("Iniciando mineração de subárvores frequentes.")
            logging.info("STEP: Frequent subtree mining.")
            try:
                miner = SubtreeMiner(**self.subtree_miner_configs)
                data = miner.miner(data=self.raw_data)
                json_result = json.dumps(data, indent=2)
                df = pd.DataFrame(data)
                csv_path = os.path.join(self.output_path, 'outputs', 'metadata.csv')
                df.to_csv(csv_path)
                logging.info(f"Metadados minerados salvos em CSV: {csv_path}")
                
                json_path = os.path.join(self.output_path, 'outputs', 'metadata.json')
                with open(json_path, 'w') as output_file:
                    output_file.write(json_result)
                logging.info(f"Metadados minerados salvos em JSON: {json_path}")
            except Exception as e:
                logging.error(f"Erro na mineração de subárvores: {e}", exc_info=True)
                raise
        
        try:
            plot_path = os.path.join(self.output_path, 'outputs', 'Plots')
            process_histogram_frequence(data, plot_path)
            logging.info(f"Histograma de frequência processado e salvo em: {plot_path}")
        except Exception as e:
            logging.error(f"Erro ao processar histograma de frequência: {e}", exc_info=True)
            raise
        logging.info(f"STEP: Completed successfully!")
        
        return data


    def builder(self, path: str, name: str) -> Dict[str, Any]:
        """
        Constrói as subárvores para um arquivo de entrada específico.

        Parameters
        ----------
        path : str
            Caminho completo do arquivo de entrada contendo a árvore.
        name : str
            Nome do arquivo de entrada.

        Returns
        -------
        Dict[str, Any]
            Dados brutos das subárvores construídas.
        """
        logging.info(f"Iniciando a construção da subárvore para: {name}")
        self.count_subtrees = 0
        try:
            builder = SubtreeBuilder(**self.subtree_kwargs)
            rawdata = builder.subtree_constructor(path, name)
            self.count_subtrees += builder.count_subtrees
            logging.info(f"Subárvore para {name} construída com sucesso. Total de subárvores construídas: {self.count_subtrees}")
        except Exception as e:
            logging.error(f"Erro na construção da subárvore para {name}: {e}", exc_info=True)
            raise
        return rawdata