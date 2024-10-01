import os, sys, time, json
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
            
        self.start = time.time()
        self.msg = Messages(logPath=os.path.join(self.output_path,'outputs'))
        print(self.msg.init_message())
        print('             - CONTRUÇÃO DE SUBÁRVORES -             \n')
        print('------------------------------------------------------')

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
        for name in tqdm(self.files, desc="Gerando subárvores....", ascii="░▒█"):
            path = os.path.join(self.input_path, name)
            self.raw_data.append(self.builder(path, name))

        json_result = json.dumps(self.raw_data, indent=2)

        if self.save_metadata:
            with open(os.path.join(self.output_path, 'outputs', 'metadata.json'), 'w') as output_file:
                output_file.write(json_result)

        self.msg.resume_subtree(start=self.start,
                                num_trees=self.count_trees,
                                output_format=self.output_format,
                                num_subtrees=self.count_subtrees)
        
        if self.subtree_miner:
            miner = SubtreeMiner(**self.subtree_miner_configs)
            data = miner.miner(data=self.raw_data)
            json_result = json.dumps(data, indent=2)

            df = pd.DataFrame(data)
            df.to_csv(os.path.join(self.output_path, 'outputs', 'metadata.csv'))

            path = os.path.join(self.output_path, 'outputs', 'metadata.json')
            with open(path, 'w') as output_file:
                output_file.write(json_result)
        
        process_histogram_frequence(data,os.path.join(self.output_path, 'outputs','Plots'))
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
        builder = SubtreeBuilder(**self.subtree_kwargs)
        self.rawdata = builder.subtree_constructor(path, name)
        self.count_subtrees += builder.count_subtrees

        return self.rawdata