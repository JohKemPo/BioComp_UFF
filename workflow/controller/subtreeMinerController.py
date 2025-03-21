import time, logging, datetime
import os, sys
import pandas as pd
import numpy as np
from typing import Any, Dict, List

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../..'))

from workflow.utils.messages import Messages
from workflow.utils.neo4jProcessing import parse_tree
from workflow.subtree_mining.miner import SubtreeMiner

class SubtreeMinerController:
    """
    Controlador responsável pela mineração de subárvores frequentes em árvores filogenéticas.

    A classe coordena a mineração de subárvores, agrupando os dados e processando as subárvores frequentes.
    
    Attributes
    ----------
    mode : str
        Modo de execução da mineração (e.g., "OFST" para mineração dentro do mesmo grupo de árvores).
    output_path : str
        Caminho para o diretório de saída.
    matriz_subtree : List
        Matriz de subárvores processadas.
    msg : Messages
        Instância da classe Messages para exibir mensagens de status.
    start : float
        Tempo de início da execução para monitoramento do tempo de execução.
    """

    def __init__(self, **kwargs: Any) -> None:
        """
        Inicializa a classe SubtreeMinerController, configurando os atributos com base nas opções fornecidas via kwargs.

        Parameters
        ----------
        **kwargs : dict
            Argumentos passados para configurar os atributos da classe, como o modo de execução e caminho de saída.
        """
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        date = datetime.datetime.now()
        logging.basicConfig(level=logging.INFO, 
                    filename=os.path.join(self.output_path,'outputs',f"log_setup_{date.year}_{date.month}_{date.day}.log"),
                    format='%(asctime)s - %(levelname)s - %(message)s')

        self.msg = Messages(logPath=os.path.join(self.output_path,'outputs'))
        self.start = time.time()
        self.matriz_subtree = []

        print(self.msg.init_message())
        print('       - MINERAÇÃO DE SUBÁRVORES FREQUENTES -       \n')
        print('------------------------------------------------------')
        logging.info("SubtreeMinerController inicializado com sucesso.")
        
    def group_data_by_tree_base(self, data: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
        """
        Agrupa os dados de subárvores com base no nome base das árvores.

        Parameters
        ----------
        data : List[Dict[str, Any]]
            Lista de dicionários contendo informações sobre árvores e subárvores.

        Returns
        -------
        Dict[str, List[Dict[str, Any]]]
            Dicionário com as subárvores agrupadas por nome base da árvore.
        """
        grouped_data = {}
        logging.debug("Iniciando agrupamento de dados por base de árvore.")
        for item in data:
            for tree_name in item.keys():
                base_name = '_'.join(tree_name.split('_')[:2])
                if base_name not in grouped_data:
                    grouped_data[base_name] = []
                    logging.debug(f"Criando novo grupo para a base: {base_name}")
                grouped_data[base_name].append(item)
                logging.debug(f"Adicionado item ao grupo {base_name}: {item}")
        logging.info("Agrupamento de dados concluído.")
        return grouped_data
    
    def miner(self, data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Executa o processo de mineração de subárvores frequentes.

        Se o modo de execução for "OFST", os dados serão agrupados por base de árvores e processados individualmente.
        Caso contrário, todo o conjunto de dados será processado de uma vez.

        Parameters
        ----------
        data : List[Dict[str, Any]]
            Lista de dicionários contendo as árvores e subárvores que serão processadas.

        Returns
        -------
        List[Dict[str, Any]]
            Lista contendo as subárvores processadas após a mineração.
        """
        processed_group = []
        logging.info("Iniciando mineração de subárvores frequentes.")
        try:
            if self.mode == "OFST":  # Only of the same tree
                logging.info("Modo OFST detectado: agrupando dados por base de árvore.")
                # Utiliza o método de agrupamento da classe SubtreeMiner
                grouped_data = SubtreeMiner.group_data_by_tree_base(data)
                logging.debug(f"Dados agrupados: {list(grouped_data.keys())}")
                for base_name, group in grouped_data.items():
                    self.matriz_subtree = []
                    logging.info(f'Analisando grupo de árvores com base "{base_name}"')
                    result = SubtreeMiner.process_group(group, base_name)
                    processed_group.append(result)
                    out_path = os.path.join(self.output_path, 'outputs')
                    parse_tree(trees=result, path=out_path, mode=base_name)
                    logging.info(f"Grupo {base_name} minerado e árvore salva em {out_path}.")
                return processed_group
            else:
                logging.info("Modo padrão detectado: processando conjunto completo de dados.")
                result = self.process_group(data)
                out_path = os.path.join(self.output_path, 'outputs')
                parse_tree(trees=result, path=out_path)
                logging.info(f"Mineração concluída e árvore salva em {out_path}.")
                return result
        except Exception as e:
            logging.error(f"Erro durante a mineração de subárvores: {e}", exc_info=True)
            raise