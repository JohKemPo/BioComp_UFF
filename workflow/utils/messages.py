import time
import resource
from Bio import Phylo
import logging
import datetime

class Messages:
    """
    Classe responsável por exibir mensagens e informações durante a execução do script de identificação de subárvores frequentes.
    
    Fornece funcionalidades para resumir a execução, exibir informações sobre árvores e subárvores, e monitorar o uso de recursos do sistema.
    """
    
    def __init__(self,
                 logPath: str):
        """
        Inicializa uma instância da classe Messages.
        """
        super().__init__()
        self.logPath = logPath
        date = datetime.datetime.now()

        logging.basicConfig(level=logging.INFO, 
                    filename=f"{logPath}/log_setup_{date.year}_{date.month}_{date.day}.log",
                    format='%(asctime)s - %(levelname)s - %(message)s')

    def init_message(self) -> str:
        """
        Retorna a mensagem inicial exibida ao iniciar o sistema.

        Returns
        -------
        str
            A mensagem de introdução do sistema.
        """
        message = """
------------------------------------------------------
                     MFSt.P                                                   
                    -------                                                   

         Identificação em Paralelo de Subárvores                           
    Frequentes em Conjuntos de Árvores Filogenéticas                              
------------------------------------------------------
            """
        return message

    def resume_tree(self,
                    start: float,
                    sum_time: list,
                    num_trees: int,
                    output_format: str,
                    n_nodes: list,
                    method: str) -> None:
        """
        Exibe o resumo da execução da construção das árvores.

        Parameters
        ----------
        start : float
            O tempo de início da execução.
        sum_time : list
            Lista contendo os tempos de execução de cada árvore.
        num_trees : int
            Número de árvores geradas.
        output_format : str
            O formato de saída das árvores.
        n_nodes : list
            Lista com o número de terminais em cada árvore.
        method : str
            O método de construção das árvores utilizado.

        Returns
        -------
        None
        """
        print('\n------------------------------------------------------\n')
        print(f'Tempo de execução: {"{:.2f}".format(time.time() - start)}s')
        print(f'Tempo médio de execução por árvore: {"{:.2f}".format(sum(sum_time) / num_trees)}s')
        print(f'Número de árvores do tipo {output_format} geradas: {num_trees}')
        print(f'Número médio de terminais por árvore: {int(sum(n_nodes) / len(n_nodes))}')
        print(f'Método utilizado: {method}')

        logging.info('CONSTRUÇÃO DE ÁRVORES')
        logging.info(f'    Tempo de execução: {"{:.2f}".format(time.time() - start)}s')
        logging.info(f'    Tempo médio de execução por árvore: {"{:.2f}".format(sum(sum_time) / num_trees)}s')
        logging.info(f'    Número de árvores do tipo {output_format} geradas: {num_trees}')
        logging.info(f'    Número médio de terminais por árvore: {int(sum(n_nodes) / len(n_nodes))}')
        logging.info(f'    Método utilizado: {method}')

    def resume_subtree(self,
                       start: float,
                       num_trees: int,
                       output_format: str,
                       num_subtrees: int) -> None:
        """
        Exibe o resumo da execução de análise das subárvores.

        Parameters
        ----------
        start : float
            O tempo de início da execução.
        num_trees : int
            Número de árvores analisadas.
        output_format : str
            O formato de saída das árvores.
        num_subtrees : int
            Número de subárvores geradas e analisadas.

        Returns
        -------
        None
        """
        print('\n------------------------------------------------------\n')
        print(f'Tempo de execução: {"{:.2f}".format(time.time() - start)}s')
        print(f'Árvores {output_format} analisadas com sucesso: {num_trees}')
        print(f'Número de subárvores {output_format} produzidas e analisadas: {num_subtrees}\n')

        logging.info('CONSTRUÇÃO DE SUBÁRVORES')
        logging.info(f'    Tempo de execução: {"{:.2f}".format(time.time() - start)}s')
        logging.info(f'    Árvores {output_format} analisadas com sucesso: {num_trees}')
        logging.info(f'    Número de subárvores {output_format} produzidas e analisadas: {num_subtrees}\n')

    def get_resource_usage(self) -> dict:
        """
        Retorna um dicionário contendo o uso de recursos do sistema durante a execução do script.

        Returns
        -------
        dict
            Dicionário contendo o tempo de CPU e a memória residente utilizada pelo processo.
        """
        rusage = resource.getrusage(resource.RUSAGE_SELF)
        return {
            'Tempo de CPU (segundos)': rusage.ru_utime + rusage.ru_stime,
            'Memória Residente (kB)': rusage.ru_maxrss
        }

    def save_resource_usage_summary(self,
                                    resource_usage: dict) -> None:
        """
        Exibe o resumo do uso de recursos no console.

        Parameters
        ----------
        resource_usage : dict
            Dicionário contendo o uso de recursos do sistema.

        Returns
        -------
        None
        """
        for key, value in resource_usage.items():
            print(f"{key}: {value}")

    def print_tree_info(self, name: str, hash_value: int) -> None:
        """
        Exibe informações sobre uma árvore e seu valor hash.

        Parameters
        ----------
        name : str
            O nome da árvore.
        hash_value : int
            O valor hash da árvore.

        Returns
        -------
        None
        """
        print(f'\nTREE: HASH VALUE\n{name}: {hash_value}\n')

    def print_subtree_info(self, subtree: Phylo.BaseTree.Tree, name: str = None) -> None:
        """
        Exibe informações sobre uma subárvore e imprime sua representação ASCII.

        Parameters
        ----------
        subtree : Phylo.BaseTree.Tree
            A subárvore a ser exibida.
        name : str, optional
            O nome da subárvore.

        Returns
        -------
        None
        """
        print(f'\nSUBTREE:\n{name}\n')
        Phylo.draw_ascii(subtree)
        print('\n')

    def print_subtree_clade_info(self, clade: Phylo.BaseTree.Clade, decode: str, hash_dict: dict) -> None:
        """
        Exibe informações sobre um clado de uma subárvore e seu valor hash.

        Parameters
        ----------
        clade : Phylo.BaseTree.Clade
            O clado da subárvore.
        decode : str
            A representação decodificada do clado.
        hash_dict : dict
            O dicionário contendo o hash da subárvore.

        Returns
        -------
        None
        """
        print(f'subtree clade: {clade}: {hash_dict["terminal_hash"]}')
        print(f'Decode: {decode}')

    def print_list_clade(self, list_clades: list, hash_value: int) -> None:
        """
        Exibe a lista de terminais de um clado e seu valor hash.

        Parameters
        ----------
        list_clades : list
            Lista contendo os terminais do clado.
        hash_value : int
            O valor hash do clado.

        Returns
        -------
        None
        """
        print(f'\nTerminais: {list_clades}:{hash_value}\n')

    def print_tree_clade_info(self, dict_terminals: dict, name: str) -> None:
        """
        Exibe informações sobre os clados de uma árvore.

        Parameters
        ----------
        dict_terminals : dict
            Dicionário contendo as informações sobre os clados da árvore.
        name : str
            O nome da árvore.

        Returns
        -------
        None
        """
        print(f'CLADES DA ARVORE {name}:\n')       
        for key, value in dict_terminals.items():
            print(f'{key}: {value}')