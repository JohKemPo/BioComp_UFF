from typing import List, Dict, Any
import os, sys
import pandas as pd
import numpy as np

from Bio import Phylo

from mlxtend.frequent_patterns import fpmax
from mlxtend.preprocessing import TransactionEncoder

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../..'))

from workflow.utils.treeUtils import (extract_terminais_list,
                                      extract_name_tree,
                                      extract_subtree_info,
                                      dict_to_tree)
from workflow.utils.neo4jProcessing import parse_tree

class SubtreeMiner:
    """
    Classe responsável pela mineração de subárvores frequentes em árvores filogenéticas.

    A classe implementa métodos para agrupar dados, processar subconjuntos exatos de subárvores frequentes, e calcular as maiores subárvores frequentes.
    
    Attributes
    ----------
    support_fpmax : str or float
        Suporte mínimo para mineração de subárvores frequentes.
    mode : str
        Modo de execução do algoritmo.
    matriz_subtree : list
        Matriz de subárvores.
    output_path : str
        Caminho para o diretório de saída.
    max_columns : int
        Número máximo de colunas na matriz de subárvores.
    max_rows : int
        Número máximo de linhas na matriz de subárvores.

    """
    
    def __init__(self, **kwargs) -> None:
        """
        Inicializa a classe SubtreeMiner e define os atributos a partir de kwargs.

        Parameters
        ----------
        **kwargs : dict
            Argumentos passados para definir os atributos da instância.
        """
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        self.matriz_subtree = list()

        

    def miner(self, data):
        processed_group = list()

        if self.mode == "OFST":# Only of the same tree
            grouped_data = self.group_data_by_tree_base(data)
            for base_name, group in grouped_data.items():
                print(f'Analisando grupo de árvores com base "{base_name}"')
                result = self.process_group(group,base_name)
                processed_group.append(result)
                parse_tree(trees=result, path=os.path.join(self.output_path,'outputs'),mode=base_name)

            return processed_group
        else:
            result = self.process_group(data)
            parse_tree(trees=result, path=os.path.join(self.output_path,'outputs'))
            return result
        
    def group_data_by_tree_base(self, data: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Agrupa os dados de subárvores com base no nome base das árvores.

        Parameters
        ----------
        data : List[Dict]
            Lista de dicionários contendo informações sobre árvores e subárvores.

        Returns
        -------
        Dict[str, List[Dict]]
            Dicionário com as subárvores agrupadas por nome base da árvore.
        """
        grouped_data = {}
        for item in data:
            for tree_name in item.keys():
                base_name = '_'.join(tree_name.split('_')[:2]) 
                if base_name not in grouped_data:
                    grouped_data[base_name] = []
                grouped_data[base_name].append(item)
        return grouped_data

    def process_group(self, data: List[Dict], base_name: str = None, sufix: str = None) -> List[Dict]:
        """
        Processa grupos de subárvores e executa a mineração de subárvores frequentes utilizando o algoritmo FPMax.

        Parameters
        ----------
        data : List[Dict]
            Lista de dicionários contendo informações sobre árvores e subárvores.
        base_name : str, optional
            Nome base do grupo de árvores.
        sufix : str, optional
            Sufixo para o nome dos arquivos de saída.

        Returns
        -------
        List[Dict]
            Lista de dicionários com subárvores processadas e atualizadas.
        """
        for metadata in data:
            self.matriz_subtree.append(extract_terminais_list(metadata))
        
        self.max_columns = max(len(row) for row in self.matriz_subtree)
        self.max_rows = len(self.matriz_subtree)

        print(f'\nMATRIZ DE SUBÁRVORES\n')

        for linha in self.matriz_subtree:
            print(linha)
        print(f'\nDIM: "{self.max_rows}x{self.max_columns}" ', '-' * 20)
        print()

        te = TransactionEncoder()
        te_ary = te.fit(self.matriz_subtree).transform(self.matriz_subtree)
        df = pd.DataFrame(te_ary, columns=te.columns_)
        
        all_results_fpmax = pd.DataFrame()
        data_aux = data

        if self.support_fpmax == "auto":
            print(f"Iniciando FPMAX no modo: Variável (0.1 a 0.9)")
            for support in np.arange(0.1, 1.1, 0.1):
                result_fpmax = fpmax(df, min_support=support, use_colnames=True)
                print(f'Resultado FPMAX com suporte {support}:\n{result_fpmax}\n')
                result_fpmax['support'] = support

                if self.mode == "OFST":
                    data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows, base_name)
                else:
                    data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows)
                
                all_results_fpmax = pd.concat([all_results_fpmax, result_fpmax], ignore_index=True)

            all_results_fpmax.to_csv(os.path.join(self.output_path, 'outputs', f'all_results_fpmax.csv'))
        else:
            print(f"Iniciando FPMAX no modo: Fixo em {self.support_fpmax}")
            result_fpmax = fpmax(df, min_support=self.support_fpmax, use_colnames=True)
            print(f'Resultado FPMAX com suporte {self.support_fpmax}:\n{result_fpmax}\n')

            if self.mode == "OFST" and base_name:
                data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows, base_name)
            else:
                data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows)
            
            result_fpmax.to_csv(os.path.join(self.output_path, 'outputs', f'resul_of_fpmax_{self.support_fpmax}.csv'))
        
        return data_aux

    def find_exact_subsets(self, data: List[Dict], result_fpmax: pd.DataFrame, rows: int, base_name: str = None) -> List[Dict]:
        """
        Encontra subconjuntos exatos de subárvores frequentes e exibe informações detalhadas sobre elas.

        Parameters
        ----------
        data : List[Dict]
            Lista de dicionários contendo informações sobre árvores e subárvores.
        result_fpmax : pd.DataFrame
            DataFrame contendo os resultados da mineração de subárvores frequentes.
        rows : int
            Número de linhas (árvores) no dataset.
        base_name : str, optional
            Nome base do grupo de árvores.

        Returns
        -------
        List[Dict]
            Lista de dicionários com subárvores atualizadas com os valores de suporte.
        """
        print('')
        print('- -'*33)
        print(f'                 Subárvores frequentes:                 \n')
        data_dict = result_fpmax.to_dict(orient='records')
        print(f'Número de Árvores que contém Subárvores frequentes no dataset: {len(data_dict)} árvores\n')

        for key in data_dict:
            itemset_list = list(key['itemsets'])
            support = key['support']
            for hash_code in itemset_list:        
                tree_name = extract_name_tree(data, hash_code)
                if not tree_name:
                    break
                print(f'\nÁrvore: {tree_name}\n')

                subtree_names, terminals_lists, metadatas = extract_subtree_info(data, hash_code)
                if not subtree_names: 
                    print(f'\n\nERROR [{hash_code}][{data[0].keys()}]\n\n')
                    break

                for tree_dict in data:
                    for i, subtree_name in enumerate(subtree_names):
                        for tree_name, tree_data in tree_dict.items():
                            if subtree_name in tree_data:
                                subtree_data = tree_data[subtree_name]

                                if 'supports' not in subtree_data:
                                    subtree_data['supports'] = []
                                
                                if support not in subtree_data['supports']:
                                    subtree_data['supports'].append(support)
                                
                                print(f'     - Subárvore: {subtree_name}')
                                print(f'     - hash: {hash_code}')
                                print(f'     - Número de aparições: {int(rows * support)}')
                                print(f'     - Suporte: {"{:.2f}".format(support * 100)}%')
                                print(f'     - Terminais: {terminals_lists[i]}')
                                print(f'     - Subárvore ASCII: \n')
                                subtree_metadata = dict_to_tree(metadatas[i])
                                decode_subtree = Phylo.BaseTree.Tree(subtree_metadata)
                                #Phylo.draw_ascii(decode_subtree)
            
        print('- -'*33)
        return data

    def get_largest_frequent_subtree_by_tree(self, data: List[Dict], tree_name: str) -> Dict[str, Any]:
        """
        Encontra a maior subárvore frequente em uma árvore específica.

        Parameters
        ----------
        data : List[Dict]
            Lista de dicionários contendo informações sobre árvores e subárvores.
        tree_name : str
            Nome da árvore onde a maior subárvore frequente será encontrada.

        Returns
        -------
        Dict[str, Any]
            Dicionário contendo informações sobre a maior subárvore frequente encontrada.
        """
        max_size = 0
        largest_subtree = None

        for tree_dict in data:
            if tree_name in tree_dict:
                subtrees = tree_dict[tree_name]
                for subtree_name, subtree_data in subtrees.items():
                    subtree_size = len(subtree_data['Terminals'])
                    if subtree_size > max_size:
                        max_size = subtree_size
                        largest_subtree = {
                            'tree_name': tree_name,
                            'subtree_name': subtree_name,
                            'subtree_data': subtree_data
                        }

        if largest_subtree:
            print(f"\nA maior subárvore frequente na árvore '{largest_subtree['tree_name']}' é '{largest_subtree['subtree_name']}' com tamanho {max_size}.\n\n")
        else:
            print(f"\nNenhuma subárvore frequente encontrada na árvore '{tree_name}'.\n\n")

        return largest_subtree

    def get_largest_frequent_subtree(self, data: List[Dict]) -> Dict[str, Any]:
        """
        Encontra a maior subárvore frequente em todas as árvores.

        Parameters
        ----------
        data : List[Dict]
            Lista de dicionários contendo informações sobre árvores e subárvores.

        Returns
        -------
        Dict[str, Any]
            Dicionário contendo informações sobre a maior subárvore frequente encontrada em todas as árvores.
        """
        max_size = 0
        largest_subtree = None

        for tree_dict in data:
            for tree_name, subtrees in tree_dict.items():
                for subtree_name, subtree_data in subtrees.items():
                    subtree_size = len(subtree_data['Terminals'])
                    if subtree_size > max_size:
                        max_size = subtree_size
                        largest_subtree = {
                            'tree_name': tree_name,
                            'subtree_name': subtree_name,
                            'subtree_data': subtree_data
                        }

        if largest_subtree:
            print(f"A maior subárvore frequente é '{largest_subtree['subtree_name']}' na árvore '{largest_subtree['tree_name']}' com tamanho {max_size}.")
        else:
            print("Nenhuma subárvore frequente encontrada.")

        return largest_subtree