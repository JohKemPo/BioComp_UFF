from typing import List, Dict, Any
import os, sys, logging, datetime
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
            
        date = datetime.datetime.now()
        logging.basicConfig(level=logging.INFO, 
                    filename=os.path.join(self.output_path,'outputs',f"log_setup_{date.year}_{date.month}_{date.day}.log"),
                    format='%(asctime)s - %(levelname)s - %(message)s')
        
        self.matriz_subtree = list()

        

    def miner(self, data):
        processed_group = list()

        if self.mode == "OFST":# Only of the same tree
            grouped_data = self.group_data_by_tree_base(data)
            for base_name, group in grouped_data.items():
                logging.info(f'Analisando grupo de árvores com base "{base_name}"')
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
            logging.info(f"Iniciando FPMAX no modo: Variável (0.1 a 0.9)")
            for threshold in np.arange(0.1, 1.1, 0.1):
                threshold = round(float(threshold), 2)
                result_fpmax = fpmax(df, min_support=threshold, use_colnames=True)
                print(f'Resultado FPMAX com limiar {threshold}:\n{result_fpmax}\n')
                # D4 — o suporte devolvido pelo mlxtend é a fração de árvores que
                # contém o itemset. O limiar da varredura é outro número e vai em
                # coluna própria: sobrescrever um com o outro fazia o mesmo padrão
                # aparecer como frágil e como robusto nas duas tabelas da UI.
                result_fpmax['min_support_threshold'] = threshold

                if self.mode == "OFST":
                    data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows, base_name)
                else:
                    data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows)

                all_results_fpmax = pd.concat([all_results_fpmax, result_fpmax], ignore_index=True)

            all_results_fpmax = self.consolidate_fpmax_results(all_results_fpmax, self.max_rows)
            all_results_fpmax.to_csv(os.path.join(self.output_path, 'outputs', f'all_results_fpmax.csv'))
        else:
            logging.info(f"Iniciando FPMAX no modo: Fixo em {self.support_fpmax}")
            result_fpmax = fpmax(df, min_support=self.support_fpmax, use_colnames=True)
            print(f'Resultado FPMAX com limiar {self.support_fpmax}:\n{result_fpmax}\n')
            result_fpmax['min_support_threshold'] = float(self.support_fpmax)

            if self.mode == "OFST" and base_name:
                data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows, base_name)
            else:
                data_aux = self.find_exact_subsets(data_aux, result_fpmax, self.max_rows)

            result_fpmax = self.consolidate_fpmax_results(result_fpmax, self.max_rows)
            result_fpmax.to_csv(os.path.join(self.output_path, 'outputs', f'resul_of_fpmax_{self.support_fpmax}.csv'))

        return data_aux

    @staticmethod
    def consolidate_fpmax_results(results: pd.DataFrame, rows: int) -> pd.DataFrame:
        """Uma linha por itemset, com o suporte real e a faixa de limiares (D4).

        A varredura por limiar devolve o mesmo itemset em vários limiares. Antes,
        cada repetição virava uma linha do CSV com um "suporte" diferente — em
        VARV-49, 7 de 7 itemsets distintos apareciam com mais de um valor, e 2
        deles nas **duas** tabelas da Deep Analysis ao mesmo tempo. O suporte
        real é propriedade do itemset e não do limiar: é o mesmo em toda
        repetição, então deduplicar não perde informação.

        Colunas:
          `support`                 fração de árvores que contêm o itemset (mlxtend)
          `min_support_threshold`   menor limiar da varredura que o devolveu
          `max_support_threshold`   maior limiar da varredura que o devolveu
          `n_trees`                 número de árvores que o contêm
        """
        if results.empty:
            return pd.DataFrame(columns=['itemsets', 'support', 'min_support_threshold',
                                         'max_support_threshold', 'n_trees'])

        # `sort=False`: a chave é um frozenset, cuja ordenação natural é a de
        # subconjunto (parcial) e não serve para ordenar linhas. A ordem final é
        # imposta logo abaixo, explicitamente.
        consolidado = (results
                       .groupby('itemsets', as_index=False, sort=False)
                       .agg(support=('support', 'max'),
                            support_minimo_visto=('support', 'min'),
                            min_support_threshold=('min_support_threshold', 'min'),
                            max_support_threshold=('min_support_threshold', 'max')))

        divergentes = consolidado[
            (consolidado['support'] - consolidado['support_minimo_visto']).abs() > 1e-9]
        if not divergentes.empty:
            # Não deveria acontecer: o suporte não depende do limiar. Se acontecer,
            # é sintoma de matriz alterada entre limiares — vale ruído no log.
            logging.warning("FPMAX: suporte divergente para o mesmo itemset em "
                            f"limiares diferentes: {divergentes['itemsets'].tolist()}")
        consolidado = consolidado.drop(columns=['support_minimo_visto'])

        consolidado['n_trees'] = (consolidado['support'] * rows).round().astype(int)
        consolidado['tamanho'] = consolidado['itemsets'].map(len)
        consolidado = (consolidado
                       .sort_values(by=['support', 'tamanho', 'min_support_threshold'],
                                    ascending=[False, False, True], kind='stable')
                       .drop(columns=['tamanho'])
                       .reset_index(drop=True))
        return consolidado

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
        logging.info(f'Número de Árvores que contém Subárvores frequentes no dataset: {len(data_dict)} árvores\n')

        for key in data_dict:
            itemset_list = list(key['itemsets'])
            # Suporte real do itemset (fração de árvores que o contêm), não o
            # limiar da varredura — ver D4 e `consolidate_fpmax_results`.
            support = key['support']
            for hash_code in itemset_list:        
                tree_name = extract_name_tree(data, hash_code)
                if not tree_name:
                    break
                print(f'\nÁrvore: {tree_name}\n')

                subtree_names, terminals_lists, metadatas = extract_subtree_info(data, hash_code)
                if not subtree_names: 
                    logging.error(f'\n\nERROR [{hash_code}][{data[0].keys()}]\n\n')
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