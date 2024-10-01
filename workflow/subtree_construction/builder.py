from Bio import Phylo

import os, sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../..'))

from workflow.utils.messages import Messages
from workflow.utils.treeUtils import (tree_to_dict, 
                                      calculate_tree_hash,
                                      encode_list_to_int)

class SubtreeBuilder:
    """
    Classe responsável pela construção e manipulação de subárvores filogenéticas.

    Esta classe fornece métodos para a construção de subárvores a partir de uma árvore filogenética,
    cálculo e manipulação de hashes de subárvores, e exportação dessas subárvores para arquivos.

    Attributes
    ----------
    count_subtrees : int
        Contador para acompanhar o número de subárvores construídas.
    """
    
    def __init__(self, **kwargs):
        """
        Inicializa a instância da classe SubtreeBuilder.
        ----------
        Esta função atribui os valores dos argumentos fornecidos via kwargs como atributos
        da instância. Também inicializa uma instância de Messages e o contador de subárvores.

        Parameters
        ----------
        **kwargs : dict
            Dicionário contendo pares chave-valor que serão configurados como atributos da instância.

        Return
        ------
        None
        """
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        self.count_subtrees = 0

    def subtree_constructor(self, path: str , name: str ) -> dict:
        """
        Constrói subárvores a partir de uma árvore filogenética e exporta os resultados.

        A função percorre a árvore principal e, para cada clado, constrói subárvores.
        Cada subárvore com mais de um terminal é exportada em um arquivo no formato especificado.
        Informações adicionais sobre as subárvores e seus hashes são impressas se resume_infos for True.

        Parameters
        ----------
        path : str
            Caminho para o arquivo contendo a árvore filogenética.
        name : str
            Nome do arquivo de saída para as subárvores.

        Return
        ------
        dict
            Um dicionário contendo informações sobre todas as subárvores geradas, incluindo seus terminais e metadados.
        """
        tree = Phylo.read(path, self.input_format)
        name_tree = str(name.rsplit(".", 1)[0])
        tree_hash = calculate_tree_hash(tree)['terminal_hash']

        if self.resume_infos:
            print('======================================================')        
            Messages.print_tree_info(name=name_tree, hash_value=tree_hash)

        dict_tree_terminals_hash = dict()
        dict_aux = dict()
        result_dict = dict()

        for clade in tree.find_clades():
            hash_list = list()
            subtree_list_termials = list()

            subtree = Phylo.BaseTree.Tree(clade)
            name_subtree = f'{name_tree}_{clade.name}'

            if clade.is_terminal():
                clade_hash = calculate_tree_hash(clade.name)['terminal_hash']
                dict_tree_terminals_hash[clade.name] = clade_hash
            
            if subtree.count_terminals() > 1:
                self.count_subtrees += 1 
                if self.resume_infos:
                    Messages.print_subtree_info(subtree=subtree, name=name_subtree)

                filepath_out = os.path.join(self.output_path, 'Subtrees', f'{name_tree}_{clade.name}.{self.output_format}')
                Phylo.write(subtree, filepath_out, self.output_format)
                            
                for subtree_clade in subtree.find_clades():
                    name_terminal = subtree_clade.name
                    if subtree_clade.is_terminal():
                        hash_result = calculate_tree_hash(name_terminal)
                        hash_list.append(hash_result)
                        subtree_list_termials.append(hash_result['terminal_hash'])

                        if self.resume_infos:
                            Messages.print_subtree_clade_info(clade=subtree_clade, decode=self._decode_tree_hash(hash_result), hash_dict=hash_result)

                dict_aux[name_subtree] = {
                    'Terminals': subtree_list_termials,
                    'List_terminals_hash': encode_list_to_int(subtree_list_termials),
                    'data_terminals': hash_list,
                    'metadata': tree_to_dict(clade)
                }

                if self.resume_infos:
                    Messages.print_list_clade(list_clades=subtree_list_termials, hash_value=encode_list_to_int(subtree_list_termials))

        result_dict[name_tree] = dict_aux 

        if self.resume_infos:
            Messages.print_tree_clade_info(name=name_tree, dict_terminals=dict_tree_terminals_hash)
        
        return result_dict





