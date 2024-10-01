from Bio import Phylo

#Types
from Bio.Phylo.BaseTree import Clade, Tree
from typing import Dict, List, Any, Optional, Tuple

Matrix = List[List[Any]]

import hashlib


def tree_to_dict(clade: Clade) -> Dict:
    """
        Converte um clado de uma árvore filogenética em um dicionário.

        A função percorre recursivamente os clados filhos, criando uma estrutura
        hierárquica representada como um dicionário.

        Parameters
        ----------
        clade : Bio.Phylo.BaseTree.Clade
            O clado a ser convertido em dicionário.

        Return
        ------
        Dict
            Dicionário representando o clado e seus filhos.
    """
    node_dict = {
        "name": clade.name,
        "branch_length": clade.branch_length if clade.branch_length is not None else 0,
        "children": [tree_to_dict(child) for child in clade.clades],
    }
    return node_dict

def dict_to_tree(node_dict: Dict) -> Clade:
    """
        Converte um dicionário de clado de volta em um objeto Clade da Bio.Phylo.

        A função percorre o dicionário de forma recursiva, criando objetos Clade
        e seus filhos a partir da estrutura hierárquica do dicionário.

        Parameters
        ----------
        node_dict : dict
            O dicionário representando a estrutura do clado.

        Return
        ------
        Bio.Phylo.BaseTree.Clade
            Objeto Clade reconstruído a partir do dicionário.
    """
    clade = Phylo.BaseTree.Clade()
    clade.name = node_dict["name"]
    clade.branch_length = node_dict["branch_length"]

    for child_dict in node_dict["children"]:
        child_clade = dict_to_tree(child_dict)
        clade.clades.append(child_clade)

    return clade

def extract_terminais_list(data: Dict) -> List:
    """
    Extrai uma lista de hashes terminais de subárvores a partir de uma estrutura de dados hierárquica.

    A função percorre o dicionário fornecido, que contém informações sobre árvores e subárvores,
    e agrega os valores de 'List_terminals_hash' de todas as subárvores em uma lista.

    Parameters
    ----------
    data : dict
        Dicionário contendo informações sobre as árvores e subárvores, onde cada subárvore
        possui uma lista de hashes terminais armazenada na chave 'List_terminals_hash'.

    Return
    ------
    list
        Lista contendo todos os valores de 'List_terminals_hash' extraídos das subárvores.
    """
    all_terminais = []
    
    for _, tree_info in data.items():
        for _, subtree_info in tree_info.items():
            all_terminais.append(subtree_info['List_terminals_hash'])
                

    return all_terminais

def preencher_matriz(matriz: Matrix, valor_preenchimento: Any, max_columns: int) -> Matrix:
    """
    Preenche uma matriz com um valor específico até que todas as linhas tenham o número máximo de colunas.

    Esta função percorre cada linha da matriz fornecida e adiciona o valor de preenchimento às linhas
    que contêm menos colunas do que o valor especificado por `max_columns`.

    Parameters
    ----------
    matriz : Matrix
        A matriz que será preenchida. Cada elemento da matriz é uma lista representando uma linha.
    valor_preenchimento : Any
        O valor que será utilizado para preencher as células vazias nas linhas da matriz.
    max_columns : int
        O número máximo de colunas que cada linha da matriz deve ter após o preenchimento.

    Return
    ------
    Matrix
        A matriz original, com todas as suas linhas preenchidas até o número máximo de colunas.
    """
    for row in matriz:
        while len(row) < max_columns:
            row.append(valor_preenchimento)

    return matriz


def extract_name_tree(data: List[Dict], hash) -> str:
    """
    Extrai o nome da árvore correspondente a um hash específico.

    Esta função percorre uma lista de dicionários contendo informações sobre árvores e subárvores,
    buscando o nome da árvore que contém uma subárvore cujo valor de 'List_terminals_hash' corresponde
    ao hash fornecido.

    Parameters
    ----------
    data : List[Dict]
        Lista de dicionários contendo informações sobre árvores e subárvores.
    hash : int
        O valor do hash que será utilizado para identificar a árvore correspondente.

    Return
    ------
    str
        O nome da árvore que contém a subárvore correspondente ao hash fornecido. Se a árvore não for
        encontrada, retorna uma string vazia.
    """
    tree_name = str()
    
    for info in data:
        for tree_name, subtree_info in info.items():
            for _, subtree_data in subtree_info.items():
                if "List_terminals_hash" in subtree_data and subtree_data["List_terminals_hash"] == hash:
                    return tree_name



def find_exact_subsets(data, result_fpmax, rows, count_trees) -> None:
    """
    Encontra e exibe subconjuntos exatos de subárvores frequentes a partir dos resultados de mineração.

    Esta função analisa os resultados da mineração de subárvores frequentes, identificando e imprimindo
    as subárvores que aparecem em múltiplas árvores no dataset. Para cada subárvore frequente, são exibidos
    o nome da árvore original, o número de aparições, o suporte percentual, a lista de terminais e uma
    representação ASCII da subárvore.

    Parameters
    ----------
    data : list
        Estrutura contendo as árvores e suas subárvores, utilizadas para recuperar as informações correspondentes.
    result_fpmax : pandas.DataFrame
        DataFrame contendo os resultados da mineração de subárvores frequentes com FPMax.
    rows : int
        O número total de linhas no dataset original, usado para calcular o número de aparições das subárvores.
    count_trees : int
        O número total de árvores no dataset, usado para calcular o suporte das subárvores.

    Return
    ------
    None
        A função não retorna nada, apenas imprime informações detalhadas sobre as subárvores frequentes.
    """
    print('')
    print('- -'*33)
    print(f'                                       Subárvores frequentes:                                      \n')
    data_dict = result_fpmax.to_dict(orient='records')
    
    print(f'Número de Árvores que contém Subárvores frequêntes no dataset: {len(data_dict)} árvores\n')
    
    for key in data_dict:
        itemset_list = list(key['itemsets'])
        suport = key['support']
        for tree in itemset_list:
            tree_name = extract_name_tree(data, tree)
            print(f'\nÁrvore: { tree_name }\n')
            
            for subtree in itemset_list:
                subtree_name, terminals_list, metadata = extract_subtree_info(data, subtree)

                print(f'     - Subárvore: { subtree_name }')
                print(f'     - hash: { subtree }')
                print(f'     - Número de aparições: {int(rows * suport)} de {count_trees} ')
                print(f'     - Suporte: {"{:.2f}".format(suport*100)}%')
                print(f'     - Terminais: { terminals_list }')
                print(f'     - Subárvore ASCII: \n')
                subtree_metadata = dict_to_tree(metadata)
                decode_subtree = Phylo.BaseTree.Tree(subtree_metadata)
                Phylo.draw_ascii(decode_subtree)
            break
    print('- -'*33)


def extract_subtree_info(data: List[Dict], hash: int) -> Tuple[List, List, List]:
    """
    Extrai informações de subárvores correspondentes a um hash específico.

    Esta função percorre uma lista de dicionários que contém informações de subárvores, 
    buscando por subárvores cujo valor de 'List_terminals_hash' corresponde ao hash fornecido.
    Para cada subárvore correspondente, a função extrai o nome da subárvore, a lista de terminais 
    (em formato Newick) e os metadados associados.

    Parameters
    ----------
    data : List[Dict]
        Lista de dicionários contendo informações sobre árvores e subárvores.
    hash : int
        O valor do hash que será utilizado para identificar as subárvores desejadas.

    Return
    ------
    Tuple[str, List, List]
        Uma tupla contendo três elementos:
        - lista de nomes das subárvores correspondentes ao hash.
        - lista de listas contendo as sequências terminais em formato Newick de cada subárvore.
        - lista de metadados associados a cada subárvore.
    """
    subtree_names = list()
    terminals_lists = list()
    metadatas = list()
    for info in data:
        for _, subtree_info in info.items():
            for subtree_name, subtree_data in subtree_info.items():
                if subtree_data["List_terminals_hash"] == hash:
                    terminals_list = [item["newick"] for item in subtree_data["data_terminals"]]
                    metadata = subtree_data['metadata']
                    subtree_names.append(subtree_name)
                    terminals_lists.append(terminals_list)
                    metadatas.append(metadata)
    return subtree_names, terminals_lists, metadatas

def decode_int_to_list(encoded_int: int, original_list: List) -> Optional[list]:
    """
    Decodifica um valor inteiro para validar sua correspondência com a lista original.

    A função verifica se o valor inteiro codificado corresponde ao hash da lista
    original, garantindo sua integridade.

    Parameters
    ----------
    encoded_int : int
        O valor inteiro codificado para validação.
    original_list : list
        A lista original a ser comparada com o valor inteiro.

    Return
    ------
    list or None
        A lista original se o hash for válido, caso contrário, None.
    """
    lst_str = ",".join(map(str, original_list))
    hash_object = hashlib.md5(lst_str.encode())
    if int(hash_object.hexdigest()[:4], 16) == encoded_int:
        return original_list
    else:
        return None

def encode_list_to_int(lst: List) -> int:
    """
    Codifica uma lista de valores em um único inteiro usando hash MD5.

    A função gera uma string da lista, calcula seu hash MD5, e retorna
    os primeiros 4 caracteres como um número inteiro.

    Parameters
    ----------
    lst : list
        A lista de valores para codificação.

    Return
    ------
    int
        Valor inteiro resultante da codificação da lista.
    """
    lst_str = str(lst)
    hash_object = hashlib.md5(lst_str.encode())
    return int(hash_object.hexdigest()[:4], 16)

def decode_tree_hash(encoded_data: Dict) -> Optional[str]:
    """
    Decodifica e valida o hash de uma subárvore.

    A função verifica se o hash armazenado no dicionário corresponde ao hash
    gerado pela string Newick, garantindo a integridade dos dados.

    Parameters
    ----------
    encoded_data : dict
        Dicionário contendo a string Newick e o hash terminal.

    Return
    ------
    str or None
        A string Newick original se o hash for válido, caso contrário, None.
    """
    newick = encoded_data['newick']
    original_hash = encoded_data['terminal_hash']
    hash_object = hashlib.md5(newick.encode())
    if int(hash_object.hexdigest()[:4], 16) == original_hash:
        return newick
    else:
        return None
    
def calculate_tree_hash(data: Tree) -> Dict:
    """
    Calcula o hash terminal de uma subárvore.

    A função gera uma representação Newick da subárvore e, em seguida,
    calcula o hash MD5 dessa string, retornando os valores hash e Newick.

    Parameters
    ----------
    data : Bio.Phylo.BaseTree.Tree
        A subárvore para a qual o hash será calculado.

    Return
    ------
    dict
        Dicionário contendo a string Newick da subárvore e o valor do hash terminal.
    """
    newick = data.format("newick")
    hash_object = hashlib.md5(newick.encode())
    return {
        'newick': newick,
        'terminal_hash': int(hash_object.hexdigest()[:4], 16)
    }
