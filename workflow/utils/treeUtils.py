from Bio import Phylo
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, Reference, FeatureLocation
from Bio.Seq import Seq
from collections import defaultdict
import json

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
                # Phylo.draw_ascii(decode_subtree)
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

def seqrecord_to_serializable_dict(record: SeqRecord) -> dict:
    def convert(obj):
        if isinstance(obj, Seq):
            return str(obj)
        elif isinstance(obj, SeqFeature):
            return {
                "type": obj.type,
                "location": str(obj.location),
                "strand": obj.strand,
                "qualifiers": {
                    k: convert(v) for k, v in obj.qualifiers.items()
                    if k not in {"translation", "seq"}
                }
            }
        elif isinstance(obj, FeatureLocation):
            return str(obj)
        elif isinstance(obj, defaultdict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, Reference):
            return {
                "title": obj.title,
                "authors": obj.authors,
                "journal": obj.journal,
                "pubmed_id": obj.pubmed_id,
                "comment": obj.comment
            }
        elif isinstance(obj, dict):
            return {
                k: convert(v) for k, v in obj.items()
                if k not in {"translation", "seq"}
            }
        elif isinstance(obj, list):
            return [convert(i) for i in obj]
        elif callable(obj):  # Ignorar métodos
            return None
        else:
            try:
                json.dumps(obj)  # Testa se é serializável
                return obj
            except:
                return str(obj)

    return {
        attr: convert(getattr(record, attr))
        for attr in dir(record)
        if not attr.startswith("_")
        and not callable(getattr(record, attr))
        and attr not in {"seq"}  # <- ESSA LINHA IGNORA 'seq' da raiz
    }



def fetch_local_record(accession: str, gbk_file: str) -> dict:
    """
    Busca informações de uma sequência em um arquivo local GenBank (.gb/.gbk)

    Parâmetros:
    - accession: str - o identificador da sequência (terminal da árvore)
    - gbk_file: str - caminho para o arquivo local no formato GenBank

    Retorna:
    - dict com informações principais da sequência ou erro
    """
    try:
        for record in SeqIO.parse(gbk_file, "genbank"):
            if record.id == accession or accession in record.annotations.get("accessions", []):
                return seqrecord_to_serializable_dict(record)
        return {"error": f"Acesso {accession} não encontrado no arquivo {gbk_file}"}
    except Exception as e:
        return {"error": str(e)}

def calculate_tree_hash(data: Tree, is_terminal: bool = False, gbk_file: str = None) -> Dict:
    """
    Calcula o hash terminal de uma subárvore e busca metadados localmente.

    Parameters
    ----------
    data : Tree or str
        Subárvore ou nome terminal (str) se is_terminal=True.
    is_terminal : bool
        Indica se é um terminal (folha) da árvore.
    gbk_file : str
        Caminho para o arquivo GenBank local (obrigatório se is_terminal=True)

    Return
    ------
    dict
        Contém a string Newick, hash MD5 e metadados (se aplicável)
    """
    if is_terminal:
        if not isinstance(data, str):
            raise ValueError("Para terminais, 'data' deve ser uma string com o nome do terminal.")
        newick = data
        hash_object = hashlib.md5(newick.encode())

        if gbk_file is None:
            raise ValueError("Arquivo GenBank local (gbk_file) deve ser fornecido para terminais.")
        
        metadata = fetch_local_record(data, gbk_file)
    else:
        newick = data.format("newick")
        hash_object = hashlib.md5(newick.encode())
        metadata = None

    return {
        'newick': newick.split('.')[0],
        'terminal_hash': int(hash_object.hexdigest()[:4], 16),
        'metadata': metadata
    }
    
def download_sequences(self, queries, output_file):
    """
    Passo 1: Baixa sequências do GenBank usando o Biopython.
    
    Parameters
    ----------

    queries: str
        String de consultas para o GenBank.
    output_file: str
        Caminho para salvar as sequências baixadas (formato GenBank).
    """
    try:
        for query in queries:
            if not query.name:
                continue  # Pula terminais sem nome
            Entrez.email = "email@email.com"
            handle = Entrez.esearch(db="nucleotide", term=query.name, retmax=1000)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
            with open(output_file, "+a") as f:
                f.write(handle.read())
            handle.close()
    except Exception as e:
        return {"error": str(e)}
