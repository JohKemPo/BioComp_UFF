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

import hashlib, io, os, time

# D5 — a identidade de clado usada em produção é a canônica, definida uma única
# vez em `workflow.stability.clade_identity`. As funções `legacy_*` deste módulo
# ficam só para auditoria e para reler artefatos antigos.
from workflow.stability.clade_identity import (canonical_bipartition,
                                               canonical_item_id,
                                               strip_accession_version)


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
    Extrai a transação de uma árvore: os itens (bipartições canônicas) de suas subárvores.

    A função percorre o dicionário fornecido, que contém informações sobre árvores e subárvores,
    e agrega os valores de 'List_terminals_hash' de todas as subárvores em uma lista.

    D29 — o item é a bipartição canônica, não o clado como escrito. Subárvore cujo
    clado não é bipartição não trivial (o clado-raiz, o clado de n-1 táxons)
    tem 'List_terminals_hash' ``None`` e fica de fora; a mesma bipartição vista
    por dois nós da mesma árvore (os dois filhos de uma raiz bifurcante, um nó
    unário) conta uma vez. Transação é conjunto: a ordem estável só mantém a
    saída reprodutível.

    Parameters
    ----------
    data : dict
        Dicionário contendo informações sobre as árvores e subárvores, onde cada subárvore
        possui o item armazenado na chave 'List_terminals_hash'.

    Return
    ------
    list
        Itens distintos da transação, na ordem em que aparecem; vazia se a árvore
        não tem bipartição não trivial (n < 4, ou só arestas triviais).
    """
    all_terminais = []
    vistos = set()
    
    for _, tree_info in data.items():
        for _, subtree_info in tree_info.items():
            item = subtree_info['List_terminals_hash']
            if item is None or item in vistos:
                continue
            vistos.add(item)
            all_terminais.append(item)

    return all_terminais

def verificar_mesmo_conjunto_de_folhas(data: List[Dict]) -> None:
    """
    Recusa a mineração de árvores que não compartilham o mesmo conjunto de táxons (D29).

    Uma bipartição só é comparável entre árvores sobre o mesmo ``X``; com
    conjuntos diferentes, "o mesmo item" não tem significado e o suporte seria
    ruído. A recusa é a de `TreeSet`. O conjunto vem de 'Tree_terminals_hash',
    gravado em toda subárvore pelo `SubtreeBuilder`; registro sem o campo (metadado
    anterior ao D29) não pode ser verificado e é ignorado.

    Parameters
    ----------
    data : list of dict
        Uma entrada por árvore, ``{nome_da_árvore: {nome_da_subárvore: registro}}``.

    Raises
    ------
    ValueError
        Se duas árvores têm conjuntos de táxons diferentes.
    """
    por_arvore = {}
    for info in data:
        for tree_name, subtree_info in info.items():
            ids = {sub.get('Tree_terminals_hash') for sub in subtree_info.values()} - {None}
            if ids:
                por_arvore[tree_name] = ids
    distintos = {i for ids in por_arvore.values() for i in ids}
    if len(distintos) > 1:
        referencia = next(iter(por_arvore))
        divergentes = [n for n, ids in por_arvore.items() if ids != por_arvore[referencia]]
        raise ValueError(
            f"Árvore(s) {sorted(divergentes)[:5]} não compartilham o mesmo conjunto de "
            f"táxons de '{referencia}': a bipartição só é comparável sobre o mesmo X.")

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
    Identidade **legada** de clado: MD5 de 16 bits sobre a lista em ordem de travessia.

    .. deprecated::
        Mantida para auditoria e para reler `metadata.json` gerado antes de M1.2.
        A identidade de produção é `encode_clade_to_int` — ver D5.

    Dois defeitos somados: o hash é calculado sobre a *representação textual* da
    lista, o que o torna dependente da ordem de travessia (o mesmo clado recebe
    identificadores diferentes em árvores que ordenam filhos de modo distinto, e
    o suporte é **subestimado**); e o espaço de 16 bits faz clados distintos
    colidirem, o que **fabrica** suporte.

    Parameters
    ----------
    lst : list
        A lista de valores para codificação.

    Return
    ------
    int
        Valor inteiro de 16 bits, dependente da ordem.
    """
    lst_str = str(lst)
    hash_object = hashlib.md5(lst_str.encode())
    return int(hash_object.hexdigest()[:4], 16)


def encode_clade_to_int(terminal_names: List[str]) -> int:
    """
    Identidade de um clado **como escrito**, a partir dos nomes de seus terminais (D5).

    Canônica quanto à ordem de travessia e ao rótulo (60 bits de espaço; os nomes
    são normalizados por `strip_accession_version`, de modo que o mesmo clado tenha
    a mesma identidade numa árvore de IQ-TREE, com rótulo truncado em 10
    caracteres, e numa de FastTree, com rótulo íntegro — ver D13). **Não é
    canônica quanto à raiz**: o mesmo ramo escrito pelo outro lado é outro clado
    e outro inteiro. Para o item de mineração use `encode_bipartition_to_int` (D29).

    Parameters
    ----------
    terminal_names : list of str
        Nomes dos terminais do clado, em qualquer ordem.

    Return
    ------
    int
        Identificador do clado como escrito.
    """
    return canonical_item_id(strip_accession_version(name) for name in terminal_names)

def encode_bipartition_to_int(terminal_names: List[str], all_terminal_names: List[str]) -> Optional[int]:
    """
    Item de mineração: identidade da bipartição canônica do clado (D29).

    O clado lido do Newick depende de onde cada programa escreveu a raiz; a
    bipartição não (`03-metricas §2.2`). Usa a mesma `canonical_bipartition` do
    relatório de estabilidade e a mesma `canonical_item_id`, de modo que um ramo
    tem um só inteiro, qualquer que seja o lado em que o arquivo o escreveu.

    Parameters
    ----------
    terminal_names : list of str
        Nomes dos terminais do clado, em qualquer ordem.
    all_terminal_names : list of str
        Nomes de todos os terminais da árvore.

    Return
    ------
    int or None
        Identificador da bipartição, ou ``None`` quando o clado não é uma
        bipartição não trivial (o clado-raiz, o de n-1 táxons, um terminal).
        ``None`` é "não é item", nunca um número.
    """
    todos = frozenset(strip_accession_version(name) for name in all_terminal_names)
    lado = canonical_bipartition((strip_accession_version(name) for name in terminal_names), todos)
    return canonical_item_id(lado) if lado is not None else None

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
    if canonical_item_id([strip_accession_version(newick)]) == original_hash:
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
        # O hash é do rótulo NORMALIZADO, que é também o que se grava em 'newick'.
        # Antes hashava-se o rótulo bruto e gravava-se o truncado: `decode_tree_hash`
        # nunca validava, e `NC_008030.` e `NC_008030.1` viravam dois terminais
        # distintos (D5 + D13).
        newick = strip_accession_version(data)

        if gbk_file is None:
            raise ValueError("Arquivo GenBank local (gbk_file) deve ser fornecido para terminais.")

        metadata = fetch_local_record(data, gbk_file)
    else:
        newick = data.format("newick")
        metadata = None

    return {
        'newick': newick,
        'terminal_hash': canonical_item_id([newick]),
        'metadata': metadata
    }
    
_arquivos_iniciados = set()


def _ordenar_por_acesso(texto, ids):
    """
    Reordena os registros GenBank de um lote na ordem dos acessos pedidos.

    O efetch não garante a ordem da resposta; sem isso o arquivo muda de
    uma execução para outra. Acesso que não casa com nenhum pedido vai ao fim,
    na ordem em que chegou.
    """
    chave = lambda acesso: acesso.split(".")[0]
    ordem = {}
    for pos, acesso in enumerate(ids):
        ordem.setdefault(chave(acesso), pos)

    registros = []
    for bloco in texto.split("//\n"):
        if not bloco.strip():
            continue
        registro = bloco.strip("\n") + "\n//\n"
        try:
            acesso = SeqIO.read(io.StringIO(registro), "genbank").id
        except Exception:
            acesso = ""
        registros.append((ordem.get(chave(acesso), len(ids)), len(registros), registro))

    return "".join(registro for _, _, registro in sorted(registros))


def download_sequences(self, queries, output_file, log, batch_size: int = 10):
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
        Entrez.email = "email@email.com"
        Entrez.max_tries = 3
        Entrez.sleep_between_tries = 2
        ncbi_api_key = os.getenv("NCBI_API_KEY")
        if ncbi_api_key:
            Entrez.api_key = ncbi_api_key
        # O builder chama esta função uma vez por árvore; só a primeira chamada de cada processo zera o arquivo.
        caminho = os.path.abspath(output_file)
        if caminho not in _arquivos_iniciados:
            with open(output_file, "w") as f:
                f.write('')
            _arquivos_iniciados.add(caminho)
                
                
        id_list = [q.name for q in queries if q.name]
        total_ids = len(id_list)
        log.info(f"     Número de queiries detectadas {total_ids}. Inidicando processamento em lotes...")

        for i in range(0, total_ids, batch_size):
            batch_ids = id_list[i:i+batch_size]
            lote_num = i//batch_size + 1
            total_lotes = (total_ids - 1)//batch_size + 1
            log.info(f"     Baixando lote {lote_num}/{total_lotes} ({len(batch_ids)} sequências)")
            
            for tentativa in range(3):
                log.info(f"     Executando queries: {'| '.join(batch_ids)}...")
                try:                   
                    handle = Entrez.efetch(
                        db="nucleotide", 
                        id=batch_ids,
                        rettype="gb",
                        retmode="text"
                    )

                    data = handle.read()
                    handle.close()
                    
                    if data:
                        log.info(f"         Resultados salvos com sucesso!")
                        with open(output_file, "a") as f:
                            f.write(_ordenar_por_acesso(data, batch_ids))
                        break
                    else:
                        log.error(f"        Dados vazios!")
                        raise Exception("Dados vazios")
                except Exception as e:
                    log.warning(f"      Tentativa {tentativa+1}/3 falhou: {e}")
                    if tentativa < 2:  
                        time.sleep(5 * (tentativa + 1))  
                    else:
                        log.error(f"        Falha definitiva no lote {lote_num}")
                    

            if lote_num < total_lotes:
                time.sleep(3)
                
    except Exception as e:
        log.info(f"     {e}")
        return {"error": str(e)}
    
    