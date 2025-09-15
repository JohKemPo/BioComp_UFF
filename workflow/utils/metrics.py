from typing import List, Dict, Tuple, Set, Any

from Bio import Phylo
from itertools import combinations

from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd

def get_bipartitions(tree: Phylo.BaseTree.Tree):
    """
    Gera as bipartições da árvore fornecida.

    Parameters
    ----------
    tree : Phylo.BaseTree.Tree
        A árvore filogenética da qual as bipartições serão extraídas.

    Yields
    ------
    Set[frozenset]
        Um conjunto imutável representando as bipartições dos terminais da árvore.
    """
    terminals = tree.get_terminals()
    for clade in tree.find_clades(order='level'):
        if clade.is_terminal():
            continue
        yield frozenset(terminals.index(t) for t in clade.get_terminals())

def rf_distance(tree1: Phylo.BaseTree.Tree, tree2: Phylo.BaseTree.Tree) -> int:
    """
    Calcula a distância de Robinson-Foulds (RF) entre duas árvores filogenéticas.

    A distância RF mede a diferença topológica entre duas árvores. Quanto maior a distância RF, 
    maior a diferença entre as árvores. Uma distância RF de 0 indica que as árvores são idênticas em termos de topologia.

    Parameters
    ----------
    tree1 : Phylo.BaseTree.Tree
        A primeira árvore filogenética.
    tree2 : Phylo.BaseTree.Tree
        A segunda árvore filogenética.

    Returns
    -------
    int
        A distância RF entre as duas árvores.
    """
    bp_tree1 = set(get_bipartitions(tree1))
    bp_tree2 = set(get_bipartitions(tree2))
    
    unique_to_tree1 = bp_tree1 - bp_tree2
    unique_to_tree2 = bp_tree2 - bp_tree1
    
    return len(unique_to_tree1) + len(unique_to_tree2)

def get_trees_data_list(data_dict: Dict) -> List[Dict[str, Phylo.BaseTree.Tree]]:
    """
    Converte um dicionário de dados de árvores em uma lista de dicionários contendo os nomes das árvores e suas estruturas.

    Parameters
    ----------
    data_dict : Dict
        Dicionário contendo informações sobre árvores, seus métodos de construção e algoritmos.

    Returns
    -------
    List[Dict[str, Phylo.BaseTree.Tree]]
        Lista de dicionários contendo os nomes das árvores e suas respectivas estruturas.
    """
    trees_data = []
    for alg, trs in data_dict.items():
        for method, trees in trs.items():
            for mtd, tr in trees.items():
                if len(tr) > 0:
                    trees_data.append({f"tree_{alg}_{method}_{mtd}": tr[0]})
    return trees_data

def process_rf_distance(data_dict: Dict) -> Dict[Tuple[str, str], int]:
    """
    Processa e calcula a distância RF entre todas as combinações de árvores em um conjunto de dados.

    Parameters
    ----------
    data_dict : Dict
        Dicionário contendo as árvores filogenéticas a serem comparadas.

    Returns
    -------
    Dict[Tuple[str, str], int]
        Dicionário contendo as distâncias RF entre todas as combinações de árvores.
    """
    rf_distances = {}
    trees_data = get_trees_data_list(data_dict)

    for tree_dict1, tree_dict2 in combinations(trees_data, 2):
        name1, tree1 = list(tree_dict1.items())[0]
        name2, tree2 = list(tree_dict2.items())[0]
        
        distance = rf_distance(tree1, tree2)
        rf_distances[(name1, name2)] = distance
    
    return rf_distances

def flatten_multi_trees(multi_trees: Dict) -> Dict[str, Any]:
    """
    Converte a estrutura aninhada de multi_trees em um dicionário plano.
    
    Parameters
    ----------
    multi_trees : Dict
        Dicionário com estrutura aninhada de árvores.
        
    Returns
    -------
    Dict[str, Any]
        Dicionário plano com nomes únicos para cada árvore.
    """
    flat_trees = {}
    
    for alignment_method, methods_dict in multi_trees.items():
        for method_type, trees in methods_dict.items():
            if method_type in ["iqtree", "fasttree", "raxml", "mrbayes"]:
                for i, tree in enumerate(trees):
                    tree_name = f"{alignment_method}_{method_type}_{i}"
                    flat_trees[tree_name] = tree
            elif method_type in ["distance", "parsimony"]:
                for submethod, sub_trees in trees.items():
                    for i, tree in enumerate(sub_trees):
                        tree_name = f"{alignment_method}_{method_type}_{submethod}_{i}"
                        flat_trees[tree_name] = tree
    
    return flat_trees

def process_rf_distance(multi_trees: Dict) -> Dict[Tuple[str, str], int]:
    """
    Processa e calcula a distância RF entre todas as combinações de árvores em multi_trees.
    
    Parameters
    ----------
    multi_trees : Dict
        Dicionário contendo as árvores filogenéticas a serem comparadas no formato multi_trees.
        
    Returns
    -------
    Dict[Tuple[str, str], int]
        Dicionário contendo as distâncias RF entre todas as combinações de árvores.
    """
    rf_distances = {}
    
    flat_trees = flatten_multi_trees(multi_trees)
    
    tree_names = list(flat_trees.keys())
    
    for (name1, name2) in combinations(tree_names, 2):
        tree1 = flat_trees[name1]
        tree2 = flat_trees[name2]
        
        try:
            distance = rf_distance(tree1, tree2)
            rf_distances[(name1, name2)] = distance
        except Exception as e:
            print(f"Erro ao calcular distância RF entre {name1} e {name2}: {e}")
            rf_distances[(name1, name2)] = -1  
    
    return rf_distances

def plot_heatmap_distances(data_dict: Dict, 
                          base_name: str, 
                          path: str, 
                          scores: Dict[Tuple[str, str], int] = None, 
                          distance_matrix: np.ndarray = None) -> np.ndarray:
    """
    Plota um mapa de calor das distâncias RF entre árvores do formato data_dict.
    
    Parameters
    ----------
    data_dict : Dict
        Dicionário contendo as árvores no formato data_dict.
    base_name : str
        Nome base para o arquivo de saída.
    path : str
        Caminho onde o arquivo de saída será salvo.
    scores : Dict[Tuple[str, str], int], optional
        Dicionário de distâncias RF calculadas entre as árvores.
    distance_matrix : np.ndarray, optional
        Matriz de distâncias para ser usada no mapa de calor.
        
    Returns
    -------
    np.ndarray
        Matriz de distâncias utilizada para o mapa de calor.
    """
    flat_trees = flatten_multi_trees(data_dict)
    tree_names = list(flat_trees.keys())
    n_trees = len(tree_names)
    
    if distance_matrix is None and scores:
        distance_matrix = np.zeros((n_trees, n_trees))
        
        for (name1, name2), distance in scores.items():
            try:
                i = tree_names.index(name1)
                j = tree_names.index(name2)
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
            except ValueError:
                print(f"Árvore não encontrada: {name1} ou {name2}")
    
    df = pd.DataFrame(distance_matrix, index=tree_names, columns=tree_names)
    
    plt.figure(figsize=(max(12, n_trees), max(10, n_trees)))
    sns.heatmap(df, annot=True, cmap='coolwarm', linewidths=.5, fmt=".0f",
                cbar_kws={'label': 'Robinson-Foulds Distance'})
    
    plt.title('Robinson-Foulds Distance Heatmap')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    plt.savefig(f"{path}/tree_{base_name}_heatmap_of_distances.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return distance_matrix

def analyze_rf_distances(multi_trees: Dict, scores: Dict[Tuple[str, str], int]) -> Dict:
    """
    Analisa estatísticas das distâncias RF.
    
    Parameters
    ----------
    multi_trees : Dict
        Dicionário com as árvores.
    scores : Dict[Tuple[str, str], int]
        Distâncias RF calculadas.
        
    Returns
    -------
    Dict
        Estatísticas das distâncias.
    """
    flat_trees = flatten_multi_trees(multi_trees)
    tree_names = list(flat_trees.keys())
    
    method_groups = {}
    for name in tree_names:
        parts = name.split('_')
        method_key = '_'.join(parts[:2])  
        if method_key not in method_groups:
            method_groups[method_key] = []
        method_groups[method_key].append(name)
    
    distances = list(scores.values())
    stats = {
        'total_trees': len(tree_names),
        'total_comparisons': len(scores),
        'mean_distance': np.mean(distances),
        'median_distance': np.median(distances),
        'min_distance': np.min(distances),
        'max_distance': np.max(distances),
        'method_groups': method_groups
    }
    
    return stats

def get_clade_leaves(clade: Phylo.BaseTree.Clade) -> List[str]:
    """
    Obtém a lista de folhas (terminais) de um clado.

    Parameters
    ----------
    clade : Phylo.BaseTree.Clade
        O clado do qual as folhas serão extraídas.

    Returns
    -------
    List[str]
        Lista contendo os nomes das folhas do clado.
    """
    return sorted([leaf.name for leaf in clade.get_terminals()])

def get_quartet_relations(tree: Phylo.BaseTree.Tree) -> Set[Tuple[Tuple[str, ...], ...]]:
    """
    Retorna um conjunto de quartetos de uma árvore.

    Parameters
    ----------
    tree : Phylo.BaseTree.Tree
        A árvore filogenética da qual os quartetos serão extraídos.

    Returns
    -------
    Set[Tuple[Tuple[str, ...], ...]]
        Conjunto de quartetos representando as relações entre quatro folhas.
    """
    terminals = tree.get_terminals()
    quartets = set()
    
    for quartet in combinations(terminals, 4):
        subclades = [tree.common_ancestor(a, b) for a, b in combinations(quartet, 2)]
        sorted_subclades = sorted([tuple(get_clade_leaves(subclade)) for subclade in subclades])
        quartet_relation = tuple(sorted_subclades)
        quartets.add(quartet_relation)
    
    return quartets

def quartet_distance(tree1: Phylo.BaseTree.Tree, tree2: Phylo.BaseTree.Tree) -> int:
    """
    Calcula a distância de quarteto entre duas árvores.

    Parameters
    ----------
    tree1 : Phylo.BaseTree.Tree
        A primeira árvore filogenética.
    tree2 : Phylo.BaseTree.Tree
        A segunda árvore filogenética.

    Returns
    -------
    int
        A distância de quarteto entre as duas árvores.
    """
    quartet_set1 = get_quartet_relations(tree1)
    quartet_set2 = get_quartet_relations(tree2)
    
    unique_to_tree1 = quartet_set1 - quartet_set2
    unique_to_tree2 = quartet_set2 - quartet_set1
    
    return len(unique_to_tree1) + len(unique_to_tree2)

def process_quartet_distances(data_dict: Dict) -> np.ndarray:
    """
    Processa e calcula a distância de quarteto entre todas as combinações de árvores em um conjunto de dados.

    Parameters
    ----------
    data_dict : Dict
        Dicionário contendo as árvores filogenéticas a serem comparadas.

    Returns
    -------
    np.ndarray
        Matriz de distâncias de quarteto entre as árvores.
    """
    trees = get_trees_data_list(data_dict)
    quartet_distances = np.zeros((len(trees), len(trees)))

    for i, j in combinations(range(len(trees)), 2):
        name1, tree1 = list(trees[i].items())[0]
        name2, tree2 = list(trees[j].items())[0]
        
        distance = quartet_distance(tree1, tree2)
        quartet_distances[i, j] = distance
        quartet_distances[j, i] = distance

    return quartet_distances

def plot_quartet_distance_dendrogram(data_dict: Dict, scores: np.ndarray, base_name: str, path: str) -> None:
    """
    Plota um dendrograma com base na distância de quarteto entre as árvores.

    Parameters
    ----------
    data_dict : Dict
        Dicionário contendo as árvores filogenéticas.
    scores : np.ndarray
        Matriz de distâncias de quarteto.
    base_name : str
        Nome base para o arquivo de saída.
    path : str
        Caminho onde o arquivo de saída será salvo.

    Returns
    -------
    None
    """
    trees = get_trees_data_list(data_dict)
    tree_names = [list(tree.keys())[0].replace('_', " ") for tree in trees]
    linked = linkage(scores, method='single')

    dendrogram(linked, labels=tree_names, orientation='top')
    plt.title('Dendrograma das Árvores com base na Distância de Quarteto')
    plt.xlabel('Árvores Filogenéticas')
    plt.ylabel('Distância de Quarteto')
    plt.xticks(rotation=90)
    plt.savefig(f"{path}/tree_{base_name}_quartet_distance_dendrogram.png")
    plt.close()

def extract_support_values(tree: Dict, support_values: List) -> List:
    """
    Extrai os valores de suporte de uma árvore filogenética.

    Parameters
    ----------
    tree : Dict
        Dicionário que representa a árvore filogenética.
    support_values : List
        Lista que acumula os valores de suporte extraídos.

    Returns
    -------
    List
        Lista contendo os valores de suporte extraídos da árvore.
    """
    if 'supports' not in tree:
        return support_values
    support_values.extend(tree['supports'])
    for child in tree['metadata']['children']:
        extract_support_values(child, support_values)
    return support_values

def process_histogram_frequence(data: List[Dict], path: str) -> None:
    """
    Processa e plota um histograma das frequências de valores de suporte de subárvores.

    Parameters
    ----------
    data : List[Dict]
        Lista de dicionários contendo as árvores filogenéticas e suas subárvores.
    path : str
        Caminho onde será salvo o plot

    Returns
    -------
    None
    """
    support_values = list()
    for tree_key in data:
        for sub_tree_key in tree_key:
            for _, subtrees in sub_tree_key.items():
                for subtree in subtrees.values():
                    support_values = extract_support_values(subtree,support_values)
    
    plt.figure(figsize=(10, 6))
    plt.hist(support_values, bins=10, edgecolor='black')
    plt.title('Histograma de Frequências das Subárvores')
    plt.xlabel('Valor de Suporte')
    plt.ylabel('Frequência')
    plt.savefig(f'{path}/support_frequence_bar_plot.png')
    plt.close()

    counts, bin_edges = np.histogram(support_values, bins=10)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    plt.figure(figsize=(10, 6))
    plt.plot(bin_centers, counts, marker='o', linestyle='-', color='blue')
    plt.title('Frequências das Subárvores')
    plt.xlabel('Valor de Suporte')
    plt.ylabel('Frequência')
    plt.grid(True)
    plt.savefig(f'{path}/support_frequence_frequency_polygon.png')
    plt.close()
