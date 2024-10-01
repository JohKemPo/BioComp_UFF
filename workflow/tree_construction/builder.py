

from Bio import Phylo
from Bio.Phylo.TreeConstruction import (DistanceCalculator, 
                                        DistanceTreeConstructor, 
                                        ParsimonyTreeConstructor, 
                                        ParsimonyScorer, 
                                        NNITreeSearcher)

class TreeBuilder:
    """
    Classe responsável pela construção de árvores filogenéticas a partir de alinhamentos de sequências.

    Esta classe oferece métodos para alinhar sequências utilizando ClustalW e MAFFT, calcular a matriz de distâncias,
    e construir árvores filogenéticas utilizando métodos baseados em distâncias e parcimônia.

    Attributes
    ----------
    count_noudes : int
        Contador que acompanha o número de nós nas árvores geradas.
    """

    def __init__(self, **kwargs):
        """
        Inicializa a instância da classe TreeBuilder.
        ----------

        Esta função atribui os valores dos argumentos fornecidos via kwargs como atributos
        da instância e inicializa o contador de nós das árvores.

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

        self.count_noudes = 0
    
    def distance_matrix(self, alignment):
        """
        Calcula a matriz de distâncias a partir de um alinhamento de sequências.

        Usa o método 'identity' para calcular a distância entre as sequências alinhadas, 
        retornando uma matriz de distâncias.

        Parameters
        ----------
        alignment : Bio.Align.MultipleSeqAlignment
            Objeto de alinhamento do Biopython.

        Return
        ------
        Bio.Phylo.TreeConstruction._DistanceMatrix
            Matriz de distâncias gerada a partir do alinhamento.
        """
        calculator = DistanceCalculator('identity')        
        return calculator.get_distance(alignment) 
    
    def distance_constructor(self, distance_matrix, construct_tree_method):
        """
        Constrói uma árvore filogenética usando métodos baseados em distâncias.

        Dependendo do método especificado (Neighbor-Joining ou UPGMA), esta função constrói
        e retorna uma árvore filogenética a partir da matriz de distâncias fornecida.

        Parameters
        ----------
        distance_matrix : Bio.Phylo.TreeConstruction._DistanceMatrix
            A matriz de distâncias gerada a partir do alinhamento.
        construct_tree_method : str
            Método de construção da árvore: 'nj' para Neighbor-Joining ou 'upgma' para UPGMA.

        Return
        ------
        Bio.Phylo.BaseTree.Tree
            A árvore filogenética construída utilizando o método especificado.
        """
        constructor = DistanceTreeConstructor()
        if construct_tree_method.lower() == 'nj':
            tree = constructor.nj(distance_matrix)
        if construct_tree_method.lower() == 'upgma':
            tree = constructor.upgma(distance_matrix)
        return tree
    
    def parsimony_constructor(self, distance_matrix, alignment, construct_tree_method):
        """
        Constrói uma árvore filogenética usando parcimônia.

        Este método utiliza uma árvore inicial gerada pelo método Neighbor-Joining ou UPGMA
        e, em seguida, aplica parcimônia para otimizar a árvore, retornando a árvore resultante.

        Parameters
        ----------
        distance_matrix : Bio.Phylo.TreeConstruction._DistanceMatrix
            A matriz de distâncias gerada a partir do alinhamento.
        alignment : Bio.Align.MultipleSeqAlignment
            O alinhamento das sequências que será utilizado na construção da árvore.
        construct_tree_method : str
            Método de construção da árvore inicial: 'nj' para Neighbor-Joining ou 'upgma' para UPGMA.

        Return
        ------
        Bio.Phylo.BaseTree.Tree
            A árvore filogenética otimizada por parcimônia.
        """
        constructor = DistanceTreeConstructor()
        scorer = ParsimonyScorer()
        searcher = NNITreeSearcher(scorer)
        if construct_tree_method.lower() == 'nj':
            starting_tree = constructor.nj(distance_matrix)
        if construct_tree_method.lower() == 'upgma':
            starting_tree = constructor.upgma(distance_matrix)
        constructor = ParsimonyTreeConstructor(searcher, starting_tree)
        tree = constructor.build_tree(alignment)

        return tree

    def save_tree(self, tree, path, format):
        """
        Salva a árvore filogenética em um arquivo.

        Exporta a árvore filogenética gerada para um arquivo no formato especificado (e.g., 'newick', 'nexus').

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            A árvore filogenética a ser salva.
        path : str
            O caminho do arquivo onde a árvore será salva.
        format : str
            O formato de saída da árvore, como 'newick' ou 'nexus'.

        Return
        ------
        None
        """
        Phylo.write(tree, path, format)
