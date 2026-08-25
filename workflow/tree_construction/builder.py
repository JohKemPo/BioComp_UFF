

from Bio import Phylo, AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import (DistanceCalculator, 
                                        DistanceTreeConstructor, 
                                        ParsimonyTreeConstructor, 
                                        ParsimonyScorer, 
                                        NNITreeSearcher)
import subprocess
import os
import re
import logging
from dendropy import Tree, DataSet
from io import StringIO


#TODO: Melhorias de parametros dos novos metodos ( esta é somente a versão estavel )
#: Padrões de reprodutibilidade. Ficam aqui, e não espalhados pelas chamadas das
#: ferramentas, para que o manifesto de execução possa declarar exatamente o
#: mesmo valor que o pipeline vai usar — uma única fonte da verdade.
#:
#: Não confundir com `num_threads` do `tree_config`, que governa apenas o
#: **alinhamento** (`mafft --thread`, `clustalo --threads`). Inferência e
#: alinhamento têm perfis de paralelismo diferentes, e os projetos existentes
#: usam valores bem distintos ali — de 1 em VARV a 16 em ZIKV-480.
REPRODUCIBILITY_DEFAULTS = {
    'random_seed': 12345,
    'raxml_threads': 4,
    'iqtree_threads': 4,
}


def reproducibility_settings(config: dict) -> dict:
    """
    Resolve semente e paralelização a partir da configuração do projeto.

    D11 — sem semente fixa, reexecutar não reproduz a árvore. D17 — mesmo com a
    semente fixa, deixar a paralelização a cargo da ferramenta muda a topologia
    entre máquinas (medido: RF = 8 no mesmo alinhamento, mesma semente).

    Parameters
    ----------
    config : dict
        `tree_config` do projeto; chaves ausentes caem no padrão.

    Return
    ------
    dict
        ``random_seed``, ``raxml_threads`` e ``iqtree_threads`` já como inteiros.
    """
    return {chave: int(config.get(chave, padrao))
            for chave, padrao in REPRODUCIBILITY_DEFAULTS.items()}


from workflow.utils.external_tools import require_tool

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

        # Reprodutibilidade (D11, D17): semente e paralelização são parâmetros
        # do experimento, não detalhe de implementação.
        for chave, valor in reproducibility_settings(kwargs).items():
            setattr(self, chave, valor)

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

    def iqtree_constructor(self, alignment, output_path_tree):
        """
        Constrói árvore usando IQ-TREE 2.
        Arquivos extras são salvos na pasta tmp.
        """
        try:
            base_name = os.path.basename(output_path_tree).replace('.nexus', '').replace('.nwk', '')
            tmp_dir = os.path.join(os.path.dirname(output_path_tree).split('/Trees')[0], 'tmp', f'iqtree_{base_name}')
            os.makedirs(tmp_dir, exist_ok=True)
            
            align_path = os.path.join(tmp_dir, f'{base_name}.phylip')
            AlignIO.write(alignment, align_path, 'phylip')
            
            prefix = os.path.join(tmp_dir, base_name)
            
            # D11 — sem `-seed`, o IQ-TREE gera a própria semente (nos logs de
            # VARV aparece `97376`) e reexecutar não reproduz a árvore. `-nt` é
            # fixado pelo mesmo motivo de D17 no RAxML: número de threads
            # decidido pela máquina torna a execução incomparável entre elas.
            cmd = [
                require_tool('iqtree'), '-s', align_path,
                '-m', 'GTR+G', '-bb', '1000',
                '-seed', str(self.random_seed),
                '-pre', prefix,
                '-nt', str(self.iqtree_threads)
            ]

            logging.info(f"IQ-TREE: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            possible_tree_files = [
                prefix + '.treefile',
                prefix + '.contree',
                prefix + '.tre'
            ]
            
            tree_file_found = None
            for tree_file in possible_tree_files:
                if os.path.exists(tree_file):
                    tree_file_found = tree_file
                    break
            
            if tree_file_found:
                tree = Phylo.read(tree_file_found, 'newick')
                Phylo.write(tree, output_path_tree, 'nexus')
                
                logging.info(f"Arquivos IQ-TREE salvos em: {tmp_dir}")
                return tree
            else:
                raise FileNotFoundError(f"Nenhum arquivo de árvore encontrado em: {tmp_dir}")
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Erro no IQ-TREE: {e.stderr}")
            logging.info(f"Output do IQ-TREE: {e.stdout}")
            raise
    
    def fasttree_constructor(self, alignment, output_path_tree):
        """
        Constrói árvore usando FastTree.
        Arquivos extras são salvos na pasta tmp.
        """
        try:
            base_name = os.path.basename(output_path_tree).replace('.nexus', '').replace('.nwk', '')
            tmp_dir = os.path.join(os.path.dirname(output_path_tree).split('/Trees')[0], 'tmp', f'fasttree_{base_name}')
            os.makedirs(tmp_dir, exist_ok=True)
            
            align_path = os.path.join(tmp_dir, f'{base_name}.fasta')
            AlignIO.write(alignment, align_path, 'fasta')
            
            cmd = [require_tool('fasttree'), '-nt', '-gtr', align_path]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            tree = Phylo.read(StringIO(result.stdout), 'newick')
            
            Phylo.write(tree, output_path_tree, 'nexus')
            
            log_path = os.path.join(tmp_dir, f'{base_name}.log')
            with open(log_path, 'w') as log_file:
                log_file.write(result.stderr)
            
            logging.info(f"Arquivos FastTree salvos em: {tmp_dir}")
            return Phylo.read(output_path_tree, 'nexus')
                    
        except subprocess.CalledProcessError as e:
            logging.error(f"Erro no FastTree: {e.stderr}")
            logging.info(f"Output do FastTree: {e.stdout}")
            raise
    
    def raxml_ng_constructor(self, alignment, output_path_tree):
        """
        Constrói árvore usando RAxML-NG.
        Arquivos extras são salvos na pasta tmp.
        """
        try:
            base_name = os.path.basename(output_path_tree).replace('.nexus', '').replace('.nwk', '')
            tmp_dir = os.path.join(os.path.dirname(output_path_tree).split('/Trees')[0], 'tmp', f'raxml_{base_name}')
            os.makedirs(tmp_dir, exist_ok=True)
            
            align_path = os.path.join(tmp_dir, f'{base_name}.phylip')
            AlignIO.write(alignment, align_path, 'phylip')
            
            prefix = os.path.join(tmp_dir, base_name)
            
            # D17 — `--threads auto` escolhe o esquema de paralelização a partir
            # do número de núcleos da máquina, e o esquema **muda a topologia**:
            # medido RF = 8 entre duas execuções com a MESMA semente, variando só
            # a paralelização (verossimilhanças −591486,234 e −591486,233, dois
            # ótimos quase equivalentes). Além disso, `auto` já escolheu
            # `5 workers x 3 threads` e derrubou o processo com SIGSEGV.
            # Um worker só, com número de threads declarado, torna a execução
            # comparável entre máquinas. Custo medido: ~10% de tempo.
            cmd = [
                require_tool('raxml-ng'), '--msa', align_path,
                '--model', 'GTR+G',
                '--threads', str(self.raxml_threads), '--workers', '1',
                '--seed', str(self.random_seed), '--tree', 'rand{10}',
                '--prefix', prefix
            ]

            logging.info(f"RAxML-NG: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            tree_file = prefix + '.raxml.bestTree'
            if os.path.exists(tree_file):
                tree = Phylo.read(tree_file, 'newick')
                
                Phylo.write(tree, output_path_tree, 'nexus')
                
                logging.info(f"Arquivos RAxML-NG salvos em: {tmp_dir}")
                return tree
            else:
                raise FileNotFoundError(f"Arquivo de árvore não encontrado: {tree_file}")
                    
        except subprocess.CalledProcessError as e:
            logging.error(f"Erro no RAxML-NG: {e.stderr}")
            logging.info(f"Output do RAxML-NG: {e.stdout}")
            raise
    
    
    
    def _clean_mrbayes_tree(self, tree_path, output_newick):
        with open(tree_path, "r") as f:
            content = f.read()

        translate_dict = {}
        translate_match = re.search(r"translate\s*((?:.|\n)+?);", content, re.IGNORECASE)
        if translate_match:
            translate_text = translate_match.group(1)
            for match in re.finditer(r'(\d+)\s+([^,\n]+)', translate_text):
                num, label = match.groups()
                translate_dict[num] = label.strip().rstrip(',')

        tree_match = re.search(r"tree.*?=.*?\((.*?)\);", content, re.DOTALL)
        if not tree_match:
            raise ValueError("Árvore não encontrada")
        
        tree_str = f"({tree_match.group(1)})"
        
        tree_str = re.sub(r"\[.*?\]", "", tree_str)
        
        for num, label in translate_dict.items():
            tree_str = re.sub(rf'(?<=[\(,]){num}(?=[:\),])', label, tree_str)
        
        #tree_str = re.sub(r":[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", "", tree_str)
        
        inner_count = 1
        result = []
        for char in tree_str:
           if char == ')':
                result.append(f')inner{inner_count}')
                inner_count += 1
           else:
                result.append(char)
        
        #final_tree = ''.join(result)
        
        with open(output_newick, "w") as f:
            f.write(tree_str + ";\n")
        
        return output_newick
    

    def mrbayes_constructor(self, alignment, output_path_tree, generations=1000000):
        """
        Constrói árvore usando MrBayes.
        Arquivos extras são salvos na pasta tmp.
        """
        try:
            base_name = os.path.basename(output_path_tree).replace('.nexus', '').replace('.nwk', '')
            tmp_dir = os.path.join((os.path.dirname(output_path_tree).split('/PhyloTreeMiner/')[-1]).split('/Trees')[0],'tmp', f'mrbayes_{base_name}')
            os.makedirs(tmp_dir, exist_ok=True)
            
            nexus_path = os.path.join(tmp_dir, 'alignment.nexus')
        
            for record in alignment:
                if not hasattr(record, 'annotations'):
                    record.annotations = {}
                record.annotations['molecule_type'] = 'DNA'
            
            AlignIO.write(alignment, nexus_path, 'nexus')
            
            mrbayes_script = f"""set autoclose=yes nowarn=yes
    execute {os.path.basename(nexus_path)}
    lset nst=6 rates=gamma
    mcmc ngen={generations} printfreq=1000 samplefreq=100
    sump
    sumt burnin=250
    quit
    """
            
            script_path = os.path.join(tmp_dir, 'run_mb.txt')
            with open(script_path, 'w') as f:
                f.write(mrbayes_script)
            
            result = subprocess.run(
                [require_tool('mrbayes')],
                stdin=open(script_path, 'r'),
                capture_output=True,
                timeout=3600,
                cwd=tmp_dir
            )
            
            stdout_text = result.stdout.decode('utf-8', errors='ignore') if result.stdout else ''
            stderr_text = result.stderr.decode('utf-8', errors='ignore') if result.stderr else ''
            
            if result.returncode != 0:
                logging.error(f"MrBayes exit code: {result.returncode}")
                logging.error(f"MrBayes stdout: {stdout_text[:1000]}")  
                logging.error(f"MrBayes stderr: {stderr_text[:1000]}")
                raise subprocess.CalledProcessError(result.returncode, [require_tool('mrbayes')], stdout_text, stderr_text)
            
            tree_files = [
                'alignment.nexus.con.tre',
                'alignment.con.tre',
                'alignment.t',
                'alignment.nex.con.tre'
            ]
            
            for tree_file in tree_files:
                tree_path = os.path.join(tmp_dir, tree_file)
                if os.path.exists(tree_path):
                    try:
                        nwk_file = self._clean_mrbayes_tree(tree_path,os.path.join(tmp_dir,"tree_clean.nwk") )
                        tree = Phylo.read(nwk_file, "newick")

                        Phylo.write(tree, output_path_tree, "nexus")
                        logging.info(f"Árvore MrBayes salva em: {output_path_tree}")
                        return tree
                    except Exception as e:
                        logging.warning(f"Erro ao ler árvore {tree_file}: {e}")
                        continue
            
            
            
            raise FileNotFoundError(f"Arquivo de árvore não encontrado em: {tmp_dir}. Arquivos: {os.listdir(tmp_dir)}")
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Erro no MrBayes - exit code: {e.returncode}")
            if hasattr(e, 'stdout') and e.stdout:
                logging.error(f"Stdout: {e.stdout[:1000]}")
            if hasattr(e, 'stderr') and e.stderr:
                logging.error(f"Stderr: {e.stderr[:1000]}")
            raise
        except subprocess.TimeoutExpired:
            logging.error("MrBayes excedeu o tempo limite de execução")
            raise
        except Exception as e:
            logging.error(f"Erro inesperado: {e}")
            raise

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
        if isinstance(tree, list):
            Phylo.write(tree, path, format)
        else:
            Phylo.write([tree], path, format)
