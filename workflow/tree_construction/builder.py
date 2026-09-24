

from Bio import Phylo, AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import (DistanceTreeConstructor,
                                        ParsimonyTreeConstructor,
                                        NNITreeSearcher)
import subprocess
import os
import re
import logging
from dendropy import Tree, DataSet
from io import StringIO
from workflow.tree_construction.distancia_p import matriz_distancia_p, FitchAmbiguidadeIupac


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


#: M7.4 — parâmetros do MrBayes, vindos do `tree_config` (chaves `mrbayes_*`).
#:
#: Os padrões reproduzem o que o pipeline **já executava de fato**, conferido
#: na saída do MrBayes 3.2.7 (caracterização de 2026-09-23, D20):
#:
#: - `ngen`, `samplefreq` e `printfreq` eram os literais do script.
#: - `nruns=2`, `nchains=4` não eram declarados; são os padrões da ferramenta
#:   (`help mcmc`), e passam a ser escritos para não depender deles.
#: - `burninfrac=0.25`: o script dizia `sumt burnin=250`, e o MrBayes **ignorava**
#:   — com `relburnin=yes` (padrão), `burnin` absoluto não é lido, e a saída
#:   dizia `discarding the first 25 % of sampled trees`. A fração efetiva sempre
#:   foi 25 %, não os 2,5 % que D20 e a ficha M7.1 deduziram lendo o script.
#: - `timeout_s=3600` era o literal do `subprocess.run`.
#:
#: Só `asdsf_max` é novo: é o gate de convergência, aditivo — ver
#: `avaliar_convergencia_mrbayes`.
MRBAYES_DEFAULTS = {
    'mrbayes_ngen': 1_000_000,
    'mrbayes_samplefreq': 100,
    'mrbayes_printfreq': 1000,
    'mrbayes_nruns': 2,
    'mrbayes_nchains': 4,
    'mrbayes_burninfrac': 0.25,
    'mrbayes_asdsf_max': 0.01,
    'mrbayes_timeout_s': 3600,
}


def mrbayes_settings(config: dict) -> dict:
    """
    Resolve e valida os parâmetros do MrBayes.

    Uma configuração inválida é recusada **aqui**, antes de a cadeia rodar — o
    MCMC é o método mais caro do pipeline, e descobrir depois de uma hora que o
    diagnóstico não era computável é desperdício que se evita de graça.

    `mrbayes_asdsf_max` ausente cai no padrão (0,01); **presente e `None`**
    (`null` no JSON) desliga o gate. A diferença é deliberada: desligar a
    verificação de convergência precisa ser um ato explícito do experimento, e
    o manifesto o declara.

    Raises
    ------
    ValueError
        Parâmetro fora do domínio, ou gate ligado com `nruns < 2` — o ASDSF
        compara frequências de bipartição **entre corridas independentes**, e
        com uma corrida só ele não existe (medido: o MrBayes 3.2.7 termina com
        código 0 e simplesmente não imprime a linha).
    """
    cfg = {}
    for chave, padrao in MRBAYES_DEFAULTS.items():
        valor = config[chave] if chave in config else padrao
        if chave == 'mrbayes_asdsf_max':
            cfg[chave] = None if valor is None else float(valor)
        elif chave == 'mrbayes_burninfrac':
            cfg[chave] = float(valor)
        else:
            cfg[chave] = int(valor)

    for chave in ('mrbayes_ngen', 'mrbayes_samplefreq', 'mrbayes_printfreq',
                  'mrbayes_nruns', 'mrbayes_nchains', 'mrbayes_timeout_s'):
        if cfg[chave] < 1:
            raise ValueError(f"{chave} precisa ser >= 1, veio {cfg[chave]}")
    if not 0.0 <= cfg['mrbayes_burninfrac'] < 1.0:
        raise ValueError(f"mrbayes_burninfrac precisa estar em [0, 1), "
                         f"veio {cfg['mrbayes_burninfrac']}")
    if cfg['mrbayes_asdsf_max'] is not None:
        if not cfg['mrbayes_asdsf_max'] > 0:
            raise ValueError(f"mrbayes_asdsf_max precisa ser > 0 ou null, "
                             f"veio {cfg['mrbayes_asdsf_max']}")
        if cfg['mrbayes_nruns'] < 2:
            raise ValueError(
                "mrbayes_nruns=1 torna o ASDSF incomputável (ele compara corridas "
                "independentes), e o gate de convergência está ligado. Use "
                "mrbayes_nruns >= 2, ou desligue o gate explicitamente com "
                "mrbayes_asdsf_max: null — o manifesto registrará a árvore como "
                "de convergência não verificada.")
    amostras = cfg['mrbayes_ngen'] // cfg['mrbayes_samplefreq'] + 1
    if int(amostras * (1 - cfg['mrbayes_burninfrac'])) < 2:
        raise ValueError(
            f"ngen/samplefreq = {amostras} amostras por corrida; depois de "
            f"descartar {cfg['mrbayes_burninfrac']:.0%} não sobra amostra para sumarizar")
    return cfg


class ConvergenciaNaoAtingida(RuntimeError):
    """O MCMC terminou, mas o diagnóstico não autoriza usar a árvore (D20)."""


#: Linhas do resumo de `sumt`. Com `=`, não `:` — a forma com `:` é a do
#: progresso do `mcmc`, calculada durante a cadeia; a do `sumt` é sobre as
#: mesmas amostras que formam o consenso, e é a que conta.
_RE_ASDSF = re.compile(r"Average standard deviation of split frequencies\s*=\s*(\S+)")
_RE_SDSF_MAX = re.compile(r"Maximum standard deviation of split frequencies\s*=\s*(\S+)")
_RE_PSRF_MAX = re.compile(r"Maximum PSRF for parameter values\s*=\s*(\S+)")


def _numero(texto):
    """`float` ou `None` — `NA`/`nan` do MrBayes não viram número (regra 5)."""
    try:
        valor = float(texto)
    except (TypeError, ValueError):
        return None
    return None if valor != valor else valor  # NaN


def ler_diagnosticos_mrbayes(stdout: str, pstat_path: str = None) -> dict:
    """
    Diagnósticos de convergência que o MrBayes imprime e ninguém lia (D20).

    Return
    ------
    dict
        ``asdsf``, ``sdsf_max``, ``psrf_max``, ``ess_min``. Cada um é ``None``
        quando a ferramenta não o produziu — por exemplo, ASDSF e PSRF com
        `nruns=1` —, nunca um número inventado.
    """
    def ultimo(regex):
        achados = regex.findall(stdout or "")
        return _numero(achados[-1]) if achados else None

    ess_min = None
    if pstat_path and os.path.exists(pstat_path):
        with open(pstat_path) as fh:
            linhas = [l.rstrip("\n").split("\t") for l in fh if "\t" in l]
        if linhas and "minESS" in linhas[0]:
            col = linhas[0].index("minESS")
            valores = [_numero(l[col]) for l in linhas[1:] if len(l) > col]
            valores = [v for v in valores if v is not None]
            ess_min = min(valores) if valores else None

    return {"asdsf": ultimo(_RE_ASDSF), "sdsf_max": ultimo(_RE_SDSF_MAX),
            "psrf_max": ultimo(_RE_PSRF_MAX), "ess_min": ess_min}


def avaliar_convergencia_mrbayes(diagnosticos: dict, asdsf_max) -> tuple:
    """
    Decide se a árvore consenso pode ser usada.

    Critério: **ASDSF < `asdsf_max`**, 0,01 por padrão. É o limiar que o manual
    do MrBayes 3.2 dá como boa indicação de convergência topológica entre
    corridas independentes (Ronquist et al., *MrBayes 3.2 manual*; ver
    `docs/science/11-auditoria-m7-inferencia.md §4` para a citação e o que ela
    não cobre). Estritamente menor: no empate, recusa.

    ESS e PSRF são registrados, mas **não** decidem — o pedido de M7.4 é o
    ASDSF, e um gate em ESS precisaria de limiar próprio e decisão própria.

    Return
    ------
    tuple
        ``(aceita, estado, motivo)``, com ``estado`` em ``atingida``,
        ``nao_atingida`` ou ``nao_verificada``.
    """
    asdsf = diagnosticos.get("asdsf")
    if asdsf_max is None:
        return (True, "nao_verificada",
                "gate desligado por configuração (mrbayes_asdsf_max=null); "
                f"ASDSF={asdsf if asdsf is not None else 'indisponível'}")
    if asdsf is None:
        return (False, "nao_atingida",
                "ASDSF ausente da saída do sumt: convergência não pôde ser "
                "verificada, e o gate está ligado")
    if asdsf < asdsf_max:
        return (True, "atingida", f"ASDSF={asdsf:g} < {asdsf_max:g}")
    return (False, "nao_atingida",
            f"ASDSF={asdsf:g} >= {asdsf_max:g}: as corridas independentes não "
            f"concordam na topologia; aumente mrbayes_ngen")


#: D20, item 7 — `prob` do comentário `[&…]` que segue um `)` (nó interno) no
#: `.con.tre` do MrBayes. `prob=` precedido de `&` ou `,`: não casa
#: `prob_stddev=`, `prob(percent)=` nem `prob+-sd=`.
_RE_PROB_NO_INTERNO = re.compile(r"\)\[&(?:[^\]]*?,)?prob=([^,\]]+)[^\]]*\]")

#: Formato da probabilidade posterior no Nexus gravado em `Trees/`.
FORMATO_POSTERIOR = "%.8f"


class ProbabilidadePosteriorInvalida(ValueError):
    """`prob` do consenso que não é número em [0, 1] — D20, item 7."""


def _rotulo_de_probabilidade(achado) -> str:
    """`)[&prob=p,…]` → `)p` (com o comentário restante apagado depois)."""
    texto = achado.group(1).strip()
    valor = _numero(texto)
    if valor is None or not 0.0 <= valor <= 1.0:
        raise ProbabilidadePosteriorInvalida(f"probabilidade posterior fora de [0, 1] no consenso "
                         f"do MrBayes: prob={texto!r}")
    return f"){texto}"


def modelo_declarado(caminho_log: str, padrao: str):
    """
    Modelo que a **própria ferramenta** declara ter usado, lido do log dela — M7.3.

    O `model` do manifesto é o literal da linha de comando (`GTR+G`); o que a
    ferramenta faz com ele depende da versão — `GTR+G` virou `GTR+F+G4` no
    IQ-TREE 3.1.3 e `GTR+FO+G4m` no RAxML-NG 2.0.2 (frequências empíricas num,
    estimadas por ML no outro). Um rótulo fixo no código repetiria o defeito de
    D26: o manifesto afirmaria o que foi pedido, não o que rodou.

    Return
    ------
    str or None
        ``None`` se o log não existe ou não traz a linha — nunca um palpite.
    """
    try:
        with open(caminho_log, errors="ignore") as fh:
            achado = re.search(padrao, fh.read(), re.MULTILINE)
    except OSError:
        return None
    return achado.group(1) if achado else None


from workflow.utils.external_tools import require_tool
from workflow.utils import tool_runs
from workflow.tree_construction.modelo_substituicao import MODELO_LEGADO

#: M7.3 — `model` no nível da ferramenta, em `tools_invoked`, é um valor por
#: ferramenta; com seleção por alinhamento ele deixa de ser invariante (dois
#: braços de alinhador podem escolher modelos diferentes). Nesse caso o nível
#: da ferramenta diz onde olhar, e o modelo de cada chamada vai em `runs[i]`.
_MODELO_POR_CHAMADA = "selecionado por alinhamento (M7.3, BIC): ver runs[i].modelo"

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

    def _modelo(self):
        """
        Modelo de substituição desta instância — M7.3.

        `modelo_substituicao` vem do controlador: a tradução do modelo que a
        seleção C′ escolheu para o alinhamento (`modelo_substituicao.traduzir`)
        ou `MODELO_LEGADO` com `model_selection: "nenhuma"`. Um `TreeBuilder`
        criado sem ele (scripts, testes) roda os literais de antes de M7.3 —
        quem escolhe o modelo é o controlador, não o construtor.
        """
        return getattr(self, 'modelo_substituicao', None) or MODELO_LEGADO

    def _modelo_selecionado(self):
        return 'modelfinder' in self._modelo()

    def _model_nivel_ferramenta(self, legado):
        """`model` de `tools_invoked.<ferramenta>`: o literal de sempre sem
        seleção (manifestos antigos e novos continuam comparáveis); com
        seleção, um apontador para o valor por chamada."""
        return _MODELO_POR_CHAMADA if self._modelo_selecionado() else legado

    def _anotar_modelo(self, ferramenta, saida, literal):
        """`modelo` (o que foi passado à ferramenta) e `modelo_origem` na chamada."""
        modelo = self._modelo()
        tool_runs.anotar(ferramenta, saida, modelo=literal,
                         modelo_origem=(f"seleção BIC (C′): ModelFinder escolheu "
                                        f"{modelo['modelfinder']}"
                                        if self._modelo_selecionado()
                                        else "fixo: literal anterior a M7.3"))
    
    def distance_matrix(self, alignment):
        """
        Calcula a matriz de distâncias a partir de um alinhamento de sequências.

        p-distância nos sítios em que ambas as sequências têm base
        determinada — **não** `DistanceCalculator('identity')` do Biopython,
        cujo denominador é o comprimento do alinhamento inteiro mesmo com
        `skip_letters` preenchido (D28, `docs/science/02-defeitos-que-
        alteram-resultado.md#d28`): um par lacuna×genoma completo media
        "quanto do alinhamento inteiro é igual", não "quanto das posições
        comparáveis é igual", e errava por até 15× no caso medido.

        Parameters
        ----------
        alignment : Bio.Align.MultipleSeqAlignment
            Objeto de alinhamento do Biopython.

        Return
        ------
        Bio.Phylo.TreeConstruction.DistanceMatrix
            Matriz de distâncias gerada a partir do alinhamento.

        Raises
        ------
        workflow.tree_construction.distancia_p.SemSitioComparavel
            Se algum par do alinhamento não tiver nenhum sítio em que ambas
            as sequências tenham base determinada — regra 5 do projeto,
            "não aplicável" não vira `0.0` nem `1.0`.
        """
        return matriz_distancia_p(alignment)
    
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
        # D28: lacuna/ambiguidade IUPAC como conjunto de bases possíveis, não
        # como um estado de caráter literal — `ParsimonyScorer` original
        # cobrava 1 passo de parcimônia por lacuna mesmo quando as demais
        # sequências concordavam (`docs/science/02-defeitos-que-alteram-
        # resultado.md#d28`).
        scorer = FitchAmbiguidadeIupac()
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
            # VARV aparece `97376`) e reexecutar não reproduz a árvore.
            #
            # D21 — e a semente **não basta**. Medido em 2026-08-26: com
            # `-nt 4`, três repetições da mesma semente, entrada, máquina e
            # versão devolveram **três topologias** (RF = 2); com `-nt 1`, uma
            # só. A ordem em que as reduções de ponto flutuante chegam decide
            # entre ótimos quase empatados, e o IQ-TREE não tem equivalente ao
            # `--workers 1` do RAxML-NG. Decisão do usuário: comprar
            # reprodutibilidade com tempo.
            #
            # `iqtree_threads` continua governando o BOOTSTRAP, que é
            # embaraçosamente paralelo e não decide topologia — só a busca de
            # ML roda em uma thread.
            cmd = [
                require_tool('iqtree'), '-s', align_path,
                '-m', self._modelo()['iqtree'], '-bb', '1000',
                '-seed', str(self.random_seed),
                '-pre', prefix,
                '-nt', '1'
            ]

            logging.info(f"IQ-TREE: {' '.join(cmd)}")
            tool_runs.registrar('iqtree', cmd, saida=output_path_tree,
                                seed=self.random_seed, threads=1,
                                threads_configurados=self.iqtree_threads,
                                model=self._model_nivel_ferramenta('GTR+G'),
                                bootstrap='UFBoot 1000',
                                nota='-nt 1 fixo por D21: com -nt N a mesma semente dá topologias diferentes')
            self._anotar_modelo('iqtree', output_path_tree, self._modelo()['iqtree'])
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            tool_runs.anotar('iqtree', output_path_tree, modelo_efetivo=modelo_declarado(
                prefix + '.iqtree', r"Model of substitution:\s*(\S+)"))
            
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

            # Sem semente e sem paralelização declarável: o FastTree não aceita
            # nem uma nem outra nesta chamada. Registrar assim mesmo é o que
            # diferencia "não se aplica" de "ninguém registrou" — a regra 5 do
            # projeto, aplicada ao manifesto.
            # M7.3 — era declarado `GTR (-nt -gtr)`, e lido como "sem variação
            # de taxa". O próprio log do FastTree diz `CAT approximation with
            # 20 rate categories`: CAT é o padrão, e o modelo efetivo é
            # GTR+CAT20.
            tool_runs.registrar('fasttree', cmd, saida=output_path_tree,
                                model='GTR+CAT20 (-nt -gtr; CAT 20 categorias é o padrão)')
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            achado = re.search(r"ML Model:\s*(.+)", result.stderr or '')
            tool_runs.anotar('fasttree', output_path_tree,
                             modelo_efetivo=achado.group(1).strip() if achado else None)
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
            #
            # M3.2/M7.2 — até aqui só busca de ML: sem `--bootstrap`/`--all`, o
            # RAxML-NG não calcula suporte nenhum (confirmado em
            # `docs/science/08-ficha-de-chamada-por-metodo.md §3`). `--all`
            # combina busca de ML + bootstrap + mapeamento de suporte num só
            # comando (verificado com `raxml-ng --help` e execução real numa
            # entrada sintética em 2026-09-02) e grava o suporte já mapeado em
            # `<prefix>.raxml.support` — o `.raxml.bestTree` continua existindo,
            # mas sem confidence. `--bs-trees 1000` casa a contagem de réplicas
            # com o `-bb 1000` do IQ-TREE, mas **não é a mesma métrica**: o
            # `--all` do RAxML-NG por padrão calcula FBP (Felsenstein bootstrap
            # proportion, bootstrap não-paramétrico clássico), e o `-bb` do
            # IQ-TREE é UFBoot (bootstrap ultrarrápido, aproximado). Ambos saem
            # em escala 0-100, mas não são o mesmo suporte e não devem ser lidos
            # com o mesmo limiar — ver `08-ficha-de-chamada-por-metodo.md §5`.
            cmd = [
                require_tool('raxml-ng'), '--all', '--msa', align_path,
                '--model', self._modelo()['raxml-ng'],
                '--threads', str(self.raxml_threads), '--workers', '1',
                '--seed', str(self.random_seed), '--tree', 'rand{10}',
                '--bs-trees', '1000',
                '--prefix', prefix
            ]

            logging.info(f"RAxML-NG: {' '.join(cmd)}")
            # `workers=1` é fixo por D17 e vai ao manifesto junto com `threads`:
            # medido que o esquema de paralelização muda a topologia com a
            # mesma semente, então ele é metadado do resultado, não da máquina.
            tool_runs.registrar('raxml-ng', cmd, saida=output_path_tree,
                                seed=self.random_seed, threads=self.raxml_threads,
                                workers=1, model=self._model_nivel_ferramenta('GTR+G'),
                                bootstrap='FBP (Felsenstein) 1000 réplicas via --all',
                                nota='FBP não é a mesma escala/interpretação do UFBoot do IQ-TREE')
            self._anotar_modelo('raxml-ng', output_path_tree, self._modelo()['raxml-ng'])
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            tool_runs.anotar('raxml-ng', output_path_tree, modelo_efetivo=modelo_declarado(
                prefix + '.raxml.log', r"^\s*Model:\s*(\S+)"))

            # `.raxml.support` só existe quando o bootstrap terminou e o
            # suporte foi mapeado na árvore de ML; `.raxml.bestTree` é o
            # fallback sem confidence, para não derrubar o pipeline inteiro se
            # o bootstrap falhar por algum motivo e a busca de ML tiver ido bem.
            possible_tree_files = [
                prefix + '.raxml.support',
                prefix + '.raxml.bestTree',
            ]

            tree_file = next((f for f in possible_tree_files if os.path.exists(f)), None)

            if tree_file:
                # M7.6 — o fallback existe, mas não pode ser silencioso: uma
                # árvore de RAxML-NG sem suporte é outro produto, e quem lê o
                # manifesto precisa saber qual das duas saiu.
                if tree_file.endswith('.raxml.bestTree'):
                    logging.warning(
                        "RAxML-NG: .raxml.support ausente; usando .raxml.bestTree, "
                        "SEM valores de suporte (bootstrap não concluiu).")
                    tool_runs.anotar('raxml-ng', output_path_tree,
                                     arvore_lida='.raxml.bestTree',
                                     suporte='ausente: .raxml.support não foi gravado')
                else:
                    tool_runs.anotar('raxml-ng', output_path_tree,
                                     arvore_lida='.raxml.support')
                tree = Phylo.read(tree_file, 'newick')

                Phylo.write(tree, output_path_tree, 'nexus')

                logging.info(f"Arquivos RAxML-NG salvos em: {tmp_dir}")
                return tree
            else:
                raise FileNotFoundError(f"Arquivo de árvore não encontrado em: {tmp_dir}")
                    
        except subprocess.CalledProcessError as e:
            logging.error(f"Erro no RAxML-NG: {e.stderr}")
            logging.info(f"Output do RAxML-NG: {e.stdout}")
            raise
    
    
    
    def _clean_mrbayes_tree(self, tree_path, output_newick):
        """
        Consenso do MrBayes (`.con.tre`) em Newick simples, **com** a
        probabilidade posterior de cada clado — D20, item 7.

        Formato caracterizado no MrBayes 3.2.7 (`sumt`, `nruns=2`): cada nó
        interno sai como
        ``)[&prob=9.84042553e-01,prob_stddev=…,prob_range={…,…},prob(percent)="98",prob+-sd="98+-0"]:0.01[&length_mean=…]``.
        Antes, todo comentário `[...]` era apagado e a probabilidade ia junto
        (`confidence = None` em todo nó interno). Agora o `prob` do nó interno
        vira o rótulo numérico `)p`, que o `Phylo.read` promove a
        `.confidence` — o mesmo mecanismo que já leva o suporte do FastTree e
        do IQ-TREE ao Nexus (`08-ficha-de-chamada-por-metodo.md §1`).

        **Escala 0-1, probabilidade posterior.** Não é comparável a UFBoot nem a
        FBP (0-100) pelo mesmo limiar, nem por reescala — é outra quantidade
        (`Backend/src/suporte_de_ramo.py`, métrica ``posterior``, sem limiar
        decidido). O manifesto declara isso em `tools_invoked.mrbayes`.

        Só o `prob` agregado das corridas é propagado; `prob_stddev` e
        `prob_range` (a discordância entre corridas) ficam no `.con.tre` em
        `tmp/`. Nó interno sem `prob` fica sem rótulo → `confidence = None`,
        nunca 0 (regra 5). A raiz do consenso não é anotada pelo MrBayes (não
        é clado de uma árvore não enraizada) e continua sem suporte.

        Raises
        ------
        ProbabilidadePosteriorInvalida
            `prob` que não é número em [0, 1] — um valor fora da escala
            gravado como suporte seria pior do que a ausência.
        """
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

        # D20, item 7 — antes de apagar os comentários, o `prob` de cada nó
        # interno vira rótulo. Só depois de `)`: folhas também trazem
        # `[&prob=1…]`, e rótulo de folha é o nome do táxon.
        tree_str = _RE_PROB_NO_INTERNO.sub(_rotulo_de_probabilidade, tree_str)
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
    

    @staticmethod
    def mrbayes_script(cfg: dict, seed: int, nexus_name: str, modelo: dict = None) -> str:
        """
        Script de comandos do MrBayes para `cfg` (saída de `mrbayes_settings`).

        Tudo o que decide o resultado é escrito, inclusive o que coincide com o
        padrão da ferramenta (`nruns`, `nchains`, `relburnin`): um padrão que
        muda de versão para versão não pode ser a fonte do número publicado.
        A mesma fração de *burn-in* vale para o diagnóstico do `mcmc`, para o
        `sump` e para o `sumt` — o ASDSF que o gate lê é calculado sobre as
        mesmas amostras que formam o consenso.

        `modelo` (M7.3) é a tradução do modelo selecionado para o alinhamento
        (`modelo_substituicao.traduzir`): `lset` e, com frequências iguais,
        `prset statefreqpr=fixed(equal)`. Sem ele, o literal de antes de M7.3.
        """
        bf = cfg['mrbayes_burninfrac']
        modelo = modelo or MODELO_LEGADO
        linhas_modelo = modelo['mrbayes_lset'] + "\n"
        if modelo.get('mrbayes_prset'):
            linhas_modelo += modelo['mrbayes_prset'] + "\n"
        return (
            f"set autoclose=yes nowarn=yes seed={seed} swapseed={seed}\n"
            f"execute {nexus_name}\n"
            f"{linhas_modelo}"
            f"mcmc ngen={cfg['mrbayes_ngen']} nruns={cfg['mrbayes_nruns']} "
            f"nchains={cfg['mrbayes_nchains']} printfreq={cfg['mrbayes_printfreq']} "
            f"samplefreq={cfg['mrbayes_samplefreq']} relburnin=yes burninfrac={bf}\n"
            f"sump relburnin=yes burninfrac={bf}\n"
            f"sumt relburnin=yes burninfrac={bf}\n"
            f"quit\n"
        )

    def mrbayes_constructor(self, alignment, output_path_tree, generations=None):
        """
        Constrói árvore usando MrBayes, e só a devolve se a cadeia convergiu.

        Parameters
        ----------
        alignment : Bio.Align.MultipleSeqAlignment
        output_path_tree : str
        generations : int, optional
            Sobrepõe `mrbayes_ngen`. Mantido por compatibilidade com a
            assinatura anterior; o caminho do pipeline usa a configuração.

        Raises
        ------
        ConvergenciaNaoAtingida
            O ASDSF não ficou abaixo de `mrbayes_asdsf_max`, ou não pôde ser
            lido. A árvore consenso **não** é gravada em `Trees/` — senão a
            próxima execução a reaproveitaria como se fosse válida; fica só em
            `tmp/` para inspeção.
        """
        config = dict(vars(self))
        if generations is not None:
            config['mrbayes_ngen'] = generations
        cfg = mrbayes_settings(config)
        seed = int(self.random_seed)
        try:
            base_name = os.path.basename(output_path_tree).replace('.nexus', '').replace('.nwk', '')
            # D20, item 1 — era `split('/PhyloTreeMiner/')[-1]`: um caminho
            # RELATIVO que dependia do nome do repositório (já renomeado uma
            # vez, de FPM-Tree). Agora a mesma regra dos outros construtores,
            # sobre o caminho absoluto: `cwd=tmp_dir` abaixo torna um caminho
            # relativo ambíguo.
            tmp_dir = os.path.join(
                os.path.dirname(os.path.abspath(output_path_tree)).split('/Trees')[0],
                'tmp', f'mrbayes_{base_name}')
            os.makedirs(tmp_dir, exist_ok=True)
            
            nexus_path = os.path.join(tmp_dir, 'alignment.nexus')
        
            for record in alignment:
                if not hasattr(record, 'annotations'):
                    record.annotations = {}
                record.annotations['molecule_type'] = 'DNA'
            
            AlignIO.write(alignment, nexus_path, 'nexus')
            
            script_path = os.path.join(tmp_dir, 'run_mb.txt')
            with open(script_path, 'w') as f:
                f.write(self.mrbayes_script(cfg, seed, os.path.basename(nexus_path),
                                            self._modelo()))
            
            # A linha de comando do MrBayes é só o binário: o que decide o
            # resultado está no script lido pela entrada padrão, então os
            # parâmetros do script vão ao manifesto. D20, item 2: a semente é a
            # mesma `random_seed` de IQ-TREE/RAxML-NG, em `seed` e `swapseed`
            # (a segunda governa a troca de estados entre cadeias aquecidas).
            cmd_mb = [require_tool('mrbayes')]
            tool_runs.registrar(
                'mrbayes', cmd_mb, saida=output_path_tree,
                model=self._model_nivel_ferramenta(
                    'GTR+G4 (lset nst=6 rates=gamma; Ngammacat=4 padrão)'),
                suporte='probabilidade posterior do clado (sumt), escala 0-1',
                nota_suporte=('probabilidade posterior NÃO é comparável a UFBoot '
                              '(IQ-TREE) nem a FBP (RAxML-NG) pelo mesmo limiar, '
                              'nem por reescala — D20 item 7'),
                seed=seed, swapseed=seed,
                ngen=cfg['mrbayes_ngen'], samplefreq=cfg['mrbayes_samplefreq'],
                nruns=cfg['mrbayes_nruns'], nchains=cfg['mrbayes_nchains'],
                burninfrac=cfg['mrbayes_burninfrac'],
                asdsf_max=cfg['mrbayes_asdsf_max'],
                gate_convergencia=('ligado' if cfg['mrbayes_asdsf_max'] is not None
                                   else 'desligado por configuração'),
                script=os.path.basename(script_path))
            modelo = self._modelo()
            self._anotar_modelo('mrbayes', output_path_tree, "; ".join(
                linha for linha in (modelo['mrbayes_lset'], modelo.get('mrbayes_prset')) if linha))
            # D20, item 5 — o descritor era aberto sem `with` e nunca fechado.
            with open(script_path, 'r') as entrada:
                result = subprocess.run(
                    cmd_mb,
                    stdin=entrada,
                    capture_output=True,
                    timeout=cfg['mrbayes_timeout_s'],
                    cwd=tmp_dir
                )
            
            stdout_text = result.stdout.decode('utf-8', errors='ignore') if result.stdout else ''
            stderr_text = result.stderr.decode('utf-8', errors='ignore') if result.stderr else ''
            with open(os.path.join(tmp_dir, 'mb_stdout.log'), 'w') as fh:
                fh.write(stdout_text)
            
            if result.returncode != 0:
                logging.error(f"MrBayes exit code: {result.returncode}")
                logging.error(f"MrBayes stdout: {stdout_text[:1000]}")  
                logging.error(f"MrBayes stderr: {stderr_text[:1000]}")
                raise subprocess.CalledProcessError(result.returncode, cmd_mb, stdout_text, stderr_text)
            
            # D20, item 4 — o gate. Lido e anotado ANTES de gravar a árvore em
            # `Trees/`: uma árvore recusada que chegasse ao disco seria
            # reaproveitada pela próxima execução (`_load_existing_tree`).
            diagnosticos = ler_diagnosticos_mrbayes(
                stdout_text, os.path.join(tmp_dir, 'alignment.nexus.pstat'))
            aceita, estado, motivo = avaliar_convergencia_mrbayes(
                diagnosticos, cfg['mrbayes_asdsf_max'])
            tool_runs.anotar('mrbayes', output_path_tree, **diagnosticos,
                             convergencia=estado, convergencia_motivo=motivo)
            if not aceita:
                logging.error(f"MrBayes: árvore RECUSADA — {motivo}")
                raise ConvergenciaNaoAtingida(motivo)
            if estado == 'nao_verificada':
                logging.warning(f"MrBayes: convergência NÃO verificada — {motivo}")
            else:
                logging.info(f"MrBayes: convergência atingida — {motivo}")

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

                        # D20, item 7 — `%1.2f` (padrão do escritor do
                        # Biopython) arredondaria a posterior a 2 casas:
                        # 0,949 viraria "0.95" e cruzaria um limiar de 0,95.
                        # Oito casas recuperam #obs/N exatamente para até
                        # ~10^7 amostras.
                        Phylo.write(tree, output_path_tree, "nexus",
                                    format_confidence=FORMATO_POSTERIOR)
                        logging.info(f"Árvore MrBayes salva em: {output_path_tree}")
                        return tree
                    except ProbabilidadePosteriorInvalida:
                        # Não é "arquivo ilegível, tente o próximo": o
                        # consenso existe e diz algo impossível.
                        raise
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
            logging.error(f"MrBayes excedeu o tempo limite de execução "
                          f"({cfg['mrbayes_timeout_s']} s, mrbayes_timeout_s)")
            raise
        except ConvergenciaNaoAtingida:
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
