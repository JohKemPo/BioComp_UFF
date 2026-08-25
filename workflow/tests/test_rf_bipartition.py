"""
Testes da unidade de comparação entre árvores — D3.

O defeito: `clade_sets` guardava clados **enraizados** e `rf_matrix` os comparava
entre pipelines. FastTree, IQ-TREE, RAxML e NJ emitem topologia **não enraizada**,
escrita em Newick com raiz trifurcante; o clado que se lê do arquivo depende de
onde a ferramenta pôs a raiz, que é convenção de escrita. Comparar clados entre
essas árvores mede a convenção, não a topologia — em VARV-6, três métodos que
produzem a mesma topologia eram reportados como 75% discordantes.

Segundo erro, somado ao primeiro: a normalização dividia por ``2(n-2)``, o
máximo para árvore **enraizada**. Sobre bipartições o máximo é ``2(n-3)``.

Executar com:

    python -m unittest workflow.tests.test_rf_bipartition
"""

from __future__ import annotations

import io
import unittest

from Bio import Phylo

from workflow.stability.clade_identity import canonical_bipartition
from workflow.stability.stability import PipelineLabel, StabilityAnalyzer, TreeSet


TODOS = frozenset("ABCDEF")


def conjunto_de_arvores(**newicks):
    """TreeSet montado em memória, sem passar por disco."""
    arvores = {nome: Phylo.read(io.StringIO(texto), "newick")
               for nome, texto in newicks.items()}
    rotulos = {nome: PipelineLabel.parse(nome) for nome in arvores}
    return TreeSet(arvores, rotulos)


class TestBipartioCanonica(unittest.TestCase):

    def test_um_clado_e_seu_complemento_sao_a_mesma_bipartição(self):
        """É o ponto todo: onde a raiz caiu não pode mudar o objeto comparado."""
        self.assertEqual(canonical_bipartition({"A", "B"}, TODOS),
                         canonical_bipartition({"C", "D", "E", "F"}, TODOS))

    def test_representante_e_o_lado_menor(self):
        self.assertEqual(canonical_bipartition({"C", "D", "E", "F"}, TODOS),
                         frozenset({"A", "B"}))

    def test_empate_resolve_por_ordem_lexicografica(self):
        """Com |A| = |X∖A|, o critério de tamanho não decide; sem desempate
        estável, a mesma bipartição teria dois representantes."""
        esquerda = canonical_bipartition({"A", "B", "C"}, TODOS)
        direita = canonical_bipartition({"D", "E", "F"}, TODOS)
        self.assertEqual(esquerda, direita)
        self.assertEqual(esquerda, frozenset({"A", "B", "C"}))

    def test_bipartição_trivial_e_none(self):
        """Aresta externa: um lado com menos de 2 táxons não carrega informação
        topológica e não pode entrar na contagem de RF."""
        self.assertIsNone(canonical_bipartition({"A"}, TODOS))
        self.assertIsNone(canonical_bipartition({"B", "C", "D", "E", "F"}, TODOS))

    def test_e_idempotente(self):
        uma_vez = canonical_bipartition({"A", "B"}, TODOS)
        self.assertEqual(canonical_bipartition(uma_vez, TODOS), uma_vez)


class TestRFNaoEnraizada(unittest.TestCase):

    def setUp(self):
        # Mesma topologia não enraizada, escrita com raízes diferentes — é
        # exatamente o que FastTree e IQ-TREE fazem com a mesma árvore.
        self.mesma_topologia = conjunto_de_arvores(
            mafft_fasttree="((A,B),(C,D),(E,F));",
            mafft_iqtree="(((A,B),(C,D)),(E,F));",
        )
        self.topologias_diferentes = conjunto_de_arvores(
            mafft_fasttree="((A,B),(C,D),(E,F));",
            mafft_nj_distance="((A,C),(B,D),(E,F));",
        )

    def test_mesma_topologia_com_raizes_diferentes_tem_rf_zero(self):
        """O caso de VARV-6: 0,75 de discordância entre árvores idênticas."""
        analisador = StabilityAnalyzer(self.mesma_topologia)
        matriz = analisador.rf_matrix(normalized=False)
        self.assertEqual(matriz["mafft_fasttree"]["mafft_iqtree"], 0)

    def test_analise_enraizada_ainda_ve_diferenca_no_mesmo_par(self):
        """Caracteriza o comportamento anterior, que continua disponível via
        `rooted=True` e continua legítimo quando a raiz é real."""
        analisador = StabilityAnalyzer(self.mesma_topologia, rooted=True)
        matriz = analisador.rf_matrix(normalized=False)
        self.assertGreater(matriz["mafft_fasttree"]["mafft_iqtree"], 0)

    def test_topologias_diferentes_continuam_diferentes(self):
        """A correção não pode zerar tudo: RF tem de seguir detectando conflito."""
        analisador = StabilityAnalyzer(self.topologias_diferentes)
        matriz = analisador.rf_matrix(normalized=False)
        self.assertGreater(matriz["mafft_fasttree"]["mafft_nj_distance"], 0)

    def test_denominador_e_2_n_menos_3(self):
        analisador = StabilityAnalyzer(self.topologias_diferentes)
        n = analisador.tree_set.n_taxa
        bruta = analisador.rf_matrix(normalized=False)["mafft_fasttree"]["mafft_nj_distance"]
        normalizada = analisador.rf_matrix()["mafft_fasttree"]["mafft_nj_distance"]
        self.assertAlmostEqual(normalizada, bruta / (2 * (n - 3)))

    def test_denominador_enraizado_continua_2_n_menos_2(self):
        analisador = StabilityAnalyzer(self.topologias_diferentes, rooted=True)
        n = analisador.tree_set.n_taxa
        bruta = analisador.rf_matrix(normalized=False)["mafft_fasttree"]["mafft_nj_distance"]
        normalizada = analisador.rf_matrix()["mafft_fasttree"]["mafft_nj_distance"]
        self.assertAlmostEqual(normalizada, bruta / (2 * (n - 2)))

    def test_diagonal_e_zero(self):
        analisador = StabilityAnalyzer(self.mesma_topologia)
        matriz = analisador.rf_matrix()
        self.assertEqual(matriz["mafft_fasttree"]["mafft_fasttree"], 0.0)

    def test_arvore_pequena_demais_devolve_none_e_nao_zero(self):
        """`04-rigor-cientifico §3`: "não aplicável" nunca é um número. Com
        n < 4 a RF não enraizada é indefinida, e devolver 0 faria topologias
        incomparáveis passarem por idênticas."""
        pequeno = conjunto_de_arvores(mafft_fasttree="(A,B,C);", mafft_iqtree="((A,B),C);")
        matriz = StabilityAnalyzer(pequeno).rf_matrix()
        self.assertIsNone(matriz["mafft_fasttree"]["mafft_iqtree"])
        self.assertEqual(matriz["mafft_fasttree"]["mafft_fasttree"], 0.0)

    def test_contagem_de_bipartições_acompanha_a_distancia(self):
        """Exigência de `03-metricas §3`: sem |B(T)| ao lado, um valor baixo é
        ambíguo entre topologias parecidas e árvore malresolvida."""
        analisador = StabilityAnalyzer(self.mesma_topologia)
        contagens = analisador.bipartition_counts()
        n = analisador.tree_set.n_taxa
        for pipeline, quantas in contagens.items():
            self.assertLessEqual(quantas, n - 3, pipeline)

    def test_politomia_reduz_a_contagem_de_bipartições(self):
        com_politomia = conjunto_de_arvores(
            mafft_fasttree="((A,B),(C,D),(E,F));",
            mafft_nj_distance="((A,B),C,D,E,F);",
        )
        contagens = StabilityAnalyzer(com_politomia).bipartition_counts()
        self.assertLess(contagens["mafft_nj_distance"], contagens["mafft_fasttree"])

    def test_fatores_ignoram_pares_sem_distancia_definida(self):
        pequeno = conjunto_de_arvores(mafft_fasttree="(A,B,C);", mafft_nj_distance="((A,B),C);")
        efeitos = StabilityAnalyzer(pequeno).factor_effects()
        self.assertEqual(efeitos["inference"]["n"], 0)
        self.assertIsNone(efeitos["inference"]["mean"])



class TestRotuloDePipeline(unittest.TestCase):
    """D19 — dois arquivos não podem designar o mesmo pipeline.

    `INFERENCE_METHODS` não conhecia `nj_parsimony` nem `upgma_parsimony`, e
    ambos casavam com o sufixo `parsimony`. `TreeSet.trees` é um dicionário
    indexado pelo rótulo: a segunda árvore sobrescrevia a primeira sem uma linha
    de log, e `M` — o denominador de todo suporte metodológico — vinha menor que
    o número de árvores em disco. No conjunto de validação eram 12 pipelines
    para 14 árvores.
    """

    def test_nj_e_upgma_parcimonia_sao_pipelines_distintos(self):
        from workflow.stability.stability import PipelineLabel
        nj = PipelineLabel.parse("tree_dataset_final_clustalo_nj_parsimony.nexus",
                                 prefix="tree_dataset_final_")
        upgma = PipelineLabel.parse("tree_dataset_final_clustalo_upgma_parsimony.nexus",
                                    prefix="tree_dataset_final_")
        self.assertNotEqual(nj.name, upgma.name)
        self.assertEqual(nj.inference, "nj_parsimony")
        self.assertEqual(upgma.inference, "upgma_parsimony")

    def test_sufixo_mais_longo_vence(self):
        """A regra não pode depender da ordem em que INFERENCE_METHODS foi
        escrita — `parsimony` e `nj_parsimony` casam os dois."""
        from workflow.stability.stability import PipelineLabel
        rotulo = PipelineLabel.parse("tree_dataset_final_mafft_upgma_distance.nexus",
                                     prefix="tree_dataset_final_")
        self.assertEqual(rotulo.inference, "upgma_distance")
        self.assertEqual(rotulo.aligner, "mafft")

    def test_colisao_de_rotulo_e_recusada_e_nao_silenciada(self):
        import os
        import shutil
        import tempfile
        from workflow.stability.stability import TreeSet

        pasta = tempfile.mkdtemp()
        try:
            # Dois nomes que o rótulo não distingue: nenhum sufixo conhecido.
            for nome in ("tree_dataset_final_mafft_metodoX.nexus",
                         "tree_dataset_final_mafft_metodoY.nexus"):
                with open(os.path.join(pasta, nome), "w") as handle:
                    handle.write("#NEXUS\nBegin Trees;\n Tree t=((A,B),(C,D));\nEnd;\n")
            with self.assertRaises(ValueError) as erro:
                TreeSet.from_directory(pasta)
            self.assertIn("mesmo pipeline", str(erro.exception))
        finally:
            shutil.rmtree(pasta, ignore_errors=True)

if __name__ == "__main__":
    unittest.main()
