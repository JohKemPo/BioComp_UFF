"""
Testes do módulo `workflow.subtree_mining` — D4.

O defeito: a varredura de limiares do FPMax sobrescrevia `result_fpmax['support']`
com o próprio limiar. A coluna passava a guardar o parâmetro da varredura, não a
fração de árvores que contém o itemset, e o mesmo padrão aparecia em vários
limiares — em VARV-49, 7 de 7 itemsets distintos com mais de um "suporte", e 2
deles exibidos ao mesmo tempo como *method-sensitive signature* e como
*topologically robust* na Deep Analysis.

Executar com:

    python -m unittest workflow.tests.test_subtree_mining
"""

from __future__ import annotations

import unittest

import numpy as np
import pandas as pd
from mlxtend.frequent_patterns import fpmax
from mlxtend.preprocessing import TransactionEncoder

from workflow.subtree_mining.miner import SubtreeMiner


def varredura(matriz, grade=np.arange(0.1, 1.1, 0.1)):
    """Reproduz a varredura de limiares do `process_group` no modo 'auto'."""
    te = TransactionEncoder()
    df = pd.DataFrame(te.fit(matriz).transform(matriz), columns=te.columns_)
    bruto = pd.DataFrame()
    for limiar in grade:
        limiar = round(float(limiar), 2)
        resultado = fpmax(df, min_support=limiar, use_colnames=True)
        resultado['min_support_threshold'] = limiar
        bruto = pd.concat([bruto, resultado], ignore_index=True)
    return bruto


class TestConsolidacaoFPMax(unittest.TestCase):
    """`consolidate_fpmax_results` — uma linha por itemset, suporte real."""

    def setUp(self):
        # 4 "árvores"; o clado A ocorre em todas, B em duas, C em uma.
        self.matriz = [
            ['A', 'B', 'C'],
            ['A', 'B'],
            ['A', 'D'],
            ['A', 'D'],
        ]
        self.bruto = varredura(self.matriz)

    def test_suporte_nao_e_o_limiar_da_varredura(self):
        """O núcleo de D4. O itemset {A} está em 4 de 4 árvores: suporte 1,0,
        seja qual for o limiar em que a varredura o devolveu."""
        consolidado = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        linha_a = consolidado[consolidado['itemsets'] == frozenset({'A'})]
        self.assertEqual(len(linha_a), 1)
        self.assertAlmostEqual(float(linha_a['support'].iloc[0]), 1.0)
        self.assertEqual(int(linha_a['n_trees'].iloc[0]), 4)

    def test_uma_linha_por_itemset(self):
        consolidado = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        self.assertEqual(len(consolidado), consolidado['itemsets'].nunique())
        self.assertGreater(len(self.bruto), len(consolidado),
                           "a varredura precisa mesmo repetir itemsets, senão o teste não prova nada")

    def test_faixa_de_limiares_e_preservada(self):
        """Deduplicar não pode apagar em que limiares o itemset foi devolvido —
        é o que permite auditar a varredura depois.

        {A} tem suporte 1,0 mas só é devolvido a partir do limiar 0,6: abaixo
        disso ele não é **maximal**, porque {A,B} e {A,D} ainda são frequentes.
        É a distinção que D4 apagava — limiar e suporte são números diferentes e
        nem sequer variam juntos."""
        consolidado = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        linha_a = consolidado[consolidado['itemsets'] == frozenset({'A'})]
        self.assertAlmostEqual(float(linha_a['support'].iloc[0]), 1.0)
        self.assertAlmostEqual(float(linha_a['min_support_threshold'].iloc[0]), 0.6)
        self.assertAlmostEqual(float(linha_a['max_support_threshold'].iloc[0]), 1.0)

        linha_ab = consolidado[consolidado['itemsets'] == frozenset({'A', 'B'})]
        self.assertAlmostEqual(float(linha_ab['support'].iloc[0]), 0.5)
        self.assertAlmostEqual(float(linha_ab['min_support_threshold'].iloc[0]), 0.3)
        self.assertAlmostEqual(float(linha_ab['max_support_threshold'].iloc[0]), 0.5)

    def test_limiar_nunca_excede_o_suporte_real(self):
        """Invariante do FPMax: só devolve itemset cujo suporte alcança o
        limiar pedido. Se algum limiar ultrapassar o suporte, a coluna errada
        voltou a ser gravada."""
        consolidado = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        for _, linha in consolidado.iterrows():
            self.assertLessEqual(linha['max_support_threshold'], linha['support'] + 1e-9,
                                 f"limiar acima do suporte em {set(linha['itemsets'])}")

    def test_n_trees_bate_com_a_contagem_direta(self):
        consolidado = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        for _, linha in consolidado.iterrows():
            itens = set(linha['itemsets'])
            direta = sum(1 for transacao in self.matriz if itens <= set(transacao))
            self.assertEqual(int(linha['n_trees']), direta,
                             f"n_trees diverge da contagem direta em {itens}")

    def test_nenhum_itemset_em_duas_categorias(self):
        """A contradição que D4 produzia na UI: o mesmo padrão listado como
        frágil (suporte <= 0,3) e como robusto (>= 0,6). Com um suporte por
        itemset, isso passa a ser impossível por construção."""
        consolidado = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        frageis = {i for i, s in zip(consolidado['itemsets'], consolidado['support']) if s <= 0.3}
        robustos = {i for i, s in zip(consolidado['itemsets'], consolidado['support']) if s >= 0.6}
        self.assertEqual(frageis & robustos, set())

    def test_entrada_vazia_devolve_quadro_com_as_colunas(self):
        vazio = SubtreeMiner.consolidate_fpmax_results(pd.DataFrame(), 4)
        self.assertTrue(vazio.empty)
        for coluna in ('itemsets', 'support', 'min_support_threshold',
                       'max_support_threshold', 'n_trees'):
            self.assertIn(coluna, vazio.columns)

    def test_ordem_e_deterministica(self):
        """Sem ordem estável, o CSV muda entre execuções e nenhuma figura é
        reproduzível a partir do commit (D14)."""
        primeira = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        for _ in range(3):
            outra = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
            self.assertEqual(list(primeira['itemsets']), list(outra['itemsets']))
            self.assertEqual(list(primeira['support']), list(outra['support']))

    def test_ordenado_por_suporte_decrescente(self):
        consolidado = SubtreeMiner.consolidate_fpmax_results(self.bruto.copy(), len(self.matriz))
        suportes = list(consolidado['support'])
        self.assertEqual(suportes, sorted(suportes, reverse=True))


class TestSemanticaAntigaEraOLimiar(unittest.TestCase):
    """Caracteriza o que a lógica antiga produzia, para que a diferença fique
    escrita em teste e não só no ledger."""

    def test_logica_antiga_repetia_o_itemset_com_suportes_diferentes(self):
        matriz = [['A', 'B'], ['A', 'B'], ['A'], ['A']]
        bruto = varredura(matriz)
        antiga = bruto.copy()
        antiga['support'] = antiga['min_support_threshold']   # a linha que D4 descreve

        por_itemset = antiga.groupby('itemsets', sort=False)['support'].agg(set)
        ambiguos = [i for i, valores in por_itemset.items() if len(valores) > 1]
        self.assertTrue(ambiguos, "a varredura tem de repetir itemsets — é a premissa de D4")

        nova = SubtreeMiner.consolidate_fpmax_results(bruto.copy(), len(matriz))
        self.assertEqual(len(nova), len(por_itemset))
        self.assertEqual(nova['itemsets'].nunique(), len(nova))


if __name__ == '__main__':
    unittest.main()
