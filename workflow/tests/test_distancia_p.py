"""
D28 — testes de `workflow.tree_construction.distancia_p`.

Caracterização, casos-limite e regressão contra `ParsimonyScorer` original
do Biopython (sem lacuna/ambiguidade, os dois devem concordar byte a byte).
"""

import unittest
from io import StringIO

from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, ParsimonyScorer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Phylo as Phylo

from workflow.tree_construction.distancia_p import (
    FitchAmbiguidadeIupac,
    SemSitioComparavel,
    matriz_distancia_p,
    p_distancia_par,
)

#: Árvore com comprimento de ramo explícito — sem isso, `root_at_midpoint`
#: do Biopython quebra em árvores sem nenhum comprimento (bug de terceiros,
#: não deste módulo; contornado aqui só no teste).
NEWICK_4_TAXONS = "(((A:1,B:1):1,C:1):1,D:1);"


def _alinhamento(seqs: dict) -> MultipleSeqAlignment:
    return MultipleSeqAlignment([SeqRecord(Seq(s), id=nome) for nome, s in seqs.items()])


def _arvore(newick: str = NEWICK_4_TAXONS):
    return Phylo.read(StringIO(newick), "newick")


class TestPDistanciaPar(unittest.TestCase):
    def test_sequencias_identicas_distancia_zero(self):
        d, n = p_distancia_par("acgtacgt", "acgtacgt")
        self.assertEqual(d, 0.0)
        self.assertEqual(n, 8)

    def test_oraculo_independente_sem_lacuna(self):
        import random

        random.seed(7)
        alfabeto = "acgt"
        s1 = "".join(random.choice(alfabeto) for _ in range(200))
        s2 = list(s1)
        posicoes = random.sample(range(200), 30)
        for p in posicoes:
            s2[p] = random.choice([c for c in alfabeto if c != s1[p]])
        s2 = "".join(s2)
        d, n = p_distancia_par(s1, s2)
        dif_oraculo = sum(1 for a, b in zip(s1, s2) if a != b)
        self.assertEqual(n, 200)
        self.assertEqual(dif_oraculo, 30)
        self.assertAlmostEqual(d, 30 / 200, places=12)

    def test_denominador_e_sitios_comparaveis_nao_alinhamento_inteiro(self):
        """O bug D28: Biopython usa `len(seq1)` como denominador, sempre."""
        genoma = "acgt" * 250  # 1000 sítios
        fragmento = "-" * 300 + genoma[300:600] + "-" * 400
        d, n = p_distancia_par(genoma, fragmento)
        self.assertEqual(n, 300)  # só os sítios com base em ambos os lados
        self.assertEqual(d, 0.0)  # concordam em todos

    def test_diferencas_contam_so_nos_sitios_comparaveis(self):
        genoma = "acgt" * 250
        fragmento = list("-" * 300 + genoma[300:600] + "-" * 400)
        for p in (310, 320, 330):
            fragmento[p] = "g" if genoma[p] != "g" else "c"
        d, n = p_distancia_par(genoma, "".join(fragmento))
        self.assertEqual(n, 300)
        self.assertAlmostEqual(d, 3 / 300, places=12)

    def test_zero_sitios_comparaveis_levanta_excecao_nao_numero(self):
        """Regra 5 do projeto: 'não aplicável' nunca vira 0.0 nem 1.0."""
        with self.assertRaises(SemSitioComparavel):
            p_distancia_par("-" * 10, "acgtacgtac")

    def test_n_conta_como_sem_informacao_igual_a_lacuna(self):
        d, n = p_distancia_par("acgtnnnnac", "acgtacgtac")
        self.assertEqual(n, 6)

    def test_comprimentos_diferentes_e_erro(self):
        with self.assertRaises(ValueError):
            p_distancia_par("acgt", "acgtac")


class TestMatrizDistanciaP(unittest.TestCase):
    def test_matriz_compativel_com_distance_tree_constructor(self):
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

        aln = _alinhamento({"A": "acgtacgt", "B": "acgtacgc", "C": "acggacgt", "D": "acgtccgt"})
        dm = matriz_distancia_p(aln)
        # não deve levantar; e o resultado deve servir de entrada normal ao NJ/UPGMA
        tree_nj = DistanceTreeConstructor().nj(dm)
        tree_upgma = DistanceTreeConstructor().upgma(dm)
        self.assertEqual(tree_nj.count_terminals(), 4)
        self.assertEqual(tree_upgma.count_terminals(), 4)

    def test_matriz_levanta_sem_sitio_comparavel_em_vez_de_construir_arvore_errada(self):
        aln = _alinhamento({"A": "acgtacgt", "B": "----acgt", "C": "acgt----", "D": "acgtacgt"})
        with self.assertRaises(SemSitioComparavel):
            matriz_distancia_p(aln)


class TestFitchAmbiguidadeIupac(unittest.TestCase):
    def test_regressao_sem_lacuna_identico_ao_biopython_original(self):
        aln = _alinhamento({"A": "acgtacgt", "B": "acgtacgc", "C": "acgtcggt", "D": "acggacgt"})
        original = ParsimonyScorer().get_score(_arvore(), aln)
        novo = FitchAmbiguidadeIupac().get_score(_arvore(), aln)
        self.assertEqual(original, novo)

    def test_lacuna_nao_custa_passo_quando_demais_concordam(self):
        aln = _alinhamento({"A": "acgt", "B": "acgt", "C": "acgt", "D": "ac-t"})
        score = FitchAmbiguidadeIupac().get_score(_arvore(), aln)
        self.assertEqual(score, 0)

    def test_biopython_original_cobra_passo_pela_lacuna_mostrando_o_defeito(self):
        """Documenta o defeito D28 diretamente: o scorer original erra aqui."""
        aln = _alinhamento({"A": "acgt", "B": "acgt", "C": "acgt", "D": "ac-t"})
        original = ParsimonyScorer().get_score(_arvore(), aln)
        self.assertGreater(original, 0, "o Biopython original deveria cobrar >0 pela lacuna")

    def test_divergencia_real_nao_e_mascarada_pela_lacuna(self):
        aln = _alinhamento({"A": "acgt", "B": "acgt", "C": "tcgt", "D": "ac-t"})
        score = FitchAmbiguidadeIupac().get_score(_arvore(), aln)
        self.assertEqual(score, 1)

    def test_coluna_so_de_lacuna_e_n_nao_custa_nada(self):
        aln = _alinhamento({"A": "a-", "B": "an", "C": "a-", "D": "an"})
        score = FitchAmbiguidadeIupac().get_score(_arvore(), aln)
        self.assertEqual(score, 0)

    def test_ambiguidade_parcial_r_intersecta_so_purinas(self):
        # r = A ou G. Coluna: A, G, R, A -> A/G intersectam bem com R, sem custo
        aln = _alinhamento({"A": "a", "B": "g", "C": "r", "D": "a"})
        score = FitchAmbiguidadeIupac().get_score(_arvore(), aln)
        # A,B ja divergem (A vs G) independente de R: pelo menos 1 passo
        self.assertGreaterEqual(score, 1)

    def test_sankoff_com_matriz_delega_ao_biopython_original(self):
        """Fora do escopo de D28 — o pipeline nunca usa Sankoff hoje."""
        import numpy as np
        from Bio.Phylo.TreeConstruction import _Matrix

        aln = _alinhamento({"A": "ac", "B": "ac", "C": "ac", "D": "ac"})
        nomes = ["a", "c", "g", "t"]
        matriz = _Matrix(nomes, [[0], [1, 0], [1, 1, 0], [1, 1, 1, 0]])
        scorer = FitchAmbiguidadeIupac(matriz)
        score = scorer.get_score(_arvore(), aln)
        self.assertEqual(score, 0)


if __name__ == "__main__":
    unittest.main()
