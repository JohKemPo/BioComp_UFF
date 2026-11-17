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
    limiar_sitios,
    pares_sem_sitio_comparavel,
)

#: Sem limiar de U4: só a sobreposição zero é indefinida (comportamento de D28/DEC-152),
#: para os testes de sequências curtas que não tratam do limiar.
SEM_LIMIAR = dict(fracao_minima=0, piso_sitios=0)

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
        d, n = p_distancia_par("acgtacgt", "acgtacgt", **SEM_LIMIAR)
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
            p_distancia_par("-" * 10, "acgtacgtac", **SEM_LIMIAR)

    def test_n_conta_como_sem_informacao_igual_a_lacuna(self):
        d, n = p_distancia_par("acgtnnnnac", "acgtacgtac", **SEM_LIMIAR)
        self.assertEqual(n, 6)

    def test_comprimentos_diferentes_e_erro(self):
        with self.assertRaises(ValueError):
            p_distancia_par("acgt", "acgtac")


class TestMatrizDistanciaP(unittest.TestCase):
    def test_matriz_compativel_com_distance_tree_constructor(self):
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

        aln = _alinhamento({"A": "acgtacgt", "B": "acgtacgc", "C": "acggacgt", "D": "acgtccgt"})
        dm = matriz_distancia_p(aln, **SEM_LIMIAR)
        # não deve levantar; e o resultado deve servir de entrada normal ao NJ/UPGMA
        tree_nj = DistanceTreeConstructor().nj(dm)
        tree_upgma = DistanceTreeConstructor().upgma(dm)
        self.assertEqual(tree_nj.count_terminals(), 4)
        self.assertEqual(tree_upgma.count_terminals(), 4)

    def test_matriz_levanta_sem_sitio_comparavel_em_vez_de_construir_arvore_errada(self):
        aln = _alinhamento({"A": "acgtacgt", "B": "----acgt", "C": "acgt----", "D": "acgtacgt"})
        with self.assertRaises(SemSitioComparavel):
            matriz_distancia_p(aln, **SEM_LIMIAR)

    def test_excecao_da_matriz_nomeia_todos_os_pares(self):
        """DEC-152: antes a mensagem dizia 'este par' sem dizer qual."""
        aln = _alinhamento({"A": "acgtacgt", "B": "----acgt", "C": "acgt----",
                            "D": "acnn----"})
        with self.assertRaises(SemSitioComparavel) as ctx:
            matriz_distancia_p(aln, **SEM_LIMIAR)
        self.assertEqual(ctx.exception.pares, [("C", "B"), ("D", "B")])
        self.assertIn("C × B", str(ctx.exception))
        self.assertIn("D × B", str(ctx.exception))
        self.assertIn("on_method_failure='continue'", str(ctx.exception))

    def test_mensagem_corta_em_cinco_pares_e_conta_o_resto(self):
        seqs = {"F": "acgt----"}
        seqs.update({f"G{i}": "----acgt" for i in range(7)})
        with self.assertRaises(SemSitioComparavel) as ctx:
            matriz_distancia_p(_alinhamento(seqs), **SEM_LIMIAR)
        self.assertEqual(len(ctx.exception.pares), 7)
        self.assertIn("e mais 2", str(ctx.exception))
        self.assertIn("8 sequências", str(ctx.exception))


class TestParesSemSitioComparavel(unittest.TestCase):
    def test_oraculo_independente_contra_p_distancia_par(self):
        """A versão vetorizada concorda par a par com o laço de `p_distancia_par`."""
        import random

        random.seed(11)
        for _ in range(30):
            n, L = random.randint(2, 9), random.randint(1, 25)
            seqs = {}
            for k in range(n):
                inicio = random.randint(0, L)
                fim = random.randint(inicio, L)
                miolo = "".join(random.choice("acgtnNrY-") for _ in range(fim - inicio))
                seqs[f"s{k}"] = "-" * inicio + miolo + "-" * (L - fim)
            aln = _alinhamento(seqs)
            esperado = []
            for i in range(n):
                for j in range(i):
                    try:
                        p_distancia_par(aln[i].seq, aln[j].seq, **SEM_LIMIAR)
                    except SemSitioComparavel:
                        esperado.append((aln[i].id, aln[j].id))
            self.assertEqual(pares_sem_sitio_comparavel(aln, **SEM_LIMIAR), esperado, seqs)

    def test_sem_par_vazio(self):
        aln = _alinhamento({"A": "acgt", "B": "a--t", "C": "nnnt"})
        self.assertEqual(pares_sem_sitio_comparavel(aln, **SEM_LIMIAR), [])

    def test_uma_sequencia_nao_tem_par(self):
        self.assertEqual(pares_sem_sitio_comparavel(_alinhamento({"A": "----"}), **SEM_LIMIAR), [])

    def test_maiuscula_e_minuscula_iguais(self):
        aln = _alinhamento({"A": "ACGT----", "B": "----acgt", "C": "NNNNACGT"})
        self.assertEqual(pares_sem_sitio_comparavel(aln, **SEM_LIMIAR), [("B", "A"), ("C", "A")])


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


def _janela(inicio: int, fim: int, comprimento: int) -> str:
    return "-" * inicio + "a" * (fim - inicio) + "-" * (comprimento - fim)


class TestLimiarDeSobreposicao(unittest.TestCase):
    """U4 / DM-10: sobreposição >= max(100, 0,5 × bases da mais curta)."""

    def _par(self, a, b, **kw):
        return pares_sem_sitio_comparavel(_alinhamento({"A": a, "B": b}), **kw)

    def test_limiar_sitios_e_o_maior_entre_piso_e_fracao(self):
        self.assertEqual(limiar_sitios(150, 900), 100)
        self.assertEqual(limiar_sitios(400, 900), 200)
        self.assertEqual(limiar_sitios(900, 400), 200)

    def test_piso_decide_exatamente_no_limiar_passa(self):
        # 150 bases cada: a fração pede 75, o piso pede 100
        self.assertEqual(self._par(_janela(0, 150, 300), _janela(50, 200, 300)), [])

    def test_piso_decide_um_sitio_abaixo_falha(self):
        self.assertEqual(self._par(_janela(0, 150, 300), _janela(51, 201, 300)),
                         [("B", "A")])

    def test_fracao_decide_exatamente_no_limiar_passa(self):
        # 400 bases cada: a fração pede 200, acima do piso
        self.assertEqual(self._par(_janela(0, 400, 600), _janela(200, 600, 600)), [])

    def test_fracao_decide_um_sitio_abaixo_falha(self):
        self.assertEqual(self._par(_janela(0, 400, 601), _janela(201, 601, 601)),
                         [("B", "A")])

    def test_comprimentos_muito_diferentes_a_fracao_e_da_mais_curta(self):
        genoma = _janela(0, 1000, 1000)
        # 0,5 × 1000 = 500 reprovaria; a mais curta (400) pede só 200
        self.assertEqual(self._par(genoma, _janela(300, 700, 1000)), [])
        self.assertEqual(self._par(genoma, _janela(0, 120, 1000)), [])
        self.assertEqual(self._par(genoma, _janela(0, 99, 1000)), [("B", "A")])

    def test_sobreposicao_zero_continua_indefinida_sem_limiar(self):
        a, b = _janela(0, 50, 100), _janela(50, 100, 100)
        self.assertEqual(self._par(a, b, **SEM_LIMIAR), [("B", "A")])
        self.assertEqual(self._par(a, b), [("B", "A")])

    def test_matriz_levanta_nomeando_o_par_abaixo_do_limiar(self):
        aln = _alinhamento({"A": _janela(0, 150, 300), "B": _janela(51, 201, 300),
                            "C": _janela(0, 300, 300)})
        with self.assertRaises(SemSitioComparavel) as ctx:
            matriz_distancia_p(aln)
        self.assertEqual(ctx.exception.pares, [("B", "A")])
        self.assertIn("B × A", str(ctx.exception))

    def test_p_distancia_par_abaixo_do_limiar_levanta_nao_devolve_numero(self):
        with self.assertRaises(SemSitioComparavel):
            p_distancia_par(_janela(0, 150, 300), _janela(51, 201, 300))
        d, n = p_distancia_par(_janela(0, 150, 300), _janela(50, 200, 300))
        self.assertEqual((d, n), (0.0, 100))

    def test_parametros_alteram_o_limiar(self):
        a, b = _janela(0, 150, 300), _janela(51, 201, 300)  # sobreposição 99
        self.assertEqual(self._par(a, b, piso_sitios=99, fracao_minima=0.5), [])
        self.assertEqual(self._par(a, b, piso_sitios=0, fracao_minima=0.7), [("B", "A")])
