"""
Enraizamento explícito por grupo externo — M2.3 / D3.

A regra que estes testes fixam: **enraizar errado é pior que não enraizar**.
Quando o grupo externo não permite uma raiz única, o resultado é uma recusa com
motivo declarado, e não uma árvore enraizada em algum lugar plausível.

Executar com:

    python -m unittest workflow.tests.test_rooting
"""

from __future__ import annotations

import io
import unittest

from Bio import Phylo

from workflow.stability.rooting import (root_at_outgroup, root_tree_set,
                                        outgroup_from_classifier)


def arvore(texto):
    return Phylo.read(io.StringIO(texto), "newick")


def rotulos(t):
    return sorted(x.name for x in t.get_terminals() if x.name)


class TestEnraizamentoPorGrupoExterno(unittest.TestCase):

    def setUp(self):
        # A, B, C são o grupo interno; X e Y são o grupo externo, e formam clado.
        self.texto = "(((A,B),C),(X,Y));"

    def test_enraiza_quando_o_grupo_externo_e_clado(self):
        enraizada, relatorio = root_at_outgroup(arvore(self.texto), ["X", "Y"])
        self.assertIsNotNone(enraizada)
        self.assertTrue(relatorio.rooted)
        self.assertTrue(relatorio.monophyletic)
        self.assertIsNone(relatorio.reason)
        self.assertEqual(relatorio.outgroup_found, frozenset({"X", "Y"}))

    def test_recusa_quando_o_grupo_externo_nao_e_monofiletico(self):
        """Grupo externo espalhado pela árvore dá mais de uma aresta candidata a
        raiz. Escolher uma seria arbitrário, então não se escolhe."""
        espalhado = arvore("((A,X),(B,(C,Y)));")
        enraizada, relatorio = root_at_outgroup(espalhado, ["X", "Y"])
        self.assertIsNone(enraizada)
        self.assertFalse(relatorio.rooted)
        self.assertFalse(relatorio.monophyletic)
        self.assertIn("monofilético", relatorio.reason)

    def test_permite_forcar_quando_declarado(self):
        """A recusa pode ser desligada para exploração — mas é preciso pedir."""
        espalhado = arvore("((A,X),(B,(C,Y)));")
        enraizada, relatorio = root_at_outgroup(espalhado, ["X", "Y"],
                                                require_monophyletic=False)
        self.assertIsNotNone(enraizada)
        self.assertTrue(relatorio.rooted)
        self.assertFalse(relatorio.monophyletic)

    def test_nao_modifica_a_arvore_original(self):
        """A mesma árvore precisa poder ser analisada enraizada e não enraizada
        sem que uma análise contamine a outra."""
        original = arvore(self.texto)
        antes = rotulos(original)
        raiz_antes = len(original.root.clades)
        root_at_outgroup(original, ["X", "Y"])
        self.assertEqual(rotulos(original), antes)
        self.assertEqual(len(original.root.clades), raiz_antes)

    def test_grupo_externo_vazio_e_recusado(self):
        """O enraizamento é declarado, nunca inferido do próprio dado."""
        _, relatorio = root_at_outgroup(arvore(self.texto), [])
        self.assertFalse(relatorio.rooted)
        self.assertIn("declarado", relatorio.reason)

    def test_grupo_externo_ausente_da_arvore_e_recusado(self):
        _, relatorio = root_at_outgroup(arvore(self.texto), ["Z1", "Z2"])
        self.assertFalse(relatorio.rooted)
        self.assertIn("nenhum táxon", relatorio.reason)
        self.assertEqual(relatorio.outgroup_missing, frozenset({"Z1", "Z2"}))

    def test_grupo_externo_igual_a_arvore_inteira_e_recusado(self):
        _, relatorio = root_at_outgroup(arvore("(X,Y);"), ["X", "Y"])
        self.assertFalse(relatorio.rooted)
        self.assertIn("grupo interno", relatorio.reason)

    def test_ausencia_parcial_e_registrada_sem_impedir(self):
        """Um táxon do grupo externo que não entrou naquele alinhamento é fato a
        registrar, não motivo para recusar — desde que sobre grupo externo."""
        enraizada, relatorio = root_at_outgroup(arvore(self.texto), ["X", "Y", "Z"])
        self.assertIsNotNone(enraizada)
        self.assertEqual(relatorio.outgroup_missing, frozenset({"Z"}))

    def test_rotulo_truncado_casa_com_o_integro(self):
        """D13: IQ-TREE e RAxML truncam o acesso em 10 caracteres. O grupo
        externo é declarado uma vez e precisa valer nas duas grafias."""
        truncada = arvore("(((A.1,B.1),C.1),(NC_008291.,NC_008030.));")
        enraizada, relatorio = root_at_outgroup(
            truncada, ["NC_008291.1", "NC_008030.1"])
        self.assertIsNotNone(enraizada)
        self.assertEqual(relatorio.outgroup_found,
                         frozenset({"NC_008291", "NC_008030"}))
        self.assertEqual(relatorio.outgroup_missing, frozenset())

    def test_um_unico_taxon_externo_enraiza(self):
        enraizada, relatorio = root_at_outgroup(arvore("((A,B),(C,X));"), ["X"])
        self.assertIsNotNone(enraizada)
        self.assertIsNone(relatorio.monophyletic,
                          "monofilia não faz sentido para um táxon só")


class TestConjuntoInteiro(unittest.TestCase):

    def test_enraizamento_comum_produz_a_mesma_raiz(self):
        """O ponto de M2.3: a mesma raiz em todos os métodos. Duas escritas da
        mesma topologia, com raízes diferentes no arquivo, têm de convergir para
        o mesmo clado externo depois de enraizadas."""
        conjunto = {
            "mafft_fasttree": arvore("(((A,B),C),(X,Y));"),
            "mafft_iqtree": arvore("((X,Y),((A,B),C));"),
        }
        enraizadas, relatorios = root_tree_set(conjunto, ["X", "Y"])
        self.assertEqual(len(enraizadas), 2)
        self.assertTrue(all(r.rooted for r in relatorios.values()))

        externos = []
        for t in enraizadas.values():
            filhos = [sorted(n.name for n in c.get_terminals() if n.name)
                      for c in t.root.clades]
            externos.append(sorted(filhos))
        self.assertEqual(externos[0], externos[1])

    def test_relatorio_cobre_todas_inclusive_as_que_falharam(self):
        """Se um método não enraizou, a análise enraizada daquele conjunto não é
        comparável — e quem decide isso precisa ver o quadro inteiro."""
        conjunto = {
            "mafft_fasttree": arvore("(((A,B),C),(X,Y));"),
            "mafft_nj_distance": arvore("((A,X),(B,(C,Y)));"),
        }
        enraizadas, relatorios = root_tree_set(conjunto, ["X", "Y"])
        self.assertEqual(set(enraizadas), {"mafft_fasttree"})
        self.assertEqual(set(relatorios), {"mafft_fasttree", "mafft_nj_distance"})
        self.assertFalse(relatorios["mafft_nj_distance"].rooted)

    def test_resumo_e_serializavel(self):
        _, relatorios = root_tree_set({"x": arvore("(((A,B),C),(X,Y));")}, ["X", "Y"])
        import json
        json.dumps(relatorios["x"].summary())   # não levanta


class TestGrupoExternoPorClassificador(unittest.TestCase):

    def test_complemento_do_grupo_interno(self):
        """Nos experimentos de Variola declara-se o grupo INTERNO (VARV) e o
        externo é o complemento — camelpox, cowpox, taterapox."""
        t = arvore("(((v1,v2),v3),(cmlv1,tatv1));")
        classificador = lambda nome: "VARV" if nome.startswith("v") else "OUTRO"
        externo = outgroup_from_classifier(t, classificador, ["VARV"])
        self.assertEqual(externo, frozenset({"cmlv1", "tatv1"}))

    def test_enraiza_com_o_grupo_derivado(self):
        t = arvore("(((v1,v2),v3),(cmlv1,tatv1));")
        classificador = lambda nome: "VARV" if nome.startswith("v") else "OUTRO"
        externo = outgroup_from_classifier(t, classificador, ["VARV"])
        enraizada, relatorio = root_at_outgroup(t, externo)
        self.assertIsNotNone(enraizada)
        self.assertTrue(relatorio.monophyletic)


if __name__ == "__main__":
    unittest.main()
