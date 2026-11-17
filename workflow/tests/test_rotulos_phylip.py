"""
Testes da restauração de rótulos truncados pelo PHYLIP (D13, lado do pipeline).
"""
import os
import tempfile
import unittest
from unittest import mock
from io import StringIO

from Bio import Phylo

from workflow.subtree_construction.builder import SubtreeBuilder
from workflow.tree_construction.builder import restaurar_rotulos_phylip


def _arvore(newick):
    return Phylo.read(StringIO(newick), "newick")


def _pontas(arv):
    return sorted(t.name for t in arv.get_terminals())


class RestaurarRotulosPhylipTest(unittest.TestCase):

    def test_rotulo_completo_nao_muda(self):
        arv = _arvore("(NC_008030.1:1,AY009089.1:1,(L22579.1:1,Y16780.1:1):1);")
        mudou = restaurar_rotulos_phylip(
            arv, ["NC_008030.1", "AY009089.1", "L22579.1", "Y16780.1"])
        self.assertEqual(mudou, {})
        self.assertEqual(_pontas(arv), ["AY009089.1", "L22579.1", "NC_008030.1", "Y16780.1"])

    def test_rotulo_truncado_e_restaurado_inclusive_refseq(self):
        arv = _arvore("(NC_001611.:1,NC_008030.:1,(L22579.1:1,AY009089.1:1):1);")
        aln = ["NC_001611.1", "NC_008030.1", "L22579.1", "AY009089.1"]
        mudou = restaurar_rotulos_phylip(arv, aln)
        self.assertEqual(mudou, {"NC_001611.": "NC_001611.1", "NC_008030.": "NC_008030.1"})
        self.assertEqual(_pontas(arv), sorted(aln))

    def test_ponta_sem_correspondente_falha_alto(self):
        arv = _arvore("(NC_001611.:1,XX_999999.:1);")
        with self.assertRaises(ValueError) as ctx:
            restaurar_rotulos_phylip(arv, ["NC_001611.1", "NC_008030.1"])
        self.assertIn("XX_999999.", str(ctx.exception))

    def test_prefixo_apenas_parcial_nao_casa(self):
        arv = _arvore("(NC_0016:1,NC_008030.:1);")
        with self.assertRaises(ValueError):
            restaurar_rotulos_phylip(arv, ["NC_001611.1", "NC_008030.1"])

    def test_colisao_de_prefixo_nomeia_os_rotulos(self):
        arv = _arvore("(NC_001611.:1,NC_008030.:1);")
        with self.assertRaises(ValueError) as ctx:
            restaurar_rotulos_phylip(arv, ["NC_001611.1", "NC_001611.2", "NC_008030.1"])
        msg = str(ctx.exception)
        self.assertIn("NC_001611.1", msg)
        self.assertIn("NC_001611.2", msg)

    def test_topologia_preservada(self):
        arv = _arvore("((NC_001611.:1,NC_008030.:2):3,(L22579.1:1,AY009089.1:1):1);")
        antes = [(t.name, arv.distance(t)) for t in arv.get_terminals()]
        restaurar_rotulos_phylip(arv, ["NC_001611.1", "NC_008030.1", "L22579.1", "AY009089.1"])
        depois = [(t.name, arv.distance(t)) for t in arv.get_terminals()]
        self.assertEqual([d for _, d in antes], [d for _, d in depois])


class CacheDeMetadadoTest(unittest.TestCase):

    def _construir(self, cache, newick):
        with tempfile.TemporaryDirectory() as d:
            caminho = os.path.join(d, "t1.nwk")
            with open(caminho, "w") as f:
                f.write(newick)
            os.makedirs(os.path.join(d, "out", "Subtrees"))
            b = SubtreeBuilder(cache, input_format="newick", output_format="nexus",
                               output_path=os.path.join(d, "out"), resume_infos=False)
            with mock.patch("workflow.subtree_construction.builder.download_sequences") as dl:
                b.subtree_constructor(caminho, "t1.nwk")
        return dl

    def test_rotulo_truncado_ja_baixado_nao_redispara_download(self):
        cache = {"NC_008030", "AY009089"}
        dl = self._construir(cache, "(NC_008030.:1,AY009089.1:1,L22579.1:1);")
        novos = [q.name for q in dl.call_args[0][1]]
        self.assertEqual(novos, ["L22579.1"])
        self.assertIn("L22579", cache)

    def test_tudo_ja_baixado_nao_chama_download(self):
        dl = self._construir({"NC_008030", "AY009089"}, "(NC_008030.:1,AY009089.1:1);")
        self.assertEqual(dl.call_count, 0)


if __name__ == "__main__":
    unittest.main()
