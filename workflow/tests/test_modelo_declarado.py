"""
M7.3 — o manifesto declara o modelo que a ferramenta **rodou**, lido do log dela.

`model` continua sendo o literal da linha de comando (`GTR+G`). `modelo_efetivo`
é o que o IQ-TREE (`.iqtree`), o RAxML-NG (`.raxml.log`) e o FastTree (stderr)
dizem ter usado — que não é o mesmo: `GTR+G` vira `GTR+F+G4` num e
`GTR+FO+G4m` no outro, e o FastTree `-gtr` é GTR+CAT20, não "GTR sem variação
de taxa" como o projeto registrava.

Executar com:

    python -m unittest workflow.tests.test_modelo_declarado
"""

from __future__ import annotations

import os
import shutil
import tempfile
import unittest

from workflow.tree_construction.builder import TreeBuilder, modelo_declarado
from workflow.utils import tool_runs
from workflow.utils.external_tools import resolve_tool


class TestLeitura(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_log_ausente_ou_sem_linha_e_none(self):
        """Nunca um palpite: sem a linha, sem modelo declarado (regra 5)."""
        self.assertIsNone(modelo_declarado(os.path.join(self.dir, "nao_existe"), r"Model:\s*(\S+)"))
        caminho = os.path.join(self.dir, "x.log")
        with open(caminho, "w") as fh:
            fh.write("nada aqui\n")
        self.assertIsNone(modelo_declarado(caminho, r"Model:\s*(\S+)"))

    def test_le_a_linha_do_log(self):
        caminho = os.path.join(self.dir, "x.raxml.log")
        with open(caminho, "w") as fh:
            fh.write("Analysis options:\n  Model: GTR+FO+G4m\n")
        self.assertEqual(modelo_declarado(caminho, r"^\s*Model:\s*(\S+)"), "GTR+FO+G4m")


def _aln(ntax=6):
    from Bio import AlignIO
    return AlignIO.read(os.path.join("projects", "Zika_21seq_validacao", "out", "Align",
                                     "dataset_final_mafft.aln"), "fasta")[:ntax]


class TestFerramentasDeVerdade(unittest.TestCase):
    """Binários reais, 6 táxons de Zika. Pulados se a ferramenta faltar."""

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.dir, "Trees"))

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def _run(self, ferramenta):
        return tool_runs.execucoes()[ferramenta]["runs"][0]

    @unittest.skipUnless(resolve_tool("iqtree"), "IQ-TREE ausente")
    def test_iqtree_declara_gtr_f_g4(self):
        saida = os.path.join(self.dir, "Trees", "t_iqtree.nexus")
        TreeBuilder().iqtree_constructor(_aln(), saida)
        self.assertEqual(tool_runs.execucoes()["iqtree"]["model"], "GTR+G")
        self.assertEqual(self._run("iqtree")["modelo_efetivo"], "GTR+F+G4")

    @unittest.skipUnless(resolve_tool("fasttree"), "FastTree ausente")
    def test_fasttree_declara_cat(self):
        saida = os.path.join(self.dir, "Trees", "t_fasttree.nexus")
        TreeBuilder().fasttree_constructor(_aln(), saida)
        self.assertIn("CAT approximation with 20 rate categories",
                      self._run("fasttree")["modelo_efetivo"])

    @unittest.skipUnless(resolve_tool("raxml-ng"), "RAxML-NG ausente")
    def test_raxml_declara_frequencias_por_ml(self):
        saida = os.path.join(self.dir, "Trees", "t_raxml.nexus")
        TreeBuilder(raxml_threads=1).raxml_ng_constructor(_aln(), saida)
        self.assertEqual(self._run("raxml-ng")["modelo_efetivo"], "GTR+FO+G4m")


if __name__ == "__main__":
    unittest.main()
