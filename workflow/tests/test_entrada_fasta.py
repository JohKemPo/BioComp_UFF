"""
DEC-146 — o diretório de entrada só entrega FASTA ao pipeline.

Em 2026-10-16 a reexecução de VARV-52 morreu em menos de 1s: o diretório
`*-clean` tinha um `PROVENIENCIA.md` ao lado do FASTA, o controlador o tratou
como conjunto de sequências e o MAFFT recebeu um Markdown.

Executar com:

    python -m unittest workflow.tests.test_entrada_fasta
"""

from __future__ import annotations

import os
import shutil
import tempfile
import unittest

from workflow.controller.treeBuilderController import TreeBuilderController, eh_fasta
from workflow.utils import tool_runs

FASTA = ">t1\nACGT\n>t2\nACGA\n"


class TestEhFasta(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _arquivo(self, nome, conteudo):
        caminho = os.path.join(self.dir, nome)
        with open(caminho, "w", encoding="utf-8") as fh:
            fh.write(conteudo)
        return caminho

    def test_fasta_com_e_sem_extensao(self):
        self.assertTrue(eh_fasta(self._arquivo("ds.fasta", FASTA)))
        self.assertTrue(eh_fasta(self._arquivo("ORTHOMCL1", FASTA)))

    def test_linhas_em_branco_antes_do_primeiro_registro(self):
        self.assertTrue(eh_fasta(self._arquivo("ds.fa", "\n\n  \n" + FASTA)))

    def test_markdown_vazio_e_diretorio_nao_sao_fasta(self):
        self.assertFalse(eh_fasta(self._arquivo("PROVENIENCIA.md", "# Proveniência\n\n> citação\n")))
        self.assertFalse(eh_fasta(self._arquivo("vazio.fasta", "")))
        self.assertFalse(eh_fasta(self._arquivo("arvore.dnd", "(t1,t2);\n")))
        os.mkdir(os.path.join(self.dir, "sub"))
        self.assertFalse(eh_fasta(os.path.join(self.dir, "sub")))
        self.assertFalse(eh_fasta(os.path.join(self.dir, "nao_existe")))


class TestControladorSoLeFasta(unittest.TestCase):
    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        self.entrada = os.path.join(self.dir, "data")
        self.saida = os.path.join(self.dir, "out")
        os.makedirs(self.entrada)
        for nome, conteudo in (("dataset_final.fasta", FASTA), ("ORTHOMCL7", FASTA),
                               ("PROVENIENCIA.md", "# Proveniência\nbaixado do NCBI\n"),
                               ("metadata.json", "{}")):
            with open(os.path.join(self.entrada, nome), "w", encoding="utf-8") as fh:
                fh.write(conteudo)
        os.mkdir(os.path.join(self.entrada, "subdir"))

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def test_lista_de_entrada_so_tem_fasta(self):
        c = TreeBuilderController(input_path=self.entrada, output_path=self.saida,
                                  num_threads=1, mode="advanced", output_format="nexus",
                                  aligners=["mafft"], align_method="mafft",
                                  ignore_mode=["mrbayes"])
        self.assertEqual(c.files, ["ORTHOMCL7", "dataset_final.fasta"])


if __name__ == "__main__":
    unittest.main()
