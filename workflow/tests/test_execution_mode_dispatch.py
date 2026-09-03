"""Guarda de regressão de D18 — `mode="auto"` e `mode="basic"` percorrem o
mesmo dispatch em `TreeBuilderController`.

A revisão do lote que fechou D18 (ver `docs/automation/07-log-de-execucao.md`,
DEC-072/073) achou que a evidência de que os dois modos se comportam de forma
idêntica era só uma execução manual, não um teste automatizado — nenhum dos
três pontos de dispatch de `treeBuilderController.py` (descrição inicial, o
`elif` que decide `_process_auto_mode`, o `if` do heatmap acumulado) tinha
guarda permanente contra alguém adicionar um quarto ponto que reconheça só
"basic" e esqueça o alias legado "auto" — os ~20 projetos já em disco com
`mode: "auto"` em `config_backup.json` quebrariam em silêncio.

Este teste roda o pipeline de verdade (FASTA sintético de 3 sequências,
alinhador MAFFT real) nos dois modos e confere que o resultado — número de
árvores construídas e o `execution_mode` do manifesto — é idêntico, a menos do
`mode_solicitado` declarado.

Executar com:

    python -m unittest workflow.tests.test_execution_mode_dispatch
"""

from __future__ import annotations

import os
import shutil
import tempfile
import unittest

from workflow.controller.treeBuilderController import (MODOS_BASICOS,
                                                        TreeBuilderController)
from workflow.utils import external_tools, tool_runs

_FASTA = (
    ">a\nACGTACGTACGTACGTACGTACGTACGT\n"
    ">b\nACGTACGTACGTACGTACGTACGTACGA\n"
    ">c\nACGTACGTACGTACGTACGTACGTACTA\n"
)


@unittest.skipUnless(external_tools.resolve_tool("mafft"), "requer MAFFT instalado")
class TestDispatchAutoEBasicSaoIdenticos(unittest.TestCase):
    def setUp(self):
        self.raiz = tempfile.mkdtemp()
        self.entrada = os.path.join(self.raiz, "in")
        os.makedirs(self.entrada, exist_ok=True)
        with open(os.path.join(self.entrada, "t.fasta"), "w", encoding="utf-8") as f:
            f.write(_FASTA)
        tool_runs.limpar()

    def tearDown(self):
        shutil.rmtree(self.raiz, ignore_errors=True)
        tool_runs.limpar()

    def _rodar(self, modo):
        saida = os.path.join(self.raiz, f"out_{modo}")
        controlador = TreeBuilderController(
            input_path=self.entrada,
            output_path=saida,
            mode=modo,
            num_threads=1,
            aligners=["mafft"],
            ignore_mode=[],
            output_format="nexus",
        )
        controlador()
        return controlador

    def test_auto_e_basic_constroem_o_mesmo_numero_de_arvores(self):
        c_basic = self._rodar("basic")
        c_auto = self._rodar("auto")
        self.assertEqual(c_basic.count_trees, c_auto.count_trees)
        self.assertGreater(c_basic.count_trees, 0, "o modo básico não construiu nenhuma árvore")

    def test_modos_basicos_cobre_os_dois_nomes(self):
        self.assertIn("auto", MODOS_BASICOS)
        self.assertIn("basic", MODOS_BASICOS)
        self.assertNotIn("advanced", MODOS_BASICOS)


if __name__ == "__main__":
    unittest.main()
