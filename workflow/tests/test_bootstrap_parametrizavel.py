"""
DEC-179 — contagem de réplicas de bootstrap parametrizável no `tree_config`.

Nenhum binário roda: `subprocess.run` é substituído e só a linha de comando
montada é inspecionada.

Executar com:

    python -m unittest workflow.tests.test_bootstrap_parametrizavel
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from workflow.tests.test_tool_runs import NEWICK, _alinhamento
from workflow.tree_construction.builder import TreeBuilder, bootstrap_settings
from workflow.utils import tool_runs


class TestValidacao(unittest.TestCase):
    def test_padrao_sao_strings_1000(self):
        self.assertEqual(bootstrap_settings({}),
                         {'raxml_bs_trees': '1000', 'iqtree_bs_trees': '1000'})

    def test_raxml_aceita_convergencia(self):
        for valor in ('50', 'autoMRE', 'autoMRE{1000}', 'autoMRC{500}', 'autoFC'):
            self.assertEqual(bootstrap_settings({'raxml_bs_trees': valor})['raxml_bs_trees'], valor)

    def test_raxml_recusa_forma_desconhecida(self):
        for valor in ('autoXYZ', '0', '-5', 'autoMRE{}'):
            with self.assertRaisesRegex(ValueError, 'raxml_bs_trees'):
                bootstrap_settings({'raxml_bs_trees': valor})

    def test_iqtree_recusa_convergencia_e_diz_que_e_do_raxml(self):
        with self.assertRaisesRegex(ValueError, 'RAxML-NG'):
            bootstrap_settings({'iqtree_bs_trees': 'autoMRE{1000}'})

    def test_iqtree_recusa_menos_de_1000(self):
        with self.assertRaisesRegex(ValueError, "iqtree_bs_trees.*'500'"):
            bootstrap_settings({'iqtree_bs_trees': '500'})

    def test_treebuilder_recusa_na_construcao(self):
        with self.assertRaises(ValueError):
            TreeBuilder(raxml_bs_trees='autoXYZ')


class _BaseComBinarioFalso(unittest.TestCase):
    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        self.trees = os.path.join(self.dir, "Trees")
        os.makedirs(self.trees)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    @staticmethod
    def _fake_run(sufixo_arvore):
        def executar(cmd, *args, **kwargs):
            for flag in ("--prefix", "-pre"):
                if flag in cmd:
                    with open(cmd[cmd.index(flag) + 1] + sufixo_arvore, "w") as handle:
                        handle.write(NEWICK + "\n")
            return subprocess.CompletedProcess(cmd, 0, stdout=NEWICK + "\n", stderr="")
        return executar

class TestLinhaDeComando(_BaseComBinarioFalso):
    def _comando(self, metodo, ferramenta, sufixo, **config):
        saida = os.path.join(self.trees, f"t_{metodo}.nexus")
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value=f"/opt/env/bin/{ferramenta}"), \
             mock.patch("subprocess.run", self._fake_run(sufixo)):
            builder = TreeBuilder(**config)
            getattr(builder, f"{metodo}_constructor")(_alinhamento(), saida)
        return tool_runs.execucoes()[ferramenta]["runs"][0]["command"]

    def test_sem_as_chaves_as_duas_saem_com_1000(self):
        raxml = self._comando("raxml_ng", "raxml-ng", ".raxml.bestTree")
        self.assertEqual(raxml[raxml.index("--bs-trees") + 1], "1000")
        iqtree = self._comando("iqtree", "iqtree", ".treefile")
        self.assertEqual(iqtree[iqtree.index("-bb") + 1], "1000")

    def test_raxml_recebe_automre(self):
        cmd = self._comando("raxml_ng", "raxml-ng", ".raxml.bestTree",
                            raxml_bs_trees="autoMRE{1000}")
        self.assertEqual(cmd[cmd.index("--bs-trees") + 1], "autoMRE{1000}")

    def test_iqtree_recebe_2000(self):
        cmd = self._comando("iqtree", "iqtree", ".treefile", iqtree_bs_trees="2000")
        self.assertEqual(cmd[cmd.index("-bb") + 1], "2000")

    def test_valor_efetivo_vai_para_effective_params(self):
        self._comando("raxml_ng", "raxml-ng", ".raxml.bestTree", raxml_bs_trees="autoMRE{1000}")
        self._comando("iqtree", "iqtree", ".treefile", iqtree_bs_trees="2000")
        self.assertEqual(tool_runs.parametros_efetivos()["bootstrap"],
                         {"raxml_bs_trees": "autoMRE{1000}", "iqtree_bs_trees": "2000"})


class TestControladorAteOBinario(_BaseComBinarioFalso):
    """`tree_config` -> `setattr` do controlador -> `TreeBuilder` -> linha de comando.

    Simulado: o `__init__` do controlador (efeitos de log/diretório), contornado
    com `object.__new__` e o mesmo laço de `setattr`; `_get_alignment` e
    `_modelo_do_alinhamento`; `require_tool` e `subprocess.run`.
    """

    def _controlador(self, tree_config):
        from workflow.controller.treeBuilderController import TreeBuilderController
        controlador = object.__new__(TreeBuilderController)
        for chave, valor in tree_config.items():
            setattr(controlador, chave, valor)
        controlador.count_nodes = []
        return controlador

    def _via_controlador(self, metodo, ferramenta, sufixo, tree_config):
        controlador = self._controlador(tree_config)
        saida = os.path.join(self.trees, f"t_{metodo}.nexus")
        with mock.patch.object(controlador, "_get_alignment", return_value=_alinhamento()), \
             mock.patch.object(controlador, "_modelo_do_alinhamento", return_value=None), \
             mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value=f"/opt/env/bin/{ferramenta}"), \
             mock.patch("subprocess.run", self._fake_run(sufixo)):
            getattr(controlador, f"build_tree_{metodo}")("x.fasta", "x.aln", saida, "mafft", "x.html")
        return tool_runs.execucoes()[ferramenta]["runs"][0]["command"]

    def test_raxml_recebe_do_tree_config(self):
        cmd = self._via_controlador("raxml", "raxml-ng", ".raxml.bestTree",
                                    {"raxml_bs_trees": "autoMRE{1000}"})
        self.assertEqual(cmd[cmd.index("--bs-trees") + 1], "autoMRE{1000}")

    def test_iqtree_recebe_do_tree_config(self):
        cmd = self._via_controlador("iqtree", "iqtree", ".treefile", {"iqtree_bs_trees": "2000"})
        self.assertEqual(cmd[cmd.index("-bb") + 1], "2000")

    def test_controlador_sem_as_chaves_mantem_1000(self):
        cmd = self._via_controlador("raxml", "raxml-ng", ".raxml.bestTree", {})
        self.assertEqual(cmd[cmd.index("--bs-trees") + 1], "1000")


if __name__ == "__main__":
    unittest.main()
