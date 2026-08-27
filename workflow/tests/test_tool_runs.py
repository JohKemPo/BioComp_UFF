"""
Testes do registro de chamadas de ferramenta externa — DEC-045.

O defeito que estes testes existem para impedir não é de cálculo: é de
**ligação**. `ExecutionManifest.register_tool_run` existia desde M2.5, com
docstring justificando-se por D17 e com teste de unidade próprio — e nenhum
ponto do pipeline a chamava. O teste passava porque chamava o método
diretamente; o artefato em disco saía com `tools_invoked: {}` de toda execução,
nas duas máquinas do projeto.

A lição é a de sempre neste repositório: **um teste que exercita a função e não
o caminho não prova que o caminho existe.** Por isso os testes abaixo entram
pelo construtor de árvore, com a ferramenta externa substituída, e conferem que
o registro aconteceu — não que ele funciona quando chamado.

Executar com:

    python -m unittest workflow.tests.test_tool_runs
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from workflow.utils import tool_runs


def _alinhamento():
    """Quatro táxons, o mínimo para uma árvore ter topologia interna."""
    return MultipleSeqAlignment([
        SeqRecord(Seq("ACGTACGTAA"), id="t1"),
        SeqRecord(Seq("ACGTACGTAG"), id="t2"),
        SeqRecord(Seq("ACGAACGTTA"), id="t3"),
        SeqRecord(Seq("ACGAACGTTG"), id="t4"),
    ])


NEWICK = "((t1:0.1,t2:0.1):0.1,(t3:0.1,t4:0.1):0.1);"


class TestRegistro(unittest.TestCase):
    """O coletor em si: acumula, não sobrescreve, e é limpável."""

    def setUp(self):
        tool_runs.limpar()

    def tearDown(self):
        tool_runs.limpar()

    def test_registro_vazio_por_padrao(self):
        self.assertEqual(tool_runs.execucoes(), {})

    def test_acumula_chamadas_da_mesma_ferramenta(self):
        tool_runs.registrar("raxml-ng", ["raxml-ng", "--msa", "a"], saida="a.nexus")
        tool_runs.registrar("raxml-ng", ["raxml-ng", "--msa", "b"], saida="b.nexus")
        runs = tool_runs.execucoes()["raxml-ng"]["runs"]
        self.assertEqual([r["saida"] for r in runs], ["a.nexus", "b.nexus"])

    def test_parametros_nulos_nao_entram(self):
        """`seed=None` é o caso do MrBayes: não há semente. Gravar `None` seria
        indistinguível de "ninguém registrou" — e a regra 5 do projeto diz que
        'não aplicável' nunca é um valor inventado."""
        tool_runs.registrar("mrbayes", ["mb"], seed=None, ngen=100000)
        entrada = tool_runs.execucoes()["mrbayes"]
        self.assertNotIn("seed", entrada)
        self.assertEqual(entrada["ngen"], 100000)

    def test_execucoes_devolve_copia(self):
        """Quem lê não pode alterar o registro por acidente."""
        tool_runs.registrar("mafft", ["mafft", "x"])
        copia = tool_runs.execucoes()
        copia["mafft"]["runs"].append({"command": ["intruso"]})
        self.assertEqual(len(tool_runs.execucoes()["mafft"]["runs"]), 1)


class TestPipelineRegistra(unittest.TestCase):
    """
    O caminho de verdade: construir uma árvore registra a chamada.

    A ferramenta externa é substituída — o que está sob teste é a instrumentação
    do pipeline, não o RAxML. Mas a substituição **grava o arquivo de árvore**
    que o método real espera encontrar, de modo que o construtor percorre o
    mesmo caminho que percorreria numa execução de verdade.
    """

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        self.trees = os.path.join(self.dir, "Trees")
        os.makedirs(self.trees)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def _builder(self, **config):
        from workflow.tree_construction.builder import TreeBuilder
        return TreeBuilder(**config)

    def _saida(self, nome):
        return os.path.join(self.trees, f"{nome}.nexus")

    @staticmethod
    def _fake_run(sufixo_arvore):
        """Substitui o processo externo: grava a árvore que o método vai ler."""
        def executar(cmd, *args, **kwargs):
            prefixo = None
            for flag in ("--prefix", "-pre"):
                if flag in cmd:
                    prefixo = cmd[cmd.index(flag) + 1]
            if prefixo:
                with open(prefixo + sufixo_arvore, "w") as handle:
                    handle.write(NEWICK + "\n")
            return subprocess.CompletedProcess(cmd, 0, stdout=NEWICK + "\n", stderr="")
        return executar

    def test_raxml_registra_semente_e_paralelizacao(self):
        """D17: a paralelização muda a topologia com a mesma semente, então
        `threads` e `workers` são metadado do resultado."""
        saida = self._saida("t_raxml")
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value="/opt/env/bin/raxml-ng"), \
             mock.patch("subprocess.run", self._fake_run(".raxml.bestTree")):
            self._builder(random_seed=777, raxml_threads=3).raxml_ng_constructor(
                _alinhamento(), saida)

        entrada = tool_runs.execucoes()["raxml-ng"]
        self.assertEqual(entrada["seed"], 777)
        self.assertEqual(entrada["threads"], 3)
        self.assertEqual(entrada["workers"], 1)
        comando = entrada["runs"][0]["command"]
        self.assertIn("--workers", comando)
        self.assertNotIn("auto", comando)
        self.assertEqual(entrada["runs"][0]["saida"], saida)

    def test_iqtree_registra_semente_e_uma_thread(self):
        """D21 — a busca de ML roda em `-nt 1`. Medido: com `-nt N`, três
        repetições da mesma semente devolvem três topologias, porque a ordem
        das reduções de ponto flutuante decide entre ótimos empatados e o
        IQ-TREE não tem equivalente ao `--workers 1` do RAxML-NG."""
        saida = self._saida("t_iqtree")
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value="/opt/env/bin/iqtree3"), \
             mock.patch("subprocess.run", self._fake_run(".treefile")):
            self._builder(random_seed=777, iqtree_threads=3).iqtree_constructor(
                _alinhamento(), saida)

        entrada = tool_runs.execucoes()["iqtree"]
        comando = entrada["runs"][0]["command"]
        self.assertEqual(entrada["seed"], 777)
        self.assertEqual(entrada["threads"], 1)
        self.assertEqual(comando[comando.index("-nt") + 1], "1")
        # O valor configurado não some: ele explica a diferença para quem
        # comparar com uma execução anterior.
        self.assertEqual(entrada["threads_configurados"], 3)

    def test_fasttree_registra_ainda_sem_semente(self):
        """O FastTree não aceita semente nesta chamada. Registrar assim mesmo é
        o que distingue 'não se aplica' de 'ninguém registrou'."""
        saida = self._saida("t_fasttree")
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value="/opt/env/bin/FastTree"), \
             mock.patch("subprocess.run", self._fake_run(".nwk")):
            self._builder().fasttree_constructor(_alinhamento(), saida)

        entrada = tool_runs.execucoes()["fasttree"]
        self.assertNotIn("seed", entrada)
        self.assertEqual(len(entrada["runs"]), 1)

    def test_dois_bracos_do_delineamento_geram_dois_registros(self):
        """É o caso real: dois alinhadores, o mesmo inferidor, duas árvores.
        Um registro por ferramenta perderia metade delas."""
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value="/opt/env/bin/raxml-ng"), \
             mock.patch("subprocess.run", self._fake_run(".raxml.bestTree")):
            for braco in ("mafft", "clustalo"):
                self._builder(random_seed=777).raxml_ng_constructor(
                    _alinhamento(), self._saida(f"t_{braco}_raxml"))

        runs = tool_runs.execucoes()["raxml-ng"]["runs"]
        self.assertEqual(len(runs), 2)
        self.assertEqual([os.path.basename(r["saida"]) for r in runs],
                         ["t_mafft_raxml.nexus", "t_clustalo_raxml.nexus"])


if __name__ == "__main__":
    unittest.main()
