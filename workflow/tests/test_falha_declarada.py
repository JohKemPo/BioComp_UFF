"""
M7.6 — falha nunca é silenciosa.

Antes: um método que quebrava derrubava a execução, e a saída de sempre era
tirá-lo à mão via `ignore_mode` e rodar de novo. O artefato da segunda execução
dizia "ignorado" — igual a um método excluído de propósito. `execution_mode`
(D18) é calculado **antes** de rodar, então tampouco sabia de falha: listava
como executado todo método disponível e não ignorado.

Depois: cada pipeline tem um desfecho no manifesto (`inference_methods`) —
``executado``, ``reaproveitado``, ``ignorado_por_configuracao`` ou
``tentado_e_falhou`` —, os três últimos com motivo obrigatório.

Os testes de caminho entram pelo **controlador real** (`__init__` incluído),
com a ferramenta externa substituída. A falha forçada é a mais comum na
prática: **binário ausente** (`require_tool` levanta `FileNotFoundError`).

Executar com:

    python -m unittest workflow.tests.test_falha_declarada
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from workflow.utils import tool_runs
from workflow.utils.manifest import ExecutionManifest

NEWICK = "((t1:0.1,t2:0.1):0.1,(t3:0.1,t4:0.1):0.1);"
ALINHAMENTO = ">t1\nACGTACGTAA\n>t2\nACGTACGTAG\n>t3\nACGAACGTTA\n>t4\nACGAACGTTG\n"

#: Caminho absoluto com nome de usuário, como o `FileNotFoundError` real do
#: `require_tool` cita. Não pode chegar ao manifesto (D15).
CAMINHO_PRIVADO = "/home/fulano/miniconda3/envs/Phylotreeminer/bin"


def _require_tool_sem(ausente):
    """`require_tool` em que a ferramenta `ausente` não existe."""
    def resolver(nome, *a, **k):
        if nome == ausente:
            raise FileNotFoundError(
                f"'{nome}' não encontrado. Procurado em {CAMINHO_PRIVADO}/{nome}")
        return f"/opt/env/bin/{nome}"
    return resolver


#: Relatório mínimo de uma seleção `-m MF` (M7.3): a linha que o controlador
#: lê. Sem ela, a seleção falha — e é outro teste (`test_selecao_modelo`).
RELATORIO_MF = "Best-fit model according to BIC: GTR+F+G4\n"


def _fake_run(cmd, *args, **kwargs):
    """Grava a árvore que o IQ-TREE/RAxML-NG gravariam e devolve sucesso."""
    if "MF" in cmd:  # M7.3 — seleção de modelo, uma por alinhamento
        with open(cmd[cmd.index("-pre") + 1] + ".iqtree", "w") as fh:
            fh.write(RELATORIO_MF)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
    for flag, sufixo in (("-pre", ".treefile"), ("--prefix", ".raxml.support")):
        if flag in cmd:
            with open(cmd[cmd.index(flag) + 1] + sufixo, "w") as fh:
                fh.write(NEWICK + "\n")
    return subprocess.CompletedProcess(cmd, 0, stdout=NEWICK + "\n", stderr="")


class TestRegistroDeDesfecho(unittest.TestCase):
    """O coletor: enumeração fechada, motivo obrigatório, primeiro prevalece."""

    def setUp(self):
        tool_runs.limpar()

    def tearDown(self):
        tool_runs.limpar()

    def test_estado_fora_da_enumeracao_e_erro(self):
        with self.assertRaises(ValueError):
            tool_runs.registrar_metodo("iqtree", "pulado", "a.nexus")

    def test_falha_sem_motivo_e_erro(self):
        """Sem motivo, `tentado_e_falhou` volta a ser só uma ausência."""
        for estado in ("tentado_e_falhou", "ignorado_por_configuracao", "reaproveitado"):
            with self.subTest(estado=estado), self.assertRaises(ValueError):
                tool_runs.registrar_metodo("iqtree", estado, "a.nexus", motivo="  ")

    def test_primeiro_desfecho_prevalece(self):
        """O modo `advanced` visita cada método avançado duas vezes (uma por
        NJ, outra por UPGMA); a segunda visita lê a árvore da primeira e não
        pode rebaixar "executado" para "reaproveitado"."""
        self.assertTrue(tool_runs.registrar_metodo("iqtree", "executado", "a.nexus"))
        self.assertFalse(tool_runs.registrar_metodo(
            "iqtree", "reaproveitado", "a.nexus", motivo="já em disco"))
        self.assertEqual([m["estado"] for m in tool_runs.metodos()], ["executado"])

    def test_motivo_longo_e_truncado(self):
        tool_runs.registrar_metodo("raxml", "tentado_e_falhou", "a.nexus",
                                   motivo="x" * 5000)
        self.assertLess(len(tool_runs.metodos()[0]["motivo"]), 700)

    def test_limpar_esvazia_desfechos(self):
        tool_runs.registrar_metodo("iqtree", "executado", "a.nexus")
        tool_runs.limpar()
        self.assertEqual(tool_runs.metodos(), [])


class TestCaminhoDoControlador(unittest.TestCase):
    """Binário ausente num método, sucesso em outro, e o resto ignorado."""

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        self.entrada = os.path.join(self.dir, "entrada")
        self.saida = os.path.join(self.dir, "out")
        os.makedirs(self.entrada)
        with open(os.path.join(self.entrada, "ds.fasta"), "w") as fh:
            fh.write(ALINHAMENTO)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def _controlador(self, **extra):
        from workflow.controller.treeBuilderController import TreeBuilderController
        config = dict(input_path=self.entrada, output_path=self.saida, num_threads=1,
                      mode="advanced", output_format="nexus", aligners=["mafft"],
                      align_method="mafft",
                      ignore_mode=["distance", "parsimony", "fasttree", "mrbayes"])
        config.update(extra)
        c = TreeBuilderController(**config)
        # Alinhamento já em disco: o que está sob teste é a inferência, e o
        # controlador reaproveita um `.aln` não vazio sem chamar o MAFFT.
        with open(os.path.join(self.saida, "tmp", "ds_mafft.aln"), "w") as fh:
            fh.write(ALINHAMENTO)
        return c

    def _rodar_avancado(self, controlador):
        caminhos = controlador._prepare_output_paths("ds")
        fasta = os.path.join(self.entrada, "ds.fasta")
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        side_effect=_require_tool_sem("raxml-ng")), \
             mock.patch("workflow.tree_construction.modelo_substituicao.require_tool",
                        side_effect=_require_tool_sem("raxml-ng")), \
             mock.patch("subprocess.run", side_effect=_fake_run), \
             mock.patch.object(controlador, "save_tree_image"):
            return controlador._process_advanced_mode("ds", fasta, caminhos)

    @staticmethod
    def _por_metodo():
        return {m["metodo"]: m for m in tool_runs.metodos()}

    def test_politica_desconhecida_e_recusada_no_inicio(self):
        with self.assertRaises(ValueError):
            self._controlador(on_method_failure="ignorar")

    def test_padrao_fail_registra_a_falha_e_para(self):
        """Comportamento de sempre (a exceção sobe), agora com a falha
        declarada antes de subir — o `finally` do workflow drena isso."""
        controlador = self._controlador()
        self.assertEqual(controlador.on_method_failure, "fail")
        with self.assertRaises(FileNotFoundError):
            self._rodar_avancado(controlador)

        estados = self._por_metodo()
        self.assertEqual(estados["iqtree"]["estado"], "executado")
        self.assertEqual(estados["fasttree"]["estado"], "ignorado_por_configuracao")
        self.assertIn("ignore_mode", estados["fasttree"]["motivo"])
        self.assertEqual(estados["raxml"]["estado"], "tentado_e_falhou")
        self.assertIn("FileNotFoundError", estados["raxml"]["motivo"])
        self.assertIn("raxml-ng", estados["raxml"]["motivo"])
        # A execução parou no RAxML-NG: o MrBayes vem depois e não foi alcançado.
        self.assertNotIn("mrbayes", estados)

    def test_continue_segue_e_distingue_os_quatro_casos(self):
        controlador = self._controlador(on_method_failure="continue")
        construidas, arvores = self._rodar_avancado(controlador)

        estados = self._por_metodo()
        self.assertEqual(estados["iqtree"]["estado"], "executado")
        self.assertEqual(estados["raxml"]["estado"], "tentado_e_falhou")
        self.assertEqual(estados["mrbayes"]["estado"], "ignorado_por_configuracao")
        self.assertEqual(estados["nj_distance"]["estado"], "ignorado_por_configuracao")
        # Uma árvore só (IQ-TREE); o RAxML-NG não entra no conjunto nem é
        # tentado de novo na segunda passada (UPGMA) do laço.
        self.assertEqual(construidas, 1)
        self.assertEqual(arvores["mafft"]["raxml"], [])
        self.assertEqual(len(arvores["mafft"]["iqtree"]), 2)  # 2ª passada lê do disco
        raxml = [m for m in tool_runs.metodos() if m["metodo"] == "raxml"]
        self.assertEqual(len(raxml), 1)

    def test_arvore_em_disco_e_reaproveitada_nao_executada(self):
        controlador = self._controlador(on_method_failure="continue")
        self._rodar_avancado(controlador)
        tool_runs.limpar()

        # Segunda execução no mesmo diretório: a árvore do IQ-TREE já existe.
        controlador = self._controlador(on_method_failure="continue")
        self._rodar_avancado(controlador)
        self.assertEqual(self._por_metodo()["iqtree"]["estado"], "reaproveitado")

    def test_modo_basico_declara_os_avancados_como_ignorados(self):
        controlador = self._controlador(mode="basic", ignore_mode=["parsimony"])
        caminhos = controlador._prepare_output_paths("ds")
        with mock.patch.object(controlador, "save_tree_image"):
            controlador._process_auto_mode("ds", os.path.join(self.entrada, "ds.fasta"),
                                           caminhos)
        estados = self._por_metodo()
        for metodo in ("iqtree", "fasttree", "raxml", "mrbayes"):
            self.assertEqual(estados[metodo]["estado"], "ignorado_por_configuracao")
            self.assertIn("D18", estados[metodo]["motivo"])
        self.assertEqual(estados["nj_distance"]["estado"], "executado")
        self.assertEqual(estados["nj_parsimony"]["estado"], "ignorado_por_configuracao")


class TestManifestoDeclara(unittest.TestCase):
    """O que o controlador registrou chega ao `manifest.json`, higienizado."""

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def test_desfecho_chega_ao_manifesto_sem_caminho_privado(self):
        arvore = os.path.join(self.dir, "out", "Trees", "tree_ds_mafft_raxml.nexus")
        tool_runs.registrar_metodo(
            "raxml", "tentado_e_falhou", arvore, alinhador="mafft",
            motivo=f"FileNotFoundError: 'raxml-ng' não encontrado. "
                   f"Procurado em {CAMINHO_PRIVADO}/raxml-ng")
        manifesto = ExecutionManifest(project_root=self.dir, params={})
        manifesto.drain_tool_runs()
        manifesto.register_outcome("falhou", f"SystemExit: 1 ({CAMINHO_PRIVADO})")

        d = manifesto.to_dict()
        (registro,) = d["inference_methods"]
        self.assertEqual(registro["estado"], "tentado_e_falhou")
        self.assertEqual(registro["saida"], os.path.join("out", "Trees",
                                                         "tree_ds_mafft_raxml.nexus"))
        self.assertIn("raxml-ng", registro["motivo"])
        self.assertEqual(d["outcome"]["status"], "falhou")
        texto = repr(d)
        self.assertNotIn("/home/fulano", texto)
        self.assertNotIn("fulano", texto)

    def test_manifesto_parcial_diz_nulo_nao_vazio(self):
        """Antes de drenar, `None` — "ainda não se sabe", não "nenhum método"."""
        d = ExecutionManifest(project_root=self.dir, params={}).to_dict()
        self.assertIsNone(d["inference_methods"])
        self.assertIsNone(d["outcome"])

    def test_desfecho_desconhecido_e_erro(self):
        with self.assertRaises(ValueError):
            ExecutionManifest(project_root=self.dir, params={}).register_outcome("ok")


class TestRaxmlSemSuporteDeclarado(unittest.TestCase):
    """O fallback para `.raxml.bestTree` existe, mas agora é declarado."""

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.dir, "Trees"))

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def _rodar(self, sufixo):
        from Bio import AlignIO
        from io import StringIO
        from workflow.tree_construction.builder import TreeBuilder

        def executar(cmd, *a, **k):
            with open(cmd[cmd.index("--prefix") + 1] + sufixo, "w") as fh:
                fh.write(NEWICK + "\n")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        saida = os.path.join(self.dir, "Trees", "t_raxml.nexus")
        aln = AlignIO.read(StringIO(ALINHAMENTO), "fasta")
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value="/opt/env/bin/raxml-ng"), \
             mock.patch("subprocess.run", side_effect=executar):
            TreeBuilder().raxml_ng_constructor(aln, saida)
        return tool_runs.execucoes()["raxml-ng"]["runs"][0]

    def test_sem_support_declara_suporte_ausente(self):
        chamada = self._rodar(".raxml.bestTree")
        self.assertEqual(chamada["arvore_lida"], ".raxml.bestTree")
        self.assertIn("ausente", chamada["suporte"])

    def test_com_support_nao_declara_ausencia(self):
        chamada = self._rodar(".raxml.support")
        self.assertEqual(chamada["arvore_lida"], ".raxml.support")
        self.assertNotIn("suporte", chamada)

    def test_anotacao_chega_ao_manifesto_por_chamada(self):
        self._rodar(".raxml.bestTree")
        manifesto = ExecutionManifest(project_root=self.dir, params={})
        manifesto.drain_tool_runs()
        run = manifesto.to_dict()["tools_invoked"]["raxml-ng"]["runs"][0]
        self.assertEqual(run["arvore_lida"], ".raxml.bestTree")


if __name__ == "__main__":
    unittest.main()
