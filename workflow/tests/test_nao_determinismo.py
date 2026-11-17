"""
O manifesto declara as etapas não determinísticas e o paralelismo efetivo.

Nenhum teste roda ferramenta: os comandos são montados aqui, como o pipeline
os registraria. Declarar não é corrigir — o comando do MAFFT não muda.
"""
import json
import os
import unittest

from workflow.utils.manifest import (ExecutionManifest, etapas_nao_deterministicas,
                                     reproducibility_efetivo)

RAIZ = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _mafft(*args, saida="out/tmp/a.aln"):
    return {"mafft": {"runs": [{"command": ["mafft", *args, "X"], "saida": saida}]}}


class TestEtapasNaoDeterministicas(unittest.TestCase):

    def test_iterativo_com_varias_threads_e_declarado(self):
        itens = etapas_nao_deterministicas(_mafft(
            "--thread", "16", "--retree", "2", "--maxiterate", "2"))
        self.assertEqual(len(itens), 1)
        self.assertEqual(itens[0]["etapa"], "alinhamento")
        self.assertEqual(itens[0]["saidas_afetadas"], ["out/tmp/a.aln"])
        self.assertNotIn("alinhador", itens[0])

    def test_maxiterate_zero_nao_e_declarado(self):
        self.assertEqual(etapas_nao_deterministicas(_mafft(
            "--thread", "16", "--retree", "2", "--maxiterate", "0")), [])

    def test_uma_thread_nao_e_declarado(self):
        self.assertEqual(etapas_nao_deterministicas(_mafft(
            "--thread", "1", "--retree", "2", "--maxiterate", "2")), [])

    def test_threadit_zero_nao_e_declarado(self):
        self.assertEqual(etapas_nao_deterministicas(_mafft(
            "--thread", "16", "--threadit", "0", "--retree", "2",
            "--maxiterate", "1000")), [])

    def test_sem_mafft_e_lista_vazia(self):
        self.assertEqual(etapas_nao_deterministicas({"iqtree": {"runs": []}}), [])
        self.assertEqual(etapas_nao_deterministicas(None), [])


class TestManifestoDeclara(unittest.TestCase):

    def _manifesto(self, tools):
        import tempfile
        m = ExecutionManifest(tempfile.mkdtemp(), {})
        for nome, dados in tools.items():
            for r in dados["runs"]:
                m.register_tool_run(nome, r["command"], saida=r.get("saida"))
        return m

    def test_iqtree_efetivo_e_1_mesmo_com_16_configurado(self):
        m = self._manifesto({"iqtree": {"runs": [{"command": [
            "iqtree3", "-s", "a.phy", "-m", "GTR", "-bb", "1000", "-seed", "7", "-nt", "1"]}]}})
        m.register_reproducibility({"iqtree_threads": 16, "random_seed": 7, "raxml_threads": 8})
        r = m.to_dict()["reproducibility"]
        self.assertEqual(r["iqtree_threads"], 16)
        self.assertEqual(r["efetivo"]["iqtree_threads"], 1)
        self.assertEqual(r["efetivo"]["random_seed"], 7)
        self.assertIsNone(r["efetivo"]["raxml_threads"])
        self.assertIsNone(r["efetivo"]["mafft_threadit"])

    def test_estrategia_por_execucao_e_nota_legada(self):
        m = self._manifesto({"mafft": {"runs": [
            {"command": ["mafft", "--thread", "4", "--retree", "2", "--maxiterate", "0", "X"]},
            {"command": ["mafft", "--thread", "4", "--retree", "2", "--maxiterate", "2", "X"]}]}})
        mafft = m.to_dict()["tools_invoked"]["mafft"]
        self.assertEqual([r["estrategia"] for r in mafft["runs"]],
                         ["--retree 2 --maxiterate 0", "--retree 2 --maxiterate 2"])
        self.assertIn("estrategia_nota", mafft)
        self.assertEqual(len(m.to_dict()["nondeterminism"]), 1)

    def test_sem_estrategias_diferentes_nao_tem_nota(self):
        m = self._manifesto(_mafft("--thread", "4", "--maxiterate", "0"))
        self.assertNotIn("estrategia_nota", m.to_dict()["tools_invoked"]["mafft"])

    def test_chaves_anteriores_continuam(self):
        d = self._manifesto({}).to_dict()
        for chave in ("manifest_version", "tools_invoked", "reproducibility", "params",
                      "outputs_sha256", "inference_methods", "execution_mode"):
            self.assertIn(chave, d)
        self.assertEqual(d["nondeterminism"], [])


class TestCaracterizacaoComDadoReal(unittest.TestCase):
    """Somente leitura: manifestos de execuções já feitas."""

    def _tools(self, projeto):
        caminho = os.path.join(RAIZ, "projects", projeto, "out", "outputs", "manifest.json")
        if not os.path.exists(caminho):
            self.skipTest(f"{projeto} não existe nesta máquina")
        with open(caminho, encoding="utf-8") as h:
            return json.load(h)["tools_invoked"]

    def test_varv52(self):
        itens = etapas_nao_deterministicas(self._tools("Variola_VARV52_reexec_20261029b"))
        self.assertEqual(len(itens), 1)
        self.assertIn("--maxiterate", itens[0]["comando"])
        self.assertEqual(itens[0]["comando"][itens[0]["comando"].index("--maxiterate") + 1], "2")

    def test_zika21(self):
        itens = etapas_nao_deterministicas(self._tools("Zika_21seq_validacao_20261029b"))
        self.assertEqual(len(itens), 1)
        self.assertEqual(itens[0]["comando"][itens[0]["comando"].index("--maxiterate") + 1], "1000")


if __name__ == "__main__":
    unittest.main()
