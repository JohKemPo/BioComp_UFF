"""
P2 / DM-10 — testes de `workflow.tree_construction.qualidade_alinhamento`.

Executar com:

    python -m unittest workflow.tests.test_qualidade_alinhamento
"""

import unittest

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from workflow.controller.treeBuilderController import TreeBuilderController
from workflow.tree_construction.builder import TreeBuilder
from workflow.tree_construction.distancia_p import FRACAO_MINIMA, PISO_SITIOS
from workflow.tree_construction.qualidade_alinhamento import (
    LACUNA_MAX_COLUNA, LACUNA_MAX_TAXON, LIMITE_APLICACAO, avaliar_alinhamento, medir_perfil)
from workflow.utils import tool_runs
from workflow.utils.manifest import ExecutionManifest


def aln(*seqs):
    return MultipleSeqAlignment(
        [SeqRecord(Seq(s), id=f"t{i}", description="") for i, s in enumerate(seqs)])


class TestPadroes(unittest.TestCase):

    def test_padroes_do_protocolo(self):
        self.assertEqual((LACUNA_MAX_COLUNA, LACUNA_MAX_TAXON), (0.50, 0.90))


class TestColuna(unittest.TestCase):

    def test_coluna_exatamente_em_50_nao_sai(self):
        # coluna 1: 2 de 4 com lacuna = 0,5, regra é "mais de"
        saida, rel = avaliar_alinhamento(aln("A-A", "A-A", "ACA", "ACA"), ligado=True, teto_aplicacao=1.0)
        self.assertEqual(rel["colunas_antes"], 3)
        self.assertEqual(rel["colunas_depois"], 3)
        self.assertEqual(saida.get_alignment_length(), 3)

    def test_coluna_um_pouco_acima_de_50_sai(self):
        saida, rel = avaliar_alinhamento(aln("A-A", "A-A", "A-A", "ACA", "ACA"), ligado=True, teto_aplicacao=1.0)
        self.assertEqual(rel["colunas_depois"], 2)
        self.assertEqual([str(r.seq) for r in saida], ["AA", "AA", "AA", "AA", "AA"])
        self.assertAlmostEqual(rel["fracao_removida"], 1 / 3)

    def test_sem_lacuna_nada_sai(self):
        saida, rel = avaliar_alinhamento(aln("ACGT", "ACGA"), ligado=True, teto_aplicacao=1.0)
        self.assertEqual(rel["colunas_depois"], 4)
        self.assertEqual(rel["fracao_removida"], 0.0)
        self.assertEqual(rel["taxons_acima_limiar"], [])

    def test_tudo_lacuna_remove_todas_as_colunas_e_mantem_os_taxons(self):
        saida, rel = avaliar_alinhamento(aln("----", "----"), ligado=True, teto_aplicacao=1.0)
        self.assertEqual(rel["colunas_depois"], 0)
        self.assertEqual(rel["fracao_removida"], 1.0)
        self.assertEqual(len(saida), 2)

    def test_n_nao_e_lacuna(self):
        _, rel = avaliar_alinhamento(aln("NNNN", "NNNN"), ligado=True, teto_aplicacao=1.0)
        self.assertEqual(rel["colunas_depois"], 4)

    def test_sem_colunas_fracao_e_none(self):
        _, rel = avaliar_alinhamento(aln("", ""), ligado=True, teto_aplicacao=1.0)
        self.assertIsNone(rel["fracao_removida"])

    def test_limiar_fora_de_0_1_e_recusado(self):
        with self.assertRaises(ValueError):
            medir_perfil(aln("AC"), lacuna_max_coluna=1.5)


class TestTaxon(unittest.TestCase):

    def _dez_colunas(self, lacunas_no_taxon):
        # 3 táxons, 10 colunas todas conservadas (<= 50% de lacuna por coluna)
        return aln("-" * lacunas_no_taxon + "A" * (10 - lacunas_no_taxon),
                   "A" * 10, "A" * 10)

    def test_taxon_exatamente_em_90_nao_e_reportado(self):
        _, rel = avaliar_alinhamento(self._dez_colunas(9), ligado=True, teto_aplicacao=1.0)
        self.assertEqual(rel["taxons_acima_limiar"], [])

    def test_taxon_acima_de_90_e_reportado_e_mantido(self):
        saida, rel = avaliar_alinhamento(self._dez_colunas(10), ligado=True, teto_aplicacao=1.0)
        # 10 de 10 lacunas: o táxon cruza 90%, as colunas ficam (1 de 3 <= 50%)
        self.assertEqual([t["id"] for t in rel["taxons_acima_limiar"]], ["t0"])
        self.assertEqual([r.id for r in saida], ["t0", "t1", "t2"])
        self.assertEqual(rel["colunas_depois"], 10)

    def test_filtro_nunca_remove_taxon(self):
        original = aln("----------", "AAAAAAAAAA", "----------", "----------")
        saida, rel = avaliar_alinhamento(original, ligado=True, teto_aplicacao=1.0)
        self.assertEqual([r.id for r in saida], [r.id for r in original])
        self.assertEqual(rel["colunas_depois"], 0)


class TestLigadoDesligado(unittest.TestCase):

    def test_desligado_devolve_o_mesmo_alinhamento_e_mede(self):
        original = aln("A-A", "A-A", "A-A", "ACA")
        saida, rel = avaliar_alinhamento(original, ligado=False)
        self.assertIs(saida, original)
        self.assertFalse(rel["aplicado"])
        self.assertEqual(rel["colunas_removiveis"], 1)
        self.assertEqual(rel["colunas_depois"], rel["colunas_antes"])
        self.assertEqual(rel["fracao_removida"], 0.0)

    def test_ligado_marca_aplicado(self):
        _, rel = avaliar_alinhamento(aln("A-A", "A-A", "A-A", "ACA"), ligado=True, teto_aplicacao=1.0)
        self.assertTrue(rel["aplicado"])


class TestLigacao(unittest.TestCase):

    def setUp(self):
        import os, tempfile
        tool_runs.limpar()
        self.addCleanup(tool_runs.limpar)
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.caminho = os.path.join(tmp.name, "x.aln")

    def _controlador(self, ligado, teto=1.0):
        c = object.__new__(TreeBuilderController)
        c.filtro_lacuna, c.lacuna_max_coluna, c.lacuna_max_taxon = ligado, 0.5, 0.9
        c.lacuna_teto_aplicacao = teto
        c._qualidade_alinhamento = {}
        return c

    def test_controlador_desligado_mede_e_registra_sem_filtrar(self):
        original = aln("A-A", "A-A", "A-A", "ACA")
        saida = self._controlador(False)._qualidade_do_alinhamento(original, self.caminho)
        self.assertIs(saida, original)
        registro, = tool_runs.qualidades_alinhamento()
        self.assertEqual(registro["colunas_removiveis"], 1)
        self.assertFalse(registro["aplicado"])

    def test_controlador_ligado_entrega_o_alinhamento_filtrado(self):
        saida = self._controlador(True)._qualidade_do_alinhamento(
            aln("A-A", "A-A", "A-A", "ACA"), self.caminho)
        self.assertEqual(saida.get_alignment_length(), 2)

    def test_avaliacao_e_uma_por_alinhamento(self):
        c = self._controlador(True)
        a = c._qualidade_do_alinhamento(aln("A-A", "A-A", "A-A", "ACA"), self.caminho)
        b = c._qualidade_do_alinhamento(aln("A-A", "A-A", "A-A", "ACA"), self.caminho)
        self.assertIs(a, b)
        self.assertEqual(len(tool_runs.qualidades_alinhamento()), 1)


class TestTeto(unittest.TestCase):

    def test_acima_do_teto_nao_aplica_e_diz_o_motivo(self):
        original = aln("A-A", "A-A", "A-A", "ACA")  # 1 de 3 colunas = 33%
        saida, rel = avaliar_alinhamento(original, ligado=True)
        self.assertIs(saida, original)
        self.assertFalse(rel["aplicado"])
        self.assertEqual(rel["motivo"], "acima_do_teto")
        self.assertEqual(rel["colunas_removiveis"], 1)
        self.assertEqual(rel["colunas_depois"], 3)

    def test_exatamente_no_teto_aplica(self):
        # 1 de 4 colunas = 25%: a regra é "mais de"
        saida, rel = avaliar_alinhamento(aln("A-AA", "A-AA", "A-AA", "ACAA"), ligado=True)
        self.assertTrue(rel["aplicado"])
        self.assertIsNone(rel["motivo"])
        self.assertEqual(saida.get_alignment_length(), 3)

    def test_desligado_tem_motivo_distinto_do_teto(self):
        _, rel = avaliar_alinhamento(aln("A-A", "A-A", "A-A", "ACA"), ligado=False)
        self.assertEqual(rel["motivo"], "desligado_por_configuracao")

    def test_padrao_do_teto_e_25(self):
        self.assertEqual(LIMITE_APLICACAO, 0.25)

    def test_teto_fora_de_0_1_e_recusado(self):
        with self.assertRaises(ValueError):
            avaliar_alinhamento(aln("AC"), ligado=True, teto_aplicacao=2)

    def test_controlador_grava_o_filtrado_e_os_dois_sha256(self):
        import hashlib, os, tempfile
        from Bio import AlignIO
        tool_runs.limpar()
        self.addCleanup(tool_runs.limpar)
        with tempfile.TemporaryDirectory() as d:
            origem = os.path.join(d, "x.aln")
            original = aln("A-A", "A-A", "A-A", "ACA")
            AlignIO.write(original, origem, "fasta")
            c = object.__new__(TreeBuilderController)
            c.filtro_lacuna, c.lacuna_max_coluna, c.lacuna_max_taxon = True, 0.5, 0.9
            c.lacuna_teto_aplicacao, c._qualidade_alinhamento = 1.0, {}
            c._qualidade_do_alinhamento(original, origem)
            reg, = tool_runs.qualidades_alinhamento()
            filtrado = os.path.join(d, "x_lacuna.aln")
            self.assertEqual(reg["alinhamento_filtrado"], filtrado)
            self.assertEqual(AlignIO.read(filtrado, "fasta").get_alignment_length(), 2)
            for chave, caminho in (("alinhamento_sha256", origem),
                                   ("alinhamento_filtrado_sha256", filtrado)):
                with open(caminho, "rb") as f:
                    self.assertEqual(reg[chave], hashlib.sha256(f.read()).hexdigest())

    def test_caminhos_do_manifesto_existem_ao_fim_da_execucao(self):
        import hashlib, os, tempfile
        from Bio import AlignIO
        from workflow.utils.dataCleaning import copiar_arquivos
        tool_runs.limpar()
        self.addCleanup(tool_runs.limpar)
        with tempfile.TemporaryDirectory() as raiz:
            for d in ("tmp", "Align"):
                os.makedirs(os.path.join(raiz, d))
            origem = os.path.join(raiz, "tmp", "x.aln")
            original = aln("A-A", "A-A", "A-A", "ACA")
            AlignIO.write(original, origem, "fasta")
            c = object.__new__(TreeBuilderController)
            c.output_path = raiz
            c.filtro_lacuna, c.lacuna_max_coluna, c.lacuna_max_taxon = True, 0.5, 0.9
            c.lacuna_teto_aplicacao, c._qualidade_alinhamento = 1.0, {}
            c._qualidade_do_alinhamento(original, origem)
            # o mesmo passo que fecha `TreeBuilderController.__call__`
            copiar_arquivos(os.path.join(raiz, "tmp"), os.path.join(raiz, "Align"))
            reg, = tool_runs.qualidades_alinhamento()
            for chave, sha in (("alinhamento", "alinhamento_sha256"),
                               ("alinhamento_filtrado", "alinhamento_filtrado_sha256")):
                self.assertTrue(os.path.isfile(reg[chave]), reg[chave])
                self.assertEqual(os.path.dirname(reg[chave]), os.path.join(raiz, "Align"))
                with open(reg[chave], "rb") as f:
                    self.assertEqual(hashlib.sha256(f.read()).hexdigest(), reg[sha])

    def test_acima_do_teto_nao_grava_arquivo(self):
        import os, tempfile
        tool_runs.limpar()
        self.addCleanup(tool_runs.limpar)
        with tempfile.TemporaryDirectory() as d:
            origem = os.path.join(d, "x.aln")
            c = object.__new__(TreeBuilderController)
            c.filtro_lacuna, c.lacuna_max_coluna, c.lacuna_max_taxon = True, 0.5, 0.9
            c.lacuna_teto_aplicacao, c._qualidade_alinhamento = 0.25, {}
            c._qualidade_do_alinhamento(aln("A-A", "A-A", "A-A", "ACA"), origem)
            self.assertEqual(os.listdir(d), [])
            self.assertIsNone(tool_runs.qualidades_alinhamento()[0]["alinhamento_filtrado"])


class TestManifesto(unittest.TestCase):

    def setUp(self):
        tool_runs.limpar()
        self.addCleanup(tool_runs.limpar)

    def test_builder_registra_os_parametros_do_u4(self):
        TreeBuilder().distance_matrix(aln("A" * 120, "A" * 119 + "C"))
        self.assertEqual(tool_runs.parametros_efetivos()["p_distancia"],
                         {"fracao_minima": FRACAO_MINIMA, "piso_sitios": PISO_SITIOS})

    def test_builder_registra_parametro_declarado(self):
        TreeBuilder(fracao_minima=0.25, piso_sitios=7).distance_matrix(aln("ACGTACGT", "ACGTACGA"))
        self.assertEqual(tool_runs.parametros_efetivos()["p_distancia"],
                         {"fracao_minima": 0.25, "piso_sitios": 7})

    def test_manifesto_leva_qualidade_e_parametros(self):
        _, rel = avaliar_alinhamento(aln("A-A", "A-A", "A-A", "ACA"), ligado=True, teto_aplicacao=1.0)
        tool_runs.registrar_qualidade_alinhamento("x.aln", **rel)
        tool_runs.registrar_parametros("p_distancia", fracao_minima=0.5, piso_sitios=100)
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            m = ExecutionManifest(d, {})
            m.drain_tool_runs()
            dados = m.to_dict()
        self.assertEqual(dados["alignment_quality"][0]["colunas_antes"], 3)
        self.assertEqual(dados["effective_params"]["p_distancia"]["piso_sitios"], 100)
        self.assertEqual(dados["manifest_version"], 2)


if __name__ == "__main__":
    unittest.main()
