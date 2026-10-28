"""
Curadoria de datas no relógio molecular (plano 16-curadoria-varv §6, I1-I3).

Invariantes testados:

1. **O GenBank vence.** `/collection_date` presente nunca é sobrescrito; a
   curadoria só preenche o que está vazio.
2. **Casamento pelo acesso sem versão** (D13), inclusive com o rótulo truncado
   de IQ-TREE/RAxML (``NC_008030.``).
3. **Ambíguo não vira número.** Linha com conflito, sem data ou incoerente é
   registrada como ignorada, com o motivo.
4. **Opcional e declarada.** Desligada, nada muda; ligada, o `protocolo` traz
   caminho, SHA-256 e a origem de cada data.

Executar com:

    python -m unittest workflow.tests.test_relogio_curadoria
"""

from __future__ import annotations

import hashlib
import os
import tempfile
import unittest
from unittest import mock

from workflow.molecular_clock import datacao
from workflow.tests import test_molecular_clock as base
from workflow.utils import external_tools

CABECALHO = ("acesso,pais,pais_iso3,data_coleta,precisao_data,origem,fonte,"
             "conflito_com_genbank,observacao,acesso_sem_versao,conjuntos")


def _linha(acesso, data="", prec="", conflito="nao", origem="curadoria_esposito2006"):
    return f"{acesso},,,{data},{prec},{origem},,{conflito},,{acesso.split('.')[0]},VARV-52"


def _csv(*linhas, cabecalho=CABECALHO) -> str:
    fd, path = tempfile.mkstemp(suffix=".csv", prefix="curadoria_")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write("\n".join([cabecalho, *linhas]) + "\n")
    return path


class _ComCsv(unittest.TestCase):
    def csv(self, *linhas, **kw) -> str:
        p = _csv(*linhas, **kw)
        self.addCleanup(os.unlink, p)
        return p


class TestLeitura(_ComCsv):
    def test_acesso_sem_versao_e_com_versao_casam_pela_mesma_chave(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("AB000003.1", "1975", "ano"), _linha("AB000004", "1980", "ano")))
        self.assertEqual(set(r["datas"]), {"AB000003", "AB000004"})

    def test_linha_com_conflito_nao_entra_e_diz_por_que(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("AB000001.1", "1975", "ano", conflito="ambiguo_persistente_na_fonte")))
        self.assertEqual(r["datas"], {})
        self.assertEqual(r["ignoradas"][0]["acesso"], "AB000001.1")
        self.assertIn("ambiguo_persistente_na_fonte", r["ignoradas"][0]["motivo"])

    def test_linhas_vazias_sao_registradas_nao_viram_data(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("AB000001.1", ""), _linha("", "1975", "ano")))
        self.assertEqual(r["datas"], {})
        self.assertEqual([i["motivo"] for i in r["ignoradas"]],
                         ["data_coleta vazia", "linha sem acesso"])

    def test_precisao_incoerente_e_data_desconhecida_sao_ignoradas(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("AB000001.1", "1975", "dia"), _linha("AB000002.1", "c. 1975", "ano")))
        self.assertEqual(r["datas"], {})
        self.assertEqual(len(r["ignoradas"]), 2)

    def test_acesso_duplicado_so_vale_se_as_datas_concordam(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("AB000001.1", "1975", "ano"), _linha("AB000001.1", "1975", "ano"),
            _linha("AB000002.1", "1975", "ano"), _linha("AB000002", "1976", "ano")))
        self.assertEqual(set(r["datas"]), {"AB000001"})
        self.assertEqual(len(r["ignoradas"]), 2)

    def test_coluna_nova_e_tolerada_e_sha256_e_do_arquivo(self):
        cab = CABECALHO + ",pais_iso3_hoje_aproximado"
        p = self.csv(_linha("AB000001.1", "1975", "ano") + ",XXX", cabecalho=cab)
        r = datacao.ler_curadoria_datas(p)
        self.assertIn("AB000001", r["datas"])
        with open(p, "rb") as fh:
            self.assertEqual(r["sha256"], hashlib.sha256(fh.read()).hexdigest())

    def test_csv_ausente(self):
        with self.assertRaises(FileNotFoundError):
            datacao.ler_curadoria_datas("/nao/existe/curadoria.csv")

    def test_csv_malformado(self):
        with self.assertRaisesRegex(ValueError, "sem a"):
            datacao.ler_curadoria_datas(self.csv("a,b", cabecalho="acesso,pais"))
        with self.assertRaisesRegex(ValueError, "malformado"):
            datacao.ler_curadoria_datas(self.csv(_linha("AB000001.1", "1975", "ano") + ",extra"))
        vazio = self.csv(cabecalho="")
        with self.assertRaises(ValueError):
            datacao.ler_curadoria_datas(vazio)


class TestAplicacao(unittest.TestCase):
    def test_genbank_nunca_e_sobrescrito(self):
        crus = {"AB000001.1": "2001-05-17", "AB000002.1": None}
        cur = {"AB000001": {"data": "1900", "precisao": "ano", "origem": "x"},
               "AB000002": {"data": "1975", "precisao": "ano", "origem": "x"}}
        novo, pre = datacao.aplicar_curadoria(crus, cur)
        self.assertEqual(novo["AB000001.1"], "2001-05-17")
        self.assertEqual(novo["AB000002.1"], "1975")
        self.assertEqual(set(pre), {"AB000002.1"})
        self.assertIsNone(crus["AB000002.1"], "a entrada não pode ser alterada")

    def test_rotulo_truncado_do_d13(self):
        novo, pre = datacao.aplicar_curadoria(
            {"NC_008030.": None}, {"NC_008030": {"data": "1970", "precisao": "ano", "origem": "x"}})
        self.assertEqual(novo["NC_008030."], "1970")

    def test_acesso_sem_registro_no_gb_tambem_e_vazio(self):
        # .gb truncado (VARV-52: 3 registros em 52 pontas): sem registro, sem campo a preservar.
        cur = {"AB000009": {"data": "1970", "precisao": "ano", "origem": "x"},
               "AB000001": {"data": "1900", "precisao": "ano", "origem": "x"}}
        novo, pre = datacao.aplicar_curadoria({"AB000001.1": "2001"}, cur,
                                              ["AB000001.1", "AB000009.1", "AB000008.1"])
        self.assertEqual(novo, {"AB000001.1": "2001", "AB000009.1": "1970"})
        self.assertTrue(pre["AB000009.1"]["sem_registro_genbank"])

    def test_sem_curadoria_nada_muda(self):
        crus = {"AB000001.1": None, "AB000002.1": "1999"}
        novo, pre = datacao.aplicar_curadoria(crus, {})
        self.assertEqual((novo, pre), (crus, {}))


class TestDatasEscritas(unittest.TestCase):
    def test_precisao_ano_mes_dia_e_origem_por_ponta(self):
        rot = {"AB000001": "AB000001.1", "AB000002": "AB000002.1",
               "AB000003": "AB000003.1", "AB000004": "AB000004.1"}
        crus = {"AB000001.1": "2001", "AB000002.1": "2002-05", "AB000003.1": "2003-06-07",
                "AB000004.1": "2010-01-01"}
        curadas = {k: {"origem": "o"} for k in ("AB000001.1", "AB000002.1", "AB000003.1")}
        with tempfile.TemporaryDirectory() as d:
            dst = os.path.join(d, "datas.txt")
            r = datacao.escrever_arquivo_datas(crus, rot, dst, frozenset(rot), False, curadas)
            linhas = dict(l.split("\t") for l in open(dst).read().split("\n") if l)
        self.assertEqual(linhas["AB000001.1"], "2001:2002")
        self.assertEqual(linhas["AB000002.1"], "2002-05")
        self.assertEqual(linhas["AB000003.1"], "2003-06-07")
        self.assertEqual(r["precisao"], {"dia": 2, "mes": 1, "ano": 1})
        self.assertEqual(r["origem_das_datas"], {"genbank_estruturado": 1, "curadoria": 3})
        self.assertEqual(sorted(r["pontas_datadas_por_curadoria"]),
                         ["AB000001.1", "AB000002.1", "AB000003.1"])

    def test_sem_curadoria_tudo_e_genbank(self):
        with tempfile.TemporaryDirectory() as d:
            r = datacao.escrever_arquivo_datas(
                {"AB000001.1": "2001"}, {"AB000001": "AB000001.1"}, os.path.join(d, "x"),
                frozenset({"AB000001"}))
        self.assertEqual(r["origem_das_datas"], {"genbank_estruturado": 1, "curadoria": 0})
        self.assertEqual(r["pontas_datadas_por_curadoria"], {})


class TestControlador(base._Projeto, _ComCsv):
    """Nenhum destes caminhos chega ao IQ-TREE (portão de datas insuficientes)."""

    def setUp(self):
        super().setUp()
        for alvo, kw in (("workflow.molecular_clock.datacao.subprocess.run",
                          {"side_effect": AssertionError("o IQ-TREE não deveria ser chamado")}),
                         ("workflow.controller.molecularClockController._versao_iqtree",
                          {"return_value": "3.1.3"})):
            p = mock.patch(alvo, **kw)
            p.start()
            self.addCleanup(p.stop)
        p = mock.patch.object(external_tools, "resolve_tool", return_value="/bin/iqtree3")
        p.start()
        self.addCleanup(p.stop)

    def _arv(self, r):
        return r["por_arvore"]["tree_dataset_final_mafft_iqtree.nexus"]

    def test_desligada_por_padrao_e_declarada(self):
        self.escrever(datas={r: None for r in base.TODAS})
        r = self.controlador()()
        self.assertFalse(r["protocolo"]["curadoria"]["ligada"])
        self.assertEqual(r["protocolo"]["curadoria"]["estado"], "desligada")
        self.assertEqual(self._arv(r)["datas"]["pontas_datadas"], 0)

    def test_ligada_registra_origem_caminho_e_hash(self):
        self.escrever(datas={r: None for r in base.TODAS})
        p = self.csv(_linha("AB000001.1", "1975", "ano"), _linha("AB000002", "1980-03", "mes"),
                     _linha("AB000003.1", "1990", "ano", conflito="ambiguo_marcado_para_decisao"))
        r = self.controlador(molecular_clock_curadoria=p)()
        cur = r["protocolo"]["curadoria"]
        self.assertEqual((cur["ligada"], cur["estado"], cur["csv"]), (True, "aplicada", p))
        with open(p, "rb") as fh:
            self.assertEqual(cur["sha256"], hashlib.sha256(fh.read()).hexdigest())
        self.assertEqual(cur["registros_do_gb_sem_collection_date"], len(base.TODAS))
        self.assertEqual(cur["linhas_ignoradas"][0]["acesso"], "AB000003.1")
        d = self._arv(r)["datas"]
        self.assertEqual(d["origem_das_datas"], {"genbank_estruturado": 0, "curadoria": 2})
        self.assertEqual(d["pontas_datadas"], 2)
        self.assertEqual(sorted(d["pontas_datadas_por_curadoria"]), ["AB000001.1", "AB000002.1"])
        self.assertIn("nunca sobrescreve", r["protocolo"]["fonte_das_datas"])

    def test_genbank_prevalece_no_arquivo_de_datas(self):
        datas = {r: None for r in base.TODAS}
        datas["AB000001.1"] = "2001-05-17"
        self.escrever(datas=datas)
        p = self.csv(_linha("AB000001.1", "1900", "ano"), _linha("AB000002.1", "1975", "ano"))
        r = self.controlador(molecular_clock_curadoria=p)()
        stem = os.path.join(self.out, "outputs", "molecular_clock",
                            "tree_dataset_final_mafft_iqtree", "tree_dataset_final_mafft_iqtree")
        escritas = dict(l.split("\t") for l in open(stem + ".datas.txt").read().split("\n") if l)
        self.assertEqual(escritas["AB000001.1"], "2001-05-17")
        self.assertEqual(escritas["AB000002.1"], "1975:1976")
        self.assertEqual(self._arv(r)["datas"]["origem_das_datas"],
                         {"genbank_estruturado": 1, "curadoria": 1})

    def test_csv_ausente_roda_como_sem_curadoria_e_declara(self):
        self.escrever(datas={r: None for r in base.TODAS})
        r = self.controlador(molecular_clock_curadoria="/nao/existe.csv")()
        cur = r["protocolo"]["curadoria"]
        self.assertEqual(cur["estado"], "indisponivel")
        self.assertTrue(cur["ligada"])
        self.assertIsNone(cur["sha256"])
        self.assertEqual(self._arv(r)["datas"]["pontas_datadas"], 0)

    def test_csv_malformado_roda_como_sem_curadoria_e_declara(self):
        self.escrever(datas={r: None for r in base.TODAS})
        r = self.controlador(molecular_clock_curadoria=self.csv("x", cabecalho="acesso"))()
        self.assertEqual(r["protocolo"]["curadoria"]["estado"], "indisponivel")
        self.assertIn("sem a", r["protocolo"]["curadoria"]["motivo"])


if __name__ == "__main__":
    unittest.main()
