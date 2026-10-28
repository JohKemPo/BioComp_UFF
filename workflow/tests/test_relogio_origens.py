"""
Origens de datação na curadoria do relógio molecular (decisão de 2026-10-27).

Só linha com `origem` numa lista de origens confiáveis data a árvore. As
provisórias (`genbank_strain_texto`: DQ441438, DQ441439, DQ441441) e as de
origem desconhecida ou vazia ficam de fora, com o motivo na lista de ignoradas.
A lista é parâmetro, para uma fonte primária nova entrar sem mexer no código.

Executar com:

    python -m unittest workflow.tests.test_relogio_origens
"""

from __future__ import annotations

import os
import unittest

from workflow.molecular_clock import datacao
from workflow.tests.test_relogio_curadoria import _ComCsv, _linha

CSV_REAL = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "docs", "science", "dados_curadoria", "varv_curadoria_pais_data.csv")


class TestOrigemConfiavel(_ComCsv):
    def test_origens_confiaveis_datam(self):
        for origem in ("genbank_estruturado", "curadoria_esposito2006",
                       "genbank_estruturado+curadoria_esposito2006"):
            with self.subTest(origem=origem):
                r = datacao.ler_curadoria_datas(self.csv(
                    _linha("AB000001.1", "1975", "ano", origem=origem)))
                self.assertIn("AB000001", r["datas"])
                self.assertEqual(r["datas"]["AB000001"]["origem"], origem)
                self.assertEqual(r["ignoradas"], [])

    def test_provisoria_e_ignorada_com_motivo(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("DQ441438.1", "2004", "ano", origem="genbank_strain_texto")))
        self.assertEqual(r["datas"], {})
        self.assertEqual(len(r["ignoradas"]), 1)
        self.assertEqual(r["ignoradas"][0]["acesso"], "DQ441438.1")
        self.assertIn("genbank_strain_texto", r["ignoradas"][0]["motivo"])

    def test_origem_desconhecida_e_vazia_sao_ignoradas(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("AB000001.1", "1975", "ano", origem="fonte_nova_nao_aprovada"),
            _linha("AB000002.1", "1975", "ano", origem="")))
        self.assertEqual(r["datas"], {})
        self.assertEqual([i["acesso"] for i in r["ignoradas"]], ["AB000001.1", "AB000002.1"])
        self.assertTrue(all("origem não confiável" in i["motivo"] for i in r["ignoradas"]))

    def test_conflito_continua_com_o_motivo_proprio(self):
        r = datacao.ler_curadoria_datas(self.csv(
            _linha("AB000001.1", "1975", "ano", conflito="ambiguo")))
        self.assertIn("conflito_com_genbank=ambiguo", r["ignoradas"][0]["motivo"])


class TestOrigensConfiguraveis(_ComCsv):
    def test_lista_passada_substitui_o_padrao(self):
        p = self.csv(_linha("AB000001.1", "1975", "ano", origem="genbank_strain_texto"),
                     _linha("AB000002.1", "1975", "ano", origem="curadoria_esposito2006"))
        r = datacao.ler_curadoria_datas(p, origens_confiaveis=("genbank_strain_texto",))
        self.assertEqual(set(r["datas"]), {"AB000001"})
        self.assertEqual([i["acesso"] for i in r["ignoradas"]], ["AB000002.1"])

    def test_fonte_primaria_nova_entra_sem_mudar_o_padrao(self):
        p = self.csv(_linha("AB000001.1", "1975", "ano", origem="fonte_primaria_nova"))
        r = datacao.ler_curadoria_datas(p, origens_confiaveis=datacao.ORIGENS_CONFIAVEIS_CURADORIA
                                        + ("fonte_primaria_nova",))
        self.assertIn("AB000001", r["datas"])
        self.assertNotIn("fonte_primaria_nova", datacao.ORIGENS_CONFIAVEIS_CURADORIA)

    def test_string_unica_nao_vira_busca_por_caractere(self):
        p = self.csv(_linha("AB000001.1", "1975", "ano", origem="curadoria_esposito2006"))
        r = datacao.ler_curadoria_datas(p, origens_confiaveis="curadoria_esposito2006")
        self.assertIn("AB000001", r["datas"])
        r = datacao.ler_curadoria_datas(p, origens_confiaveis="curadoria")
        self.assertEqual(r["datas"], {})

    def test_padrao_e_o_conjunto_decidido(self):
        self.assertEqual(set(datacao.ORIGENS_CONFIAVEIS_CURADORIA), {
            "genbank_estruturado", "curadoria_esposito2006",
            "genbank_estruturado+curadoria_esposito2006"})


class TestCsvRealDaCuradoria(unittest.TestCase):
    def test_tres_provisorias_ficam_fora_do_csv_real(self):
        r = datacao.ler_curadoria_datas(os.path.abspath(CSV_REAL))
        for acc in ("DQ441438", "DQ441439", "DQ441441"):
            self.assertNotIn(acc, r["datas"], acc)
        motivos = {i["acesso"]: i["motivo"] for i in r["ignoradas"]}
        for acc in ("DQ441438.1", "DQ441439.1", "DQ441441.1"):
            self.assertIn("genbank_strain_texto", motivos[acc], acc)
        self.assertTrue(all(o in datacao.ORIGENS_CONFIAVEIS_CURADORIA
                            for o in {d["origem"] for d in r["datas"].values()}))


if __name__ == "__main__":
    unittest.main()
