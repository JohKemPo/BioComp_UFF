"""
Testes do download de GenBank — F2.5.

O defeito: `download_sequences` zerava `raw_data_sequences.gb` a cada chamada.
O builder chama a função uma vez por árvore com acessos novos, então o arquivo
ficava só com os registros da última árvore (VARV-52: 3 de 51 terminais).
Estes testes fixam que as chamadas acumulam, que uma reexecução não duplica,
e que a ordem do arquivo segue a lista de acessos pedida.
"""
import io
import os
import tempfile
import types
import unittest
from unittest import mock

from workflow.utils import treeUtils


def _registro(acesso):
    return (
        f"LOCUS       {acesso:<16}10 bp    DNA     linear   VRL 01-JAN-2000\n"
        f"DEFINITION  teste {acesso}.\n"
        f"ACCESSION   {acesso}\n"
        f"VERSION     {acesso}.1\n"
        "FEATURES             Location/Qualifiers\n"
        "ORIGIN\n"
        "        1 acgtacgtac\n"
        "//\n"
    )


class _Servidor:
    """Imita o efetch: devolve os registros pedidos, na ordem indicada por `inverter`."""

    def __init__(self, inverter=False):
        self.inverter = inverter
        self.chamadas = []

    def __call__(self, db, id, rettype, retmode):
        ids = list(id)
        self.chamadas.append(ids)
        ordem = list(reversed(ids)) if self.inverter else ids
        handle = mock.Mock()
        handle.read.return_value = "".join(_registro(a) for a in ordem)
        return handle


def _consultas(*acessos):
    return [types.SimpleNamespace(name=a) for a in acessos]


def _log():
    return mock.Mock()


class DownloadSequencesAcumulaTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.saida = os.path.join(self.tmp.name, "raw_data_sequences.gb")
        self.estado = mock.patch.object(treeUtils, "_arquivos_iniciados", set(), create=True)
        self.estado.start()
        self.sono = mock.patch.object(treeUtils.time, "sleep")
        self.sono.start()

    def tearDown(self):
        self.sono.stop()
        self.estado.stop()
        self.tmp.cleanup()

    def _ler(self):
        with open(self.saida) as f:
            return f.read()

    def _contar(self):
        return self._ler().count("LOCUS       ")

    def _baixar(self, servidor, acessos, batch_size=2):
        with mock.patch.object(treeUtils.Entrez, "efetch", side_effect=servidor):
            return treeUtils.download_sequences(None, _consultas(*acessos), self.saida, _log(), batch_size=batch_size)

    def test_dois_lotes_numa_chamada_somam_registros(self):
        servidor = _Servidor()
        self._baixar(servidor, ["A1", "B2", "C3"], batch_size=2)
        self.assertEqual(len(servidor.chamadas), 2)
        self.assertEqual(self._contar(), 3)

    def test_um_lote_so_nao_muda_nada(self):
        servidor = _Servidor()
        self._baixar(servidor, ["A1", "B2"], batch_size=10)
        self.assertEqual(len(servidor.chamadas), 1)
        self.assertEqual(self._ler(), _registro("A1") + _registro("B2"))

    def test_chamadas_sucessivas_acumulam_como_no_builder(self):
        servidor = _Servidor()
        self._baixar(servidor, ["A1", "B2"])
        self._baixar(servidor, ["C3"])
        self._baixar(servidor, ["D4", "E5", "F6"])
        self.assertEqual(self._contar(), 6)
        self.assertIn("LOCUS       A1", self._ler())

    def test_reexecucao_sobre_arquivo_existente_nao_duplica(self):
        servidor = _Servidor()
        self._baixar(servidor, ["A1", "B2", "C3"])
        with open(self.saida, "w") as f:
            f.write(_registro("A1") + _registro("B2") + _registro("C3"))
        treeUtils._arquivos_iniciados.clear()
        self._baixar(servidor, ["A1", "B2", "C3"])
        self.assertEqual(self._contar(), 3)

    def test_ordem_segue_a_lista_mesmo_com_resposta_invertida(self):
        servidor = _Servidor(inverter=True)
        self._baixar(servidor, ["C3", "A1", "B2"], batch_size=10)
        self.assertEqual(self._ler(), _registro("C3") + _registro("A1") + _registro("B2"))

    def test_lista_vazia_nao_chama_rede_e_zera_o_arquivo(self):
        servidor = _Servidor()
        resultado = self._baixar(servidor, [])
        self.assertEqual(servidor.chamadas, [])
        self.assertIsNone(resultado)
        self.assertEqual(self._ler(), "")


if __name__ == "__main__":
    unittest.main()
