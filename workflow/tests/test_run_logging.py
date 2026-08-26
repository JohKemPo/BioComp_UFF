"""
Um arquivo de log por execução — item 4 de D22.

O que estes testes fixam não é como o log é formatado: é que **duas execuções
nunca compartilham arquivo**. Enquanto compartilhavam, nenhuma leitura
conseguia separá-las, e a duração reportada somava as duas mais o intervalo
ocioso entre elas — 1 960 s onde a última execução levara 396 s.

Executar com:

    python -m unittest workflow.tests.test_run_logging
"""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
import unittest

from workflow.utils import run_logging


class TestNomeDoArquivo(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self._limpar_raiz()

    def tearDown(self):
        self._limpar_raiz()
        run_logging._ATUAL = None
        shutil.rmtree(self.dir, ignore_errors=True)

    @staticmethod
    def _limpar_raiz():
        raiz = logging.getLogger()
        for h in list(raiz.handlers):
            raiz.removeHandler(h)
            h.close()

    def test_run_id_entra_no_nome(self):
        caminho = run_logging.configurar(self.dir, "4a8ad78f90d246acaf36b98d")
        self.assertIn("4a8ad78f90d2", os.path.basename(caminho))

    def test_duas_execucoes_nao_compartilham_arquivo(self):
        """É o defeito inteiro: com o nome por dia e o modo `append`, duas
        execuções do mesmo dia escreviam no mesmo arquivo."""
        a = run_logging.configurar(self.dir, "aaaaaaaaaaaa1111")
        b = run_logging.configurar(self.dir, "bbbbbbbbbbbb2222")
        self.assertNotEqual(a, b)

    def test_sem_run_id_mantem_o_nome_antigo(self):
        """Módulo usado fora do workflow continua tendo para onde escrever."""
        caminho = run_logging.configurar(self.dir)
        self.assertTrue(os.path.basename(caminho).startswith("log_setup_"))


class TestDestinoEfetivo(unittest.TestCase):
    """Não basta escolher o nome: o `logging` tem de acabar naquele arquivo."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        TestNomeDoArquivo._limpar_raiz()

    def tearDown(self):
        TestNomeDoArquivo._limpar_raiz()
        run_logging._ATUAL = None
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_mensagem_vai_para_o_arquivo_da_execucao(self):
        caminho = run_logging.configurar(self.dir, "cccccccccccc3333")
        logging.info("STEP: Aligning seqs...")
        logging.shutdown()
        with open(caminho, encoding="utf-8") as handle:
            self.assertIn("STEP: Aligning seqs...", handle.read())

    def test_configurar_de_novo_redireciona(self):
        """A segunda execução não pode continuar escrevendo no arquivo da
        primeira. `logging.basicConfig` é no-op quando a raiz já tem handler, e
        era por isso que a configuração dependia da ordem de importação."""
        primeiro = run_logging.configurar(self.dir, "dddddddddddd4444")
        logging.info("da primeira execução")
        segundo = run_logging.configurar(self.dir, "eeeeeeeeeeee5555")
        logging.info("da segunda execução")
        logging.shutdown()

        with open(primeiro, encoding="utf-8") as handle:
            texto_a = handle.read()
        with open(segundo, encoding="utf-8") as handle:
            texto_b = handle.read()

        self.assertIn("da primeira execução", texto_a)
        self.assertNotIn("da segunda execução", texto_a)
        self.assertIn("da segunda execução", texto_b)

    def test_garantir_nao_rouba_o_log_da_execucao(self):
        """Sete módulos chamavam `basicConfig` no construtor. Se qualquer um
        deles pudesse trocar o destino no meio da execução, o log da execução
        se partiria em dois arquivos."""
        da_execucao = run_logging.configurar(self.dir, "ffffffffffff6666")
        outro = os.path.join(self.dir, "outro")
        os.makedirs(outro)
        run_logging.garantir(outro)
        logging.info("continua no log da execução")
        logging.shutdown()

        with open(da_execucao, encoding="utf-8") as handle:
            self.assertIn("continua no log da execução", handle.read())
        self.assertEqual(os.listdir(outro), [])

    def test_garantir_configura_quando_ninguem_configurou(self):
        caminho = run_logging.garantir(self.dir)
        self.assertTrue(os.path.exists(caminho))


if __name__ == "__main__":
    unittest.main()
