"""
Testes do AC-09 (resíduo `_NoPipe` reprocessado como conjunto de dados) e da
soma do heatmap acumulado.

Os dois defeitos apareceram juntos na execução do ZIKV-479 de 2026-11-18: o
resíduo de uma execução morta entrou na lista de entradas, e a soma do heatmap
contra um `8x8` fixo devolveu None no primeiro arquivo e derrubou o processo no
segundo, antes da mineração e do relógio.
"""

import os
import shutil
import tempfile
import unittest

from workflow.controller.treeBuilderController import TreeBuilderController


FASTA = ">A\nACGT\n>B\nACGA\n"


class TestEntradaIgnoraNoPipe(unittest.TestCase):
    """O resíduo da deduplicação não é conjunto de dados."""

    def setUp(self):
        self.base = tempfile.mkdtemp()
        self.entrada = os.path.join(self.base, 'entrada')
        self.saida = os.path.join(self.base, 'saida')
        os.makedirs(self.entrada)
        os.makedirs(self.saida)

    def tearDown(self):
        shutil.rmtree(self.base, ignore_errors=True)

    def _controlador(self):
        return TreeBuilderController(
            input_path=self.entrada, output_path=self.saida,
            mode='advanced', num_threads=1)

    def _escrever(self, nome):
        with open(os.path.join(self.entrada, nome), 'w') as fh:
            fh.write(FASTA)

    def test_residuo_nopipe_nao_entra_na_lista(self):
        self._escrever('dataset_final.fasta')
        self._escrever('dataset_final_NoPipe')

        c = self._controlador()

        self.assertEqual(c.files, ['dataset_final.fasta'])

    def test_residuo_duplo_tambem_fica_fora(self):
        self._escrever('dataset_final.fasta')
        self._escrever('dataset_final_NoPipe')
        self._escrever('dataset_final_NoPipe_NoPipe')

        c = self._controlador()

        self.assertEqual(c.files, ['dataset_final.fasta'])

    def test_limpeza_acontece_antes_da_enumeracao(self):
        self._escrever('dataset_final.fasta')
        self._escrever('dataset_final_NoPipe')

        self._controlador()

        # A limpeza roda na largada: o diretório de entrada fica só com a entrada.
        self.assertEqual(sorted(os.listdir(self.entrada)), ['dataset_final.fasta'])

    def test_entrada_legitima_sobrevive(self):
        self._escrever('dataset_final.fasta')
        self._escrever('outro_conjunto.fasta')

        c = self._controlador()

        self.assertEqual(c.files, ['dataset_final.fasta', 'outro_conjunto.fasta'])


class TestSomarMatrizes(unittest.TestCase):
    """A soma do acumulado não devolve None nem exige uma ordem fixa."""

    def setUp(self):
        self.base = tempfile.mkdtemp()
        entrada = os.path.join(self.base, 'entrada')
        saida = os.path.join(self.base, 'saida')
        os.makedirs(entrada)
        os.makedirs(saida)
        with open(os.path.join(entrada, 'd.fasta'), 'w') as fh:
            fh.write(FASTA)
        self.c = TreeBuilderController(
            input_path=entrada, output_path=saida, mode='advanced', num_threads=1)

    def tearDown(self):
        shutil.rmtree(self.base, ignore_errors=True)

    def test_acumulado_vazio_assume_a_ordem_do_primeiro(self):
        seis = [[1] * 6 for _ in range(6)]

        resultado = self.c.somarMatrizes(None, seis)

        self.assertEqual(resultado, seis)

    def test_soma_elemento_por_elemento(self):
        a = [[1, 2], [3, 4]]
        b = [[10, 20], [30, 40]]

        self.assertEqual(self.c.somarMatrizes(a, b), [[11, 22], [33, 44]])

    def test_ordem_incompativel_preserva_o_acumulado(self):
        seis = [[1] * 6 for _ in range(6)]
        oito = [[9] * 8 for _ in range(8)]

        resultado = self.c.somarMatrizes(seis, oito)

        self.assertEqual(resultado, seis)

    def test_segunda_matriz_ausente_preserva_o_acumulado(self):
        seis = [[1] * 6 for _ in range(6)]

        self.assertEqual(self.c.somarMatrizes(seis, None), seis)

    def test_duas_ausentes_devolvem_none(self):
        self.assertIsNone(self.c.somarMatrizes(None, None))

    def test_acumulado_sobrevive_a_sequencia_que_derrubou_o_zikv479(self):
        # Era esta a sequência: matriz de 6 (U4 derrubou NJ e UPGMA) somada ao
        # 8x8 inicial, e o arquivo seguinte somando contra o None resultante.
        seis = [[1] * 6 for _ in range(6)]

        acumulado = self.c.somarMatrizes(None, seis)
        acumulado = self.c.somarMatrizes(acumulado, seis)

        self.assertEqual(acumulado, [[2] * 6 for _ in range(6)])


if __name__ == '__main__':
    unittest.main()
