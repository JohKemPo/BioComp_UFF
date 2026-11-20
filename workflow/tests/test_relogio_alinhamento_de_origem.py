"""Testes da resolução do alinhamento que o relógio molecular data.

O defeito, medido em 2026-11-19: o nome da árvore não carrega o sufixo do filtro
de coluna (P2), então o controlador reconstruía sempre o `.aln` **não** filtrado
e datava uma matriz diferente da que gerou a árvore — no VARV-48, 191 248 colunas
no lugar de 185 988; no ZIKV-479, 10 817 no lugar de 10 698.
"""

import hashlib
import json
import os
import shutil
import tempfile
import unittest

from workflow.controller.molecularClockController import MolecularClockController


def _sha256(caminho):
    with open(caminho, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


class _Out(unittest.TestCase):
    """Um `out/` mínimo com as duas matrizes e um manifesto configurável."""

    def setUp(self):
        self.base = tempfile.mkdtemp()
        self.out = os.path.join(self.base, 'out')
        self.align = os.path.join(self.out, 'Align')
        os.makedirs(self.align)
        os.makedirs(os.path.join(self.out, 'outputs'))
        self.simples = os.path.join(self.align, 'dataset_final_mafft.aln')
        self.filtrado = os.path.join(self.align, 'dataset_final_mafft_lacuna.aln')
        self._escrever(self.simples, 'ACGT' * 10)
        self._escrever(self.filtrado, 'ACGT' * 3)
        self.c = MolecularClockController(output_path=self.out)

    def tearDown(self):
        shutil.rmtree(self.base, ignore_errors=True)

    def _escrever(self, caminho, seq):
        with open(caminho, 'w') as fh:
            fh.write(f">A\n{seq}\n>B\n{seq}\n")

    def _manifesto(self, aplicado, **sobrepor):
        bloco = {
            'alinhamento': 'out/tmp/dataset_final_mafft.aln',
            'alinhamento_filtrado': 'out/tmp/dataset_final_mafft_lacuna.aln',
            'alinhamento_sha256': _sha256(self.simples),
            'alinhamento_filtrado_sha256': _sha256(self.filtrado),
            'aplicado': aplicado,
            'colunas_antes': 40,
            'colunas_depois': 12,
        }
        bloco.update(sobrepor)
        caminho = os.path.join(self.out, 'outputs', 'manifest.json')
        with open(caminho, 'w') as fh:
            json.dump({'alignment_quality': [bloco]}, fh)
        return caminho


class TestEscolhaPelaProveniencia(_Out):

    def test_filtro_aplicado_data_a_matriz_filtrada(self):
        self._manifesto(aplicado=True)

        aln, info = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertEqual(aln, self.filtrado)
        self.assertTrue(info['filtro_de_coluna_aplicado'])
        self.assertEqual(info['colunas'], 12)

    def test_filtro_nao_aplicado_data_a_matriz_simples(self):
        self._manifesto(aplicado=False)

        aln, info = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertEqual(aln, self.simples)
        self.assertFalse(info['filtro_de_coluna_aplicado'])
        self.assertEqual(info['colunas'], 40)

    def test_a_proveniencia_vai_ao_relatorio(self):
        # A lacuna de proveniência era a segunda metade do defeito: o relatório
        # não dizia qual matriz tinha sido datada.
        self._manifesto(aplicado=True)

        _, info = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertEqual(info['alinhamento'], 'dataset_final_mafft_lacuna.aln')
        self.assertEqual(info['sha256'], _sha256(self.filtrado))
        self.assertIn('aplicado=True', info['proveniencia'])


class TestCasosLimite(_Out):

    def test_sem_manifesto_e_com_as_duas_matrizes_recusa(self):
        # Não se adivinha pelo nome: qualquer escolha pode datar a errada.
        aln, motivo = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertIsNone(aln)
        self.assertIn('proveniência do alinhamento indisponível', motivo)

    def test_sem_manifesto_e_so_a_simples_em_disco_segue(self):
        os.remove(self.filtrado)

        aln, info = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertEqual(aln, self.simples)
        self.assertFalse(info['filtro_de_coluna_aplicado'])

    def test_manifesto_declara_aplicado_e_o_arquivo_falta(self):
        self._manifesto(aplicado=True)
        os.remove(self.filtrado)

        aln, motivo = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertIsNone(aln)
        self.assertIn('não está em Align', motivo)

    def test_matriz_alterada_depois_da_execucao_recusa(self):
        self._manifesto(aplicado=True)
        self._escrever(self.filtrado, 'TTTT' * 3)

        aln, motivo = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertIsNone(aln)
        self.assertIn('mudou depois da execução', motivo)

    def test_manifesto_sem_o_bloco_do_alinhador_recusa(self):
        self._manifesto(aplicado=True,
                        alinhamento='out/tmp/dataset_final_outro.aln')

        aln, motivo = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertIsNone(aln)
        self.assertIn('não registra', motivo)

    def test_manifesto_ilegivel_cai_na_regra_sem_proveniencia(self):
        with open(os.path.join(self.out, 'outputs', 'manifest.json'), 'w') as fh:
            fh.write('{ isto nao e json')

        aln, motivo = self.c._alinhamento_da_arvore('dataset_final_mafft')

        self.assertIsNone(aln)
        self.assertIn('proveniência do alinhamento indisponível', motivo)


class TestRegressaoDosCasosMedidos(_Out):
    """Os dois conjuntos do artigo, com as contagens medidas em 2026-11-19."""

    def _caso(self, colunas_antes, colunas_depois):
        self._manifesto(aplicado=True, colunas_antes=colunas_antes,
                        colunas_depois=colunas_depois)
        aln, info = self.c._alinhamento_da_arvore('dataset_final_mafft')
        return os.path.basename(aln), info['colunas']

    def test_varv48_data_as_185988_colunas_e_nao_as_191248(self):
        self.assertEqual(self._caso(191248, 185988),
                         ('dataset_final_mafft_lacuna.aln', 185988))

    def test_zikv479_data_as_10698_colunas_e_nao_as_10817(self):
        self.assertEqual(self._caso(10817, 10698),
                         ('dataset_final_mafft_lacuna.aln', 10698))


if __name__ == '__main__':
    unittest.main()
