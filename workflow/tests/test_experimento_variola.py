"""
M2.1 — o baseline de Li *et al.* (2007) deixa de ser um bloco comentado.

Um experimento comentado no rodapé de um módulo não é reprodutível: não roda,
não é testado, e ninguém consegue afirmar que ainda corresponde ao que gerou os
artefatos em disco. Estes testes fixam o que o bloco antigo não garantia.

Executar com:

    python -m unittest workflow.tests.test_experimento_variola
"""

from __future__ import annotations

import os
import unittest
from unittest import mock

from workflow.experimentos import variola_li_2007 as exp


class TestAcessos(unittest.TestCase):

    def test_as_48_accessions_do_estudo(self):
        """O número vem de `08-ficha-de-fatos §4`, e é o que define o conjunto."""
        acessos = exp.carregar_acessos()
        self.assertEqual(len(acessos), 48)
        self.assertEqual(len(set(acessos)), 48, "há acesso repetido na lista")

    def test_faixas_do_estudo(self):
        acessos = exp.carregar_acessos()
        self.assertIn("DQ437580", acessos)
        self.assertIn("DQ441448", acessos)
        self.assertTrue(all(a.startswith("DQ4") for a in acessos))

    def test_comentario_e_linha_vazia_sao_ignorados(self):
        """O arquivo de acessos carrega o aviso de D23 no cabeçalho; se os
        comentários entrassem na consulta, o Entrez receberia lixo."""
        acessos = exp.carregar_acessos()
        self.assertFalse([a for a in acessos if a.startswith("#") or not a])


class TestEmail(unittest.TestCase):
    """O e-mail do Entrez é dado pessoal e não entra em repositório (regra 7)."""

    def test_ausencia_falha_com_instrucao(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(RuntimeError) as ctx:
                exp.email_do_ambiente()
        self.assertIn("NCBI_EMAIL", str(ctx.exception))

    def test_lido_do_ambiente(self):
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "a@b.br"}):
            self.assertEqual(exp.email_do_ambiente(), "a@b.br")

    def test_sem_email_no_ambiente_nao_ha_recurso_embutido(self):
        """A propriedade que importa é comportamental, não textual: sem
        `NCBI_EMAIL`, **não existe caminho** que monte o workflow. Um teste que
        procurasse a string `@` no fonte erraria na própria documentação — e,
        pior, passaria se alguém pusesse o endereço numa variável."""
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(RuntimeError):
                exp.montar_workflow("/tmp/qualquer")


class TestParametros(unittest.TestCase):

    def test_parametros_do_estudo(self):
        self.assertEqual(exp.PARAMETROS["initial_min_length"], 180_000)
        self.assertEqual(exp.PARAMETROS["refined_min_length"], 183_000)
        self.assertEqual(exp.PARAMETROS["similarity_threshold"], 0.999)
        self.assertEqual(exp.PARAMETROS["retmax"], 200)

    def test_utr_e_none_para_virus(self):
        """`None` é 'não se aplica', e não um corte em zero: vírus não têm UTR
        no sentido eucarioto. É a regra 5 do projeto."""
        self.assertIsNone(exp.PARAMETROS["utr5_end"])
        self.assertIsNone(exp.PARAMETROS["utr3_start"])

    def test_grupo_externo_declarado_por_organismo(self):
        self.assertIn("Taterapox virus", exp.OUTGROUP_QUERY)
        self.assertIn("Camelpox virus", exp.OUTGROUP_QUERY)


class TestD23(unittest.TestCase):
    """
    A consulta traz RefSeq e GenBank do mesmo genoma. Não é corrigido aqui — a
    decisão do usuário foi declarar e corrigir depois —, mas tem de ser
    **nomeado**, e não descoberto de novo daqui a seis meses.
    """

    def test_os_pares_conhecidos_estao_declarados(self):
        pares = dict(exp.ACESSOS_DUPLICADOS_CONHECIDOS)
        self.assertEqual(pares["NC_008291"], "DQ437594")   # Taterapox
        self.assertEqual(pares["NC_003391"], "AF438165")   # Camelpox

    def test_par_presente_e_apontado(self):
        encontrados = exp.conferir_duplicatas(["NC_008291.1", "DQ437594.1", "DQ437580"])
        self.assertEqual(encontrados, [("NC_008291", "DQ437594")])

    def test_versao_do_acesso_nao_esconde_o_par(self):
        """`NC_008291.1` e `NC_008291` são o mesmo acesso — é a armadilha de
        D13 noutro lugar."""
        self.assertTrue(exp.conferir_duplicatas(["NC_003391.2", "AF438165.1"]))

    def test_lista_sem_par_nao_acusa(self):
        self.assertEqual(exp.conferir_duplicatas(exp.carregar_acessos()), [])


if __name__ == "__main__":
    unittest.main()
