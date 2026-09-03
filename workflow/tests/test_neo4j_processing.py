"""
Testes de `neo4jProcessing.py` — identidade de projeto no `.cql` gerado.

Antes deste lote, o nó `Tree` gravado no grafo não carregava nenhuma
propriedade que dissesse de qual projeto/experimento os dados vieram — um
grafo com mais de um projeto carregado não dava para separar por origem.
`project` vem da própria estrutura de diretórios do pipeline
(`projects/<project_name>/out/outputs`), não de um campo novo de config.

Executar com:

    python -m unittest workflow.tests.test_neo4j_processing
"""

from __future__ import annotations

import os
import tempfile
import shutil
import unittest

from workflow.utils.neo4jProcessing import (
    generate_cypher,
    parse_tree,
    project_name_from_output_path,
)


class TestProjectNameFromOutputPath(unittest.TestCase):
    def test_deriva_do_diretorio_do_projeto(self):
        path = "/home/geomesh/.../projects/Variola_VARV49_reexec_20260901/out/outputs"
        self.assertEqual(
            project_name_from_output_path(path), "Variola_VARV49_reexec_20260901"
        )

    def test_funciona_com_nome_de_projeto_simples(self):
        path = "/x/projects/teste52/out/outputs"
        self.assertEqual(project_name_from_output_path(path), "teste52")


class TestGenerateCypherCarregaProjeto(unittest.TestCase):
    def test_no_tree_grava_a_propriedade_project(self):
        statements = generate_cypher("nj_distance_mafft", {}, "Variola_VARV49_reexec_20260901")
        tree_stmt = statements[0]
        self.assertIn("project: 'Variola_VARV49_reexec_20260901'", tree_stmt)
        self.assertIn("name: 'nj_distance_mafft'", tree_stmt)


class TestParseTreeEscreveProjetoNoArquivo(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out_path = os.path.join(
            self.tmp, "projects", "Meu_Experimento_X", "out", "outputs"
        )
        os.makedirs(self.out_path, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_cql_gerado_identifica_o_projeto(self):
        trees = [{"nj_distance": {}}]
        parse_tree(trees=trees, path=self.out_path)

        gerado = os.path.join(self.out_path, "neo4j_commands_completo.cql")
        self.assertTrue(os.path.exists(gerado))
        with open(gerado) as f:
            conteudo = f.read()
        self.assertIn("project: 'Meu_Experimento_X'", conteudo)


if __name__ == "__main__":
    unittest.main()
