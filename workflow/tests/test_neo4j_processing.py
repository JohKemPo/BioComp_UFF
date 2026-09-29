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
from workflow.utils.deduplicar_metadata_cql import (
    blocos,
    contar_original,
    contar_reparado,
    reparar_texto,
    segmentos,
    verificar,
)
from workflow.tests.dados_neo4j_g1 import PROJETO, arvores_fixture, ocorrencias_esperadas


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


class TestG1MetadataUmNoPorTerminal(unittest.TestCase):
    """G1 (`docs/science/05-grafo-neo4j.md`): `Metadata` por terminal, não por ocorrência.

    Oráculo: o golden `cql_legado_pre_G1.cql`, gerado pelo código ANTES da
    correção sobre a mesma fixture. A contagem de arestas HAS_METADATA e a
    sequência (subárvore, terminal) têm de ser iguais; só os nós deixam de
    se repetir, e o registro GenBank de cada terminal tem de sair idêntico.
    """

    GOLDEN = os.path.join(os.path.dirname(__file__), "cql_legado_pre_G1.cql")

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out_path = os.path.join(self.tmp, "projects", PROJETO, "out", "outputs")
        os.makedirs(self.out_path, exist_ok=True)
        parse_tree(trees=arvores_fixture(), path=self.out_path)
        with open(os.path.join(self.out_path, "neo4j_commands_completo.cql")) as f:
            self.novo = f.read()
        with open(self.GOLDEN) as f:
            self.legado = f.read()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_golden_legado_tem_a_duplicacao_caracterizada(self):
        antes = contar_original(self.legado)
        self.assertEqual(
            (antes["Metadata"], antes["Feature"], antes["Qualifier"], antes["HAS_METADATA"]),
            (11, 8, 12, 11))

    def test_um_metadata_por_terminal_e_mesmas_ligacoes(self):
        depois = contar_reparado(self.novo, PROJETO)
        self.assertEqual(depois["Metadata"], 3)
        self.assertEqual(depois["Feature"], 4)       # T1: 2, T2: 2, T3: 0
        self.assertEqual(depois["Qualifier"], 6)     # T1: 3+1, T2: 2+0
        self.assertEqual(depois["HAS_METADATA"], 11)
        self.assertEqual(depois["ligacoes"],
                         [(s, str(h)) for s, h in ocorrencias_esperadas()])

    def test_conteudo_de_cada_terminal_identico_ao_legado(self):
        erros = verificar(contar_original(self.legado), contar_reparado(self.novo, PROJETO))
        self.assertEqual(erros, [])

    def test_primeira_ocorrencia_vazia_e_completada_depois(self):
        # T1 aparece primeiro sem registro (iqtree) e completo só no fasttree:
        # duas criações — a vazia e a completa —, e só a completa tem features.
        criacoes = contar_reparado(self.novo, PROJETO)["criacoes"]["1111111111111111"]
        self.assertEqual(len(criacoes), 2)
        self.assertFalse(criacoes[0][1])
        self.assertTrue(criacoes[1][0].startswith('molecule_type: "RNA"'))
        # T3 aparece completo primeiro: uma criação só, mesmo com ocorrência vazia depois.
        self.assertEqual(len(contar_reparado(self.novo, PROJETO)["criacoes"]["3333333333333333"]), 1)

    def test_chave_inclui_uid_e_projeto(self):
        self.assertIn(
            "MERGE (m:Metadata {uid: '<<USER_UID>>', project: \"Projeto_Fixture_G1\", "
            "terminal_hash: 1111111111111111})", self.novo)
        self.assertNotIn("CREATE (m:Metadata", self.novo)

    def test_features_protegidas_contra_reexecucao(self):
        for bloco in blocos(self.novo):
            if "HAS_FEATURE]->(f_" in bloco:
                self.assertIn("WITH m WHERE NOT EXISTS { (m)-[:HAS_FEATURE]->() }", bloco)

    def test_reparo_do_legado_reproduz_o_gerador_byte_a_byte(self):
        reparado = "".join(reparar_texto(self.legado, PROJETO))
        self.assertEqual(reparado, self.novo)

    def test_reparo_e_idempotente(self):
        reparado = "".join(reparar_texto(self.legado, PROJETO))
        self.assertEqual("".join(reparar_texto(reparado, PROJETO)), reparado)

    def test_blocos_do_reparo_batem_com_o_tokenizer_do_backend(self):
        # Mesmo número de blocos que o laço caractere a caractere de
        # `parse_cql_blocks` (sem a remoção de comentário, ver módulo).
        self.assertEqual(list(blocos(self.novo)), _parse_cql_blocks_referencia(self.novo))
        self.assertEqual(list(blocos(self.legado)), _parse_cql_blocks_referencia(self.legado))

    def test_segmentos_preservam_o_texto(self):
        self.assertEqual("".join(segmentos(self.legado)), self.legado)

    def test_terminal_sem_hash_falha_declarado(self):
        arvores = [{"t": {"t_Inner1": {"data_terminals": [
            {"newick": "X1", "terminal_hash": None, "metadata": {}}]}}}]
        with self.assertRaises(ValueError):
            parse_tree(trees=arvores, path=self.out_path)

    def test_reparo_recusa_registros_divergentes(self):
        adulterado = self.legado.replace(
            '"passage 1; Vero cells"', '"passage 2; Vero cells"', 1)
        with self.assertRaises(ValueError):
            "".join(reparar_texto(adulterado, PROJETO))


def _parse_cql_blocks_referencia(content: str) -> list:
    """Cópia do laço de `Backend/src/services/cql_batch_service.py::parse_cql_blocks`
    (C-5e), sem a remoção de comentários."""
    blocks, current = [], []
    in_single_quote = in_double_quote = escaped = False
    for char in content:
        current.append(char)
        if escaped:
            escaped = False
            continue
        if char == '\\':
            escaped = True
            continue
        if char == "'" and not in_double_quote:
            in_single_quote = not in_single_quote
        elif char == '"' and not in_single_quote:
            in_double_quote = not in_double_quote
        elif char == ';' and not in_single_quote and not in_double_quote:
            block = ''.join(current).strip()
            if block:
                blocks.append(block)
            current = []
    last = ''.join(current).strip()
    if last:
        blocks.append(last)
    return blocks


if __name__ == "__main__":
    unittest.main()
