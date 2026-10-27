"""
Testes do D29 — o item da mineração é a bipartição canônica, não o clado escrito.

O defeito: `SubtreeBuilder` gravava como `List_terminals_hash` o conjunto de
terminais do nó **como está escrito no arquivo**. A raiz que cada programa
escreve é convenção (FastTree, IQ-TREE, RAxML e NJ: raiz trifurcante; UPGMA e
NJ/parcimônia: bifurcante), então a mesma bipartição escrita pelos dois lados
virava dois itens, o suporte gravado só podia ficar menor, e o clado-raiz entrava
como item trivial. O D3 fora corrigido no relatório de estabilidade e não aqui.
Parecer: docs/manuscritos/frente-2026-10/bmc/06-parecer-d29.md, §4.2 (os 11
casos-limite, numerados abaixo).

Executar com:

    python -m unittest workflow.tests.test_item_biparticao
"""

from __future__ import annotations

import glob
import json
import os
import re
import shutil
import tempfile
import unittest
from unittest import mock

import pandas as pd

from workflow import remine
from workflow.stability.clade_identity import canonical_item_id
from workflow.subtree_construction.builder import SubtreeBuilder
from workflow.subtree_mining.miner import SubtreeMiner
from workflow.utils import manifest
from workflow.utils.treeUtils import (encode_bipartition_to_int, encode_clade_to_int,
                                      extract_terminais_list,
                                      verificar_mesmo_conjunto_de_folhas)

T = [f"T{i}" for i in range(1, 9)]

#: Quatro pares (T1T2 | T3T4 | T5T6 | T7T8) e a aresta central {T1..T4 | T5..T8}:
#: 5 bipartições não triviais, escritas de formas diferentes abaixo.
RAIZ_TRIFURCANTE = "((T1,T2),(T3,T4),((T5,T6),(T7,T8)));"
RAIZ_NO_T1 = "(T1,(T2,((T3,T4),((T5,T6),(T7,T8)))));"
RAIZ_BIFURCANTE = "(((T1,T2),(T3,T4)),((T5,T6),(T7,T8)));"


class _TodosBaixados:
    """Cache que diz já ter toda sequência: o teste não vai ao NCBI."""
    def __contains__(self, _):
        return True

    def add(self, _):
        pass


class _Base(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)
        os.makedirs(os.path.join(self.dir, "out", "Subtrees"))
        os.makedirs(os.path.join(self.dir, "out", "outputs"))
        patcher = mock.patch("workflow.utils.run_logging.garantir")
        patcher.start()
        self.addCleanup(patcher.stop)

    def construir(self, newick, nome="t1"):
        """Roda o `SubtreeBuilder` real sobre um Newick; devolve {subárvore: registro}."""
        caminho = os.path.join(self.dir, f"{nome}.nwk")
        with open(caminho, "w") as f:
            f.write(newick)
        builder = SubtreeBuilder(_TodosBaixados(), input_format="newick", output_format="nexus",
                                 output_path=os.path.join(self.dir, "out"), resume_infos=False)
        return builder.subtree_constructor(caminho, f"{nome}.nwk")[nome]

    def transacao(self, newick, nome="t1"):
        return extract_terminais_list({nome: self.construir(newick, nome)})


def item(*taxa, todos=T):
    return encode_bipartition_to_int(list(taxa), todos)


class TestCasosLimite(_Base):
    """Os 11 casos-limite do parecer G9 §4.2."""

    def test_01_mesma_topologia_com_raizes_diferentes_da_a_mesma_transacao(self):
        """O núcleo do D29: onde o arquivo pôs a raiz não pode mudar a transação."""
        a = self.transacao(RAIZ_TRIFURCANTE)
        b = self.transacao(RAIZ_NO_T1, "t2")
        c = self.transacao(RAIZ_BIFURCANTE, "t3")
        self.assertEqual(len(a), 5)
        self.assertEqual(set(a), set(b))
        self.assertEqual(set(a), set(c))

    def test_02_raiz_bifurcante_com_dois_filhos_grandes_e_um_item_so(self):
        """UPGMA: os dois filhos da raiz são a mesma aresta."""
        registros = self.construir(RAIZ_BIFURCANTE)
        filhos = [r["List_terminals_hash"] for n, r in registros.items()
                  if len(r["Terminals"]) == 4]
        self.assertEqual(len(filhos), 2)
        self.assertEqual(filhos[0], filhos[1])
        self.assertEqual(len(extract_terminais_list({"t1": registros})), 5)

    def test_03_filho_da_raiz_com_n_menos_1_taxons_nao_e_item(self):
        registros = self.construir(RAIZ_NO_T1)
        n_menos_1 = [r for r in registros.values() if len(r["Terminals"]) == 7]
        self.assertEqual(len(n_menos_1), 1)
        self.assertIsNone(n_menos_1[0]["List_terminals_hash"])
        self.assertIsNone(item(*T[1:]))

    def test_04_clado_raiz_nao_e_item(self):
        registros = self.construir(RAIZ_TRIFURCANTE)
        raiz = [r for r in registros.values() if len(r["Terminals"]) == 8]
        self.assertEqual(len(raiz), 1)
        self.assertIsNone(raiz[0]["List_terminals_hash"])
        # o clado como escrito continua gravado, só para auditoria
        self.assertEqual(raiz[0]["List_terminals_hash_rooted"], encode_clade_to_int(T))

    def test_05_politomia_interna_tem_menos_biparticoes_sem_erro(self):
        transacao = self.transacao("((T1,T2,T3,T4),(T5,T6),(T7,T8));")
        self.assertEqual(len(transacao), 3)
        self.assertLess(len(transacao), 5)

    def test_06_rotulo_truncado_da_o_mesmo_item(self):
        """D13: IQ-TREE grava `NC_008030.` onde FastTree grava `NC_008030.1`."""
        integro = "((NC_008030.1,T2),(T3,T4),((T5,T6),(T7,T8)));"
        truncado = "((NC_008030.,T2),(T3,T4),((T5,T6),(T7,T8)));"
        self.assertEqual(set(self.transacao(integro)), set(self.transacao(truncado, "t2")))
        todos = ["NC_008030.1"] + T[1:]
        todos_trunc = ["NC_008030."] + T[1:]
        self.assertEqual(item("NC_008030.1", "T2", todos=todos), item("NC_008030.", "T2", todos=todos_trunc))

    def test_07_empate_de_tamanho_desempata_pela_ordem_lexicografica(self):
        """|A| = n/2: os dois lados da aresta central dão o mesmo item, o de T1..T4."""
        lado_a = item("T1", "T2", "T3", "T4")
        lado_b = item("T5", "T6", "T7", "T8")
        self.assertEqual(lado_a, lado_b)
        self.assertEqual(lado_a, canonical_item_id(["T1", "T2", "T3", "T4"]))
        # e no construtor: a aresta central escrita só pelo lado de T1..T4 ou só pelo de T5..T8
        outro_lado = "((T5,T6),(T7,T8),((T1,T2),(T3,T4)));"
        self.assertEqual(set(self.transacao(RAIZ_TRIFURCANTE)), set(self.transacao(outro_lado, "t2")))

    def test_08_n_menor_que_4_ou_so_arestas_triviais_e_transacao_vazia(self):
        self.assertEqual(self.transacao("(T1,(T2,T3));"), [])
        self.assertEqual(self.transacao("(T1,T2,T3,T4);", "t2"), [])

    def test_08_mineracao_nao_aplicavel_nao_grava_csv_nem_chama_fpmax(self):
        """Nenhuma árvore com item: o resultado não é "0 padrões", é "não aplicável"."""
        dados = [{"t1": self.construir("(T1,(T2,T3));")}, {"t2": self.construir("(T1,T2,T3);", "t2")}]
        minerador = SubtreeMiner(mode="OFST", support_fpmax="auto",
                                 output_path=os.path.join(self.dir, "out"))
        with mock.patch("workflow.subtree_mining.miner.fpmax") as fpmax, \
                self.assertLogs(level="WARNING") as registro:
            resultado = minerador.process_group(dados)
        fpmax.assert_not_called()
        self.assertIs(resultado, dados)
        self.assertFalse(os.path.exists(os.path.join(self.dir, "out", "outputs", "all_results_fpmax.csv")))
        self.assertIn("não aplicável", registro.output[0])

    def test_09_conjuntos_de_folhas_diferentes_sao_recusados(self):
        dados = [{"t1": self.construir(RAIZ_TRIFURCANTE)},
                 {"t2": self.construir("((T1,T2),(T3,T4),((T5,T6),(T7,T9)));", "t2")}]
        with self.assertRaisesRegex(ValueError, "mesmo conjunto de táxons"):
            verificar_mesmo_conjunto_de_folhas(dados)
        minerador = SubtreeMiner(mode="OFST", support_fpmax="auto",
                                 output_path=os.path.join(self.dir, "out"))
        with self.assertRaises(ValueError):
            minerador.process_group(dados)

    def test_09_mesmo_conjunto_de_folhas_em_ordens_diferentes_passa(self):
        verificar_mesmo_conjunto_de_folhas(
            [{"t1": self.construir(RAIZ_TRIFURCANTE)}, {"t2": self.construir(RAIZ_NO_T1, "t2")}])

    def test_10_no_unario_repete_a_biparticao_e_conta_uma_vez(self):
        registros = self.construir("((T1,T2),(T3,T4),(((T5,T6)),(T7,T8)));")
        itens_dos_nos = [r["List_terminals_hash"] for r in registros.values()
                         if r["List_terminals_hash"] is not None]
        self.assertGreater(len(itens_dos_nos), len(set(itens_dos_nos)))
        self.assertEqual(len(extract_terminais_list({"t1": registros})), 5)

    def test_11_arvore_pre_rotulada_com_raiz_sem_nome_nao_perde_item(self):
        """A raiz sem nome recebe `Inner1`, que a árvore já usa (R13). Com o clado-raiz
        fora dos itens, a colisão deixa de apagar um item da transação."""
        pre_rotulada = "((T1,T2)Inner2,(T3,T4)Inner3,((T5,T6)Inner4,(T7,T8)Inner5)Inner1);"
        self.assertEqual(set(self.transacao(pre_rotulada)), set(self.transacao(RAIZ_TRIFURCANTE, "t2")))

    @unittest.expectedFailure
    def test_11_colisao_de_nome_inner1_ainda_apaga_um_registro(self):
        """R13, fora do D29 (lote próprio): `_add_internal_labels` não pula nome já
        usado, e a raiz e o `Inner1` verdadeiro dividem a chave. Deve passar quando o
        lote da colisão entrar, e então perder o `expectedFailure`."""
        pre_rotulada = "((T1,T2)Inner2,(T3,T4)Inner3,((T5,T6)Inner4,(T7,T8)Inner5)Inner1);"
        self.assertEqual(len(self.construir(pre_rotulada)), 6)


class TestItem(unittest.TestCase):
    def test_o_lado_escrito_nao_muda_o_item(self):
        self.assertEqual(item("T1", "T2"), item(*T[2:]))

    def test_terminal_e_clado_de_n_menos_1_nao_sao_item(self):
        self.assertIsNone(item("T1"))
        self.assertIsNone(item(*T[:7]))

    def test_o_clado_como_escrito_continua_nao_canonico(self):
        """O que `encode_clade_to_int` não faz, e agora diz que não faz."""
        self.assertNotEqual(encode_clade_to_int(["T1", "T2"]), encode_clade_to_int(T[2:]))

    def test_a_transacao_ignora_none_e_repeticao_e_mantem_a_ordem(self):
        dados = {"t": {"a": {"List_terminals_hash": 3}, "b": {"List_terminals_hash": None},
                       "c": {"List_terminals_hash": 1}, "d": {"List_terminals_hash": 3}}}
        self.assertEqual(extract_terminais_list(dados), [3, 1])


class TestProvenienciaDaRemineracao(unittest.TestCase):
    def setUp(self):
        self.raiz = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.raiz, ignore_errors=True)
        self.out_dir = os.path.join(self.raiz, "out")
        os.makedirs(os.path.join(self.out_dir, "outputs"))

    def _manifesto(self, **extra):
        dados = {"run_id": "abc", "git": {"BioComp_UFF": {"commit": "x"}},
                 "outputs_sha256": {"out/outputs/a.csv": "velho"}, **extra}
        with open(os.path.join(self.out_dir, "outputs", manifest.MANIFEST_FILENAME), "w") as f:
            json.dump(dados, f)

    def _lido(self):
        with open(os.path.join(self.out_dir, "outputs", manifest.MANIFEST_FILENAME)) as f:
            return json.load(f)

    def _anexar(self, motivo="D29"):
        return manifest.anexar_remineracao(
            self.out_dir, motivo, {"out/Trees/t.nwk": "e1"}, {"out/outputs/a.csv": "velho"},
            {"out/outputs/a.csv": "novo"}, repos={}, detalhes={"padroes_gravados": 49})

    def test_sem_manifesto_nao_inventa_um(self):
        self.assertFalse(self._anexar())
        self.assertFalse(os.path.exists(os.path.join(self.out_dir, "outputs", manifest.MANIFEST_FILENAME)))

    def test_anexa_sem_apagar_o_que_havia(self):
        self._manifesto()
        antes = self._lido()
        self.assertTrue(self._anexar())
        depois = self._lido()
        for chave in ("run_id", "git"):
            self.assertEqual(antes[chave], depois[chave])
        registro = depois["remineracoes"][0]
        self.assertEqual((registro["via"], registro["motivo"], registro["padroes_gravados"]), ("remine", "D29", 49))
        self.assertEqual(registro["saidas_anteriores_sha256"], {"out/outputs/a.csv": "velho"})
        self.assertEqual(depois["outputs_sha256"], {"out/outputs/a.csv": "novo"})

    def test_remineracoes_se_acumulam(self):
        self._manifesto()
        self._anexar("D29")
        self._anexar("outro")
        self.assertEqual([r["motivo"] for r in self._lido()["remineracoes"]], ["D29", "outro"])


class TestRemineracaoDePonta(_Base):
    """`remine.remineirar` num projeto pequeno: sem rede, sem tocar nas árvores,
    com proveniência, e com o resultado que o D29 exige."""

    def _projeto(self, com_manifesto=True):
        projeto = os.path.join(self.dir, "proj_teste")
        out = os.path.join(projeto, "out")
        os.makedirs(os.path.join(out, "Trees"))
        os.makedirs(os.path.join(out, "Subtrees"))
        os.makedirs(os.path.join(out, "outputs"))
        for nome, newick in (("tree_dataset_a.nwk", RAIZ_TRIFURCANTE), ("tree_dataset_b.nwk", RAIZ_NO_T1),
                             ("tree_dataset_c.nwk", RAIZ_BIFURCANTE)):
            with open(os.path.join(out, "Trees", nome), "w") as f:
                f.write(newick)
        open(os.path.join(out, "outputs", "raw_data_sequences.gb"), "w").close()
        cfg = {"subtree_config": {"construct_tree_method": "distance", "input_format": "newick",
                                  "output_format": "nexus", "resume_infos": False, "save_metadata": True,
                                  "subtree_miner": True,
                                  "subtree_miner_configs": {"mode": "OFST", "save_fpmax": True,
                                                            "support_fpmax": "auto"}}}
        with open(os.path.join(out, "outputs", "config_backup.json"), "w") as f:
            json.dump(cfg, f)
        with open(os.path.join(out, "outputs", "all_results_fpmax.csv"), "w") as f:
            f.write("antigo\n")
        if com_manifesto:
            with open(os.path.join(out, "outputs", manifest.MANIFEST_FILENAME), "w") as f:
                json.dump({"run_id": "orig", "outputs_sha256": {}}, f)
        return out

    def test_tres_raizes_da_mesma_topologia_dao_um_padrao_de_cinco_itens(self):
        out = self._projeto()
        trees_antes = {f: manifest.file_digest(os.path.join(out, "Trees", f))
                       for f in os.listdir(os.path.join(out, "Trees"))}
        resumo = remine.remineirar(out, "D29")

        self.assertEqual(resumo["padroes_gravados"], 1)
        self.assertEqual(set(resumo["itens_por_arvore"].values()), {5})
        tabela = pd.read_csv(os.path.join(out, "outputs", "all_results_fpmax.csv"))
        self.assertEqual(tabela["n_trees"].tolist(), [3])
        self.assertEqual(len(re.findall(r"\d{6,}", tabela["itemsets"][0])), 5)
        self.assertEqual({f: manifest.file_digest(os.path.join(out, "Trees", f))
                          for f in os.listdir(os.path.join(out, "Trees"))}, trees_antes)

    def test_anexa_ao_manifesto_sem_regravar_do_zero(self):
        out = self._projeto()
        remine.remineirar(out, "D29")
        with open(os.path.join(out, "outputs", manifest.MANIFEST_FILENAME)) as f:
            dados = json.load(f)
        self.assertEqual(dados["run_id"], "orig")
        registro, = dados["remineracoes"]
        self.assertEqual(registro["motivo"], "D29")
        self.assertIn("out/Trees/tree_dataset_a.nwk", registro["entradas_sha256"])
        self.assertIn("out/outputs/all_results_fpmax.csv", registro["saidas_anteriores_sha256"])
        self.assertEqual(dados["outputs_sha256"]["out/outputs/all_results_fpmax.csv"],
                         manifest.file_digest(os.path.join(out, "outputs", "all_results_fpmax.csv")))
        self.assertEqual(glob.glob(os.path.join(os.path.dirname(os.path.dirname(out)), ".remine_*")), [])

    def test_metadado_de_terminal_e_ordem_das_arvores_vem_do_metadata_json_anterior(self):
        """O `.gb` que sobra no projeto é o da última rodada de download (em VARV-52,
        3 de 51 registros): relê-lo trocaria metadado por `{"error": ...}`. O D29
        não muda `data_terminals`, então a remineração tem de preservá-lo."""
        out = self._projeto()
        antigo = [[{"tree_dataset_c": {}}, {"tree_dataset_a": {"tree_dataset_a_Inner1": {"data_terminals": [
            {"newick": "T1", "metadata": {"organism": "Variola", "nota": 1.5}}]}}}]]
        with open(os.path.join(out, "outputs", "metadata.json"), "w") as f:
            json.dump(antigo, f)

        resumo = remine.remineirar(out, "D29")

        with open(os.path.join(out, "outputs", "metadata.json")) as f:
            novo = json.load(f)
        metadados = {}
        for arvore in novo[0]:
            for nome, subarvores in arvore.items():
                for registro in subarvores.values():
                    for terminal in registro["data_terminals"]:
                        metadados.setdefault((nome, terminal["newick"]), terminal["metadata"])
        # a ordem das árvores é a da execução original (c, a), e o que ela não conhecia vem depois (b)
        self.assertEqual([nome for arvore in novo[0] for nome in arvore],
                         ["tree_dataset_c", "tree_dataset_a", "tree_dataset_b"])
        self.assertEqual(metadados[("tree_dataset_a", "T1")], {"organism": "Variola", "nota": 1.5})
        self.assertIn("error", metadados[("tree_dataset_a", "T2")])
        self.assertIn("error", metadados[("tree_dataset_b", "T1")])
        self.assertEqual(resumo["metadado_de_terminal"]["preservado_do_metadata_json_anterior"], 1)
        self.assertEqual(resumo["metadado_de_terminal"]["lido_do_gb"], 23)

    def test_projeto_sem_manifesto_e_recusado_sem_tocar_em_nada(self):
        out = self._projeto(com_manifesto=False)
        with self.assertRaises(FileNotFoundError):
            remine.remineirar(out, "D29")
        with open(os.path.join(out, "outputs", "all_results_fpmax.csv")) as f:
            self.assertEqual(f.read(), "antigo\n")

    def test_guarda_contra_download_levanta_erro_em_vez_de_truncar_o_gb(self):
        with self.assertRaisesRegex(RuntimeError, "não baixa"):
            remine._proibido()


if __name__ == "__main__":
    unittest.main()
