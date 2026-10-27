"""
Testes das CLIs do workflow granular (`workflow.akoflow`) e do comparador de equivalência.

Não precisam de MAFFT, IQ-TREE nem RAxML-NG: cobrem as atividades que são só
orquestração de funções de biblioteca (preparar, baixar, subarvores, minerar,
gerar-cql, graficos) e provam que a saída delas é a do `SubtreeBuilderController`
sobre as mesmas árvores. A equivalência com o monolito inteiro (as ferramentas
externas) é o teste do Zika-21 de `docs/akoflow/06-g4b-relatorio.md`.
"""

import collections
import contextlib
import io
import json
import os
import tempfile
import unittest
from unittest import mock

from workflow import akoflow, equivalencia
from workflow.utils import tool_runs

FASTA = """>AB000001.1
ACGTACGTACGTAAAA
>AB000002.1
ACGTACGTACGTCCCC
>AB000003.1
ACGTACGTACGTAAAA
>AB000004.1
TTGTACGTACGTAAAC
"""
ARVORES = {
    "tree_dataset_a_nj": "((AB000001.1:1,AB000002.1:1):1,(AB000003.1:1,(AB000004.1:1,AB000005.1:1):1):1);",
    "tree_dataset_a_fast": "((AB000001.1:1,AB000003.1:1):1,(AB000002.1:1,(AB000004.1:1,AB000005.1:1):1):1);",
}
CONFIG = {
    "log_file": True, "project_name": "Proj_teste", "output_log": "#",
    "tree_config": {"mode": "advanced", "ignore_mode": ["mrbayes"], "num_threads": 1,
                    "random_seed": 12345, "raxml_threads": 1, "iqtree_threads": 1,
                    "output_format": "nexus", "on_method_failure": "continue"},
    "subtree_config": {"construct_tree_method": "distance", "input_format": "nexus",
                       "output_format": "nexus", "resume_infos": False, "save_metadata": True,
                       "subtree_miner": True,
                       "subtree_miner_configs": {"mode": "OFST", "save_fpmax": True,
                                                 "support_fpmax": "auto"}},
}


def _executar(*argv):
    with contextlib.redirect_stdout(io.StringIO()):
        akoflow.main(list(argv))


class Base(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.raiz = os.path.join(self._tmp.name, "Proj_teste")
        self.out = os.path.join(self.raiz, "out")
        self.cfg = os.path.join(self._tmp.name, "cfg.json")
        with open(self.cfg, "w") as fh:
            json.dump(CONFIG, fh)
        self.entrada = os.path.join(self._tmp.name, "entrada")
        os.makedirs(self.entrada)
        with open(os.path.join(self.entrada, "dataset_a.fasta"), "w") as fh:
            fh.write(FASTA)
        self.backup = os.path.join(self.out, "outputs", "config_backup.json")

    def preparar(self):
        _executar("preparar", "--config", self.cfg, "--entrada", self.entrada, "--out", self.out)

    def prov(self, nome):
        with open(os.path.join(self.out, "outputs", "prov", f"{nome}.json")) as fh:
            return json.load(fh)


class TestPreparar(Base):
    def test_deduplica_sem_tocar_na_entrada_e_grava_acessos(self):
        antes = os.listdir(self.entrada)
        self.preparar()
        self.assertEqual(os.listdir(self.entrada), antes, "AC-09: a entrada não pode receber _NoPipe")
        self.assertTrue(os.path.isfile(os.path.join(self.out, "entrada", "dataset_a_NoPipe")))
        with open(os.path.join(self.out, "entrada", "acessos.txt")) as fh:
            self.assertEqual(fh.read().split(), ["AB000001.1", "AB000002.1", "AB000004.1"])

    def test_registra_a_deduplicacao_e_nao_vaza_caminho_absoluto(self):
        self.preparar()
        frag = self.prov("preparar-entrada")
        dedup = frag["tools_invoked"]["deduplicacao"]["descartados"]
        self.assertEqual(dedup, [{"descartado": "AB000003.1", "mantido": "AB000001.1"}])
        self.assertNotIn(self._tmp.name, json.dumps(frag), "D15: caminho absoluto no fragmento")

    def test_config_backup_tem_caminhos_relativos_e_nome_do_projeto(self):
        self.preparar()
        with open(self.backup) as fh:
            cfg = json.load(fh)
        self.assertEqual(cfg["tree_config"]["input_path"], "out/entrada")
        self.assertEqual(cfg["project_name"], "Proj_teste")
        self.assertEqual(cfg["tree_config"]["random_seed"], 12345)

    def test_sem_fasta_falha(self):
        os.remove(os.path.join(self.entrada, "dataset_a.fasta"))
        with self.assertRaises(FileNotFoundError):
            self.preparar()


class TestBaixar(Base):
    def test_gb_congelado_e_copiado_byte_a_byte(self):
        self.preparar()
        gb = os.path.join(self._tmp.name, "x.gb")
        with open(gb, "w") as fh:
            fh.write("LOCUS       X\n//\n")
        _executar("baixar-metadados", "--config", self.backup, "--out", self.out,
                  "--acessos", os.path.join(self.out, "entrada", "acessos.txt"), "--gb-congelado", gb)
        destino = os.path.join(self.out, "outputs", "raw_data_sequences.gb")
        self.assertEqual(open(destino).read(), open(gb).read())
        self.assertEqual(self.prov("baixar-metadados")["chave"]["modo"], "congelado")

    def test_download_com_erro_vira_falha_declarada(self):
        self.preparar()
        with mock.patch("workflow.utils.treeUtils.download_sequences", return_value={"error": "sem rede"}):
            with self.assertRaises(RuntimeError):
                _executar("baixar-metadados", "--config", self.backup, "--out", self.out,
                          "--acessos", os.path.join(self.out, "entrada", "acessos.txt"))
        frag = self.prov("baixar-metadados")
        self.assertEqual(frag["desfecho"]["status"], "falhou")
        self.assertIn("sem rede", frag["desfecho"]["motivo"])


class TestSubarvoresAteCql(Base):
    """Granular × `SubtreeBuilderController` sobre as mesmas árvores."""

    def montar(self, destino_out):
        from Bio import Phylo
        os.makedirs(os.path.join(destino_out, "Trees"), exist_ok=True)
        os.makedirs(os.path.join(destino_out, "outputs", "Plots"), exist_ok=True)
        for nome, newick in ARVORES.items():
            import io as _io
            arvore = Phylo.read(_io.StringIO(newick), "newick")
            Phylo.write(arvore, os.path.join(destino_out, "Trees", f"{nome}.nexus"), "nexus")
        open(os.path.join(destino_out, "outputs", "raw_data_sequences.gb"), "w").close()

    def rodar_granular(self):
        self.preparar()
        self.montar(self.out)
        for nome in ARVORES:
            _executar("subarvores", "--config", self.backup, "--out", self.out, "--stem", "dataset_a",
                      "--arvore", os.path.join(self.out, "Trees", f"{nome}.nexus"))
        _executar("minerar", "--config", self.backup, "--out", self.out, "--stem", "dataset_a")
        _executar("gerar-cql", "--config", self.backup, "--out", self.out, "--stem", "dataset_a")
        _executar("graficos-suporte", "--config", self.backup, "--out", self.out, "--stem", "dataset_a")

    def rodar_monolito(self):
        from workflow.controller.subtreeBuilderController import SubtreeBuilderController
        from workflow.subtree_construction import builder as builder_mod
        ref = os.path.join(self._tmp.name, "Ref", "out")
        self.montar(ref)
        cfg = json.loads(json.dumps(CONFIG["subtree_config"]))
        cfg.update(input_path=os.path.join(ref, "Trees"), output_path=ref)
        cfg["subtree_miner_configs"]["output_path"] = ref
        with mock.patch.object(builder_mod, "download_sequences", lambda *a, **k: None), \
                contextlib.redirect_stdout(io.StringIO()):
            controlador = SubtreeBuilderController(**cfg)
            controlador.files = sorted(controlador.files)
            controlador()
        return ref

    def setUp(self):
        super().setUp()
        # o erro de "acesso não encontrado" cita o caminho absoluto do .gb, que muda de diretório
        patcher = mock.patch("workflow.utils.treeUtils.fetch_local_record",
                             lambda acesso, gbk: {"id": acesso})
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_fpmax_metadata_e_cql_iguais_aos_do_monolito(self):
        self.rodar_granular()
        ref = self.rodar_monolito()
        resultado = equivalencia.comparar(ref, self.out)
        for classe in ("S2", "S3", "S4"):
            self.assertEqual(resultado[classe]["diferentes"], 0, (classe, resultado[classe]))
        self.assertGreater(resultado["S2"]["iguais"], 0)
        cql = [f for f in os.listdir(os.path.join(self.out, "outputs")) if f.endswith(".cql")]
        self.assertEqual(cql, ["neo4j_commands_tree_dataset.cql"])
        with open(os.path.join(self.out, "outputs", cql[0])) as fh:
            self.assertIn("Proj_teste", fh.read(), "o .cql deve levar o nome do projeto do config")

    def test_subarvores_sem_arvore_e_sem_declaracao_falha(self):
        self.preparar()
        self.montar(self.out)
        with self.assertRaises(FileNotFoundError):
            _executar("subarvores", "--config", self.backup, "--out", self.out, "--stem", "dataset_a",
                      "--arvore", os.path.join(self.out, "Trees", "tree_dataset_a_raxml.nexus"))

    def test_subarvores_de_arvore_que_falhou_e_no_op_declarado(self):
        self.preparar()
        self.montar(self.out)
        falha = os.path.join(self.out, "Trees", "tree_dataset_a_raxml.nexus")
        os.makedirs(os.path.join(self.out, "outputs", "prov"), exist_ok=True)
        with open(os.path.join(self.out, "outputs", "prov", "arvore-x-raxml.json"), "w") as fh:
            json.dump({"inference_methods": [{"metodo": "raxml", "estado": "tentado_e_falhou",
                                              "saida": falha, "motivo": "RuntimeError: x"}]}, fh)
        _executar("subarvores", "--config", self.backup, "--out", self.out, "--stem", "dataset_a",
                  "--arvore", falha)
        frag = self.prov("subarvores-raxml")
        self.assertIn("tentado_e_falhou", frag["no_op"])
        self.assertFalse(os.path.exists(os.path.join(self.out, "outputs", "subarvores", "tree_dataset_a_raxml.json")))


class TestProvenanciaEOrdem(unittest.TestCase):
    def test_higienizar_e_reidratar_sao_inversos_dentro_da_raiz(self):
        raiz = "/ws/Proj"
        frag = {"tools_invoked": {"mafft": {"runs": [{"command": ["/opt/env/bin/mafft", "--thread", "4",
                                                                 "/ws/Proj/out/entrada/x.fasta"],
                                                      "saida": "/ws/Proj/out/tmp/x.aln"}]}},
                "model_selection": [{"alinhamento": "/ws/Proj/out/tmp/x.aln", "estado": "concluida"}]}
        limpo = akoflow._higienizar(frag, raiz)
        cmd = limpo["tools_invoked"]["mafft"]["runs"][0]["command"]
        self.assertEqual(cmd, ["mafft", "--thread", "4", "out/entrada/x.fasta"])
        volta = akoflow._reidratar(limpo, raiz)
        self.assertEqual(volta["tools_invoked"]["mafft"]["runs"][0]["saida"], "/ws/Proj/out/tmp/x.aln")
        self.assertEqual(volta["model_selection"][0]["alinhamento"], "/ws/Proj/out/tmp/x.aln")

    def test_ordem_de_visita_poe_o_modelo_antes_do_primeiro_metodo_com_modelo(self):
        passos = [n for n, _ in akoflow.ordem_de_visita(["mafft"])]
        self.assertEqual(passos[:4], ["alinhar-mafft", "arvore-mafft-nj-distance",
                                      "arvore-mafft-nj-parsimony", "modelo-mafft"])
        self.assertEqual(passos[4], "arvore-mafft-iqtree")
        self.assertEqual(passos[-2:], ["arvore-mafft-upgma-distance", "arvore-mafft-upgma-parsimony"])

    def test_atividade_que_falha_grava_fragmento_com_o_motivo(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "P", "out")
            os.makedirs(out)
            with self.assertRaises(ValueError):
                with akoflow._atividade("teste", out, {}):
                    tool_runs.registrar("mafft", ["mafft", "x"], saida=os.path.join(out, "a.aln"))
                    raise ValueError("boom")
            with open(os.path.join(out, "outputs", "prov", "teste.json")) as fh:
                frag = json.load(fh)
        self.assertEqual(frag["desfecho"], {"motivo": "ValueError: boom", "status": "falhou"})
        self.assertEqual(frag["tools_invoked"]["mafft"]["runs"][0]["saida"], "out/a.aln")


class TestManifestoAbsorver(unittest.TestCase):
    def test_absorver_registros_equivale_a_drain(self):
        from workflow.utils.manifest import ExecutionManifest
        tool_runs.limpar()
        tool_runs.registrar("raxml-ng", ["raxml-ng", "--seed", "1"], saida="/p/out/a.nexus", seed=1)
        tool_runs.registrar_metodo("raxml", "executado", "/p/out/a.nexus", alinhador="mafft")
        tool_runs.registrar_selecao_modelo("/p/out/a.aln", "falhou", motivo="x")
        a = ExecutionManifest("/p", {})
        a.drain_tool_runs()
        b = ExecutionManifest("/p", {})
        b.absorver_registros(tool_runs.execucoes(), tool_runs.metodos(), tool_runs.selecoes_modelo())
        tool_runs.limpar()
        da, db = a.to_dict(), b.to_dict()
        for chave in ("tools_invoked", "inference_methods", "model_selection"):
            self.assertEqual(da[chave], db[chave])
        self.assertEqual(da["tools_invoked"]["raxml-ng"]["runs"][0]["saida"], "out/a.nexus")


class TestComparador(unittest.TestCase):
    def _out(self, d, nome, conteudo):
        out = os.path.join(d, nome, "out")
        os.makedirs(os.path.join(out, "Trees"))
        with open(os.path.join(out, "Trees", "tree_x_mafft_nj_distance.nexus"), "w") as fh:
            fh.write(conteudo)
        with open(os.path.join(out, "Trees", "tree_x_mafft_raxml.nexus"), "w") as fh:
            fh.write("r")
        return out

    def test_classes_h2_e_h3_separam_basicas_de_avancadas(self):
        with tempfile.TemporaryDirectory() as d:
            r = equivalencia.comparar(self._out(d, "a", "1"), self._out(d, "b", "2"))
        self.assertEqual((r["H2"]["iguais"], r["H2"]["diferentes"]), (0, 1))
        self.assertEqual((r["H3"]["iguais"], r["H3"]["diferentes"]), (1, 0))
        self.assertFalse(equivalencia.aprovado(r))

    def test_fpmax_ignora_ordem_de_linhas_e_do_itemset(self):
        with tempfile.TemporaryDirectory() as d:
            def csv(nome, linhas):
                p = os.path.join(d, nome)
                with open(p, "w") as fh:
                    fh.write(",itemsets,support,min_support_threshold\n" + "\n".join(linhas) + "\n")
                return p
            a = equivalencia._fpmax(csv("a.csv", ['0,"frozenset({3, 1})",0.5,0.1', '1,"frozenset({7})",1.0,0.2']))
            b = equivalencia._fpmax(csv("b.csv", ['5,"frozenset({7})",1.0,0.2', '9,"frozenset({1, 3})",0.5,0.1']))
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
