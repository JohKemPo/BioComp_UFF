"""
Testes do manifesto de execução — M2.5 / D11.

Uma figura só é reproduzível se for possível dizer qual código, qual entrada e
qual ambiente a produziram. O manifesto é o objeto que carrega esses três fatos,
e estes testes fixam as três propriedades que o tornam utilizável: ele registra
o que precisa, **não** registra o que não pode (nome de usuário, caminho
absoluto — o vazamento de D15), e sobrevive a uma execução que falha.

Executar com:

    python -m unittest workflow.tests.test_manifest
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest

from workflow.utils.manifest import (ExecutionManifest, MANIFEST_FILENAME,
                                     file_digest, tool_versions)


class TestDigest(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_digest_conhecido(self):
        caminho = os.path.join(self.dir, "entrada.fasta")
        with open(caminho, "w") as handle:
            handle.write(">a\nACGT\n")
        import hashlib
        esperado = hashlib.sha256(b">a\nACGT\n").hexdigest()
        self.assertEqual(file_digest(caminho), esperado)

    def test_arquivo_ausente_e_none(self):
        """Ausente é `None`, nunca string vazia: um arquivo que não existe é um
        fato e o manifesto tem de dizê-lo."""
        self.assertIsNone(file_digest(os.path.join(self.dir, "nao_existe")))

    def test_conteudo_diferente_muda_o_digest(self):
        a = os.path.join(self.dir, "a"); b = os.path.join(self.dir, "b")
        with open(a, "w") as handle:
            handle.write("ACGT")
        with open(b, "w") as handle:
            handle.write("ACGA")
        self.assertNotEqual(file_digest(a), file_digest(b))


class TestVersoesDeFerramenta(unittest.TestCase):

    def test_devolve_uma_entrada_por_ferramenta(self):
        versoes = tool_versions()
        for esperada in ("mafft", "iqtree2", "raxml-ng", "FastTree", "mrbayes"):
            self.assertIn(esperada, versoes)

    def test_ferramenta_ausente_e_none(self):
        versoes = tool_versions()
        for nome, versao in versoes.items():
            self.assertTrue(versao is None or isinstance(versao, str), nome)


class TestManifesto(unittest.TestCase):

    def setUp(self):
        self.raiz = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.raiz, "out", "outputs"))
        os.makedirs(os.path.join(self.raiz, "out", "Trees"))
        os.makedirs(os.path.join(self.raiz, "out", "tmp", "raxml_x"))
        self.entrada = os.path.join(self.raiz, "dataset.fasta")
        with open(self.entrada, "w") as handle:
            handle.write(">a\nACGT\n>b\nACGA\n")
        self.params = {"project_name": "teste", "tree_config": {"input_path": self.entrada}}

    def tearDown(self):
        shutil.rmtree(self.raiz, ignore_errors=True)

    def _manifesto(self):
        return ExecutionManifest(self.raiz, self.params)

    def _lido(self, manifesto):
        with open(manifesto.path(), encoding="utf-8") as handle:
            return json.load(handle)

    def test_grava_no_lugar_esperado(self):
        manifesto = self._manifesto()
        manifesto.write()
        self.assertTrue(os.path.isfile(
            os.path.join(self.raiz, "out", "outputs", MANIFEST_FILENAME)))

    def test_registra_entrada_com_digest(self):
        manifesto = self._manifesto()
        manifesto.register_input(self.entrada)
        dados = self._lido(manifesto) if manifesto.write() else None
        self.assertEqual(list(dados["inputs_sha256"]), ["dataset.fasta"])
        self.assertEqual(dados["inputs_sha256"]["dataset.fasta"], file_digest(self.entrada))

    def test_registra_entrada_que_e_diretorio(self):
        """O `input_path` do `tree_config` é um diretório de FASTA. Registrar só
        o caminho dele não diz nada sobre o conteúdo — e é o conteúdo que
        precisa ser idêntico para a execução ser a mesma."""
        pasta = os.path.join(self.raiz, "entradas")
        os.makedirs(pasta)
        for nome in ("a.fasta", "b.fasta"):
            with open(os.path.join(pasta, nome), "w") as handle:
                handle.write(f">{nome}\nACGT\n")
        with open(os.path.join(pasta, "leiame.txt"), "w") as handle:
            handle.write("ignorado")

        manifesto = self._manifesto()
        manifesto.register_input(pasta)
        registradas = manifesto.to_dict()["inputs_sha256"]
        self.assertEqual(sorted(registradas),
                         [os.path.join("entradas", "a.fasta"),
                          os.path.join("entradas", "b.fasta")])
        self.assertTrue(all(v for v in registradas.values()))

    def test_nao_grava_caminho_absoluto(self):
        """D15: o manifesto vai para o repositório e pode ir para o material
        suplementar. Caminho absoluto expõe estrutura de diretórios e nome de
        usuário de quem rodou."""
        manifesto = self._manifesto()
        manifesto.register_input(self.entrada)
        arvore = os.path.join(self.raiz, "out", "Trees", "tree_x.nexus")
        with open(arvore, "w") as handle:
            handle.write("#NEXUS\n")
        manifesto.register_outputs(os.path.join(self.raiz, "out"))
        manifesto.finish()

        with open(manifesto.path(), encoding="utf-8") as handle:
            texto = handle.read()
        for chave in manifesto.to_dict()["inputs_sha256"]:
            self.assertFalse(os.path.isabs(chave), chave)
        for chave in manifesto.to_dict()["outputs_sha256"]:
            self.assertFalse(os.path.isabs(chave), chave)
        # `params` carrega o caminho de entrada como o usuário o passou; o que
        # não pode aparecer é caminho absoluto FABRICADO pelo manifesto.
        dados = json.loads(texto)
        self.assertNotIn("hostname", dados["environment"])
        self.assertNotIn("user", dados["environment"])

    def test_ignora_o_diretorio_tmp(self):
        """`out/tmp` guarda os arquivos intermediários das ferramentas, que são
        gigabytes e não são resultado."""
        manifesto = self._manifesto()
        with open(os.path.join(self.raiz, "out", "Trees", "tree_x.nexus"), "w") as handle:
            handle.write("#NEXUS\n")
        with open(os.path.join(self.raiz, "out", "tmp", "raxml_x", "lixo.nexus"), "w") as handle:
            handle.write("x")
        manifesto.register_outputs(os.path.join(self.raiz, "out"))
        chaves = manifesto.to_dict()["outputs_sha256"]
        self.assertIn(os.path.join("out", "Trees", "tree_x.nexus"), chaves)
        self.assertTrue(all("tmp" not in k for k in chaves), chaves)

    def test_registra_linha_de_comando_efetiva(self):
        """Depois de D17 sabe-se que a paralelização muda a topologia mesmo com
        a semente fixa — então a linha de comando **é** metadado do resultado."""
        manifesto = self._manifesto()
        manifesto.register_tool_run(
            "raxml-ng",
            ["raxml-ng", "--msa", "x.phylip", "--threads", "4", "--workers", "1",
             "--seed", "12345"],
            seed=12345, threads=4, workers=1)
        invocadas = manifesto.to_dict()["tools_invoked"]
        self.assertEqual(invocadas["raxml-ng"]["seed"], 12345)
        self.assertEqual(invocadas["raxml-ng"]["workers"], 1)
        self.assertIn("--workers", invocadas["raxml-ng"]["command"])

    def test_reproducibilidade_vem_da_mesma_fonte_do_builder(self):
        """O manifesto não pode declarar uma semente e o pipeline usar outra:
        os dois leem de `builder.reproducibility_settings`."""
        from workflow.tree_construction.builder import (TreeBuilder,
                                                        reproducibility_settings)
        config = {"random_seed": 777, "raxml_threads": 2}
        manifesto = self._manifesto()
        manifesto.register_reproducibility(reproducibility_settings(config))
        declarado = manifesto.to_dict()["reproducibility"]

        usado = TreeBuilder(**config)
        self.assertEqual(declarado["random_seed"], usado.random_seed)
        self.assertEqual(declarado["raxml_threads"], usado.raxml_threads)
        self.assertEqual(declarado["iqtree_threads"], usado.iqtree_threads)

    def test_reproducibilidade_tem_padrao_explicito(self):
        from workflow.tree_construction.builder import reproducibility_settings
        manifesto = self._manifesto()
        manifesto.register_reproducibility(reproducibility_settings({}))
        declarado = manifesto.to_dict()["reproducibility"]
        self.assertEqual(declarado["random_seed"], 12345)
        self.assertGreaterEqual(declarado["raxml_threads"], 1)

    def test_horario_de_termino_so_existe_depois_de_finish(self):
        manifesto = self._manifesto()
        manifesto.write()
        self.assertIsNone(self._lido(manifesto)["finished_at_utc"])
        manifesto.finish()
        self.assertIsNotNone(self._lido(manifesto)["finished_at_utc"])

    def test_manifesto_parcial_ja_e_valido(self):
        """Se a execução morrer no meio — e D17 mostra que morre —, o manifesto
        parcial ainda diz em que ambiente ela morreu."""
        manifesto = self._manifesto()
        manifesto.write()
        dados = self._lido(manifesto)
        self.assertIsNotNone(dados["run_id"])
        self.assertIsNotNone(dados["started_at_utc"])
        self.assertIn("environment", dados)
        self.assertIn("tools_available", dados)

    def test_registra_o_commit_dos_dois_repositorios(self):
        manifesto = ExecutionManifest(self.raiz, self.params,
                                      repos={"BioComp_UFF": ".", "PhyloTreeMiner": ".."})
        estado = manifesto.to_dict()["git"]
        self.assertIn("BioComp_UFF", estado)
        self.assertIn("PhyloTreeMiner", estado)
        for repo in estado.values():
            self.assertIn("commit", repo)
            self.assertIn("dirty", repo)

    def test_repositorio_inexistente_nao_quebra(self):
        manifesto = ExecutionManifest(self.raiz, self.params,
                                      repos={"fantasma": "/nao/existe"})
        self.assertIsNone(manifesto.to_dict()["git"]["fantasma"]["commit"])

    def test_serializa_em_json_estavel(self):
        manifesto = self._manifesto()
        manifesto.register_input(self.entrada)
        manifesto.finish()
        with open(manifesto.path(), encoding="utf-8") as handle:
            texto = handle.read()
        json.loads(texto)   # não levanta
        self.assertIn('"manifest_version": 1', texto)


if __name__ == "__main__":
    unittest.main()
