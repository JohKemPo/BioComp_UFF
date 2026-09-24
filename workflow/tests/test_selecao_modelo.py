"""
M7.3, opção C′ — seleção de modelo por BIC e tradução entre ferramentas.

Quatro camadas:

1. **Tabela de tradução**, entrada por entrada: as 24 combinações de
   `--mset mrbayes -mrate E,I,G,I+G` (3 `nst` × 2 frequências × 4 taxas),
   contra uma tabela **escrita à mão aqui** (não gerada pelo mesmo código).
2. **Os três binários reais** aceitam cada entrada e declaram o modelo
   esperado: `raxml-ng --parse`, `iqtree -te` (modelo e frequências que o
   `.iqtree` declara) e `mb` (`showmodel`). Pulados se a ferramenta faltar.
3. **Leitura do relatório** do ModelFinder, sobre trecho literal do IQ-TREE
   3.1.3 (VARV-6), incluindo a armadilha de nome (`GTR` = frequências iguais).
4. **Caminho do controlador**: uma seleção por alinhamento, o modelo chega aos
   três comandos, falha da seleção vira ``tentado_e_falhou`` nos métodos que
   dependiam dela (M7.6) e o FastTree segue; `model_selection: "nenhuma"`
   reproduz os literais antigos.

Executar com:

    python -m unittest workflow.tests.test_selecao_modelo
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from workflow.tree_construction.modelo_substituicao import (
    MODELO_LEGADO, TABELA_TRADUCAO, ModeloCanonico, ModeloIntraduzivel,
    SelecaoDeModeloFalhou, interpretar_nome_modelfinder, ler_relatorio_modelfinder,
    selecionar_modelo, traduzir)
from workflow.utils import tool_runs
from workflow.utils.external_tools import resolve_tool
from workflow.utils.manifest import ExecutionManifest

#: A tabela esperada, escrita à mão: nome como o ModelFinder o imprime →
#: (IQ-TREE `-m`, RAxML-NG `--model`, MrBayes `lset`, frequências iguais?).
#: Para as bases F81/HKY/GTR, o nome sem `+F` é o que o ModelFinder imprime
#: para a variante de frequências IGUAIS — por isso sai `JC`/`K2P`/`SYM`.
ESPERADO = {
    "JC":          ("JC",          "JC",       "lset nst=1 rates=equal", True),
    "JC+I":        ("JC+I",        "JC+I",     "lset nst=1 rates=propinv", True),
    "JC+G4":       ("JC+G4",       "JC+G",     "lset nst=1 rates=gamma ngammacat=4", True),
    "JC+I+G4":     ("JC+I+G4",     "JC+I+G",   "lset nst=1 rates=invgamma ngammacat=4", True),
    "F81+F":       ("F81+F",       "F81",      "lset nst=1 rates=equal", False),
    "F81+F+I":     ("F81+F+I",     "F81+I",    "lset nst=1 rates=propinv", False),
    "F81+F+G4":    ("F81+F+G4",    "F81+G",    "lset nst=1 rates=gamma ngammacat=4", False),
    "F81+F+I+G4":  ("F81+F+I+G4",  "F81+I+G",  "lset nst=1 rates=invgamma ngammacat=4", False),
    "K2P":         ("K2P",         "K80",      "lset nst=2 rates=equal", True),
    "K2P+I":       ("K2P+I",       "K80+I",    "lset nst=2 rates=propinv", True),
    "K2P+G4":      ("K2P+G4",      "K80+G",    "lset nst=2 rates=gamma ngammacat=4", True),
    "K2P+I+G4":    ("K2P+I+G4",    "K80+I+G",  "lset nst=2 rates=invgamma ngammacat=4", True),
    "HKY+F":       ("HKY+F",       "HKY",      "lset nst=2 rates=equal", False),
    "HKY+F+I":     ("HKY+F+I",     "HKY+I",    "lset nst=2 rates=propinv", False),
    "HKY+F+G4":    ("HKY+F+G4",    "HKY+G",    "lset nst=2 rates=gamma ngammacat=4", False),
    "HKY+F+I+G4":  ("HKY+F+I+G4",  "HKY+I+G",  "lset nst=2 rates=invgamma ngammacat=4", False),
    "SYM":         ("SYM",         "SYM",      "lset nst=6 rates=equal", True),
    "SYM+I":       ("SYM+I",       "SYM+I",    "lset nst=6 rates=propinv", True),
    "SYM+G4":      ("SYM+G4",      "SYM+G",    "lset nst=6 rates=gamma ngammacat=4", True),
    "SYM+I+G4":    ("SYM+I+G4",    "SYM+I+G",  "lset nst=6 rates=invgamma ngammacat=4", True),
    "GTR+F":       ("GTR+F",       "GTR",      "lset nst=6 rates=equal", False),
    "GTR+F+I":     ("GTR+F+I",     "GTR+I",    "lset nst=6 rates=propinv", False),
    "GTR+F+G4":    ("GTR+F+G4",    "GTR+G",    "lset nst=6 rates=gamma ngammacat=4", False),
    "GTR+F+I+G4":  ("GTR+F+I+G4",  "GTR+I+G",  "lset nst=6 rates=invgamma ngammacat=4", False),
}

#: Sinônimos que o ModelFinder também imprime (medido na lista de candidatos
#: do Zika-20 com `--mset mrbayes`): a variante sem `+F` das bases de
#: frequência livre é a de frequências iguais.
SINONIMOS = {"F81": "JC", "F81+I": "JC+I", "F81+G4": "JC+G4", "F81+I+G4": "JC+I+G4",
             "HKY": "K2P", "HKY+I": "K2P+I", "HKY+G4": "K2P+G4", "HKY+I+G4": "K2P+I+G4",
             "GTR": "SYM", "GTR+I": "SYM+I", "GTR+G4": "SYM+G4", "GTR+I+G4": "SYM+I+G4",
             "K80+G4": "K2P+G4"}

#: O que o RAxML-NG 2.0.2 declara (`--parse`, linha `Model:`) para cada nome.
RAXML_DECLARA = {
    "JC": "JC", "JC+I": "JC+I", "JC+G": "JC+G4m", "JC+I+G": "JC+I+G4m",
    "F81": "F81+FO", "F81+I": "F81+FO+I", "F81+G": "F81+FO+G4m", "F81+I+G": "F81+FO+I+G4m",
    "K80": "K80", "K80+I": "K80+I", "K80+G": "K80+G4m", "K80+I+G": "K80+I+G4m",
    "HKY": "HKY+FO", "HKY+I": "HKY+FO+I", "HKY+G": "HKY+FO+G4m", "HKY+I+G": "HKY+FO+I+G4m",
    "SYM": "SYM", "SYM+I": "SYM+I", "SYM+G": "SYM+G4m", "SYM+I+G": "SYM+I+G4m",
    "GTR": "GTR+FO", "GTR+I": "GTR+FO+I", "GTR+G": "GTR+FO+G4m", "GTR+I+G": "GTR+FO+I+G4m",
}

#: Trecho literal do `.iqtree` de `iqtree -m MF --mset mrbayes -mrate E,I,G,I+G`
#: no VARV-6 (IQ-TREE 3.1.3, 2026-09-23). O IQ-TREE repete `SYM` e a escolha
#: `GTR` é a de frequências IGUAIS (lnL = a de SYM, não a de GTR+F).
RELATORIO_VARV6 = """\
ModelFinder
-----------

Best-fit model according to BIC: GTR

List of models sorted by BIC scores:

Model                  LogL         AIC      w-AIC        AICc     w-AICc         BIC      w-BIC
GTR             -541924.103 1083876.206 +    0.208 1083876.208 +    0.208 1084022.244 +    0.342
SYM             -541924.145 1083876.291 +    0.199 1083876.292 +    0.199 1084022.329 +    0.327
SYM             -541924.145 1083876.291 +    0.199 1083876.292 +    0.199 1084022.329 +    0.327
GTR+G4          -541923.783 1083877.567 +    0.105 1083877.569 +    0.105 1084034.036 -  0.00094
F81+F+I+G4      -556460.671 1112949.342 -        0 1112949.343 -        0 1113095.380 -        0

AIC, w-AIC   : Akaike information criterion scores and weights.
"""


class TestTabela(unittest.TestCase):

    def test_cada_entrada_da_tabela(self):
        self.assertEqual(len(TABELA_TRADUCAO), 24)
        self.assertEqual(len(ESPERADO), 24)
        for nome, (iq, rx, lset, iguais) in ESPERADO.items():
            with self.subTest(modelo=nome):
                t = traduzir(nome)
                self.assertEqual(t["iqtree"], iq)
                self.assertEqual(t["raxml-ng"], rx)
                self.assertEqual(t["mrbayes_lset"], lset)
                self.assertEqual(t["mrbayes_prset"],
                                 "prset statefreqpr=fixed(equal)" if iguais else None)
                self.assertEqual(t["modelfinder"], nome)
        # As 24 entradas de ESPERADO cobrem as 24 da tabela, sem sobra.
        self.assertEqual({interpretar_nome_modelfinder(n) for n in ESPERADO},
                         set(TABELA_TRADUCAO))

    def test_sinonimos_do_modelfinder_sao_frequencias_iguais(self):
        """A armadilha medida no VARV-6: `GTR` impresso pelo ModelFinder é
        GTR de frequências iguais — `-m GTR` rodaria GTR+F."""
        for nome, canonico in SINONIMOS.items():
            with self.subTest(modelo=nome):
                self.assertEqual(interpretar_nome_modelfinder(nome),
                                 interpretar_nome_modelfinder(canonico))
                self.assertEqual(interpretar_nome_modelfinder(nome).frequencias, "iguais")
        self.assertEqual(traduzir("GTR")["iqtree"], "SYM")

    def test_fora_da_tabela_e_recusado(self):
        """Nunca passar adiante por aproximação: K3Pu não existe no RAxML-NG,
        TIM2 do RAxML-NG é outro modelo, FreeRate não existe no MrBayes."""
        for nome in ("TIM2+F+G4", "K3Pu+F+G4", "TN+F+I", "TVM", "GTR+F+R3",
                     "GTR+F+I+R2", "GTR+FO+G4", "GTR+FQ", "GTR+F+F", "GTR+F+G8",
                     "GTR+F+G4+ASC", "", "LG+G4"):
            with self.subTest(modelo=nome), self.assertRaises(ModeloIntraduzivel):
                traduzir(nome)

    def test_frequencias_declaradas_por_ferramenta(self):
        """Com +F, as três ferramentas NÃO fazem o mesmo (11-auditoria §3.6,
        ponto 1, em aberto); a divergência é declarada, não escondida."""
        f = traduzir("GTR+F+G4")["frequencias"]
        self.assertIn("empíricas", f["iqtree"])
        self.assertIn("+FO", f["raxml-ng"])
        self.assertIn("Dirichlet", f["mrbayes"])
        self.assertTrue(all("iguais" in v for v in traduzir("SYM")["frequencias"].values()))

    def test_legado_e_o_literal_de_antes(self):
        self.assertEqual(MODELO_LEGADO["iqtree"], "GTR+G")
        self.assertEqual(MODELO_LEGADO["raxml-ng"], "GTR+G")
        self.assertEqual(MODELO_LEGADO["mrbayes_lset"], "lset nst=6 rates=gamma")
        self.assertIsNone(MODELO_LEGADO["mrbayes_prset"])
        self.assertNotIn("modelfinder", MODELO_LEGADO)

    def test_referencias_traduzem_para_o_mesmo_comando_do_raxml(self):
        """VARV-49 e Zika-20 escolhem GTR+F+G4: o RAxML-NG recebe exatamente
        o `--model GTR+G` de antes (Δ = 0 por construção no comando)."""
        self.assertEqual(traduzir("GTR+F+G4")["raxml-ng"], MODELO_LEGADO["raxml-ng"])


class TestRelatorio(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_le_escolha_e_candidatos_sem_repeticao(self):
        caminho = os.path.join(self.dir, "m.iqtree")
        with open(caminho, "w") as fh:
            fh.write(RELATORIO_VARV6)
        r = ler_relatorio_modelfinder(caminho)
        self.assertEqual(r["modelo_escolhido"], "GTR")
        self.assertEqual([c["modelo"] for c in r["candidatos"]],
                         ["GTR", "SYM", "GTR+G4", "F81+F+I+G4"])
        self.assertEqual(r["candidatos"][0]["BIC"], 1084022.244)
        self.assertEqual(r["candidatos"][0]["lnL"], -541924.103)
        self.assertEqual(r["candidatos"][3]["w_BIC"], 0.0)

    def test_sem_relatorio_ou_sem_linha_e_none(self):
        self.assertIsNone(ler_relatorio_modelfinder(os.path.join(self.dir, "x"))["modelo_escolhido"])
        caminho = os.path.join(self.dir, "vazio.iqtree")
        open(caminho, "w").close()
        self.assertIsNone(ler_relatorio_modelfinder(caminho)["modelo_escolhido"])


def _aln(ntax=6):
    from Bio import AlignIO
    return AlignIO.read(os.path.join("projects", "Zika_21seq_validacao", "out", "Align",
                                     "dataset_final_mafft.aln"), "fasta")[:ntax]


class _ComDiretorio(unittest.TestCase):
    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()


@unittest.skipUnless(resolve_tool("raxml-ng"), "RAxML-NG ausente")
class TestRaxmlAceitaCadaEntrada(_ComDiretorio):

    def test_parse_declara_o_modelo_esperado(self):
        from Bio import AlignIO
        msa = os.path.join(self.dir, "a.phylip")
        AlignIO.write(_aln(), msa, "phylip")
        for nome, (_, rx, _, _) in ESPERADO.items():
            with self.subTest(modelo=nome, raxml=rx):
                saida = subprocess.run(
                    [resolve_tool("raxml-ng"), "--parse", "--msa", msa, "--model", rx,
                     "--prefix", os.path.join(self.dir, rx.replace("+", "_")),
                     "--threads", "1", "--redo"],
                    capture_output=True, text=True)
                declarado = re.search(r"^\s*Model:\s*(\S+)", saida.stdout, re.MULTILINE)
                self.assertIsNotNone(declarado, saida.stdout[-400:])
                self.assertEqual(declarado.group(1), RAXML_DECLARA[rx])


@unittest.skipUnless(resolve_tool("iqtree"), "IQ-TREE ausente")
class TestIqtreeAceitaCadaEntrada(_ComDiretorio):

    def test_modelo_e_frequencias_declarados(self):
        from Bio import AlignIO
        aln = _aln()
        msa = os.path.join(self.dir, "a.phylip")
        AlignIO.write(aln, msa, "phylip")
        nomes = [r.id[:10] for r in aln]
        arvore = os.path.join(self.dir, "t.nwk")
        with open(arvore, "w") as fh:
            fh.write(f"(({nomes[0]},{nomes[1]}),({nomes[2]},{nomes[3]}),({nomes[4]},{nomes[5]}));\n")
        for nome, (iq, _, _, iguais) in ESPERADO.items():
            with self.subTest(modelo=nome, iqtree=iq):
                prefixo = os.path.join(self.dir, iq.replace("+", "_"))
                subprocess.run([resolve_tool("iqtree"), "-s", msa, "-m", iq, "-te", arvore,
                                "-nt", "1", "-seed", "1", "-pre", prefixo, "-redo", "-quiet"],
                               capture_output=True, text=True, check=True)
                with open(prefixo + ".iqtree") as fh:
                    texto = fh.read()
                self.assertIn(f"Model of substitution: {iq}\n", texto)
                self.assertIn("(equal frequencies)" if iguais
                              else "(empirical counts from alignment)", texto)


@unittest.skipUnless(resolve_tool("mrbayes"), "MrBayes (mb) ausente")
class TestMrBayesAceitaCadaEntrada(_ComDiretorio):

    TAXA = {"rates=equal": "Equal", "rates=propinv": "Propinv",
            "rates=gamma": "Gamma", "rates=invgamma": "Invgamma"}

    def test_showmodel_declara_nst_taxa_e_frequencias(self):
        from Bio import AlignIO
        aln = _aln()
        for r in aln:
            r.annotations["molecule_type"] = "DNA"
        AlignIO.write(aln, os.path.join(self.dir, "a.nexus"), "nexus")
        for nome, (_, _, lset, iguais) in ESPERADO.items():
            with self.subTest(modelo=nome, lset=lset):
                t = traduzir(nome)
                script = "set autoclose=yes nowarn=yes\nexecute a.nexus\n" + lset + "\n"
                if t["mrbayes_prset"]:
                    script += t["mrbayes_prset"] + "\n"
                script += "showmodel\nquit\n"
                saida = subprocess.run([resolve_tool("mrbayes")], input=script,
                                       capture_output=True, text=True, cwd=self.dir).stdout
                self.assertNotIn("Error", saida)
                self.assertRegex(saida, rf"Nst\s+=\s+{t['nst']}\b")
                taxa = self.TAXA[lset.split()[2]]
                self.assertRegex(saida, rf"Rates\s+=\s+{taxa}\b")
                self.assertEqual("fixed to be equal" in saida, iguais)
                if "gamma" in lset:
                    self.assertIn("approximated using 4 categories", saida)


@unittest.skipUnless(resolve_tool("iqtree"), "IQ-TREE ausente")
class TestSelecaoDeVerdade(_ComDiretorio):

    def test_escolhe_o_menor_bic_e_traduz(self):
        r = selecionar_modelo(_aln(), os.path.join(self.dir, "mf"), 12345)
        self.assertEqual(r["modelo_escolhido"], r["candidatos"][0]["modelo"])
        self.assertEqual(r["candidatos"][0]["BIC"], min(c["BIC"] for c in r["candidatos"]))
        self.assertEqual(r["traducao"]["modelfinder"], r["modelo_escolhido"])
        # Nenhum candidato fora da tabela: `-mrate E,I,G,I+G` tirou o FreeRate.
        for c in r["candidatos"]:
            interpretar_nome_modelfinder(c["modelo"])
        chamada = tool_runs.execucoes()["iqtree"]["runs"][0]
        self.assertIn("-mrate", chamada["command"])
        self.assertEqual(chamada["seed"], 12345)


# --------------------------------------------------------------------------- #
# Caminho do controlador
# --------------------------------------------------------------------------- #

NEWICK = "((t1:0.1,t2:0.1):0.1,(t3:0.1,t4:0.1):0.1);"
ALINHAMENTO = ">t1\nACGTACGTAA\n>t2\nACGTACGTAG\n>t3\nACGAACGTTA\n>t4\nACGAACGTTG\n"
CON_TRE = ("#NEXUS\nbegin trees;\n  translate\n    1 t1,\n    2 t2,\n    3 t3,\n"
           "    4 t4;\n  tree con_50_majrule = [&U] ((1[&prob=1]:0.1,2:0.1)"
           "[&prob=0.9]:0.1,3:0.1,4:0.1);\nend;\n")
SUMT_OK = "       Average standard deviation of split frequencies = 0.004\n"


class TestCaminhoDoControlador(_ComDiretorio):

    def setUp(self):
        super().setUp()
        self.entrada = os.path.join(self.dir, "entrada")
        self.saida = os.path.join(self.dir, "out")
        os.makedirs(self.entrada)
        with open(os.path.join(self.entrada, "ds.fasta"), "w") as fh:
            fh.write(ALINHAMENTO)
        self.chamadas = []
        self.scripts_mb = []

    def _controlador(self, **extra):
        from workflow.controller.treeBuilderController import TreeBuilderController
        config = dict(input_path=self.entrada, output_path=self.saida, num_threads=1,
                      mode="advanced", output_format="nexus", aligners=["mafft"],
                      align_method="mafft", ignore_mode=["distance", "parsimony"],
                      on_method_failure="continue", mrbayes_ngen=1000)
        config.update(extra)
        c = TreeBuilderController(**config)
        with open(os.path.join(self.saida, "tmp", "ds_mafft.aln"), "w") as fh:
            fh.write(ALINHAMENTO)
        return c

    def _executar(self, modelo_escolhido="GTR+F+G4", falha_mf=None):
        def run(cmd, *a, **k):
            self.chamadas.append(list(cmd))
            if "MF" in cmd:
                if falha_mf == "timeout":
                    raise subprocess.TimeoutExpired(cmd, 5)
                if falha_mf == "codigo":
                    raise subprocess.CalledProcessError(2, cmd, "", "ERROR: algo")
                if falha_mf != "sem_linha":
                    with open(cmd[cmd.index("-pre") + 1] + ".iqtree", "w") as fh:
                        fh.write(f"Best-fit model according to BIC: {modelo_escolhido}\n")
                return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
            if k.get("stdin") is not None:  # MrBayes
                self.scripts_mb.append(k["stdin"].read())
                with open(os.path.join(k["cwd"], "alignment.nexus.con.tre"), "w") as fh:
                    fh.write(CON_TRE)
                return subprocess.CompletedProcess(cmd, 0, stdout=SUMT_OK.encode(), stderr=b"")
            for flag, sufixo in (("-pre", ".treefile"), ("--prefix", ".raxml.support")):
                if flag in cmd:
                    with open(cmd[cmd.index(flag) + 1] + sufixo, "w") as fh:
                        fh.write(NEWICK + "\n")
            return subprocess.CompletedProcess(cmd, 0, stdout=NEWICK + "\n", stderr="")

        def require(nome, *a, **k):
            if falha_mf == "binario" and nome == "iqtree":
                raise FileNotFoundError("iqtree não encontrado. Procurado em "
                                        "/home/fulano/miniconda3/envs/x/bin/iqtree")
            return f"/opt/env/bin/{nome}"

        controlador = self._controlador(**getattr(self, "extra", {}))
        caminhos = controlador._prepare_output_paths("ds")
        with mock.patch("workflow.tree_construction.builder.require_tool", side_effect=require), \
             mock.patch("workflow.tree_construction.modelo_substituicao.require_tool",
                        side_effect=require), \
             mock.patch("subprocess.run", side_effect=run), \
             mock.patch.object(controlador, "save_tree_image"):
            controlador._process_advanced_mode("ds", os.path.join(self.entrada, "ds.fasta"),
                                               caminhos)

    def _mf(self):
        return [c for c in self.chamadas if "MF" in c]

    def _cmd(self, binario):
        return next(c for c in self.chamadas if c[0].endswith(binario) and "MF" not in c)

    @staticmethod
    def _estados():
        return {m["metodo"]: m for m in tool_runs.metodos()}

    def test_uma_selecao_por_alinhamento_e_modelo_nos_tres_comandos(self):
        self._executar("HKY+F+I")
        self.assertEqual(len(self._mf()), 1, "uma seleção para IQ-TREE, RAxML-NG e MrBayes")
        self.assertEqual(self._cmd("iqtree")[self._cmd("iqtree").index("-m") + 1], "HKY+F+I")
        raxml = self._cmd("raxml-ng")
        self.assertEqual(raxml[raxml.index("--model") + 1], "HKY+I")
        self.assertIn("lset nst=2 rates=propinv\n", self.scripts_mb[0])
        self.assertNotIn("prset", self.scripts_mb[0])
        # O FastTree não segue modelo escolhido: continua GTR+CAT20.
        self.assertNotIn("-m", self._cmd("fasttree"))
        self.assertIn("-gtr", self._cmd("fasttree"))
        for metodo in ("iqtree", "fasttree", "raxml", "mrbayes"):
            self.assertEqual(self._estados()[metodo]["estado"], "executado")
        (selecao,) = tool_runs.selecoes_modelo()
        self.assertEqual(selecao["estado"], "concluida")
        self.assertEqual(selecao["modelo_escolhido"], "HKY+F+I")
        self.assertEqual(selecao["traducao"]["mrbayes_lset"], "lset nst=2 rates=propinv")

    def test_frequencias_iguais_levam_prset_ao_mrbayes(self):
        self._executar("GTR+I")          # como o ModelFinder imprime: = SYM+I
        self.assertEqual(self._cmd("iqtree")[self._cmd("iqtree").index("-m") + 1], "SYM+I")
        self.assertIn("lset nst=6 rates=propinv\nprset statefreqpr=fixed(equal)\n",
                      self.scripts_mb[0])

    def test_modelo_por_chamada_no_manifesto(self):
        self._executar("GTR+F+G4")
        entrada = tool_runs.execucoes()["iqtree"]
        self.assertIn("runs[i].modelo", entrada["model"])
        inferencia = [r for r in entrada["runs"] if "etapa" not in r][0]
        self.assertEqual(inferencia["modelo"], "GTR+F+G4")
        self.assertIn("GTR+F+G4", inferencia["modelo_origem"])
        selecao = [r for r in entrada["runs"] if "etapa" in r][0]
        self.assertEqual(selecao["seed"], 12345)
        manifesto = ExecutionManifest(project_root=self.dir, params={})
        manifesto.drain_tool_runs()
        (registro,) = manifesto.to_dict()["model_selection"]
        self.assertFalse(os.path.isabs(registro["alinhamento"]))
        self.assertFalse(os.path.isabs(registro["relatorio"]))
        self.assertNotIn(self.dir, repr(manifesto.to_dict()["model_selection"]))

    def test_politica_nenhuma_reproduz_os_literais(self):
        self.extra = {"model_selection": "nenhuma"}
        self._executar()
        self.assertEqual(self._mf(), [])
        self.assertEqual(self._cmd("iqtree")[self._cmd("iqtree").index("-m") + 1], "GTR+G")
        raxml = self._cmd("raxml-ng")
        self.assertEqual(raxml[raxml.index("--model") + 1], "GTR+G")
        self.assertIn("lset nst=6 rates=gamma\nmcmc", self.scripts_mb[0])
        self.assertEqual(tool_runs.execucoes()["iqtree"]["model"], "GTR+G")
        (selecao,) = tool_runs.selecoes_modelo()
        self.assertEqual(selecao["estado"], "desligada_por_configuracao")

    def test_politica_desconhecida_e_recusada_no_inicio(self):
        with self.assertRaises(ValueError):
            self._controlador(model_selection="MFP")

    def _conferir_falha(self, trecho):
        estados = self._estados()
        for metodo in ("iqtree", "raxml", "mrbayes"):
            with self.subTest(metodo=metodo):
                self.assertEqual(estados[metodo]["estado"], "tentado_e_falhou")
                self.assertIn("SelecaoDeModeloFalhou", estados[metodo]["motivo"])
        self.assertEqual(estados["fasttree"]["estado"], "executado")
        self.assertLessEqual(len(self._mf()), 1, "a seleção não é tentada de novo")
        (selecao,) = tool_runs.selecoes_modelo()
        self.assertEqual(selecao["estado"], "falhou")
        self.assertIn(trecho, selecao["motivo"])

    def test_binario_ausente_nao_derruba_o_fasttree(self):
        self._executar(falha_mf="binario")
        self._conferir_falha("FileNotFoundError")
        manifesto = ExecutionManifest(project_root=self.dir, params={})
        manifesto.drain_tool_runs()
        self.assertNotIn("fulano", repr(manifesto.to_dict()))

    def test_tempo_esgotado(self):
        self._executar(falha_mf="timeout")
        self._conferir_falha("excedeu")

    def test_codigo_de_saida(self):
        self._executar(falha_mf="codigo")
        self._conferir_falha("código 2")

    def test_relatorio_sem_modelo(self):
        self._executar(falha_mf="sem_linha")
        self._conferir_falha("Best-fit")

    def test_modelo_sem_traducao(self):
        """Se um dia o ModelFinder devolver algo fora do conjunto (versão nova
        do IQ-TREE mudando `--mset mrbayes`), recusa — não roda GTR+G."""
        self._executar("TIM2+F+G4")
        self._conferir_falha("TIM2")

    def test_politica_fail_para_no_primeiro(self):
        self.extra = {"on_method_failure": "fail"}
        with self.assertRaises(SelecaoDeModeloFalhou):
            self._executar(falha_mf="binario")
        self.assertEqual(self._estados()["iqtree"]["estado"], "tentado_e_falhou")
        self.assertNotIn("raxml", self._estados())


if __name__ == "__main__":
    unittest.main()
