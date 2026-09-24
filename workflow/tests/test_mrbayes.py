"""
M7.4 — MrBayes correto ([D20](../../../docs/science/02-defeitos-que-alteram-resultado.md#d20)).

Três camadas:

1. **Configuração** (`mrbayes_settings`): padrões iguais ao que o pipeline já
   executava de fato; domínio validado antes de a cadeia rodar; `nruns=1` com
   o gate ligado recusado na largada, porque o ASDSF não existe com uma
   corrida só.
2. **Gate** (`ler_diagnosticos_mrbayes` + `avaliar_convergencia_mrbayes`),
   sobre texto literal da saída do MrBayes 3.2.7.
3. **Caminho real**, com o binário `mb` de verdade (pulado se ausente): mesma
   semente, mesma árvore; cadeia curta demais recusada e não gravada em
   `Trees/`; cadeia que converge aceita.

Executar com:

    python -m unittest workflow.tests.test_mrbayes
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from workflow.tree_construction.builder import (
    MRBAYES_DEFAULTS, ConvergenciaNaoAtingida, ProbabilidadePosteriorInvalida,
    TreeBuilder, avaliar_convergencia_mrbayes, ler_diagnosticos_mrbayes,
    mrbayes_settings)
from workflow.tree_construction.modelo_substituicao import traduzir
from workflow.utils import tool_runs
from workflow.utils.external_tools import resolve_tool

#: Trecho literal do `sumt` do MrBayes 3.2.7 (caracterização de 2026-09-23,
#: 8 táxons de Zika, ngen=20000). A linha com `:` é do progresso do `mcmc` e
#: NÃO é a que o gate lê.
STDOUT_SUMT = """\
   Average standard deviation of split frequencies: 0.052745
   Using relative burnin ('relburnin=yes'), discarding the first 25 % of sampled trees
       Average standard deviation of split frequencies = 0.045125
       Maximum standard deviation of split frequencies = 0.107705
       Average PSRF for parameter values (excluding NA and >10.0) = 1.003
       Maximum PSRF for parameter values = 1.034
"""

PSTAT = ("[ID: 6154880905]\n"
         "Parameter\tMean\tVariance\tLower\tUpper\tMedian\tminESS\tavgESS\tPSRF\n"
         "TL\t3.5e-01\t1.4e-03\t2.9e-01\t4.3e-01\t3.4e-01\t3.148599e+01\t4.38e+01\t1.05\n"
         "r(A<->C)\t3.1e-02\t4.0e-05\t2.0e-02\t4.3e-02\t3.0e-02\t1.590616e+01\t3.9e+01\t1.05\n")


class TestConfiguracao(unittest.TestCase):

    def test_padroes_reproduzem_o_que_ja_se_executava(self):
        """ngen/samplefreq/printfreq eram literais do script; nruns/nchains, o
        padrão da ferramenta; 25 %, o burn-in que o MrBayes aplicava de fato
        (o `burnin=250` do script era ignorado — relburnin=yes)."""
        cfg = mrbayes_settings({})
        self.assertEqual(cfg['mrbayes_ngen'], 1_000_000)
        self.assertEqual(cfg['mrbayes_samplefreq'], 100)
        self.assertEqual(cfg['mrbayes_printfreq'], 1000)
        self.assertEqual(cfg['mrbayes_nruns'], 2)
        self.assertEqual(cfg['mrbayes_nchains'], 4)
        self.assertEqual(cfg['mrbayes_burninfrac'], 0.25)
        self.assertEqual(cfg['mrbayes_timeout_s'], 3600)
        self.assertEqual(cfg['mrbayes_asdsf_max'], 0.01)

    def test_nruns_1_com_gate_ligado_e_recusado_na_largada(self):
        with self.assertRaisesRegex(ValueError, "incomputável"):
            mrbayes_settings({'mrbayes_nruns': 1})

    def test_nruns_1_com_gate_desligado_explicitamente_e_aceito(self):
        cfg = mrbayes_settings({'mrbayes_nruns': 1, 'mrbayes_asdsf_max': None})
        self.assertIsNone(cfg['mrbayes_asdsf_max'])

    def test_ausente_nao_e_o_mesmo_que_null(self):
        self.assertEqual(mrbayes_settings({})['mrbayes_asdsf_max'], 0.01)
        self.assertIsNone(mrbayes_settings({'mrbayes_asdsf_max': None})['mrbayes_asdsf_max'])

    def test_dominio(self):
        for ruim in ({'mrbayes_burninfrac': 1.0}, {'mrbayes_burninfrac': -0.1},
                     {'mrbayes_ngen': 0}, {'mrbayes_nchains': 0},
                     {'mrbayes_asdsf_max': 0}, {'mrbayes_asdsf_max': -1},
                     # 11 amostras por corrida, 99 % de burn-in: nada sobra.
                     {'mrbayes_ngen': 1000, 'mrbayes_burninfrac': 0.99}):
            with self.subTest(ruim=ruim), self.assertRaises(ValueError):
                mrbayes_settings(ruim)

    def test_script_declara_tudo_e_nao_tem_burnin_absoluto(self):
        script = TreeBuilder.mrbayes_script(mrbayes_settings({}), 777, 'alignment.nexus')
        self.assertIn("seed=777 swapseed=777", script)
        self.assertIn("nruns=2 nchains=4", script)
        self.assertIn("lset nst=6 rates=gamma", script)
        self.assertEqual(script.count("burninfrac=0.25"), 3)  # mcmc, sump, sumt
        self.assertNotIn("burnin=250", script)

    def test_controlador_repassa_mrbayes_e_preserva_null(self):
        """D26 de novo: o builder sabe ler, o controlador precisa repassar —
        e `null` (gate desligado) não pode ser filtrado como ausente."""
        from workflow.controller.treeBuilderController import TreeBuilderController
        c = TreeBuilderController.__new__(TreeBuilderController)
        c.mrbayes_ngen = 5000
        c.mrbayes_asdsf_max = None
        self.assertEqual(c._mrbayes_kwargs(),
                         {'mrbayes_ngen': 5000, 'mrbayes_asdsf_max': None})
        self.assertEqual(set(MRBAYES_DEFAULTS) >= set(c._mrbayes_kwargs()), True)


class TestGate(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_le_a_linha_do_sumt_nao_a_do_mcmc(self):
        pstat = os.path.join(self.dir, "x.pstat")
        with open(pstat, "w") as fh:
            fh.write(PSTAT)
        d = ler_diagnosticos_mrbayes(STDOUT_SUMT, pstat)
        self.assertEqual(d["asdsf"], 0.045125)      # não 0.052745
        self.assertEqual(d["sdsf_max"], 0.107705)
        self.assertEqual(d["psrf_max"], 1.034)
        self.assertAlmostEqual(d["ess_min"], 15.90616)

    def test_sem_linha_e_none_nao_zero(self):
        """nruns=1: o MrBayes termina com código 0 e não imprime o ASDSF.
        Zero seria "convergiu perfeitamente" — o oposto do fato (regra 5)."""
        d = ler_diagnosticos_mrbayes("   Using relative burnin ...\n", None)
        self.assertIsNone(d["asdsf"])
        self.assertIsNone(d["ess_min"])

    def test_nan_nao_vira_numero(self):
        d = ler_diagnosticos_mrbayes(
            "Average standard deviation of split frequencies = NA\n")
        self.assertIsNone(d["asdsf"])

    def test_decisoes(self):
        casos = [
            ({"asdsf": 0.004}, 0.01, (True, "atingida")),
            ({"asdsf": 0.045}, 0.01, (False, "nao_atingida")),
            ({"asdsf": 0.01}, 0.01, (False, "nao_atingida")),   # empate recusa
            ({"asdsf": None}, 0.01, (False, "nao_atingida")),   # não lido recusa
            ({"asdsf": 0.2}, None, (True, "nao_verificada")),   # gate desligado
            ({"asdsf": None}, None, (True, "nao_verificada")),
        ]
        for diag, limite, (aceita, estado) in casos:
            with self.subTest(diag=diag, limite=limite):
                a, e, motivo = avaliar_convergencia_mrbayes(diag, limite)
                self.assertEqual((a, e), (aceita, estado))
                self.assertTrue(motivo)


def _aln_zika(ntax):
    from Bio import AlignIO
    caminho = os.path.join("projects", "Zika_21seq_validacao", "out", "Align",
                           "dataset_final_mafft.aln")
    return AlignIO.read(caminho, "fasta")[:ntax]


class TestConstrutorComFerramentaSubstituida(unittest.TestCase):
    """Caminho do construtor com o `mb` substituído: o que se testa é a
    instrumentação — diretório, gate, manifesto —, não o MCMC."""

    NEWICK_CON = ("#NEXUS\nbegin trees;\n  translate\n    1 t1,\n    2 t2,\n    3 t3,\n"
                  "    4 t4;\n  tree con_50_majrule = [&U] ((1[&prob=1]:0.1,2:0.1)"
                  "[&prob=0.9]:0.1,3:0.1,4:0.1);\nend;\n")

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.dir, "Trees"))
        self.saida = os.path.join(self.dir, "Trees", "t_mrbayes.nexus")
        self.aln = MultipleSeqAlignment([
            SeqRecord(Seq("ACGTACGTAA"), id="t1"), SeqRecord(Seq("ACGTACGTAG"), id="t2"),
            SeqRecord(Seq("ACGAACGTTA"), id="t3"), SeqRecord(Seq("ACGAACGTTG"), id="t4")])

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def _rodar(self, stdout, **config):
        def executar(cmd, *a, cwd=None, stdin=None, **k):
            self.script = stdin.read()
            self.cwd = cwd
            with open(os.path.join(cwd, "alignment.nexus.con.tre"), "w") as fh:
                fh.write(self.NEWICK_CON)
            return subprocess.CompletedProcess(cmd, 0, stdout=stdout.encode(), stderr=b"")
        with mock.patch("workflow.tree_construction.builder.require_tool",
                        return_value="/opt/env/bin/mb"), \
             mock.patch("subprocess.run", side_effect=executar):
            return TreeBuilder(random_seed=4242, **config).mrbayes_constructor(
                self.aln, self.saida)

    def test_converge_grava_e_anota(self):
        self._rodar(STDOUT_SUMT.replace("= 0.045125", "= 0.004"))
        self.assertTrue(os.path.exists(self.saida))
        entrada = tool_runs.execucoes()["mrbayes"]
        self.assertEqual(entrada["seed"], 4242)
        run = entrada["runs"][0]
        self.assertEqual(run["convergencia"], "atingida")
        self.assertEqual(run["asdsf"], 0.004)
        self.assertIn("seed=4242 swapseed=4242", self.script)

    def test_nao_converge_recusa_e_nao_grava_em_trees(self):
        with self.assertRaises(ConvergenciaNaoAtingida):
            self._rodar(STDOUT_SUMT)
        self.assertFalse(os.path.exists(self.saida),
                         "árvore recusada em Trees/ seria reaproveitada depois")
        run = tool_runs.execucoes()["mrbayes"]["runs"][0]
        self.assertEqual(run["convergencia"], "nao_atingida")
        self.assertEqual(run["asdsf"], 0.045125)

    def test_gate_desligado_aceita_e_declara_nao_verificada(self):
        self._rodar(STDOUT_SUMT, mrbayes_asdsf_max=None)
        self.assertTrue(os.path.exists(self.saida))
        entrada = tool_runs.execucoes()["mrbayes"]
        self.assertEqual(entrada["gate_convergencia"], "desligado por configuração")
        self.assertEqual(entrada["runs"][0]["convergencia"], "nao_verificada")

    def test_diretorio_absoluto_e_independente_do_nome_do_repositorio(self):
        """D20, item 1: era `split('/PhyloTreeMiner/')[-1]` — relativo."""
        self._rodar(STDOUT_SUMT.replace("= 0.045125", "= 0.004"))
        self.assertTrue(os.path.isabs(self.cwd))
        self.assertEqual(self.cwd, os.path.join(self.dir, "tmp", "mrbayes_t_mrbayes"))

    def test_posterior_chega_a_arvore_gravada(self):
        """D20, item 7, pelo caminho do construtor: `[&prob=0.9]` do clado
        (t1,t2) vira `.confidence` na árvore devolvida e no Nexus gravado."""
        arvore = self._rodar(STDOUT_SUMT.replace("= 0.045125", "= 0.004"))
        clado = arvore.common_ancestor("t1", "t2")
        self.assertEqual(clado.confidence, 0.9)
        with open(self.saida) as fh:
            self.assertIn(")0.90000000:", fh.read())
        entrada = tool_runs.execucoes()["mrbayes"]
        self.assertIn("0-1", entrada["suporte"])
        self.assertIn("NÃO é comparável", entrada["nota_suporte"])

    def test_prob_fora_da_escala_e_recusada_e_nao_grava(self):
        for ruim in ("1.5", "-0.1", "NA", "nan"):
            with self.subTest(prob=ruim):
                self.NEWICK_CON = TestConstrutorComFerramentaSubstituida.NEWICK_CON.replace(
                    "[&prob=0.9]", f"[&prob={ruim}]")
                with self.assertRaises(ProbabilidadePosteriorInvalida):
                    self._rodar(STDOUT_SUMT.replace("= 0.045125", "= 0.004"))
                self.assertFalse(os.path.exists(self.saida))

    def test_script_com_modelo_selecionado(self):
        """M7.3: o `lset` vem da tradução; frequências iguais levam `prset`."""
        cfg = mrbayes_settings({})
        script = TreeBuilder.mrbayes_script(cfg, 1, "a.nexus", traduzir("K2P+I"))
        self.assertIn("execute a.nexus\nlset nst=2 rates=propinv\n"
                      "prset statefreqpr=fixed(equal)\nmcmc ", script)
        script = TreeBuilder.mrbayes_script(cfg, 1, "a.nexus", traduzir("GTR+F+G4"))
        self.assertIn("\nlset nst=6 rates=gamma ngammacat=4\nmcmc ", script)
        legado = TreeBuilder.mrbayes_script(cfg, 1, "a.nexus")
        self.assertIn("\nlset nst=6 rates=gamma\nmcmc ", legado)

    def test_nruns_1_falha_antes_de_chamar_a_ferramenta(self):
        with mock.patch("subprocess.run") as run:
            with self.assertRaises(ValueError):
                TreeBuilder(mrbayes_nruns=1).mrbayes_constructor(self.aln, self.saida)
            run.assert_not_called()


#: Anotação literal de nó interno do `.con.tre` do MrBayes 3.2.7 (`sumt`,
#: `nruns=2`; caracterização de 2026-09-23, 12 táxons de Zika).
ANOT_098 = ('[&prob=9.84042553e-01,prob_stddev=3.76120628e-03,prob_range='
            '{9.81382979e-01,9.86702128e-01},prob(percent)="98",prob+-sd="98+-0"]')
#: Clado presente em **uma só** das duas corridas (Nruns = 1 no `sumt`,
#: frequências 0 e 1): entra no consenso com 0,5 exato. Medido com 16 táxons e
#: 3 000 gerações, semente 12 — ASDSF 0,21, que o gate padrão recusa.
ANOT_UMA_CORRIDA = ('[&prob=5.00000000e-01,prob_stddev=7.07106781e-01,prob_range='
                    '{0.00000000e+00,1.00000000e+00},prob(percent)="50",prob+-sd="50+-71"]')
ANOT_1 = ('[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range='
          '{1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]')
COMPRIMENTO = ('[&length_mean=7.17464723e-03,length_median=7.40857200e-03,'
               'length_95%HPD={5.45463300e-03,8.20095500e-03}]')


class TestProbabilidadePosterior(unittest.TestCase):
    """D20, item 7 — `_clean_mrbayes_tree` sobre o formato real do `.con.tre`."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def _limpar(self, arvore):
        entrada = os.path.join(self.dir, "c.con.tre")
        with open(entrada, "w") as fh:
            fh.write("#NEXUS\nbegin trees;\n\ttranslate\n\t\t1\tA.1,\n\t\t2\tB.1,\n"
                     "\t\t3\tC.1,\n\t\t4\tD.1,\n\t\t5\tE.1\n\t\t;\n"
                     f"   tree con_50_majrule = [&U] {arvore};\nend;\n")
        from Bio import Phylo
        return Phylo.read(TreeBuilder()._clean_mrbayes_tree(
            entrada, os.path.join(self.dir, "c.nwk")), "newick")

    def _folha(self, n):
        return f"{n}{ANOT_1}:0.01{COMPRIMENTO}"

    def test_formato_real_prob_vira_confidence_folha_e_raiz_nao(self):
        arvore = self._limpar(
            f"({self._folha(1)},{self._folha(2)},(({self._folha(3)},{self._folha(4)})"
            f"{ANOT_098}:0.02{COMPRIMENTO},{self._folha(5)}){ANOT_UMA_CORRIDA}:0.03{COMPRIMENTO})")
        self.assertEqual(arvore.common_ancestor("C.1", "D.1").confidence, 0.984042553)
        # Presente em uma só corrida: o valor agregado (0,5) é o que o MrBayes
        # pôs no consenso; a discordância (prob_range {0,1}) fica no .con.tre.
        self.assertEqual(arvore.common_ancestor("C.1", "E.1").confidence, 0.5)
        self.assertIsNone(arvore.root.confidence)
        self.assertEqual(sorted(f.name for f in arvore.get_terminals()),
                         ["A.1", "B.1", "C.1", "D.1", "E.1"])
        self.assertTrue(all(f.confidence is None for f in arvore.get_terminals()))
        self.assertAlmostEqual(arvore.common_ancestor("C.1", "D.1").branch_length, 0.02)

    def test_no_interno_sem_prob_fica_none_nao_zero(self):
        arvore = self._limpar(
            f"({self._folha(1)},{self._folha(2)},(({self._folha(3)},{self._folha(4)})"
            f":0.02{COMPRIMENTO},{self._folha(5)}){ANOT_098}:0.03)")
        self.assertIsNone(arvore.common_ancestor("C.1", "D.1").confidence)
        self.assertEqual(arvore.common_ancestor("C.1", "E.1").confidence, 0.984042553)

    def test_prob_um_exato(self):
        arvore = self._limpar(
            f"({self._folha(1)},{self._folha(2)},(({self._folha(3)},{self._folha(4)})"
            f"{ANOT_1}:0.02,{self._folha(5)}){ANOT_1}:0.03)")
        clado = arvore.common_ancestor("C.1", "D.1")
        self.assertEqual(clado.confidence, 1.0)
        self.assertIsInstance(clado.confidence, float)

    def test_prob_fora_de_0_1_e_erro(self):
        with self.assertRaises(ProbabilidadePosteriorInvalida):
            self._limpar(f"({self._folha(1)},{self._folha(2)},({self._folha(3)},"
                         f"{self._folha(4)})[&prob=1.00000001e+00]:0.02,{self._folha(5)})")


@unittest.skipUnless(resolve_tool("mrbayes"), "MrBayes (mb) ausente")
class TestMrBayesDeVerdade(unittest.TestCase):
    """O binário real. Poucos segundos cada: 5-8 táxons, ngen pequeno."""

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.dir, "Trees"))

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)
        tool_runs.limpar()

    def _rodar(self, nome, ntax, **config):
        saida = os.path.join(self.dir, "Trees", f"{nome}.nexus")
        TreeBuilder(**config).mrbayes_constructor(_aln_zika(ntax), saida)
        return saida

    def test_mesma_semente_mesma_arvore(self):
        """D20, item 2. Antes, `Seed = <aleatório>` a cada execução."""
        cfg = dict(random_seed=31, mrbayes_ngen=5000, mrbayes_asdsf_max=None)
        a = self._rodar("a", 6, **cfg)
        b = self._rodar("b", 6, **cfg)
        with open(a) as fa, open(b) as fb:
            self.assertEqual(fa.read(), fb.read())

    def test_cadeia_curta_e_recusada(self):
        """8 táxons, 20 000 gerações: medido ASDSF ~0,045 na caracterização."""
        with self.assertRaises(ConvergenciaNaoAtingida):
            self._rodar("curta", 8, random_seed=12345, mrbayes_ngen=20000)
        run = tool_runs.execucoes()["mrbayes"]["runs"][0]
        self.assertGreaterEqual(run["asdsf"], 0.01)
        self.assertFalse(os.path.exists(os.path.join(self.dir, "Trees", "curta.nexus")))

    def test_cadeia_suficiente_e_aceita(self):
        """Mesmos 8 táxons, 100 000 gerações: medido ASDSF 0,007 em 2026-09-23.
        É o par do teste acima — o gate não recusa tudo."""
        saida = self._rodar("longa", 8, random_seed=12345, mrbayes_ngen=100000)
        self.assertTrue(os.path.exists(saida))
        run = tool_runs.execucoes()["mrbayes"]["runs"][0]
        self.assertEqual(run["convergencia"], "atingida")
        self.assertLess(run["asdsf"], 0.01)

        # D20, item 7 — o suporte gravado é a probabilidade do `sumt`
        # (`.tstat`), clado a clado. O oráculo completo, por bipartição e
        # contra #obs/N, é `docs/science/scripts/oraculo_posterior_mrbayes.py`.
        from Bio import Phylo
        arvore = Phylo.read(saida, "nexus")
        internos = [c for c in arvore.get_nonterminals() if c is not arvore.root]
        self.assertTrue(internos)
        with open(os.path.join(self.dir, "tmp", "mrbayes_longa", "alignment.nexus.tstat")) as fh:
            linhas = [l.rstrip("\n").split("\t") for l in fh][2:]
        probabilidades = [float(l[3]) for l in linhas if len(l) > 3]
        for clado in internos:
            self.assertIsNotNone(clado.confidence)
            self.assertTrue(any(abs(clado.confidence - p) <= 5e-7 for p in probabilidades),
                            (clado.confidence, probabilidades))


if __name__ == "__main__":
    unittest.main()
