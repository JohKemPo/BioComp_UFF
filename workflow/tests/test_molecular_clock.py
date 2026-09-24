"""
E12.10/E12.11 — relógio molecular leve (LSD2) como etapa do pipeline.

Três coisas são testadas aqui, e cada uma corresponde a uma parte do parecer
de zona sagrada (`docs/science/16-relogio-molecular-producao.md`):

1. **Falha declarada, nunca silenciosa.** Desligado → ``ignorado_por_configuracao``;
   sem grupo externo, sem datas, grupo externo não monofilético, sem sinal
   temporal → ``tentado_e_falhou`` com motivo, e o JSON **é** escrito dizendo
   por quê. Nenhum desses caminhos chama o IQ-TREE.
2. **Identidade de clado.** A chave de `datas_por_clado` é o `canonical_item_id`
   da bipartição contra todas as pontas da árvore de origem — o mesmo que
   `Backend/src/suporte_de_ramo.py` calcula lendo o mesmo arquivo. O teste
   refaz esse cálculo por fora (Bio.Phylo, como o Backend) e compara.
3. **Oráculo com verdade conhecida** (só se o IQ-TREE estiver no PATH): árvore
   de tempo sorteada, relógio estrito, alinhamento simulado com AliSim; a data
   recuperada de cada clado é comparada com a verdadeira, casando pelo
   `clade_id`.

Executar com:

    python -m unittest workflow.tests.test_molecular_clock
"""

from __future__ import annotations

import json
import os
import random
import shutil
import statistics
import subprocess
import tempfile
import unittest
from io import StringIO
from unittest import mock

from Bio import Phylo

from workflow.controller.molecularClockController import MolecularClockController
from workflow.molecular_clock import datacao
from workflow.stability.clade_identity import (canonical_bipartition, canonical_item_id,
                                               strip_accession_version)
from workflow.utils import external_tools, tool_runs

# Árvore de 8 pontas; o grupo externo {OG000001.1, OG000002.1} é monofilético.
NEWICK = ("((OG000001.1:0.02,OG000002.1:0.02):0.05,"
          "((AB000001.1:0.010,AB000002.1:0.012):0.004,(AB000003.1:0.011,AB000004.1:0.013):0.003):0.01,"
          "((AB000005.1:0.014,AB000006.1:0.016):0.005,AB000007.1:0.015):0.01);")
TODAS = ["OG000001.1", "OG000002.1"] + [f"AB00000{i}.1" for i in range(1, 8)]


def _nexus(newick: str) -> str:
    return f"#NEXUS\nBegin Trees;\n Tree tree1={newick}\nEnd;\n"


def _gb(datas: dict) -> str:
    """GenBank mínimo: só VERSION e, quando houver, /collection_date."""
    regs = []
    for acc, data in datas.items():
        cd = f'                     /collection_date="{data}"\n' if data else ""
        regs.append(f"LOCUS       {acc.split('.')[0]}  10 bp\nVERSION     {acc}\n"
                    f"FEATURES             Location/Qualifiers\n     source          1..10\n{cd}//")
    return "\n".join(regs) + "\n"


def _ids_como_o_backend(newick: str) -> set:
    """Refaz `suporte_de_ramo._ramos_com_identidade` por fora, com Bio.Phylo."""
    arv = Phylo.read(StringIO(newick), "newick")
    todos = frozenset(strip_accession_version(t.name) for t in arv.get_terminals())
    ids = set()
    for cl in arv.find_clades():
        if cl.is_terminal():
            continue
        bip = canonical_bipartition(
            frozenset(strip_accession_version(t.name) for t in cl.get_terminals()), todos)
        if bip is not None:
            ids.add(canonical_item_id(bip))
    return ids


class _Projeto(unittest.TestCase):
    """Monta um `out/` mínimo num diretório temporário."""

    def setUp(self):
        tool_runs.limpar()
        self.dir = tempfile.mkdtemp(prefix="relogio_")
        self.out = os.path.join(self.dir, "out")
        for d in ("Trees", "Align", "outputs"):
            os.makedirs(os.path.join(self.out, d))
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)

    def escrever(self, newick=NEWICK, datas=None, rotulos=TODAS, arvore="mafft_iqtree"):
        with open(os.path.join(self.out, "Trees", f"tree_dataset_final_{arvore}.nexus"), "w") as fh:
            fh.write(_nexus(newick))
        with open(os.path.join(self.out, "Align", "dataset_final_mafft.aln"), "w") as fh:
            fh.write("".join(f">{r}\nACGTACGTAC\n" for r in rotulos))
        if datas is None:
            datas = {r: f"{2000 + i}-06-15" for i, r in enumerate(rotulos)}
        with open(os.path.join(self.out, "outputs", "raw_data_sequences.gb"), "w") as fh:
            fh.write(_gb(datas))

    def controlador(self, **extra):
        cfg = {"output_path": self.out, "run_molecular_clock": True,
               "molecular_clock_outgroup": ["OG000001.1", "OG000002.1"], "random_seed": 7}
        cfg.update(extra)
        return MolecularClockController(**cfg)

    def json(self):
        with open(os.path.join(self.out, "outputs", "molecular_clock",
                               "relogio_molecular.json")) as fh:
            return json.load(fh)

    def desfechos(self):
        return [(m["estado"], m.get("motivo", "")) for m in tool_runs.metodos()
                if m["metodo"] == "relogio_molecular"]


# --------------------------------------------------------------------------- #
class TestConversaoDeDatas(unittest.TestCase):

    def test_formatos_do_genbank(self):
        casos = {
            "19-JUN-2016": ("2016-06-19", "dia"),
            "2016-05-24": ("2016-05-24", "dia"),
            "Jun-2007": ("2007-06", "mes"),
            "2010-06": ("2010-06", "mes"),
            "2010": ("2010:2011", "ano"),
        }
        for crua, esperado in casos.items():
            self.assertEqual(datacao.converter_data(crua), esperado, crua)
        self.assertEqual(datacao.converter_data("2010", ano_pontual=True), ("2010", "ano"))

    def test_formato_desconhecido_e_erro_nao_palpite(self):
        for crua in ("circa 1990", "1990/2000", "Xyz-2001"):
            with self.assertRaises(ValueError, msg=crua):
                datacao.converter_data(crua)


class TestLinhaDeResultado(unittest.TestCase):

    def test_decimal_com_intervalo(self):
        r = datacao.interpretar_linha_de_resultado(
            "rate 0.000196644 [5.90887e-05; 0.000325283], tMRCA 1685.49 [1047.93; 1809.59], "
            "objective function 0.0913245")
        self.assertTrue(r["solucao_unica"])
        self.assertAlmostEqual(r["taxa_subs_por_sitio_por_ano"], 0.000196644)
        self.assertEqual(r["intervalo_da_taxa"], [5.90887e-05, 0.000325283])
        self.assertAlmostEqual(r["tmrca"], 1685.49)
        self.assertEqual(r["intervalo_do_tmrca"], [1047.93, 1809.59])

    def test_data_de_calendario(self):
        r = datacao.interpretar_linha_de_resultado(
            "rate 0.000523787 [0.000434841; 0.000634798], tMRCA 1944-12-05 "
            "[1932-04-16; 1951-06-21], objective function 1.2")
        self.assertEqual(r["tmrca_lsd2"], "1944-12-05")
        self.assertAlmostEqual(r["tmrca"], 1944 + 339 / 366, places=6)
        self.assertAlmostEqual(r["intervalo_do_tmrca"][0], 1932.29, places=2)

    def test_na_nao_vira_numero(self):
        # E12 §4 e §8.1: com :NA o LSD2 imprime um número que não é estimativa.
        for linha in ("rate 0.000125:NA, tMRCA 2521-09-14 , objective function 3",
                      "rate NA:1e-10, tMRCA NA:-266334 , objective function inf"):
            r = datacao.interpretar_linha_de_resultado(linha)
            self.assertFalse(r["solucao_unica"], linha)
            self.assertIsNone(r["taxa_subs_por_sitio_por_ano"])
            self.assertIsNone(r["tmrca"])


class TestIdentidadeDeClado(unittest.TestCase):
    """A chave de `datas_por_clado` é o `clade_id` que o Backend calcula."""

    def test_clade_id_casa_com_o_da_arvore_de_origem(self):
        # Timetree como o LSD2 grava: sem o grupo externo, enraizado, rótulo
        # truncado (D13) numa ponta, CI_date entre aspas.
        timetree = (
            '#NEXUS\nBegin trees;\ntree 1 = (((AB000001.[&date=2001]:1,AB000002.1[&date=2002]:2)'
            '[&date="2000.5",CI_height={1,3},CI_date="{1999.1,2000.9}"]:3,'
            '(AB000003.1[&date=2003]:1,AB000004.1[&date=2004]:1)[&date="2002",'
            'CI_date="{2001,2002.5}"]:2)[&date="1997.2"]:1,'
            '((AB000005.1[&date=2005]:1,AB000006.1[&date=2006]:1)[&date="2004"]:1,'
            'AB000007.1[&date=2007]:3)[&date="2003"]:3)[&date="1995-07-02"];\nEnd;\n')
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d)
        caminho = os.path.join(d, "t.timetree.nex")
        with open(caminho, "w") as fh:
            fh.write(timetree)
        todos = frozenset(strip_accession_version(x) for x in TODAS)
        datas, contagem = datacao.datas_por_clado(caminho, todos)

        ids_backend = _ids_como_o_backend(NEWICK)
        self.assertTrue({v["clade_id"] for v in datas.values()} <= ids_backend)
        self.assertEqual(set(datas), {str(v["clade_id"]) for v in datas.values()})
        # 6 nós internos no timetree, e todos são bipartições não triviais
        # contra as 9 pontas da origem (o grupo externo conta).
        self.assertEqual(contagem["nos_internos_na_arvore_de_tempo"], 6)
        self.assertEqual(len(datas), 6)

        c12 = canonical_item_id(frozenset({"AB000001", "AB000002"}))
        self.assertEqual(datas[str(c12)]["data_estimada"], 2000.5)
        self.assertEqual(datas[str(c12)]["intervalo"], [1999.1, 2000.9])
        self.assertTrue(datas[str(c12)]["clado_datado_e_o_lado_canonico"])

        # A raiz do timetree é o grupo interno inteiro (7 pontas); o lado
        # canônico da aresta é o grupo externo (2 pontas). A data é a do MRCA
        # do grupo interno, e o JSON diz que não é o lado canônico.
        raiz = datas[str(canonical_item_id(frozenset({"OG000001", "OG000002"})))]
        self.assertFalse(raiz["clado_datado_e_o_lado_canonico"])
        self.assertEqual(raiz["n_taxa"], 2)
        self.assertEqual(raiz["n_taxa_do_clado_datado"], 7)
        self.assertEqual(raiz["data_estimada_lsd2"], "1995-07-02")
        self.assertAlmostEqual(raiz["data_estimada"], 1995.5, places=1)

    def test_ramo_colapsado_nao_tem_data(self):
        timetree = ('#NEXUS\nBegin trees;\ntree 1 = ((AB000001.1[&date=2001]:1,AB000002.1'
                    '[&date=2002]:1,AB000003.1[&date=2003]:1,AB000004.1[&date=2004]:1)'
                    '[&date="1999"]:1,((AB000005.1[&date=2005]:1,AB000006.1[&date=2006]:1)'
                    '[&date="2004"]:1,AB000007.1[&date=2007]:3)[&date="2003"]:3)[&date="1990"];\nEnd;\n')
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d)
        caminho = os.path.join(d, "t.timetree.nex")
        with open(caminho, "w") as fh:
            fh.write(timetree)
        todos = frozenset(strip_accession_version(x) for x in TODAS)
        datas, _ = datacao.datas_por_clado(caminho, todos)
        self.assertNotIn(str(canonical_item_id(frozenset({"AB000001", "AB000002"}))), datas)
        self.assertIn(str(canonical_item_id(frozenset({"AB000001", "AB000002", "AB000003",
                                                       "AB000004"}))), datas)


# --------------------------------------------------------------------------- #
class TestFalhaDeclarada(_Projeto):
    """Nenhum destes caminhos pode chegar ao IQ-TREE."""

    def setUp(self):
        super().setUp()
        patch = mock.patch("workflow.molecular_clock.datacao.subprocess.run",
                           side_effect=AssertionError("o IQ-TREE não deveria ser chamado"))
        patch.start()
        self.addCleanup(patch.stop)
        patch2 = mock.patch.object(external_tools, "resolve_tool", return_value="/bin/iqtree3")
        patch2.start()
        self.addCleanup(patch2.stop)
        patch3 = mock.patch("workflow.controller.molecularClockController._versao_iqtree",
                            return_value="3.1.3")
        patch3.start()
        self.addCleanup(patch3.stop)

    def test_desligado_e_ignorado_por_configuracao(self):
        self.escrever()
        self.assertIsNone(self.controlador(run_molecular_clock=False)())
        self.assertEqual(self.desfechos()[0][0], "ignorado_por_configuracao")
        self.assertFalse(os.path.exists(os.path.join(self.out, "outputs", "molecular_clock")))

    def test_sem_grupo_externo(self):
        self.escrever()
        r = self.controlador(molecular_clock_outgroup=None)()
        self.assertEqual(r["estado"], "tentado_e_falhou")
        self.assertIn("grupo externo não declarado", r["motivo"])
        self.assertEqual(self.json()["estado"], "tentado_e_falhou")
        self.assertEqual(self.desfechos()[0][0], "tentado_e_falhou")
        self.assertIsNone(self.json()["tmrca"])

    def test_datas_insuficientes_caso_varv49(self):
        # Nenhum registro com collection_date: o LSD2 não falharia sozinho
        # (grava date="0" em tudo, E12 §8.1). O portão tem de recusar antes.
        self.escrever(datas={r: None for r in TODAS})
        r = self.controlador()()
        self.assertEqual(r["estado"], "tentado_e_falhou")
        arv = r["por_arvore"]["tree_dataset_final_mafft_iqtree.nexus"]
        self.assertIn("datas de coleta insuficientes", arv["motivo"])
        self.assertEqual(arv["datas"]["pontas_datadas"], 0)
        self.assertEqual(self.desfechos(), [("tentado_e_falhou", arv["motivo"])])
        self.assertEqual(self.json()["datas_por_clado"], {})

    def test_grupo_externo_nao_monofiletico(self):
        self.escrever()
        r = self.controlador(molecular_clock_outgroup="OG000001.1,AB000005.1")()
        arv = r["por_arvore"]["tree_dataset_final_mafft_iqtree.nexus"]
        self.assertIn("não é monofilético", arv["motivo"])
        self.assertFalse(arv["grupo_externo"]["monofiletico"])

    def test_grupo_externo_como_arquivo_e_rotulo_truncado(self):
        # D13: a árvore de IQ-TREE/RAxML grava `OG000001.`; o arquivo declara
        # `OG000001.1`. Os dois têm de casar.
        self.escrever(newick=NEWICK.replace("OG000001.1", "OG000001."))
        arquivo = os.path.join(self.dir, "og.txt")
        with open(arquivo, "w") as fh:
            fh.write("OG000001.1\nOG000002.1\n")
        with mock.patch("workflow.molecular_clock.datacao.sinal_temporal",
                        return_value={"inclinacao_subs_sitio_ano": -1.0, "r2": 0.5, "n": 7}):
            r = self.controlador(molecular_clock_outgroup=arquivo)()
        arv = r["por_arvore"]["tree_dataset_final_mafft_iqtree.nexus"]
        self.assertEqual(arv["grupo_externo"]["encontrado"], ["OG000001", "OG000002"])
        self.assertTrue(arv["grupo_externo"]["monofiletico"])
        self.assertIn("sem sinal temporal", arv["motivo"])

    def test_sem_sinal_temporal(self):
        # Pontas mais distantes da raiz são as mais ANTIGAS: inclinação negativa.
        arv = Phylo.read(StringIO(NEWICK), "newick")
        arv.root_with_outgroup("OG000001.1", "OG000002.1")
        dist = {t.name: arv.distance(t) for t in arv.get_terminals()}
        datas = {n: f"{int(2020 - 1000 * d)}-01-01" for n, d in dist.items()}
        self.escrever(datas=datas)
        r = self.controlador()()
        arv_r = r["por_arvore"]["tree_dataset_final_mafft_iqtree.nexus"]
        self.assertIn("sem sinal temporal", arv_r["motivo"])
        self.assertLess(arv_r["sinal_temporal"]["inclinacao_subs_sitio_ano"], 0)

    def test_metodo_nao_datavel_fica_de_fora_declarado(self):
        self.escrever(arvore="mafft_nj_distance")
        r = self.controlador()()
        self.assertEqual(r["estado"], "tentado_e_falhou")
        self.assertIn("tree_dataset_final_mafft_nj_distance.nexus", r["arvores_nao_dataveis"])
        self.assertIn("nenhuma árvore", r["motivo"])


# --------------------------------------------------------------------------- #
def _cronograma(n_interno: int, semente: int):
    """
    Árvore de tempo com verdade conhecida: grupo interno de `n_interno` pontas
    com raiz em 1950, grupo externo de 2 pontas ligado numa raiz em 1850.

    Return
    ------
    tuple
        ``(newick em substituições, datas das pontas, {frozenset: data do nó})``.
    """
    rnd = random.Random(semente)
    taxa = 1e-3
    datas, verdade = {}, {}
    contador = iter(range(1, 10_000))

    def sub(n, t):
        if n == 1:
            nome = f"AB{next(contador):06d}.1"
            datas[nome] = min(2020.0, t + rnd.uniform(3, 30))
            return nome, datas[nome] - t, frozenset({nome.split('.')[0]})
        k = rnd.randint(1, n - 1)
        filhos = []
        folhas = frozenset()
        for m in (k, n - k):
            t_filho = t + rnd.uniform(1, 4) if m > 1 else t
            nwk, dur, fs = sub(m, t_filho)
            filhos.append((nwk, dur + (t_filho - t)))
            folhas |= fs
        verdade[folhas] = t
        corpo = "(" + ",".join(f"{c}:{d * taxa:.8f}" for c, d in filhos) + ")"
        return corpo, 0.0, folhas

    interno, _, _ = sub(n_interno, 1950.0)
    og = [f"OG{i:06d}.1" for i in (1, 2)]
    for o in og:
        datas[o] = 2000.0
    externo = f"({og[0]}:{50 * taxa:.8f},{og[1]}:{50 * taxa:.8f})"
    # raiz 1850: externo 1850 -> 1950 (+50 até as pontas) e interno 1850 -> 1950
    newick = f"({externo}:{100 * taxa:.8f},{interno}:{100 * taxa:.8f});"
    return newick, datas, verdade


@unittest.skipUnless(external_tools.resolve_tool("iqtree"), "IQ-TREE ausente do PATH")
class TestOraculoVerdadeConhecida(_Projeto):
    """Relógio estrito, alinhamento simulado: o método tem de recuperar a verdade."""

    def test_recupera_as_datas_de_no(self):
        newick, datas, verdade = _cronograma(16, semente=20260924)
        iq = external_tools.resolve_tool("iqtree")
        nwk = os.path.join(self.dir, "verdadeira.nwk")
        with open(nwk, "w") as fh:
            fh.write(newick + "\n")
        subprocess.run([iq, "--alisim", os.path.join(self.dir, "sim"), "-m", "JC", "-t", nwk,
                        "--length", "8000", "--seed", "11", "-af", "fasta", "-redo", "-quiet"],
                       check=True, capture_output=True)
        shutil.copyfile(os.path.join(self.dir, "sim.fa"),
                        os.path.join(self.out, "Align", "dataset_final_mafft.aln"))
        with open(os.path.join(self.out, "Trees", "tree_dataset_final_mafft_iqtree.nexus"), "w") as fh:
            fh.write(_nexus(newick))
        with open(os.path.join(self.out, "outputs", "raw_data_sequences.gb"), "w") as fh:
            # dia exato: a verdade é contínua, e a convenção de ano-só é testada à parte
            fh.write(_gb({n: _decimal_para_dia(d) for n, d in datas.items()}))

        r = MolecularClockController(output_path=self.out, run_molecular_clock=True,
                                     molecular_clock_outgroup=["OG000001.1", "OG000002.1"],
                                     random_seed=7, molecular_clock_ci=50)()
        self.assertEqual(r["estado"], "executado", r.get("motivo"))
        self.assertEqual(self.desfechos()[0][0], "executado")
        self.assertTrue(os.path.isfile(os.path.join(self.out, "outputs", "molecular_clock",
                                                    "timetree.nex")))

        # Toda chave é clade_id de bipartição da árvore de origem.
        self.assertTrue({v["clade_id"] for v in r["datas_por_clado"].values()}
                        <= _ids_como_o_backend(newick))

        todos = frozenset(n.split(".")[0] for n in datas)
        erros = []
        for fs, t in verdade.items():
            bip = canonical_bipartition(fs, todos)
            if bip is None:
                continue
            est = r["datas_por_clado"].get(str(canonical_item_id(bip)))
            if est is not None:
                erros.append(abs(est["data_estimada"] - t))
        self.assertGreaterEqual(len(erros), 5, "poucos nós comparáveis")
        self.assertLess(statistics.median(erros), 3.0, erros)
        self.assertLess(abs(r["tmrca"] - 1950.0), 15.0, r["tmrca"])
        self.assertGreater(r["taxa_subs_por_sitio_por_ano"], 5e-4)
        self.assertLess(r["taxa_subs_por_sitio_por_ano"], 2e-3)

    def test_reprodutivel_com_uma_thread_e_semente(self):
        newick, datas, _ = _cronograma(12, semente=3)
        iq = external_tools.resolve_tool("iqtree")
        nwk = os.path.join(self.dir, "v.nwk")
        with open(nwk, "w") as fh:
            fh.write(newick + "\n")
        subprocess.run([iq, "--alisim", os.path.join(self.dir, "sim"), "-m", "JC", "-t", nwk,
                        "--length", "5000", "--seed", "5", "-af", "fasta", "-redo", "-quiet"],
                       check=True, capture_output=True)
        shutil.copyfile(os.path.join(self.dir, "sim.fa"),
                        os.path.join(self.out, "Align", "dataset_final_mafft.aln"))
        with open(os.path.join(self.out, "Trees", "tree_dataset_final_mafft_iqtree.nexus"), "w") as fh:
            fh.write(_nexus(newick))
        with open(os.path.join(self.out, "outputs", "raw_data_sequences.gb"), "w") as fh:
            fh.write(_gb({n: _decimal_para_dia(d) for n, d in datas.items()}))
        cfg = dict(output_path=self.out, run_molecular_clock=True, random_seed=7,
                   molecular_clock_outgroup=["OG000001.1", "OG000002.1"], molecular_clock_ci=20)
        a = MolecularClockController(**cfg)()
        tool_runs.limpar()
        b = MolecularClockController(**cfg)()
        for k in ("taxa_subs_por_sitio_por_ano", "tmrca", "intervalo_do_tmrca", "datas_por_clado"):
            self.assertEqual(a[k], b[k], k)


def _decimal_para_dia(ano_decimal: float) -> str:
    import datetime as dt

    ano = int(ano_decimal)
    ini = dt.date(ano, 1, 1)
    dias = (dt.date(ano + 1, 1, 1) - ini).days
    return (ini + dt.timedelta(days=int((ano_decimal - ano) * dias))).isoformat()


if __name__ == "__main__":
    unittest.main()
