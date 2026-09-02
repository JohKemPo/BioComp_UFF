"""
M3.2/M7.2 — o RAxML-NG passa a rodar com bootstrap e a devolver suporte.

Antes desta mudança, `raxml_ng_constructor` só fazia busca de ML (sem
`--bootstrap`/`--all`) e lia `<prefix>.raxml.bestTree`, que não carrega
suporte nenhum — confirmado em `docs/science/08-ficha-de-chamada-por-metodo.md
§3` contra artefato real (`Zika_21seq_validacao`). Isto deixava o RAxML-NG
como o único dos três métodos avançados sem suporte de ramo propagável.

Este teste roda o RAxML-NG **de verdade**, numa entrada sintética pequena
fixa no próprio arquivo — não depende de nenhum genoma real nem de rede — só
para confirmar o que a leitura de código e do `--help` não bastam para provar:
que `--all --bs-trees N` grava `<prefix>.raxml.support` com valores de
suporte, e que esses valores sobrevivem à leitura do Newick pelo Biopython
como `.confidence` em nós internos (o mesmo mecanismo genérico documentado
para FastTree e IQ-TREE na ficha de chamada por método, §1/§2).

`raxml_threads=1`: com o default de 4 threads, o RAxML-NG recusa rodar em
cima de um alinhamento deste tamanho («ERROR: Too few patterns per thread!»,
medido em 2026-09-02) — comportamento pré-existente, independente desta
mudança, e fora do escopo deste lote (ver relatório do lote).

Executar com:

    python -m unittest workflow.tests.test_raxml_bootstrap
"""

from __future__ import annotations

import os
import shutil
import tempfile
import unittest

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from workflow.utils.external_tools import resolve_tool


def _raxml_ng_disponivel() -> bool:
    return resolve_tool('raxml-ng') is not None


def _alinhamento_sintetico():
    """Seis táxons, com SNPs suficientes para dar topologia interna não
    trivial. Não é dado real — sequências inventadas só para exercitar o
    inferidor."""
    bases = [
        ("tax1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ("tax2", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAACGTACGTACGT"),
        ("tax3", "ACGTACGAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
        ("tax4", "ACGTACGTACGTACCTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA"),
        ("tax5", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACCTACGTACGAACGTACGT"),
        ("tax6", "ACGAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTCGT"),
    ]
    return MultipleSeqAlignment([SeqRecord(Seq(seq), id=nome) for nome, seq in bases])


@unittest.skipUnless(_raxml_ng_disponivel(), "raxml-ng não resolvido no PATH/env")
class TestRaxmlNgProduzSuporte(unittest.TestCase):
    """Execução real, em alinhamento sintético pequeno — sem rede, sem dado
    de projeto, sem genoma de Variola/Zika."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.trees = os.path.join(self.dir, "Trees")
        os.makedirs(self.trees)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_arvore_devolvida_tem_confidence_em_nos_internos(self):
        from workflow.tree_construction.builder import TreeBuilder

        saida = os.path.join(self.trees, "t_raxml_bootstrap.nexus")
        builder = TreeBuilder(random_seed=12345, raxml_threads=1)
        tree = builder.raxml_ng_constructor(_alinhamento_sintetico(), saida)

        internos = [c for c in tree.get_nonterminals()]
        self.assertTrue(internos, "árvore sem nó interno — não dá para checar suporte")
        com_suporte = [c for c in internos if c.confidence is not None]
        self.assertTrue(
            com_suporte,
            "nenhum nó interno com .confidence — o suporte não chegou à árvore lida",
        )
        # Suporte é fração/percentual de bootstrap, nunca negativo.
        for clade in com_suporte:
            self.assertGreaterEqual(clade.confidence, 0)

    def test_suporte_sobrevive_ao_nexus_gravado(self):
        """M3.1 (parte BioComp_UFF): confirma que o suporte do RAxML-NG
        também sobrevive ao `Phylo.write` para Nexus, do mesmo jeito que já
        acontece com FastTree e IQ-TREE (ficha §1/§2)."""
        from Bio import Phylo
        from workflow.tree_construction.builder import TreeBuilder

        saida = os.path.join(self.trees, "t_raxml_bootstrap_nexus.nexus")
        builder = TreeBuilder(random_seed=12345, raxml_threads=1)
        builder.raxml_ng_constructor(_alinhamento_sintetico(), saida)

        self.assertTrue(os.path.exists(saida))
        releitura = Phylo.read(saida, "nexus")
        internos = releitura.get_nonterminals()
        com_suporte = [c for c in internos if c.confidence is not None]
        self.assertTrue(
            com_suporte,
            "o suporte não sobreviveu ao Nexus gravado em disco",
        )

    def test_comando_usa_all_e_bs_trees(self):
        """Trava a linha de comando efetiva: modo `--all`, mesma contagem de
        réplicas do UFBoot do IQ-TREE (`-bb 1000`), `--workers 1` fixo (D17)."""
        from unittest import mock

        from workflow.tree_construction import builder as builder_mod
        from workflow.utils import tool_runs

        tool_runs.limpar()
        saida = os.path.join(self.trees, "t_raxml_comando.nexus")
        original_run = builder_mod.subprocess.run
        comandos = []

        def espiao(cmd, *args, **kwargs):
            comandos.append(cmd)
            return original_run(cmd, *args, **kwargs)

        with mock.patch("workflow.tree_construction.builder.subprocess.run", espiao):
            builder_mod.TreeBuilder(random_seed=12345, raxml_threads=1).raxml_ng_constructor(
                _alinhamento_sintetico(), saida
            )

        self.assertEqual(len(comandos), 1)
        cmd = comandos[0]
        self.assertIn("--all", cmd)
        self.assertIn("--bs-trees", cmd)
        self.assertEqual(cmd[cmd.index("--bs-trees") + 1], "1000")
        self.assertIn("--workers", cmd)
        self.assertEqual(cmd[cmd.index("--workers") + 1], "1")
        self.assertNotIn("auto", cmd)
        tool_runs.limpar()


if __name__ == "__main__":
    unittest.main()
