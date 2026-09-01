"""
Testes do módulo `workflow.stability`.

Executar com:

    python -m unittest workflow.tests.test_stability
"""

from __future__ import annotations

import io
import os
import sys
import unittest

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, "../.."))

from Bio import Phylo

from workflow.stability.clade_identity import (
    audit_legacy_identity,
    canonical_digest,
    clade_identities,
    legacy_clade_hash,
    legacy_terminal_hash,
)
from workflow.stability.stability import (
    PipelineLabel,
    StabilityAnalyzer,
    TreeSet,
    strip_accession_version,
)


def newick(text: str):
    """Constrói uma árvore a partir de uma string Newick."""
    return Phylo.read(io.StringIO(text), "newick")


class TestCladeIdentity(unittest.TestCase):
    """A identidade canônica deve ser invariante à ordem; a legada, não."""

    def setUp(self):
        # Mesma topologia, filhos em ordem trocada.
        self.left = newick("(((A,B),C),(D,E));")
        self.right = newick("((C,(B,A)),(E,D));")

    def test_canonical_identity_is_order_invariant(self):
        left = {identity.taxa for identity in clade_identities(self.left)}
        right = {identity.taxa for identity in clade_identities(self.right)}
        self.assertEqual(left, right)

    def test_canonical_digest_is_order_invariant(self):
        self.assertEqual(canonical_digest(["A", "B", "C"]), canonical_digest(["C", "A", "B"]))

    def test_legacy_identity_is_order_dependent(self):
        hashes = [legacy_terminal_hash(name) for name in ("A", "B", "C")]
        self.assertNotEqual(legacy_clade_hash(hashes), legacy_clade_hash(list(reversed(hashes))))

    def test_audit_detects_fragmentation(self):
        audit = audit_legacy_identity({"left": self.left, "right": self.right})
        self.assertGreater(len(audit.fragmented), 0)
        self.assertEqual(audit.n_canonical, 3)  # (A,B), (A,B,C) e (D,E); a raiz é excluída
        self.assertGreater(audit.n_legacy, audit.n_canonical)

    def test_root_clade_is_excluded_by_default(self):
        sizes = {identity.size for identity in clade_identities(self.left)}
        self.assertNotIn(5, sizes)
        self.assertIn(5, {identity.size for identity in clade_identities(self.left, True)})


class TestPipelineLabel(unittest.TestCase):
    """O rótulo do arquivo deve decompor-se nos fatores do delineamento."""

    def test_parses_aligner_and_inference(self):
        label = PipelineLabel.parse("tree_dataset_final_mafft_iqtree.nexus", prefix="tree_")
        self.assertEqual(label.aligner, "mafft")
        self.assertEqual(label.inference, "iqtree")

    def test_parses_two_token_inference_method(self):
        label = PipelineLabel.parse("tree_dataset_final_clustalo_nj_distance.nexus", prefix="tree_")
        self.assertEqual(label.aligner, "clustalo")
        self.assertEqual(label.inference, "nj_distance")

    def test_distinguishes_mafft_from_mafft_iterative(self):
        """Achado da reexecução de 2026-09-01 (Zika-21, VARV-6, VARV-49): antes
        desta correção, `mafft_iterative` nunca era reconhecido como alinhador
        — "mafft" é um prefixo de "mafft_iterative", e o antigo `tokens = set
        (stem.split("_"))` via os dois como tokens separados. O pipeline do
        braço iterativo herdava o rótulo do braço progressivo em silêncio, e
        `TreeSet.from_directory` só não perdia a árvore por sorte: a guarda
        contra rótulo duplicado (D19) barrava com `ValueError` em vez de
        sobrescrever — o que bloqueava M1.3 em toda reexecução com os dois
        braços do MAFFT, não corrigia o rótulo."""
        progressivo = PipelineLabel.parse(
            "tree_dataset_final_mafft_fasttree.nexus", prefix="tree_dataset_final_"
        )
        iterativo = PipelineLabel.parse(
            "tree_dataset_final_mafft_iterative_fasttree.nexus", prefix="tree_dataset_final_"
        )
        self.assertEqual(progressivo.aligner, "mafft")
        self.assertEqual(iterativo.aligner, "mafft_iterative")
        self.assertEqual(progressivo.inference, iterativo.inference)
        self.assertNotEqual(progressivo.name, iterativo.name)

    def test_mafft_iterative_com_metodo_de_dois_tokens(self):
        """O braço iterativo combinado com um método de inferência que também
        tem "_" no nome — garante que a correção do alinhador não quebra a
        escolha do sufixo mais longo já corrigida por D19."""
        label = PipelineLabel.parse(
            "tree_dataset_final_mafft_iterative_nj_parsimony.nexus",
            prefix="tree_dataset_final_",
        )
        self.assertEqual(label.aligner, "mafft_iterative")
        self.assertEqual(label.inference, "nj_parsimony")

    def test_alinhador_desconhecido_nao_quebra(self):
        label = PipelineLabel.parse(
            "tree_dataset_final_desconhecido_iqtree.nexus", prefix="tree_dataset_final_"
        )
        self.assertEqual(label.aligner, "unknown")
        self.assertEqual(label.inference, "iqtree")

    def test_directory_com_os_dois_braços_do_mafft_nao_colide(self):
        """Reprodução direta do crash de `conferir_correcoes_m1.py` nas três
        reexecuções de 2026-09-01: dois arquivos do mesmo método de inferência,
        um por braço do alinhador, no mesmo diretório."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            arvore = newick("((A,B),(C,D));")
            for nome in ("mafft_fasttree", "mafft_iterative_fasttree"):
                caminho = os.path.join(tmp, f"tree_dataset_final_{nome}.nexus")
                Phylo.write(arvore, caminho, "nexus")
            tree_set = TreeSet.from_directory(tmp, normalizer=None)
            self.assertEqual(len(tree_set), 2)
            self.assertEqual(
                {label.name for label in tree_set.labels.values()},
                {"mafft_fasttree", "mafft_iterative_fasttree"},
            )


class TestAccessionNormalisation(unittest.TestCase):
    """Rótulos truncados por IQ-TREE/RAxML devem reconciliar-se."""

    def test_truncated_version_matches_full_version(self):
        self.assertEqual(strip_accession_version("NC_008030."), strip_accession_version("NC_008030.1"))

    def test_quotes_and_whitespace_are_removed(self):
        self.assertEqual(strip_accession_version(" 'AF438165.1' "), "AF438165")

    def test_tree_set_reconciles_and_reports_mangling(self):
        trees = {
            "mafft_fasttree": newick("((NC_008030.1,B),(C,D));"),
            "mafft_iqtree": newick("((NC_008030.,B),(C,D));"),
        }
        labels = {name: PipelineLabel.parse(name, prefix="") for name in trees}
        tree_set = TreeSet(trees, labels)
        self.assertEqual(tree_set.n_taxa, 4)
        self.assertEqual(tree_set.mangled_labels["mafft_iqtree"], ["NC_008030"])
        self.assertEqual(tree_set.mangled_labels["mafft_fasttree"], [])


class TestStabilityAnalyzer(unittest.TestCase):
    """Suporte, padrões maximais e efeitos de fator."""

    def setUp(self):
        # Três pipelines: dois concordam totalmente, o terceiro discorda em (C,D).
        trees = {
            "mafft_fasttree": newick("(((A,B),(C,D)),E);"),
            "clustalo_fasttree": newick("(((B,A),(D,C)),E);"),
            "mafft_nj_distance": newick("(((A,B),(C,E)),D);"),
        }
        labels = {name: PipelineLabel.parse(name, prefix="") for name in trees}
        self.analyzer = StabilityAnalyzer(TreeSet(trees, labels))

    def test_support_counts_pipelines_not_subtrees(self):
        records = {record.taxa: record for record in self.analyzer.clade_records()}
        self.assertEqual(len(records[frozenset({"A", "B"})].pipelines), 3)
        self.assertEqual(len(records[frozenset({"C", "D"})].pipelines), 2)

    def test_support_profile_is_cumulative(self):
        profile = self.analyzer.support_profile()
        self.assertGreaterEqual(profile[1]["cumulative"], profile[3]["cumulative"])
        self.assertEqual(profile[3]["cumulative"], len(self.analyzer.consensus_clades(1.0)))

    def test_maximal_patterns_are_maximal(self):
        patterns = self.analyzer.maximal_patterns(min_support=0.5)
        self.assertTrue(patterns)
        for a in patterns:
            for b in patterns:
                if a is not b:
                    self.assertFalse(set(a.clades) < set(b.clades))

    def test_identical_topologies_have_zero_rf(self):
        matrix = self.analyzer.rf_matrix(normalized=False)
        self.assertEqual(matrix["mafft_fasttree"]["clustalo_fasttree"], 0)
        self.assertGreater(matrix["mafft_fasttree"]["mafft_nj_distance"], 0)

    def test_factor_effects_separate_aligner_from_inference(self):
        effects = self.analyzer.factor_effects(normalized=False)
        self.assertEqual(effects["aligner"]["mean"], 0.0)
        self.assertGreater(effects["inference"]["mean"], 0.0)

    def test_taxonomic_purity_flags_single_label_clades(self):
        rows = self.analyzer.taxonomic_purity(lambda t: "X" if t in "AB" else "Y", min_size=2)
        pure = [row for row in rows if row["monophyletic_group"]]
        self.assertTrue(any(row["purity"] == 1.0 for row in pure))


if __name__ == "__main__":
    unittest.main()
