"""
O `akoflow/workflow.yaml` tem o formato do exemplo real do AkôFlow e é executável
pelo parser da CLI (G4b). Referência do formato: `UFFeScience/akoflow-examples`
@ d890c3171e36, `machine-learning/infrastructure-selection/workflow.yaml`.
"""

import graphlib
import json
import os
import unittest

import yaml

from workflow import akoflow

RAIZ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAMPOS_SPEC = {"activities", "image", "namespace"}          # do exemplo
CAMPOS_ATIVIDADE = {"command", "cpuLimit", "dependsOn", "memoryLimit", "name", "runtime"}


class TestWorkflowYaml(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(os.path.join(RAIZ, "akoflow", "workflow.yaml"), encoding="utf-8") as fh:
            cls.wf = yaml.safe_load(fh)
        with open(os.path.join(RAIZ, "akoflow", "config", "zika21.json"), encoding="utf-8") as fh:
            cls.cfg = json.load(fh)
        cls.ativ = {a["name"]: a for a in cls.wf["spec"]["activities"]}

    def test_campos_sao_os_do_exemplo(self):
        self.assertEqual(set(self.wf), {"name", "spec"})
        # `dataDependencies` vem do exemplo montage-050 (mesmo commit)
        self.assertLessEqual(set(self.wf["spec"]), CAMPOS_SPEC | {"dataDependencies"})
        for a in self.ativ.values():
            self.assertLessEqual(set(a), CAMPOS_ATIVIDADE, a["name"])
            self.assertEqual(set(a["command"]), {"arguments", "entrypoint"})

    def test_nomes_unicos_dag_aciclico_e_dependencias_existentes(self):
        self.assertEqual(len(self.ativ), len(self.wf["spec"]["activities"]))
        grafo = {n: set(a.get("dependsOn", [])) for n, a in self.ativ.items()}
        for n, deps in grafo.items():
            self.assertLessEqual(deps, set(self.ativ), n)
        list(graphlib.TopologicalSorter(grafo).static_order())      # levanta CycleError se houver

    def test_data_dependencies_correspondem_a_depends_on(self):
        for d in self.wf["spec"]["dataDependencies"]:
            self.assertIn(d["producerActivity"], self.ativ[d["consumerActivity"]]["dependsOn"])

    def test_cada_comando_e_aceito_pela_cli(self):
        parser = akoflow._parser()
        for nome, a in self.ativ.items():
            self.assertEqual(a["command"]["entrypoint"], "python")
            args = a["command"]["arguments"]
            self.assertEqual(args[:2], ["-m", "workflow.akoflow"], nome)
            ns = parser.parse_args(args[2:])
            self.assertEqual(ns.out, "out")

    def test_cpu_e_o_do_config_e_d17_vale(self):
        tree = self.cfg["tree_config"]
        for nome, a in self.ativ.items():
            esperado = (tree["num_threads"] if nome.startswith("alinhar-")
                        else tree["raxml_threads"] if nome.startswith("arvore-") and nome.endswith("-raxml") else 1)
            self.assertEqual(a["cpuLimit"], str(esperado), nome)
        self.assertEqual(tree["iqtree_threads"], 4)   # o IQ-TREE roda com `-nt 1` (D21) pelo builder

    def test_cada_arvore_tem_sua_extracao_e_a_mineracao_depende_de_todas(self):
        arvores = [n for n in self.ativ if n.startswith("arvore-")]
        extracoes = {n for n in self.ativ if n.startswith("subarvores-")}
        self.assertEqual(len(arvores), len(extracoes))
        self.assertEqual(set(self.ativ["minerar-fpmax"]["dependsOn"]), extracoes)
        self.assertEqual(len(self.ativ), 39)


if __name__ == "__main__":
    unittest.main()
