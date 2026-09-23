"""
`TreeBuilderController._reproducibility_kwargs` — [D26](../../../docs/science/02-defeitos-que-alteram-resultado.md#d26).

Os quatro métodos avançados (`build_tree_{iqtree,fasttree,raxml,mrbayes}`)
instanciavam `TreeBuilder` sem repassar `random_seed`/`raxml_threads`/
`iqtree_threads` — `reproducibility_settings` (`builder.py`) sempre caía nos
defaults do módulo (12345/4/4), não importa o que o `tree_config` do
experimento pedisse. `workflow.py` grava no manifesto o valor **pedido**
(`reproducibility_settings(params['tree_config'])`); a chamada real usava
outro. Confirmado em `Variola_VARV49_reexec_20260901/out/outputs/manifest.json`:
declarado `raxml_threads: 8, iqtree_threads: 16`, executado `4/4` (defaults).

A lição de `test_deduplicacao` vale aqui: um teste que exercita
`reproducibility_settings`/`TreeBuilder` direto (já existe em
`test_manifest.test_reproducibilidade_vem_da_mesma_fonte_do_builder`) não
prova que o **caminho do controlador** os alcança — é exatamente essa lacuna
que D26 documentava. Este arquivo testa `_reproducibility_kwargs`, o método
que fecha esse caminho.

Executar com:

    python -m unittest workflow.tests.test_reproducibility_kwargs
"""

from __future__ import annotations

import unittest

from workflow.controller.treeBuilderController import TreeBuilderController
from workflow.tree_construction.builder import (REPRODUCIBILITY_DEFAULTS,
                                                reproducibility_settings)


def _controlador(**tree_config):
    """`__new__`, não `__init__`: nada de I/O, diretórios ou logger — só o
    que `_reproducibility_kwargs` lê (mesmo padrão de
    `test_deduplicacao.TestCaminhoDoControlador._controlador`)."""
    controlador = TreeBuilderController.__new__(TreeBuilderController)
    for chave, valor in tree_config.items():
        setattr(controlador, chave, valor)
    return controlador


class TestReproducibilityKwargs(unittest.TestCase):
    def test_repassa_o_que_o_tree_config_pediu(self):
        """D26 — o caso medido: `raxml_threads`/`iqtree_threads` diferentes
        do default chegam ao `TreeBuilder`, não são substituídos por 4/4."""
        controlador = _controlador(random_seed=777, raxml_threads=8, iqtree_threads=16)

        kwargs = controlador._reproducibility_kwargs()

        self.assertEqual(kwargs, {"random_seed": 777, "raxml_threads": 8, "iqtree_threads": 16})

    def test_sem_tree_config_nao_inventa_chave(self):
        """Sem `random_seed`/`raxml_threads`/`iqtree_threads` no `tree_config`
        do experimento, `_reproducibility_kwargs` não inclui a chave — quem
        decide o default é `reproducibility_settings`, uma fonte só."""
        controlador = _controlador()

        kwargs = controlador._reproducibility_kwargs()

        self.assertEqual(kwargs, {})
        # A mesma função que o manifesto usa resolve os defaults do módulo
        # quando a chave está ausente — não duplicado aqui.
        self.assertEqual(reproducibility_settings(kwargs), REPRODUCIBILITY_DEFAULTS)

    def test_parcial_nao_quebra_com_none(self):
        """Caso-limite: só parte do `tree_config` declarada. `reproducibility_
        settings` faz `config.get(chave, padrao)` — uma chave presente com
        valor `None` quebraria em `int(None)`; `_reproducibility_kwargs` não
        pode introduzir esse valor."""
        controlador = _controlador(random_seed=42)

        kwargs = controlador._reproducibility_kwargs()

        self.assertEqual(kwargs, {"random_seed": 42})
        resolvido = reproducibility_settings(kwargs)
        self.assertEqual(resolvido["random_seed"], 42)
        self.assertEqual(resolvido["raxml_threads"], REPRODUCIBILITY_DEFAULTS["raxml_threads"])
        self.assertEqual(resolvido["iqtree_threads"], REPRODUCIBILITY_DEFAULTS["iqtree_threads"])


if __name__ == "__main__":
    unittest.main()
