"""
Testes da deduplicação por sequência — [D23].

O defeito que estes testes existem para impedir é de **silêncio**. A função
chamava-se `remove_pipe`, não removia pipe nenhum, deduplicava por conteúdo e
descartava registros sem dizer quais. Nos conjuntos de *Variola* isso significa
que o projeto chamado VARV-49 recebe 52 registros e produz árvores com 49
folhas — e `n` é um número de *Methods*.

A lição de `test_tool_runs` vale aqui: **um teste que exercita a função e não o
caminho não prova que o caminho existe.** Por isso há, além dos testes de
unidade, um que entra pelo `_validate_and_prepare_fasta` do controlador e
confere que o descarte chegou ao registro que alimenta o manifesto.

Os testes de sobrevivência (`test_ordem_do_arquivo_decide_o_sobrevivente`,
`test_forma_do_par_variola`) são de **caracterização**: eles fixam o
comportamento arbitrário de hoje — primeira ocorrência vence — para que a
escolha de preferência entre RefSeq e GenBank, quando o usuário a tomar, tenha
de ser feita alterando um teste deliberadamente, e não por acidente.

Executar com:

    python -m unittest workflow.tests.test_deduplicacao
"""

from __future__ import annotations

import os
import shutil
import tempfile
import unittest

from Bio import SeqIO

from workflow.utils import tool_runs
from workflow.utils.dataValidation import deduplicar_por_sequencia, remove_pipe

# Genoma de brinquedo. Curto de propósito: o que se testa é a identidade de
# conteúdo, não o alinhamento.
_SEQ_TATERAPOX = "ACGTACGTACGTAAGGCCTT"
_SEQ_CAMELPOX = "TTGGCCAAGGTTACGTACGT"
_SEQ_OUTRA = "GGGGCCCCAAAATTTTGGCC"


def _escrever_fasta(destino, registros):
    """Escreve `[(id, descricao, sequencia), ...]` como FASTA."""
    with open(destino, "w") as fh:
        for ident, descricao, seq in registros:
            fh.write(f">{ident} {descricao}\n{seq}\n")
    return destino


def _ids(caminho):
    return [r.id for r in SeqIO.parse(caminho, "fasta")]


class TestDeduplicacaoPorSequencia(unittest.TestCase):
    """Unidade: o que a função devolve, grava e avisa."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="d23_")
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)

    def _rodar(self, registros, nome="conjunto"):
        entrada = _escrever_fasta(os.path.join(self.dir, f"{nome}.fasta"), registros)
        return deduplicar_por_sequencia(nome, entrada, self.dir)

    def test_duas_sequencias_identicas_com_rotulos_diferentes(self):
        """O caso mínimo de D23: RefSeq e GenBank do mesmo genoma."""
        saida, descartados = self._rodar([
            ("DQ437594.1", "Taterapox virus strain Dahomey 1968", _SEQ_TATERAPOX),
            ("NC_008291.1", "Taterapox virus", _SEQ_TATERAPOX),
        ])

        self.assertEqual(descartados, [("NC_008291.1", "DQ437594.1")],
                         "o descarte precisa dizer quem saiu e em favor de quem")
        self.assertEqual(_ids(saida), ["DQ437594.1"],
                         "o FASTA deduplicado guarda uma entrada por sequência distinta")

    def test_descarte_e_anunciado_no_log_com_os_dois_acessos(self):
        """Descarte silencioso é o defeito; o WARNING é a correção."""
        with self.assertLogs(level="WARNING") as capturado:
            self._rodar([
                ("DQ437594.1", "Taterapox virus strain Dahomey 1968", _SEQ_TATERAPOX),
                ("NC_008291.1", "Taterapox virus", _SEQ_TATERAPOX),
            ])

        mensagem = "\n".join(capturado.output)
        self.assertIn("NC_008291.1", mensagem)
        self.assertIn("DQ437594.1", mensagem)
        self.assertIn("2 registros", mensagem)
        self.assertIn("1 sequências distintas", mensagem)

    def test_acesso_repetido_nao_e_anunciado_como_par_de_acessos(self):
        """`descartado == mantido` é registro repetido, não gêmeo RefSeq/GenBank.

        Os dois fenômenos ocorrem juntos em VARV-49 e VARV-121, e escrever o
        primeiro como "X (idêntico a X)" lê como contradição.
        """
        with self.assertLogs(level="WARNING") as capturado:
            self._rodar([
                ("DQ437594.1", "Taterapox virus strain Dahomey 1968", _SEQ_TATERAPOX),
                ("DQ437594.1", "Taterapox virus strain Dahomey 1968", _SEQ_TATERAPOX),
            ])

        mensagem = "\n".join(capturado.output)
        self.assertIn("registro repetido do mesmo acesso", mensagem)
        self.assertNotIn("acesso distinto", mensagem)

    def test_sem_duplicata_nao_descarta_nada_e_nao_avisa(self):
        """Conjunto limpo — os oito de Zika — não pode gerar ruído nem perda."""
        saida, descartados = self._rodar([
            ("A.1", "alfa", _SEQ_TATERAPOX),
            ("B.1", "beta", _SEQ_CAMELPOX),
            ("C.1", "gama", _SEQ_OUTRA),
        ])

        self.assertEqual(descartados, [])
        self.assertEqual(_ids(saida), ["A.1", "B.1", "C.1"])

    def test_ordem_do_arquivo_decide_o_sobrevivente(self):
        """CARACTERIZAÇÃO — o comportamento arbitrário de hoje, fixado.

        Os mesmos dois registros em ordem invertida produzem o **outro**
        sobrevivente. É por isso que o mesmo táxon aparece como `DQ437594.1` em
        VARV-49 e como `NC_008291.1` em VARV-121, e é a pergunta de curadoria
        que aguarda decisão do usuário. Quando ela for tomada, este teste tem de
        ser reescrito de propósito.
        """
        _, genbank_primeiro = self._rodar([
            ("DQ437594.1", "GenBank", _SEQ_TATERAPOX),
            ("NC_008291.1", "RefSeq", _SEQ_TATERAPOX),
        ], nome="ordem_a")
        _, refseq_primeiro = self._rodar([
            ("NC_008291.1", "RefSeq", _SEQ_TATERAPOX),
            ("DQ437594.1", "GenBank", _SEQ_TATERAPOX),
        ], nome="ordem_b")

        self.assertEqual(genbank_primeiro, [("NC_008291.1", "DQ437594.1")])
        self.assertEqual(refseq_primeiro, [("DQ437594.1", "NC_008291.1")])

    def test_forma_do_par_variola(self):
        """CARACTERIZAÇÃO — a forma medida em `data/replication-RetMax200-ITRs`.

        Um grupo de três (acesso repetido **e** gêmeo RefSeq) e um grupo de
        dois, sobre um conjunto que também tem sequências únicas: 4 registros
        de Taterapox/Camelpox + 1 distinta → 3 sequências distintas.
        """
        saida, descartados = self._rodar([
            ("DQ437594.1", "Taterapox GenBank", _SEQ_TATERAPOX),
            ("AF438165.1", "Camelpox GenBank", _SEQ_CAMELPOX),
            ("NC_008291.1", "Taterapox RefSeq", _SEQ_TATERAPOX),
            ("NC_003391.1", "Camelpox RefSeq", _SEQ_CAMELPOX),
            ("DQ437594.1", "Taterapox GenBank de novo", _SEQ_TATERAPOX),
            ("KP123456.1", "outro genoma", _SEQ_OUTRA),
        ], nome="variola")

        self.assertEqual(descartados, [
            ("NC_008291.1", "DQ437594.1"),
            ("NC_003391.1", "AF438165.1"),
            ("DQ437594.1", "DQ437594.1"),
        ])
        self.assertEqual(_ids(saida), ["DQ437594.1", "AF438165.1", "KP123456.1"])

    def test_chave_e_sensivel_a_caixa(self):
        """CARACTERIZAÇÃO da escolha silenciosa nº 1 da docstring.

        A chave é `str(seq)` cru. Soft-masking minúsculo faz o mesmo genoma
        entrar duas vezes na árvore. Não se corrige aqui: mudar a chave muda a
        composição do conjunto e é decisão do usuário.
        """
        _, descartados = self._rodar([
            ("X.1", "maiúscula", _SEQ_TATERAPOX),
            ("Y.1", "minúscula", _SEQ_TATERAPOX.lower()),
        ], nome="caixa")

        self.assertEqual(descartados, [],
                         "hoje a caixa distingue; se isto falhar, a composição mudou")

    def test_nome_antigo_continua_devolvendo_so_o_caminho(self):
        """`remove_pipe` é casca de compatibilidade e não pode voltar a mentir."""
        entrada = _escrever_fasta(os.path.join(self.dir, "compat.fasta"), [
            ("A.1", "alfa", _SEQ_TATERAPOX),
            ("B.1", "beta", _SEQ_TATERAPOX),
        ])
        caminho = remove_pipe("compat", entrada, self.dir)

        self.assertIsInstance(caminho, str)
        self.assertEqual(_ids(caminho), ["A.1"])


class TestCaminhoDoControlador(unittest.TestCase):
    """O caminho, não a função: o descarte chega ao registro do manifesto."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="d23_ctrl_")
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)
        tool_runs.limpar()
        self.addCleanup(tool_runs.limpar)

    def _controlador(self):
        # `__init__` abre log, cria diretórios e instancia o alinhador; nada
        # disso participa da validação de FASTA. O que se testa é o corpo de
        # `_validate_and_prepare_fasta`, que é quem chama a deduplicação.
        from workflow.controller.treeBuilderController import TreeBuilderController

        controlador = TreeBuilderController.__new__(TreeBuilderController)
        controlador.input_path = self.dir
        return controlador

    def test_descartados_entram_em_tool_runs(self):
        entrada = _escrever_fasta(os.path.join(self.dir, "dataset_final.fasta"), [
            ("DQ437594.1", "Taterapox GenBank", _SEQ_TATERAPOX),
            ("NC_008291.1", "Taterapox RefSeq", _SEQ_TATERAPOX),
            ("KP123456.1", "outro genoma", _SEQ_OUTRA),
        ])

        self._controlador()._validate_and_prepare_fasta(
            "dataset_final.fasta", entrada, "dataset_final")

        registro = tool_runs.execucoes().get("deduplicacao")
        self.assertIsNotNone(
            registro, "sem registro, o manifesto sai com o descarte invisível — é D23")
        self.assertEqual(registro["descartados"],
                         [{"descartado": "NC_008291.1", "mantido": "DQ437594.1"}])
        self.assertIn("D23", registro["nota"])
        self.assertEqual(registro["runs"][0]["command"][0], "deduplicar_por_sequencia")

    def test_conjunto_limpo_nao_registra_deduplicacao(self):
        """Sem duplicata, nenhum ruído entra no manifesto."""
        entrada = _escrever_fasta(os.path.join(self.dir, "limpo.fasta"), [
            ("A.1", "alfa", _SEQ_TATERAPOX),
            ("B.1", "beta", _SEQ_CAMELPOX),
        ])

        self._controlador()._validate_and_prepare_fasta("limpo.fasta", entrada, "limpo")

        self.assertNotIn("deduplicacao", tool_runs.execucoes())


if __name__ == "__main__":
    unittest.main()
