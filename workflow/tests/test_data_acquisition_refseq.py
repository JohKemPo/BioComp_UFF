"""
Preferência RefSeq/GenBank na aquisição — [D23](../../../docs/science/02-defeitos-que-alteram-resultado.md#d23),
implementação de [DEC-082](../../../docs/automation/07-log-de-execucao.md).

`workflow_dataAcquisition.py.filter_sequences`/`refine_dataset` deduplicavam
por conteúdo guardando a **primeira ocorrência no arquivo** — sem preferência
declarada entre RefSeq e GenBank, o mesmo genoma virava um acesso diferente
conforme a ordem de download (a mesma classe de defeito caracterizada em
`workflow.tests.test_deduplicacao` para `dataValidation.deduplicar_por_sequencia`,
que é outro caminho — o do `TreeBuilderController` — e agora recebe a mesma
preferência, ver `test_deduplicacao.test_refseq_decide_o_sobrevivente_independente_da_ordem`).

O que estes testes travam:

* **Ordem deixa de decidir a identidade.** As duas ordens de chegada do
  mesmo par produzem o MESMO sobrevivente agora — o acesso RefSeq — o que
  antes de DEC-082 não valia (`test_deduplicacao.
  test_ordem_do_arquivo_decide_o_sobrevivente` mostra o comportamento antigo
  no outro caminho).
* **A posição não muda.** O conjunto não é reordenado — só o rótulo na
  posição já ocupada pelo sobrevivente troca, quando o recém-chegado é
  RefSeq e o que já estava lá não é.
* **Regressão do `UnboundLocalError`.** `filter_sequences` tinha `del filtered`
  antes de `return filtered`, o que sempre lançava, sempre caía no `except`
  mais externo, e sempre reescrevia `output_file` como vazio — nenhuma
  filtragem por `initial_min_length` jamais sobrevivia. Achado incidental
  ao mexer nestas linhas para D23; corrigido junto.

Executar com:

    python -m unittest workflow.tests.test_data_acquisition_refseq
"""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from workflow.workflow_dataAcquisition import workflowAquisitionDatasetNCBI

# Genoma de brinquedo, longo o bastante para passar em filtros de
# comprimento mínimo pequenos usados aqui. O que se testa é identidade de
# conteúdo, não o alinhamento.
_SEQ_TATERAPOX = Seq("ACGTACGTACGTAAGGCCTT" * 40)
_SEQ_CAMELPOX = Seq("TTGGCCAAGGTTACGTACGT" * 40)
_SEQ_OUTRA = Seq("GGGGCCCCAAAATTTTGGCC" * 40)


def _registro(acesso, descricao, seq):
    return SeqRecord(seq, id=acesso, name=acesso.split(".")[0],
                      description=descricao, annotations={"molecule_type": "DNA"})


def _escrever_gb(destino, registros):
    SeqIO.write(registros, destino, "genbank")
    return destino


class _WorkflowDeTeste(unittest.TestCase):
    """Base: instancia o workflow num diretório temporário e limpa depois."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="d23_aquisicao_")
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)
        self.workflow = workflowAquisitionDatasetNCBI(
            email="mail@mail.com",
            work_dir=self.dir,
            initial_min_length=10,
            refined_min_length=10,
        )

    def _entrada(self, nome, registros):
        return _escrever_gb(os.path.join(self.dir, nome), registros)


class TestFilterSequencesRefSeq(_WorkflowDeTeste):
    def test_nao_lanca_e_nao_apaga_a_saida(self):
        """Regressão do `UnboundLocalError` — a função tinha de sobreviver ao
        próprio `return`, e o arquivo escrito não podia ser sobrescrito vazio
        pelo `except` externo."""
        entrada = self._entrada("raw.gb", [
            _registro("DQ437594.1", "Taterapox GenBank", _SEQ_TATERAPOX),
        ])
        saida = os.path.join(self.dir, "filtered.gb")

        resultado = self.workflow.filter_sequences(entrada, saida)

        self.assertEqual(len(resultado), 1)
        self.assertEqual(os.path.getsize(saida) > 0, True,
                          "output_file não pode sair vazio quando há sequência válida")
        self.assertEqual([r.id for r in SeqIO.parse(saida, "genbank")], ["DQ437594.1"])

    def test_refseq_prevalece_independente_da_ordem(self):
        """O cerne de DEC-082: as duas ordens convergem para o mesmo acesso."""
        saida_a = os.path.join(self.dir, "a.gb")
        self.workflow.filter_sequences(
            self._entrada("genbank_primeiro.gb", [
                _registro("DQ437594.1", "Taterapox GenBank", _SEQ_TATERAPOX),
                _registro("NC_008291.1", "Taterapox RefSeq", _SEQ_TATERAPOX),
            ]),
            saida_a,
        )

        saida_b = os.path.join(self.dir, "b.gb")
        self.workflow.filter_sequences(
            self._entrada("refseq_primeiro.gb", [
                _registro("NC_008291.1", "Taterapox RefSeq", _SEQ_TATERAPOX),
                _registro("DQ437594.1", "Taterapox GenBank", _SEQ_TATERAPOX),
            ]),
            saida_b,
        )

        ids_a = [r.id for r in SeqIO.parse(saida_a, "genbank")]
        ids_b = [r.id for r in SeqIO.parse(saida_b, "genbank")]
        self.assertEqual(ids_a, ["NC_008291.1"])
        self.assertEqual(ids_b, ["NC_008291.1"])
        self.assertEqual(ids_a, ids_b,
                          "a identidade do sobrevivente não pode depender da ordem de chegada")

    def test_posicao_nao_muda_so_o_rotulo(self):
        """RefSeq chega depois: a posição do grupo no conjunto (índice 1, entre
        Camelpox e a sequência sem par) não muda — só o rótulo naquela posição."""
        saida = os.path.join(self.dir, "posicao.gb")
        self.workflow.filter_sequences(
            self._entrada("posicao.gb_in", [
                _registro("AF438165.1", "Camelpox GenBank", _SEQ_CAMELPOX),
                _registro("DQ437594.1", "Taterapox GenBank", _SEQ_TATERAPOX),
                _registro("KP123456.1", "outro genoma", _SEQ_OUTRA),
                _registro("NC_008291.1", "Taterapox RefSeq", _SEQ_TATERAPOX),
            ]),
            saida,
        )

        ids = [r.id for r in SeqIO.parse(saida, "genbank")]
        self.assertEqual(ids, ["AF438165.1", "NC_008291.1", "KP123456.1"],
                          "Taterapox continua na posição 1 (entre Camelpox e o genoma "
                          "solto); só o acesso ali muda de DQ437594.1 para NC_008291.1")

    def test_sem_refseq_no_par_mantem_a_primeira_ocorrencia(self):
        """Sem RefSeq envolvido, nada para preferir — comportamento por ordem
        de chegada continua valendo, como sempre valeu para duplicata comum."""
        saida = os.path.join(self.dir, "sem_refseq.gb")
        self.workflow.filter_sequences(
            self._entrada("sem_refseq.gb_in", [
                _registro("DQ437594.1", "primeira submissão GenBank", _SEQ_TATERAPOX),
                _registro("DQ437595.1", "segunda submissão GenBank", _SEQ_TATERAPOX),
            ]),
            saida,
        )

        ids = [r.id for r in SeqIO.parse(saida, "genbank")]
        self.assertEqual(ids, ["DQ437594.1"])

    def test_log_da_conta_da_substituicao(self):
        entrada = self._entrada("log.gb", [
            _registro("DQ437594.1", "Taterapox GenBank", _SEQ_TATERAPOX),
            _registro("NC_008291.1", "Taterapox RefSeq", _SEQ_TATERAPOX),
        ])
        saida = os.path.join(self.dir, "log_saida.gb")

        with self.assertLogs("ZikaWorkflow", level="INFO") as capturado:
            self.workflow.filter_sequences(entrada, saida)

        mensagem = "\n".join(capturado.output)
        self.assertIn("D23/DEC-082", mensagem)
        self.assertIn("NC_008291.1", mensagem)
        self.assertIn("DQ437594.1", mensagem)


class TestRefineDatasetRefSeq(_WorkflowDeTeste):
    def test_refseq_prevalece_independente_da_ordem(self):
        saida_a = os.path.join(self.dir, "ra.gb")
        self.workflow.refine_dataset(
            self._entrada("refine_genbank_primeiro.gb", [
                _registro("AF438165.1", "Camelpox GenBank", _SEQ_CAMELPOX),
                _registro("NC_003391.1", "Camelpox RefSeq", _SEQ_CAMELPOX),
            ]),
            saida_a,
        )

        saida_b = os.path.join(self.dir, "rb.gb")
        self.workflow.refine_dataset(
            self._entrada("refine_refseq_primeiro.gb", [
                _registro("NC_003391.1", "Camelpox RefSeq", _SEQ_CAMELPOX),
                _registro("AF438165.1", "Camelpox GenBank", _SEQ_CAMELPOX),
            ]),
            saida_b,
        )

        ids_a = [r.id for r in SeqIO.parse(saida_a, "genbank")]
        ids_b = [r.id for r in SeqIO.parse(saida_b, "genbank")]
        self.assertEqual(ids_a, ["NC_003391.1"])
        self.assertEqual(ids_a, ids_b)

    def test_conjunto_sem_duplicata_nao_muda(self):
        saida = os.path.join(self.dir, "limpo.gb")
        self.workflow.refine_dataset(
            self._entrada("limpo.gb_in", [
                _registro("A.1", "alfa", _SEQ_TATERAPOX),
                _registro("B.1", "beta", _SEQ_CAMELPOX),
                _registro("C.1", "gama", _SEQ_OUTRA),
            ]),
            saida,
        )

        ids = [r.id for r in SeqIO.parse(saida, "genbank")]
        self.assertEqual(ids, ["A.1", "B.1", "C.1"])


if __name__ == "__main__":
    unittest.main()
