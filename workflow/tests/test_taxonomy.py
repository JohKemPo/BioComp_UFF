"""
Filtro taxonômico declarado e verificação pós-download — M2.2 / D6.

A aquisição nunca restringiu o táxon, e o conjunto foi contaminado: *Nile
crocodilepox* está em VARV-6, VARV-52 e VARV-121; *Yokapox* e dois *Saltwater
crocodilepox* estão em VARV-121. Nenhum é *Orthopoxvirus*.

O que estes testes fixam:

* o clado é **declarado**, nunca presumido — sem declaração, nada é conferido,
  e isso aparece como aviso e não como silêncio;
* a verificação é **offline** e funciona nos conjuntos que já existem;
* **ausência de linhagem não é ausência de pertencimento** — um registro sem
  anotação é indecidível, não reprovado.

Executar com:

    python -m unittest workflow.tests.test_taxonomy
"""

from __future__ import annotations

import unittest

from workflow.utils.taxonomy import (ORTHOFLAVIVIRUS, ORTHOPOXVIRUS, TaxonFilter,
                                     TaxonomyAudit, audit_records, entrez_term,
                                     record_lineage, within_taxon)


class RegistroFalso:
    """Um SeqRecord no que importa aqui: id e anotações."""

    def __init__(self, acesso, organismo=None, linhagem=None):
        self.id = acesso
        self.annotations = {}
        if organismo is not None:
            self.annotations["organism"] = organismo
        if linhagem is not None:
            self.annotations["taxonomy"] = list(linhagem)


ORTHO = ["Viruses", "Poxviridae", "Chordopoxvirinae", "Orthopoxvirus"]
CROCO = ["Viruses", "Poxviridae", "Chordopoxvirinae", "Crocodylidpoxvirus"]


class TestTermoDeBusca(unittest.TestCase):

    def test_injeta_o_filtro_declarado(self):
        termo = entrez_term("Variola virus[Organism]", ORTHOPOXVIRUS)
        self.assertIn("txid10242[Organism:exp]", termo)
        self.assertIn("Variola virus[Organism]", termo)

    def test_sem_filtro_a_consulta_passa_intacta(self):
        """`None` é uma escolha legítima — desde que registrada. O que não pode
        é o código decidir um clado por conta própria."""
        self.assertEqual(entrez_term("Variola virus[Organism]"), "Variola virus[Organism]")

    def test_nao_duplica_filtro_ja_escrito_a_mao(self):
        consulta = "Variola virus[Organism] AND txid10242[Organism:exp]"
        self.assertEqual(entrez_term(consulta, ORTHOPOXVIRUS).count("txid10242"), 1)

    def test_consulta_vazia_continua_vazia(self):
        self.assertEqual(entrez_term("", ORTHOPOXVIRUS), "")


class TestPertencimento(unittest.TestCase):

    def test_dentro_do_clado(self):
        self.assertTrue(within_taxon(RegistroFalso("X1", "Variola virus", ORTHO), ORTHOPOXVIRUS))

    def test_fora_do_clado(self):
        self.assertFalse(within_taxon(
            RegistroFalso("NC_008030", "Nile crocodilepox virus", CROCO), ORTHOPOXVIRUS))

    def test_sem_linhagem_e_indecidivel(self):
        """Nem `True` nem `False`: devolver `False` descartaria o registro por
        falta de metadado, e `True` o aceitaria pelo mesmo motivo."""
        self.assertIsNone(within_taxon(RegistroFalso("X2", "Alguma coisa"), ORTHOPOXVIRUS))

    def test_linhagem_e_lida_do_registro(self):
        self.assertEqual(record_lineage(RegistroFalso("X", "V", ORTHO))[-1], "Orthopoxvirus")
        self.assertEqual(record_lineage(RegistroFalso("X")), ())


class TestAuditoria(unittest.TestCase):

    def setUp(self):
        self.registros = [
            RegistroFalso("DQ437591.1", "Variola virus", ORTHO),
            RegistroFalso("NC_001611.1", "Variola virus", ORTHO),
            RegistroFalso("NC_008291.1", "Taterapox virus", ORTHO),
            RegistroFalso("NC_008030.1", "Nile crocodilepox virus", CROCO),
            RegistroFalso("XX000000.1", "Sem linhagem"),
        ]

    def test_separa_dentro_fora_e_indecidivel(self):
        a = audit_records(self.registros, ORTHOPOXVIRUS)
        self.assertEqual(len(a.within), 3)
        self.assertEqual(set(a.outside), {"NC_008030"})
        self.assertEqual(a.unknown, ["XX000000"])
        self.assertEqual(a.total, 5)

    def test_conjunto_com_taxon_fora_nao_e_limpo(self):
        self.assertFalse(audit_records(self.registros, ORTHOPOXVIRUS).clean)

    def test_indecidivel_tambem_impede_declarar_limpo(self):
        """Um registro sem linhagem pode ser contaminação não detectada. Chamar
        o conjunto de limpo com ele dentro seria afirmar mais do que se sabe."""
        so_indecidivel = [RegistroFalso("A.1", "V", ORTHO), RegistroFalso("B.1", "?")]
        self.assertFalse(audit_records(so_indecidivel, ORTHOPOXVIRUS).clean)

    def test_conjunto_limpo(self):
        limpos = [r for r in self.registros if "Orthopoxvirus" in (r.annotations.get("taxonomy") or [])]
        self.assertTrue(audit_records(limpos, ORTHOPOXVIRUS).clean)

    def test_acesso_repetido_conta_uma_vez(self):
        """O `raw_data_sequences.gb` guarda o mesmo registro uma vez por árvore."""
        a = audit_records(self.registros + self.registros, ORTHOPOXVIRUS)
        self.assertEqual(a.total, 5)

    def test_versao_do_acesso_e_normalizada(self):
        a = audit_records([RegistroFalso("NC_008030.1", "Nile crocodilepox virus", CROCO)],
                          ORTHOPOXVIRUS)
        self.assertIn("NC_008030", a.outside)

    def test_modo_estrito_levanta_com_o_motivo(self):
        a = audit_records(self.registros, ORTHOPOXVIRUS)
        with self.assertRaises(ValueError) as erro:
            a.raise_if_contaminated()
        self.assertIn("NC_008030", str(erro.exception))
        self.assertIn("Nile crocodilepox", str(erro.exception))

    def test_modo_estrito_nao_levanta_por_indecidivel(self):
        """Indecidível é aviso, não bloqueio: bloquear ali descartaria registro
        por falta de metadado."""
        a = audit_records([RegistroFalso("A.1", "V", ORTHO), RegistroFalso("B.1", "?")],
                          ORTHOPOXVIRUS)
        a.raise_if_contaminated()   # não levanta

    def test_resumo_e_serializavel(self):
        import json
        resumo = audit_records(self.registros, ORTHOPOXVIRUS).summary()
        json.dumps(resumo)
        self.assertEqual(resumo["taxon"]["taxid"], "txid10242")
        self.assertEqual(resumo["outside"]["NC_008030"], "Nile crocodilepox virus")


class TestOClaudeEDoExperimento(unittest.TestCase):
    """O filtro é do experimento, não do projeto."""

    def test_zika_reprovado_contra_orthopoxvirus(self):
        zika = [RegistroFalso("KF383116.1", "Zika virus",
                              ["Viruses", "Flaviviridae", "Orthoflavivirus"])]
        self.assertFalse(audit_records(zika, ORTHOPOXVIRUS).clean)

    def test_zika_aprovado_contra_o_proprio_genero(self):
        zika = [RegistroFalso("KF383116.1", "Zika virus",
                              ["Viruses", "Flaviviridae", "Orthoflavivirus"])]
        self.assertTrue(audit_records(zika, ORTHOFLAVIVIRUS).clean)

    def test_clado_arbitrario_pode_ser_declarado(self):
        proprio = TaxonFilter(taxid="txid11320", name="Alphainfluenzavirus")
        self.assertIn("txid11320[Organism:exp]", entrez_term("influenza", proprio))


if __name__ == "__main__":
    unittest.main()
