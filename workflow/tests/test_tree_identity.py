"""
Testes da identidade de clado usada em produção — D5.

O defeito: `encode_list_to_int` calculava MD5 de 16 bits sobre a *representação
textual* da lista de hashes na ordem de travessia. Isso fragmentava o mesmo clado
em vários itens quando duas árvores ordenavam os filhos de modo distinto
(**subestimando** o suporte) e fundia clados distintos por colisão em 65 536
valores (**fabricando** suporte). Em VARV-49, 155 itens legados para 100 clados
reais — 40,6% fragmentados.

A identidade correta já existia em `workflow.stability.clade_identity` e não era
usada pelo pipeline. Estes testes fixam que agora é.

Executar com:

    python -m unittest workflow.tests.test_tree_identity
"""

from __future__ import annotations

import unittest

from Bio import Phylo
import io

from workflow.stability.clade_identity import (ITEM_ID_BITS, canonical_item_id,
                                               strip_accession_version)
from workflow.utils.treeUtils import (calculate_tree_hash, decode_tree_hash,
                                      encode_clade_to_int, encode_list_to_int)


def arvore(texto):
    return Phylo.read(io.StringIO(texto), "newick")


def clados_por_identidade(texto):
    """Identidade canônica de todo clado com mais de um terminal — a mesma regra
    que `SubtreeBuilder.build` aplica ao gravar o metadata.json."""
    identidades = set()
    for clado in arvore(texto).find_clades():
        terminais = [t.name for t in clado.get_terminals() if t.name]
        if len(terminais) > 1:
            identidades.add(encode_clade_to_int(terminais))
    return identidades


class TestIdentidadeCanonica(unittest.TestCase):

    def test_invariante_a_ordem_dos_terminais(self):
        """O coração de D5. Mesmo conjunto de terminais, mesma identidade."""
        self.assertEqual(encode_clade_to_int(['A.1', 'B.1', 'C.1']),
                         encode_clade_to_int(['C.1', 'A.1', 'B.1']))

    def test_invariante_a_ordem_de_travessia_em_arvores_reais(self):
        """Duas árvores com a mesma topologia e filhos em ordem trocada: o clado
        (A,B) tem de ser o MESMO item nas duas. Sob a identidade legada não era,
        e é por isso que FastTree e IQ-TREE não compartilhavam clados enquanto NJ
        e UPGMA — que vêm do mesmo construtor do Biopython — compartilhavam."""
        uma = clados_por_identidade("((A,B),(C,D));")
        outra = clados_por_identidade("((D,C),(B,A));")
        self.assertEqual(uma, outra)

    def test_clados_distintos_recebem_identidades_distintas(self):
        self.assertNotEqual(encode_clade_to_int(['A', 'B']),
                            encode_clade_to_int(['A', 'C']))

    def test_rotulo_truncado_e_integro_sao_o_mesmo_clado(self):
        """D5 encontrando D13: IQ-TREE e RAxML gravam `NC_008030.` onde FastTree
        grava `NC_008030.1`. Sem normalizar, o mesmo clado vira dois itens e o
        suporte cai pela metade sem que nada na topologia tenha mudado."""
        self.assertEqual(encode_clade_to_int(['NC_008030.', 'NC_001611.']),
                         encode_clade_to_int(['NC_008030.1', 'NC_001611.1']))

    def test_identidade_e_exata_em_javascript(self):
        """O consumidor mais estreito da cadeia manda no tamanho: o valor viaja
        no JSON da API até o navegador, e `Number` do JavaScript só é exato até
        2^53 - 1. Passar disso trocaria a colisão de 16 bits por arredondamento
        silencioso no cliente. Cabe também no inteiro de 64 bits do Neo4j."""
        self.assertLessEqual(ITEM_ID_BITS, 53)
        limite_seguro_js = 2 ** 53 - 1
        for nomes in (['A'], ['A', 'B'], [f'TAXON_{i}' for i in range(200)]):
            identidade = encode_clade_to_int(nomes)
            self.assertGreaterEqual(identidade, 0)
            self.assertLessEqual(identidade, limite_seguro_js)

    def test_espaco_e_muito_maior_que_o_legado(self):
        """16 bits davam 65 536 valores para milhares de clados."""
        self.assertGreater(2 ** ITEM_ID_BITS, 2 ** 16 * 10 ** 10)

    def test_nome_repetido_nao_muda_a_identidade(self):
        """A identidade é o CONJUNTO de terminais. Dois rótulos que normalizam
        para o mesmo acesso não podem produzir um item diferente."""
        self.assertEqual(canonical_item_id(['A', 'B']),
                         canonical_item_id(['A', 'B', 'A']))

    def test_bits_fora_da_faixa_e_erro(self):
        with self.assertRaises(ValueError):
            canonical_item_id(['A'], bits=0)
        with self.assertRaises(ValueError):
            canonical_item_id(['A'], bits=129)


class TestIdentidadeLegada(unittest.TestCase):
    """Caracteriza o esquema antigo — ele continua existindo, só que agora
    apenas para auditoria e para reler artefatos anteriores a M1.2."""

    def test_legada_depende_da_ordem(self):
        self.assertNotEqual(encode_list_to_int([1, 2, 3]),
                            encode_list_to_int([3, 2, 1]))

    def test_legada_cabe_em_16_bits(self):
        self.assertLess(encode_list_to_int([1, 2, 3]), 2 ** 16)

    def test_canonica_resolve_o_caso_que_a_legada_erra(self):
        """Lado a lado: mesma entrada conceitual, veredito oposto."""
        self.assertNotEqual(encode_list_to_int(['A', 'B']), encode_list_to_int(['B', 'A']))
        self.assertEqual(encode_clade_to_int(['A', 'B']), encode_clade_to_int(['B', 'A']))


class TestHashDeTerminal(unittest.TestCase):

    def test_grava_o_rotulo_normalizado(self):
        resultado = calculate_tree_hash('NC_008030.1', True, gbk_file='/dev/null')
        self.assertEqual(resultado['newick'], 'NC_008030')

    def test_terminal_truncado_e_integro_tem_o_mesmo_hash(self):
        truncado = calculate_tree_hash('NC_008030.', True, gbk_file='/dev/null')
        integro = calculate_tree_hash('NC_008030.1', True, gbk_file='/dev/null')
        self.assertEqual(truncado['terminal_hash'], integro['terminal_hash'])

    def test_decode_valida_o_que_calculate_gravou(self):
        """`decode_tree_hash` **nunca** validava: `calculate_tree_hash` hashava o
        rótulo bruto (`NC_008030.1`) e gravava o truncado (`NC_008030`), então a
        conferência de integridade sempre devolvia None. Agora fecha."""
        resultado = calculate_tree_hash('NC_008030.1', True, gbk_file='/dev/null')
        self.assertEqual(decode_tree_hash(resultado), 'NC_008030')

    def test_decode_recusa_hash_adulterado(self):
        resultado = calculate_tree_hash('NC_008030.1', True, gbk_file='/dev/null')
        resultado['terminal_hash'] += 1
        self.assertIsNone(decode_tree_hash(resultado))

    def test_terminais_distintos_tem_hashes_distintos(self):
        a = calculate_tree_hash('NC_008030.1', True, gbk_file='/dev/null')
        b = calculate_tree_hash('NC_001611.1', True, gbk_file='/dev/null')
        self.assertNotEqual(a['terminal_hash'], b['terminal_hash'])


class TestNormalizacao(unittest.TestCase):

    def test_remove_sufixo_de_versao(self):
        self.assertEqual(strip_accession_version('NC_008030.1'), 'NC_008030')
        self.assertEqual(strip_accession_version('NC_008030.'), 'NC_008030')

    def test_preserva_nome_sem_versao(self):
        self.assertEqual(strip_accession_version('L22579'), 'L22579')

    def test_reexportado_por_stability(self):
        """`audit_variola.py` e `TreeSet` importam de `workflow.stability.stability`.
        Mover a definição não pode quebrar esse caminho."""
        from workflow.stability.stability import strip_accession_version as pela_stability
        self.assertIs(pela_stability, strip_accession_version)


if __name__ == '__main__':
    unittest.main()
