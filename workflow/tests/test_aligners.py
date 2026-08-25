"""
Biblioteca de alinhadores e política de substituição — D1.

D1 é o defeito mais caro do projeto, e ele **não** é "o Clustal Omega não
aguenta 250 kb". É que o controlador, ao descobrir isso, trocava para MAFFT
**mantendo o nome de arquivo `*_clustalo_*`** — metade dos "pipelines" dos
experimentos de *Variola* são cópias byte a byte, e o fator alinhador, que é
metade do delineamento, não existe ali.

O que estes testes fixam:

* o limite é **declarado com o motivo**, para poder ser mostrado antes da
  execução, quando ainda há escolha;
* substituição sem autorização é **erro**, não desvio;
* quando autorizada, a substituição devolve **o nome do que rodou**.

Executar com:

    python -m unittest workflow.tests.test_aligners
"""

from __future__ import annotations

import unittest

from workflow.alignment.aligners import (ALIGNERS, MAQUINA_DEV_BYTES, AlignerPolicy,
                                         available_aligners, resolve_aligner, viability)


class TestBiblioteca(unittest.TestCase):

    def test_os_tres_alinhadores_estao_declarados(self):
        self.assertEqual(sorted(ALIGNERS), ["clustalo", "mafft", "muscle"])

    def test_cada_limite_traz_o_motivo(self):
        """Limite sem explicação vira superstição — e este projeto já tem um
        limite de 20 kb herdado sem ninguém saber de onde veio."""
        for chave, a in ALIGNERS.items():
            if a.max_sequence_bp or a.max_sequences:
                self.assertTrue(a.note.strip(), f"{chave} tem limite sem motivo declarado")

    def test_versao_e_none_quando_nao_instalado(self):
        for chave, versao in available_aligners().items():
            self.assertTrue(versao is None or isinstance(versao, str), chave)


class TestViabilidade(unittest.TestCase):

    def test_conjunto_pequeno_aceita_os_tres(self):
        """10 788 pb é o maior comprimento **verificado funcionando** nos três
        alinhadores: é o conjunto de validação, medido em 2026-08-25."""
        v = viability(n_sequences=20, max_sequence_bp=10_788)
        instalados = {k for k, x in v.items() if x.installed}
        for chave in instalados:
            self.assertTrue(v[chave].viable, f"{chave}: {v[chave].reasons}")

    def test_variola_reprova_os_dois_alinhadores_medidos(self):
        """Medido nesta máquina: em 52 sequências de até 228 kb, Clustal Omega e
        MUSCLE são mortos pelo OOM killer. Sobra o MAFFT — e é por isso que o
        fator alinhador não existe nos conjuntos de Variola."""
        v = viability(n_sequences=52, max_sequence_bp=228_250,
                      available_bytes=MAQUINA_DEV_BYTES)
        self.assertTrue(v["mafft"].viable)
        self.assertFalse(v["clustalo"].viable)
        self.assertFalse(v["muscle"].viable)


class TestLimiteEhRelativoAMaquina(unittest.TestCase):
    """R2 — um limite de ferramenta não é propriedade da ferramenta.

    A lei de escala é do algoritmo e viaja entre arquiteturas; o ponto onde ela
    cruza o orçamento é da máquina. Declarar só o segundo, como um limite em
    pares de base, produz um número que parece universal e não é.
    """

    CONJUNTO_VARIOLA = dict(n_sequences=52, max_sequence_bp=228_250)

    def test_o_mesmo_conjunto_muda_de_veredito_com_a_maquina(self):
        """O ponto inteiro do R2: mesma ferramenta, mesmo dado, máquinas
        diferentes, vereditos diferentes — e isso é o comportamento correto."""
        aqui = viability(**self.CONJUNTO_VARIOLA, available_bytes=MAQUINA_DEV_BYTES)
        grande = viability(**self.CONJUNTO_VARIOLA, available_bytes=128 * 10**9)
        self.assertFalse(aqui["muscle"].viable)
        self.assertTrue(grande["muscle"].viable,
                        "numa máquina de 128 GB o requisito cabe; vetar ali seria "
                        "transportar um limite que é desta máquina")

    def test_falha_observada_vence_estimativa(self):
        """A estimativa vem de um pico registrado no instante da morte — é piso,
        não requisito, e subestima por construção. Sem esta regra, o MUSCLE
        apareceria como viável justamente no conjunto em que o vimos morrer."""
        v = viability(**self.CONJUNTO_VARIOLA, available_bytes=MAQUINA_DEV_BYTES)["muscle"]
        self.assertFalse(v.viable)
        self.assertIn("observado falhar", " ".join(v.reasons))

    def test_falha_observada_nao_condena_maquina_maior(self):
        v = viability(**self.CONJUNTO_VARIOLA, available_bytes=512 * 10**9)["muscle"]
        self.assertTrue(v.viable)

    def test_maquina_menor_tambem_e_condenada(self):
        v = viability(**self.CONJUNTO_VARIOLA, available_bytes=8 * 10**9)["muscle"]
        self.assertFalse(v.viable)

    def test_tolerancia_impede_erro_de_arredondamento(self):
        """A medição foi anotada com o valor exato reportado pelo sistema. Com
        `>=` estrito, a própria máquina deixaria de reconhecer a própria falha
        por alguns MB de diferença entre leituras."""
        modelo = ALIGNERS["muscle"].resources
        quase = int(MAQUINA_DEV_BYTES * 1.02)
        self.assertIsNotNone(
            modelo.blocking_failure(52, 228_250, quase),
            "diferença de 2% na leitura de memória não pode apagar a observação")

    def test_o_veredito_carrega_requisito_e_orcamento(self):
        """A mensagem útil não é 'indisponível', é 'precisa de X e há Y' — a
        primeira é um veto, a segunda é um requisito."""
        v = viability(**self.CONJUNTO_VARIOLA, available_bytes=128 * 10**9)["muscle"]
        self.assertIsNotNone(v.estimated_bytes)
        self.assertEqual(v.available_bytes, 128 * 10**9)

    def test_estimativa_nao_ajustada_e_declarada_como_tal(self):
        """Com pontos de uma máquina só, expoente e deslocamento ficam
        confundidos: qualquer curva passa por dois pontos. `fitted=False` é o
        que impede a estimativa de ser lida como previsão."""
        for chave, a in ALIGNERS.items():
            if a.resources is not None:
                self.assertFalse(a.resources.fitted,
                                 f"{chave} afirma curva ajustada; ajustar exige "
                                 f"bissectar em pelo menos duas máquinas")

    def test_medicao_sem_maquina_nao_transfere(self):
        """Uma medição sem a memória da máquina não é interpretável e, portanto,
        não pode vetar em lugar nenhum."""
        from workflow.alignment.aligners import Measurement, ResourceModel
        modelo = ResourceModel(scaling="n*L", measurements=(
            Measurement(n_sequences=10, max_sequence_bp=1_000, outcome="oom",
                        machine="desconhecida", date="2026-01-01"),
        ))
        self.assertIsNone(modelo.blocking_failure(50, 5_000, 8 * 10**9))

    def test_sequencia_longa_reprova_clustalo(self):
        """O caso de Variola: 240 kb por genoma, acima do limite de 20 kb."""
        v = viability(n_sequences=49, max_sequence_bp=240_000)
        self.assertFalse(v["clustalo"].viable)
        self.assertTrue(any("20.000" in r or "20,000" in r for r in v["clustalo"].reasons))
        self.assertTrue(v["mafft"].viable)

    def test_muitas_sequencias_reprovam_muscle(self):
        v = viability(n_sequences=5_000, max_sequence_bp=10_000)
        self.assertFalse(v["muscle"].viable)
        self.assertTrue(v["mafft"].viable)

    def test_o_veredito_diz_por_que(self):
        v = viability(n_sequences=49, max_sequence_bp=240_000)
        self.assertTrue(v["clustalo"].reasons)
        self.assertIn("OOM", " ".join(v["clustalo"].reasons))

    def test_resumo_e_serializavel(self):
        import json
        json.dumps({k: x.summary() for k, x in viability(20, 10_000).items()})


class TestPoliticaDeSubstituicao(unittest.TestCase):

    def test_viavel_passa_sem_substituicao(self):
        efetivo, motivo = resolve_aligner("mafft", 20, 10_000)
        self.assertEqual(efetivo, "mafft")
        self.assertIsNone(motivo)

    def test_inviavel_falha_por_padrao(self):
        """O padrão é falhar. Trocar em silêncio é exatamente D1."""
        with self.assertRaises(ValueError) as erro:
            resolve_aligner("clustalo", 49, 240_000)
        mensagem = str(erro.exception)
        self.assertIn("não é viável", mensagem)
        self.assertIn("fallback", mensagem, "a mensagem precisa dizer como autorizar")

    def test_substituicao_autorizada_devolve_o_nome_do_que_rodou(self):
        """O núcleo de D1: quem chama tem de nomear a saída pelo alinhador
        efetivo, e por isso ele é devolvido."""
        efetivo, motivo = resolve_aligner(
            "clustalo", 49, 240_000, AlignerPolicy(on_unavailable="fallback"))
        self.assertNotEqual(efetivo, "clustalo")
        self.assertIsNotNone(motivo)
        self.assertIn("substituído", motivo)

    def test_alinhador_desconhecido_e_erro(self):
        with self.assertRaises(ValueError) as erro:
            resolve_aligner("clustalw", 20, 10_000)
        self.assertIn("Alinhador desconhecido", str(erro.exception))

    def test_sem_alternativa_viavel_tambem_falha(self):
        """Substituir por nada não é substituir."""
        politica = AlignerPolicy(on_unavailable="fallback", fallback_order=("muscle",))
        with self.assertRaises(ValueError) as erro:
            resolve_aligner("clustalo", 5_000, 240_000, politica)
        self.assertIn("Nenhum alinhador viável", str(erro.exception))


if __name__ == "__main__":
    unittest.main()
