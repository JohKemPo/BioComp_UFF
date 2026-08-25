"""
Testes da resolução de binários das ferramentas externas.

O defeito que este módulo evita: o pipeline chamava `iqtree2` fixo, e o pacote
`iqtree` 3.x do bioconda instala `iqtree` e `iqtree3` — sem `iqtree2`. Quem
seguisse a receita do projeto instalava o IQ-TREE com sucesso e mesmo assim
não conseguia rodar.
"""

import os
import stat
import tempfile
import unittest

from workflow.utils.external_tools import (
    CANDIDATOS,
    require_tool,
    resolve_tool,
    resolved_tools,
)


def _binario_falso(diretorio, nome):
    caminho = os.path.join(diretorio, nome)
    with open(caminho, "w") as handle:
        handle.write("#!/bin/sh\nexit 0\n")
    os.chmod(caminho, os.stat(caminho).st_mode | stat.S_IEXEC)
    return caminho


class TestResolucao(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp()

    def test_encontra_no_diretorio_preferido(self):
        esperado = _binario_falso(self.dir, "iqtree3")
        self.assertEqual(resolve_tool("iqtree", prefer_dir=self.dir), esperado)

    def test_diretorio_preferido_ganha_do_path(self):
        # É o caso desta base: o env tinha FastTree 2.2.0 e o PATH resolvia
        # 2.1.11 de /usr/bin. Medir o segundo e chamá-lo de "versão do projeto"
        # produziu um registro errado que durou dias.
        esperado = _binario_falso(self.dir, "FastTree")
        self.assertEqual(resolve_tool("fasttree", prefer_dir=self.dir), esperado)

    def test_respeita_a_ordem_de_preferencia(self):
        _binario_falso(self.dir, "iqtree2")
        esperado = _binario_falso(self.dir, "iqtree3")
        self.assertEqual(resolve_tool("iqtree", prefer_dir=self.dir), esperado)

    def test_encontra_nome_alternativo(self):
        esperado = _binario_falso(self.dir, "mb")
        self.assertEqual(resolve_tool("mrbayes", prefer_dir=self.dir), esperado)

    def test_arquivo_sem_permissao_de_execucao_nao_conta(self):
        caminho = os.path.join(self.dir, "iqtree3")
        with open(caminho, "w") as handle:
            handle.write("não é executável")
        # Sem candidato utilizável no diretório, cai no PATH; o que importa é
        # que o arquivo inerte não seja devolvido como se fosse o binário.
        self.assertNotEqual(resolve_tool("iqtree", prefer_dir=self.dir), caminho)

    def test_ferramenta_inexistente_devolve_none(self):
        self.assertIsNone(resolve_tool("ferramenta-que-nao-existe"))

    def test_diretorio_preferido_inexistente_nao_quebra(self):
        resolve_tool("iqtree", prefer_dir="/caminho/que/nao/existe")


class TestRequireTool(unittest.TestCase):

    def test_levanta_com_mensagem_acionavel(self):
        with self.assertRaises(FileNotFoundError) as ctx:
            require_tool("ferramenta-que-nao-existe")
        mensagem = str(ctx.exception)
        # Um "command not found" cru não diz a quem o lê que o pacote pode
        # estar instalado sob outro nome de binário.
        self.assertIn("Procurei por", mensagem)
        self.assertIn("setup_env.sh", mensagem)

    def test_nomeia_todos_os_candidatos(self):
        # Nomes inventados: usar uma ferramenta real faria o teste passar ou
        # falhar conforme o que estivesse instalado na máquina que o roda.
        candidatos = ("ferramenta-a", "ferramenta-b", "ferramenta-c")
        CANDIDATOS["ferramenta-de-teste"] = candidatos
        try:
            with self.assertRaises(FileNotFoundError) as ctx:
                require_tool("ferramenta-de-teste")
        finally:
            del CANDIDATOS["ferramenta-de-teste"]
        for nome in candidatos:
            self.assertIn(nome, str(ctx.exception))


class TestResolvedTools(unittest.TestCase):

    def test_cobre_todas_as_ferramentas(self):
        self.assertEqual(set(resolved_tools()), set(CANDIDATOS))

    def test_valores_sao_caminho_ou_none(self):
        for nome, caminho in resolved_tools().items():
            self.assertTrue(caminho is None or os.path.isabs(caminho), nome)


if __name__ == "__main__":
    unittest.main()
