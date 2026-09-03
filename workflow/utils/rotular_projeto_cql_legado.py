"""Acrescenta a propriedade `project` aos nós `Tree` de `.cql` já gerados.

Antes desta sessão (2026-09-03), `generate_cypher` (`neo4jProcessing.py`) não
gravava de qual projeto/experimento cada árvore vinha — um grafo com mais de
um projeto carregado não dava para separar por origem. O gerador já está
corrigido; este script é o reparo retroativo dos `.cql` que **já existem em
disco**, gerados antes do fix, sem reexecutar o pipeline (que os produziu).

Deriva o nome do projeto do próprio caminho do arquivo — mesma convenção que
`neo4jProcessing.project_name_from_output_path` usa na geração:
`.../projects/<project_name>/out/outputs/*.cql`.

Uso:
    python workflow/utils/rotular_projeto_cql_legado.py [--projects-dir DIR] [--apply]

Sem `--apply`, roda em modo de conferência (mostra quantos `CREATE (t:Tree
...)` seriam corrigidos por arquivo, sem gravar). Com `--apply`, sobrescreve
(faz backup em `<arquivo>.bak-preRotuloProjeto` antes, se ainda não existir).

**Não atualiza um grafo Neo4j já em execução.** Um `.cql` que já foi
executado contra um banco fez `CREATE` sem a propriedade nova — os nós
existentes não ganham `project` por mágica. Para refletir no grafo, rode o
`.cql` corrigido de novo (é Cypher puro, segundos — não é reexecutar o
pipeline de bioinformática) ou aplique um `MATCH (t:Tree {name: ..., uid:
...}) SET t.project = ...` equivalente.
"""
import argparse
import re
import shutil
import sys
from pathlib import Path

from workflow.utils.neo4jProcessing import project_name_from_output_path

TREE_RE = re.compile(r"CREATE \(t:Tree \{([^}]*)\}\)")


def reparar_linha(linha: str, project_name: str, stats: dict) -> str:
    m = TREE_RE.search(linha)
    if not m:
        return linha

    stats["nos_tree"] += 1
    corpo = m.group(1)
    if "project:" in corpo:
        stats["ja_rotulados"] += 1
        return linha

    novo_corpo = f"{corpo}, project: '{project_name}'"
    stats["corrigidos"] += 1
    return linha[: m.start()] + f"CREATE (t:Tree {{{novo_corpo}}})" + linha[m.end() :]


def reparar_arquivo(caminho: Path, aplicar: bool) -> dict:
    project_name = project_name_from_output_path(str(caminho.parent))
    stats = {"nos_tree": 0, "corrigidos": 0, "ja_rotulados": 0, "project": project_name}

    linhas_saida = []
    with caminho.open("r", encoding="utf-8") as f:
        for linha in f:
            linhas_saida.append(reparar_linha(linha, project_name, stats))

    if aplicar and stats["corrigidos"] > 0:
        backup = caminho.with_suffix(caminho.suffix + ".bak-preRotuloProjeto")
        if not backup.exists():
            shutil.copy2(caminho, backup)
        caminho.write_text("".join(linhas_saida), encoding="utf-8")

    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--projects-dir", type=Path, default=Path("projects"))
    parser.add_argument("--apply", action="store_true", dest="aplicar")
    args = parser.parse_args()

    arquivos = sorted(args.projects_dir.glob("*/out/outputs/*.cql"))
    if not arquivos:
        print(f"Nenhum .cql encontrado em {args.projects_dir}/*/out/outputs/")
        sys.exit(0)

    modo = "APLICADO" if args.aplicar else "conferência (sem gravar)"
    total_corrigidos = 0
    for arquivo in arquivos:
        resultado = reparar_arquivo(arquivo, args.aplicar)
        print(f"{arquivo} [{modo}] — projeto: {resultado['project']}")
        print(f"  nós Tree: {resultado['nos_tree']}  "
              f"a rotular: {resultado['corrigidos']}  "
              f"já rotulados: {resultado['ja_rotulados']}")
        total_corrigidos += resultado["corrigidos"]

    print(f"\n{len(arquivos)} arquivo(s), {total_corrigidos} nó(s) Tree a rotular no total.")
