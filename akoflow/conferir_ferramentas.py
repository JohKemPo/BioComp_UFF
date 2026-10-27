"""
Confere as ferramentas externas contra a receita pinada (`environment.yml` da raiz).

Equivale, para a imagem do AkôFlow, ao `check_dependencies.sh --strict`. Sai com 1
se alguma versão divergir, e imprime as versões encontradas — é assim que se
registra o ambiente de uma máquina fora da receita.
"""
import sys

from workflow.utils.manifest import tool_versions

PINADAS = {"mafft": "v7.526", "iqtree": "3.1.3", "raxml-ng": "2.0.2", "fasttree": "2.2.0",
           "mrbayes": "3.2.7", "clustalo": "1.2.4"}

if __name__ == "__main__":
    encontradas = tool_versions()
    divergentes = {t: (v, encontradas.get(t)) for t, v in PINADAS.items() if encontradas.get(t) != v}
    for t, v in sorted(encontradas.items()):
        print(f"{t:10} {v}   (receita: {PINADAS.get(t, '-')})")
    if divergentes:
        print(f"FORA DA RECEITA: {divergentes}", file=sys.stderr)
        sys.exit(1)
