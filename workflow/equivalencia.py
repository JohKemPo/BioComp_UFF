"""
Comparador de equivalência entre duas execuções do mesmo projeto (G4b, 02 §5.4).

Compara o ``out/`` de uma execução de referência (monolito) com o de outra
(granular, ou o próprio monolito de novo, que é o controle nulo) pelas classes
fixadas antes de rodar:

* **H1..H5** — SHA-256 de alinhamentos, árvores, subárvores e do ``.gb``;
* **S1..S5** — valores canonicalizados: seleção de modelo, FPMax, metadados,
  ``.cql`` e manifesto (a ordem de agregação depende do sistema de arquivos e
  não é um resultado);
* **P** — o conjunto de nomes dos PNG;
* **X** — logs, carimbos de tempo e ambiente: excluídos, de propósito.

Uso::

    python -m workflow.equivalencia REF/out TESTE/out [--json saida.json]

Sai com 0 só se todas as classes H e S forem iguais.
"""

import argparse
import collections
import csv
import glob
import json
import os
import re
import sys
from typing import Dict, List, Optional

from workflow.utils.manifest import file_digest

_BASICA = re.compile(r".*_(nj|upgma)_(distance|parsimony)\.nexus$")


def _hashes(out: str, *padroes: str, filtro=None) -> Dict[str, Optional[str]]:
    achados = {}
    for padrao in padroes:
        for caminho in glob.glob(os.path.join(out, padrao), recursive=True):
            relativo = os.path.relpath(caminho, out)
            if os.path.isfile(caminho) and (filtro is None or filtro(relativo)):
                achados[relativo] = file_digest(caminho)
    return achados


def _comparar_mapas(ref: dict, teste: dict) -> dict:
    em_comum = sorted(set(ref) & set(teste))
    diferentes = [k for k in em_comum if ref[k] != teste[k]]
    return {"iguais": len(em_comum) - len(diferentes), "diferentes": len(diferentes),
            "so_ref": sorted(set(ref) - set(teste)), "so_teste": sorted(set(teste) - set(ref)),
            "detalhe": diferentes[:20]}


def _comparar_valores(ref, teste, rotulo: str = "") -> dict:
    # Artefato ausente nos dois lados não é igualdade: é uma comparação que não aconteceu.
    ausente = ref is None and teste is None
    igual = ref == teste and not ausente
    detalhe = [] if igual else [f"{rotulo}: " + ("ausente nas duas execuções" if ausente else "valores diferentes")]
    return {"iguais": int(igual), "diferentes": int(not igual), "so_ref": [], "so_teste": [],
            "detalhe": detalhe}


def _ler_json(caminho: str):
    if not os.path.isfile(caminho):
        return None
    with open(caminho, encoding="utf-8") as fh:
        return json.load(fh)


def _itemset(texto: str) -> tuple:
    return tuple(sorted(int(n) for n in re.findall(r"-?\d+", texto)))


def _fpmax(caminho: str):
    """Multiconjunto de linhas, sem a coluna de índice e com o itemset como tupla ordenada."""
    if not os.path.isfile(caminho):
        return None
    with open(caminho, newline="") as fh:
        leitor = csv.reader(fh)
        cab = next(leitor)
        pos = cab.index("itemsets")
        return collections.Counter(
            (_itemset(l[pos]),) + tuple(v for i, v in enumerate(l) if i not in (0, pos)) for l in leitor)


def _metadata_json(caminho: str):
    dados = _ler_json(caminho)
    if dados is None:
        return None
    arvores = {}
    def achatar(x):
        for item in x:
            if isinstance(item, list):
                achatar(item)
            else:
                arvores.update(item)
    achatar(dados)
    return arvores


def _metadata_csv(caminho: str):
    if not os.path.isfile(caminho):
        return None
    csv.field_size_limit(sys.maxsize)
    with open(caminho, newline="") as fh:
        leitor = csv.reader(fh)
        next(leitor)
        return sorted(tuple(sorted(c for c in linha[1:])) for linha in leitor)


def _cql(caminho: str):
    if not os.path.isfile(caminho):
        return None
    with open(caminho, encoding="utf-8") as fh:
        texto = re.sub(r"""project: (['"])[^'"]*\1""", "project: <PROJETO>", fh.read())
    return collections.Counter(" ".join(l.split()) for l in texto.splitlines() if l.strip())


def _selecao_modelo(manifesto: Optional[dict]):
    if not manifesto:
        return None
    campos = ("estado", "modelo_escolhido", "candidatos", "traducao", "criterio",
              "conjunto_candidatos", "seed", "politica")
    return {os.path.basename(s["alinhamento"]): {c: s.get(c) for c in campos}
            for s in manifesto.get("model_selection") or []}


def _normalizar_comando(cmd: List[str]) -> List[str]:
    # Caminho vira nome: a entrada do MAFFT é `data/...` num e `out/entrada/...` no outro (esperado).
    return [os.path.basename(t) if os.sep in t else t for t in cmd]


def _manifesto(manifesto: Optional[dict]):
    if not manifesto:
        return None
    params = json.loads(json.dumps(manifesto.get("params", {})))
    for sec in (params, params.get("tree_config", {}), params.get("subtree_config", {}),
                params.get("subtree_config", {}).get("subtree_miner_configs", {})):
        for k in ("input_path", "output_path", "output_log", "project_name"):
            sec.pop(k, None)
    return {
        "inference_methods": sorted(
            (m["metodo"], m.get("alinhador"), m["estado"], os.path.basename(m["saida"]))
            for m in manifesto.get("inference_methods") or [] if m["metodo"] != "relogio_molecular"),
        "tools_invoked": {
            t: [_normalizar_comando(r["command"]) for r in d.get("runs", [])]
            for t, d in sorted((manifesto.get("tools_invoked") or {}).items())},
        "params": params,
        "reproducibility": manifesto.get("reproducibility"),
        "execution_mode": manifesto.get("execution_mode"),
        "outcome": (manifesto.get("outcome") or {}).get("status"),
        # os CSV entram canonicalizados em S2/S3, não por hash cru
        "outputs_sha256": {k: v for k, v in (manifesto.get("outputs_sha256") or {}).items()
                           if not k.endswith(".csv")},
    }


def comparar(ref_out: str, teste_out: str) -> Dict[str, dict]:
    """
    Compara dois ``out/`` pelas classes de 02 §5.4.

    Return
    ------
    dict
        ``classe -> {iguais, diferentes, so_ref, so_teste, detalhe}``.
    """
    r, t = os.path.abspath(ref_out), os.path.abspath(teste_out)
    mr = _ler_json(os.path.join(r, "outputs", "manifest.json"))
    mt = _ler_json(os.path.join(t, "outputs", "manifest.json"))
    resultado = {}
    for classe, padroes, filtro in (
            ("H1", ("Align/*.aln", "tmp/*.aln"), None),
            ("H2", ("Trees/*.nexus",), lambda p: bool(_BASICA.match(p))),
            ("H3", ("Trees/*.nexus", "tmp/**/*.raxml.support", "tmp/**/*.raxml.bestTree",
                    "tmp/**/*.treefile", "tmp/**/*.contree"),
             lambda p: not _BASICA.match(p) if p.startswith("Trees/") else True),
            ("H4", ("Subtrees/*.nexus",), None),
            ("H5", ("outputs/raw_data_sequences.gb",), None)):
        resultado[classe] = _comparar_mapas(_hashes(r, *padroes, filtro=filtro),
                                            _hashes(t, *padroes, filtro=filtro))
    resultado["S1"] = _comparar_valores(_selecao_modelo(mr), _selecao_modelo(mt), "model_selection")
    resultado["S2"] = _comparar_valores(_fpmax(os.path.join(r, "outputs", "all_results_fpmax.csv")),
                                        _fpmax(os.path.join(t, "outputs", "all_results_fpmax.csv")),
                                        "all_results_fpmax.csv")
    jr, jt = (_metadata_json(os.path.join(x, "outputs", "metadata.json")) for x in (r, t))
    cr, ct = (_metadata_csv(os.path.join(x, "outputs", "metadata.csv")) for x in (r, t))
    resultado["S3"] = _comparar_valores((jr, cr), (jt, ct), "metadata.json/csv")
    cqr = {os.path.basename(p): _cql(p) for p in glob.glob(os.path.join(r, "outputs", "*.cql"))}
    cqt = {os.path.basename(p): _cql(p) for p in glob.glob(os.path.join(t, "outputs", "*.cql"))}
    resultado["S4"] = _comparar_valores(cqr, cqt, "*.cql")
    resultado["S5"] = _comparar_valores(_manifesto(mr), _manifesto(mt), "manifest.json")
    if mr and mt and resultado["S5"]["diferentes"]:
        a, b = _manifesto(mr), _manifesto(mt)
        resultado["S5"]["detalhe"] = [k for k in a if a[k] != b[k]]
    pngs = lambda x: {os.path.relpath(p, x): os.path.getsize(p) > 0
                      for p in glob.glob(os.path.join(x, "outputs", "Plots", "*.png"))}
    pr, pt = pngs(r), pngs(t)
    resultado["P"] = {"iguais": len(set(pr) & set(pt)), "diferentes": 0,
                      "so_ref": sorted(set(pr) - set(pt)), "so_teste": sorted(set(pt) - set(pr)),
                      "detalhe": []}
    resultado["P"]["diferentes"] = len(resultado["P"]["so_ref"]) + len(resultado["P"]["so_teste"]) + \
        sum(1 for k in set(pr) & set(pt) if not (pr[k] and pt[k]))
    return resultado


def aprovado(resultado: Dict[str, dict]) -> bool:
    return all(v["diferentes"] == 0 and not v["so_ref"] and not v["so_teste"]
               for k, v in resultado.items() if k[0] in "HS")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref_out")
    ap.add_argument("teste_out")
    ap.add_argument("--json", help="grava o resultado completo neste arquivo")
    a = ap.parse_args(argv)
    resultado = comparar(a.ref_out, a.teste_out)
    print(f"{'classe':6} {'iguais':>7} {'difer.':>7} {'só ref':>7} {'só teste':>8}")
    for classe in sorted(resultado):
        v = resultado[classe]
        print(f"{classe:6} {v['iguais']:>7} {v['diferentes']:>7} {len(v['so_ref']):>7} {len(v['so_teste']):>8}"
              + (f"   {v['detalhe'][:3]}" if v["detalhe"] else ""))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(resultado, fh, indent=2, ensure_ascii=False)
    return 0 if aprovado(resultado) else 1


if __name__ == "__main__":
    sys.exit(main())
