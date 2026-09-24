"""
Controlador da datação leve (LSD2) — E12.10 e E12.11.

Roda **depois** de `SubtreeBuilderController`, em sequência (DEC-127: o
pipeline é sequencial e fica assim), sobre as árvores que já estão em
`out/Trees`. Nenhuma busca de topologia nova: a topologia é fixada com ``-te``.

Entrada que o controlador lê do `tree_config` (via `getattr`, sem schema, como
`model_selection`/`mrbayes_*`):

``run_molecular_clock`` : bool, padrão ``False``
    Liga a etapa. Desligada, o manifesto registra ``ignorado_por_configuracao``.
``molecular_clock_outgroup`` : list of str | str
    Acessos do grupo externo — lista, texto separado por vírgula, ou caminho de
    um arquivo com um acesso por linha. **Obrigatório**: sem ele a datação é
    ``tentado_e_falhou``. O protocolo validado em E12 (§9.1) exclui o grupo
    externo da árvore de tempo; a alternativa sem grupo externo (raiz
    procurada pelo LSD2) deu tMRCA 2521 e 5774 em ZIKV. Declarado, nunca
    inferido (`stability/rooting.py`, regra 1).
``molecular_clock_tree`` : str, opcional
    Nome do arquivo em `out/Trees` cuja datação vai para o topo do JSON.
``molecular_clock_model`` : str, padrão ``GTR+G``
``molecular_clock_ci`` : int, padrão ``100`` (réplicas do ``--date-ci``)
``molecular_clock_reoptimize_branches`` : bool, padrão ``True``
``molecular_clock_threads`` : int, padrão ``1``
    **Não** herda `iqtree_threads`: com mais de uma thread o ponto estimado
    muda entre corridas (medido; ver `workflow.molecular_clock.datacao`).
``random_seed`` : int — o mesmo do resto do pipeline, repassado em ``-seed``.
``align_method`` : str — só desempata a escolha da árvore principal.

Saída (E12.11), previsível mesmo quando a etapa roda por fora do workflow::

    out/outputs/molecular_clock/
        relogio_molecular.json      contrato consumido pelo Backend
        timetree.nex                cópia do timetree da árvore principal
        <arvore>/                   uma pasta por árvore datada (datas, Newick
                                    de entrada, saída integral do IQ-TREE/LSD2)

Execução avulsa (sem o workflow)::

    python -m workflow.controller.molecularClockController \\
        --out projects/<projeto>/out --outgroup-file og.txt
"""

from __future__ import annotations

import argparse
import datetime
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from typing import Dict, List, Optional

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../..'))

from workflow.molecular_clock import datacao
from workflow.stability.clade_identity import (canonical_bipartition, canonical_item_id,
                                               strip_accession_version)
from workflow.stability.stability import PipelineLabel
from workflow.utils import external_tools, manifest, tool_runs

__all__ = ["MolecularClockController", "METODOS_DATAVEIS", "VERSAO_DO_CONTRATO", "main"]

#: Métodos cuja árvore é datada. São os três validados em E12.5 (amplitude de
#: 1,9 ano no nó mais profundo). Fora: NJ/UPGMA (RF 0,92-0,94 contra a árvore
#: publicada, E12 §7) e parcimônia (comprimento de ramo não é subs/sítio sob
#: modelo); MrBayes não foi validado — o consenso do `sumt` não foi datado em E12.
METODOS_DATAVEIS = ("iqtree", "raxml", "fasttree")

#: Muda quando o formato de `relogio_molecular.json` mudar de modo incompatível.
VERSAO_DO_CONTRATO = 1

METODO_NO_MANIFESTO = "relogio_molecular"

_ADVERTENCIAS = [
    "Estimativa pontual do LSD2 (mínimos quadrados, relógio estrito), NÃO posterior bayesiana.",
    "O intervalo é de reamostragem do LSD2 (--date-ci), condicional à topologia, aos "
    "comprimentos de ramo, às datas e ao relógio estrito; cobriu a verdade em 96,9 % dos "
    "nós sob relógio estrito e em 60,6 % sob taxa variável por ramo (E12 §5).",
    "O resultado depende do enraizamento tanto quanto do dado: grupo externo declarado, "
    "excluído da árvore de tempo (E12 §9.1).",
    "Ramo colapsado pelo LSD2 (menos de meia substituição esperada) não tem data; "
    "a bipartição não aparece em datas_por_clado.",
    "Data de divergência não é evento de transmissão: nada aqui indica cadeia, "
    "direção nem contato.",
]


def _ler_grupo_externo(valor) -> List[str]:
    """Lista, texto com vírgulas ou caminho de arquivo — nesta ordem de teste."""
    if valor is None:
        return []
    if isinstance(valor, (list, tuple, set)):
        return [str(v).strip() for v in valor if str(v).strip()]
    texto = str(valor).strip()
    if not texto:
        return []
    if os.path.isfile(texto):
        with open(texto) as fh:
            return [x for x in fh.read().split() if x]
    return [x.strip() for x in re.split(r"[,\s]+", texto) if x.strip()]


def _versao_iqtree(binario: str) -> Optional[str]:
    try:
        saida = subprocess.run([binario, "--version"], capture_output=True,
                               text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r"IQ-TREE (?:multicore )?version ([\d.]+)", saida)
    return m.group(1) if m else None


def _versao_lsd2(relatorio: str) -> Optional[str]:
    """Versão no cabeçalho do relatório ``.lsd`` (``... DATES - v.2.4.4``)."""
    try:
        with open(relatorio) as fh:
            cabecalho = fh.read(2000)
    except OSError:
        return None
    m = re.search(r"RATES AND DATES\s*-\s*v\.?\s*([\d.]+\d)", cabecalho)
    return m.group(1) if m else None


class MolecularClockController:
    """
    Datação leve (LSD2 via ``iqtree --date``) das árvores de `out/Trees`.

    Toda saída declara seu desfecho: cada árvore datável termina ``executado``
    ou ``tentado_e_falhou`` (com motivo) em `tool_runs.registrar_metodo`, e o
    JSON é escrito **também quando falha** — o Backend precisa poder dizer por
    que não há data, não só que não há.
    """

    def __init__(self, **kwargs) -> None:
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.output_path = getattr(self, 'output_path', None)
        self.via = getattr(self, 'via', 'workflow')
        if not self.output_path:
            raise ValueError("MolecularClockController exige `output_path` (o `out/` do projeto).")
        self.dir_trees = os.path.join(self.output_path, 'Trees')
        self.dir_align = os.path.join(self.output_path, 'Align')
        self.gb = os.path.join(self.output_path, 'outputs', 'raw_data_sequences.gb')
        self.destino = os.path.join(self.output_path, 'outputs', 'molecular_clock')
        self.json_path = os.path.join(self.destino, 'relogio_molecular.json')

    # ------------------------------------------------------------------ #
    @property
    def ligado(self) -> bool:
        return bool(getattr(self, 'run_molecular_clock', False))

    def _parametros(self) -> Dict:
        from workflow.tree_construction.builder import reproducibility_settings

        cfg = {k: getattr(self, k) for k in ('random_seed',) if hasattr(self, k)}
        return {
            "modelo": str(getattr(self, 'molecular_clock_model', 'GTR+G')),
            "ci": int(getattr(self, 'molecular_clock_ci', 100)),
            "reotimizar_ramos": bool(getattr(self, 'molecular_clock_reoptimize_branches', True)),
            "threads": int(getattr(self, 'molecular_clock_threads', 1)),
            "semente": reproducibility_settings(cfg)["random_seed"],
            "ano_pontual": False,
        }

    def _arvores_candidatas(self) -> List[str]:
        if not os.path.isdir(self.dir_trees):
            return []
        return sorted(n for n in os.listdir(self.dir_trees) if n.endswith('.nexus'))

    # ------------------------------------------------------------------ #
    def __call__(self) -> Optional[Dict]:
        """
        Executa a etapa. Não levanta: falha vira ``tentado_e_falhou`` declarado.

        Return
        ------
        dict or None
            O conteúdo de `relogio_molecular.json`, ou ``None`` quando a etapa
            está desligada (nada é escrito).
        """
        if not self.ligado:
            tool_runs.registrar_metodo(
                METODO_NO_MANIFESTO, "ignorado_por_configuracao", saida=self.json_path,
                motivo="run_molecular_clock ausente ou falso no tree_config")
            logging.info("Relógio molecular: desligado por configuração.")
            return None

        logging.info("STEP: Molecular clock (LSD2).")
        os.makedirs(self.destino, exist_ok=True)
        params = self._parametros()
        grupo_externo = _ler_grupo_externo(getattr(self, 'molecular_clock_outgroup', None))
        resumo = self._esqueleto(params, grupo_externo)

        motivo_global = None
        binario = external_tools.resolve_tool('iqtree')
        if not grupo_externo:
            motivo_global = ("grupo externo não declarado (molecular_clock_outgroup): o "
                             "enraizamento da datação tem de ser declarado, nunca inferido — "
                             "sem ele o LSD2 procura a raiz e dá tMRCA absurdo (E12 §9.1)")
        elif binario is None:
            motivo_global = "IQ-TREE não encontrado (iqtree3/iqtree2/iqtree) — ative o env conda"
        elif not os.path.isfile(self.gb):
            motivo_global = f"sem datas de coleta: {os.path.relpath(self.gb, self.output_path)} não existe"

        candidatas = self._arvores_candidatas()
        datavel = [n for n in candidatas
                   if PipelineLabel.parse(n, prefix="tree_").inference in METODOS_DATAVEIS]
        resumo["arvores_nao_dataveis"] = {
            n: f"método '{PipelineLabel.parse(n, prefix='tree_').inference}' fora de "
               f"{list(METODOS_DATAVEIS)} (validados em E12.5)"
            for n in candidatas if n not in datavel}
        if motivo_global is None and not datavel:
            motivo_global = (f"nenhuma árvore de {list(METODOS_DATAVEIS)} em "
                             f"{os.path.relpath(self.dir_trees, self.output_path)}")

        if motivo_global is not None:
            tool_runs.registrar_metodo(METODO_NO_MANIFESTO, "tentado_e_falhou",
                                       saida=self.json_path, motivo=motivo_global)
            resumo.update(estado="tentado_e_falhou", motivo=motivo_global)
            return self._gravar(resumo)

        resumo["ferramentas"]["iqtree"] = _versao_iqtree(binario)
        crus = datacao.extrair_collection_date(self.gb)
        for nome in datavel:
            resumo["por_arvore"][nome] = self._datar(nome, binario, crus, grupo_externo, params)
            if resumo["ferramentas"].get("lsd2") is None:
                rel = resumo["por_arvore"][nome].get("_relatorio_lsd")
                if rel:
                    resumo["ferramentas"]["lsd2"] = _versao_lsd2(rel)
        for r in resumo["por_arvore"].values():
            r.pop("_relatorio_lsd", None)

        self._escolher_principal(resumo)
        return self._gravar(resumo)

    # ------------------------------------------------------------------ #
    def _esqueleto(self, params: Dict, grupo_externo: List[str]) -> Dict:
        return {
            "versao_do_contrato": VERSAO_DO_CONTRATO,
            "estado": None,
            "motivo": None,
            "arvore_principal": None,
            "criterio_da_arvore_principal": None,
            "taxa_subs_por_sitio_por_ano": None,
            "intervalo_da_taxa": None,
            "tmrca": None,
            "tmrca_lsd2": None,
            "intervalo_do_tmrca": None,
            "solucao_unica": None,
            "percentual_de_ramos_internos_colapsados": None,
            "datas_por_clado": {},
            "timetree": None,
            "unidades": {
                "taxa": "substituições por sítio por ano",
                "datas": "ano decimal (ex.: 1944.93 = dez/1944); *_lsd2 traz o texto cru do LSD2",
                "intervalo": "[inferior, superior] do --date-ci, mesma unidade da data",
            },
            "identidade_de_clado": (
                "chave = str(canonical_item_id) da bipartição canônica (D3), calculada "
                "contra TODAS as pontas da árvore de origem em out/Trees (grupo externo "
                "incluído) — o mesmo clade_id de suporte_de_ramo.py, metadata.json e FPMax."),
            "semantica_da_data": (
                "data do nó do lado das pontas da aresta, isto é, do MRCA do clado "
                "enraizado; quando esse clado não é o lado canônico (menor) da "
                "bipartição, clado_datado_e_o_lado_canonico = false."),
            "protocolo": {
                "ferramenta": "iqtree --date (LSD2 embutido)",
                "raiz": "grupo externo declarado, excluído da árvore de tempo "
                        "(-o ... --date-no-outgroup) — E12 §9.1",
                "grupo_externo_declarado": sorted(grupo_externo),
                "topologia": "fixa (-te), árvore de out/Trees",
                "modelo": params["modelo"],
                "comprimentos_de_ramo": ("reotimizados sob o modelo na topologia fixa"
                                         if params["reotimizar_ramos"] else
                                         "os do pipeline (-blfix)"),
                "intervalo": (f"--date-ci {params['ci']}" if params["ci"] else "nenhum"),
                "threads": params["threads"],
                "semente": params["semente"],
                "keep_ident": True,
                "convencao_ano_so": "intervalo YYYY:YYYY+1",
                "minimo_de_pontas_datadas": "max(3, 10% das pontas do grupo interno)",
                "portao_de_sinal_temporal": ("inclinação raiz-ponta > 0 no grupo interno, "
                                             "raiz na aresta do grupo externo"),
            },
            "advertencias": list(_ADVERTENCIAS),
            "arvores_nao_dataveis": {},
            "por_arvore": {},
            "ferramentas": {"iqtree": None, "lsd2": None},
            "execucao": {
                "via": self.via,
                "gerado_em": datetime.datetime.now(datetime.timezone.utc)
                .isoformat(timespec="seconds"),
            },
        }

    def _falha(self, rel: Dict, saida: str, alinhador: str, motivo: str) -> Dict:
        tool_runs.registrar_metodo(METODO_NO_MANIFESTO, "tentado_e_falhou", saida=saida,
                                   alinhador=alinhador, motivo=motivo)
        logging.warning(f"Relógio molecular — {rel['arvore']}: {motivo}")
        rel.update(estado="tentado_e_falhou", motivo=motivo)
        return rel

    def _datar(self, nome: str, binario: str, crus: Dict, grupo_externo: List[str],
               params: Dict) -> Dict:
        """Data uma árvore. Toda saída é um dict com `estado`; nada levanta."""
        rotulo = PipelineLabel.parse(nome, prefix="tree_")
        stem = os.path.splitext(nome)[0]
        pasta = os.path.join(self.destino, stem)
        os.makedirs(pasta, exist_ok=True)
        prefixo = os.path.join(pasta, stem)
        timetree = f"{prefixo}.timetree.nex"
        rel: Dict = {"arvore": nome, "pipeline": rotulo.name, "alinhador": rotulo.aligner,
                     "metodo_de_inferencia": rotulo.inference, "estado": None, "motivo": None,
                     "timetree": None}
        try:
            # tree_dataset_final_mafft_iqtree -> dataset_final_mafft.aln
            base = stem[len("tree_"):] if stem.startswith("tree_") else stem
            aln = os.path.join(self.dir_align,
                               base[: -len("_" + rotulo.inference)] + ".aln")
            if not os.path.isfile(aln):
                return self._falha(rel, timetree, rotulo.aligner,
                                   f"alinhamento não encontrado: {os.path.basename(aln)}")

            arv = datacao.ler_arvore_de_origem(os.path.join(self.dir_trees, nome))
            todos = frozenset(strip_accession_version(l.taxon.label)
                              for l in arv.leaf_node_iter())
            og = datacao.verificar_grupo_externo(arv, grupo_externo)
            rel["grupo_externo"] = og
            if og.get("motivo"):
                return self._falha(rel, timetree, rotulo.aligner, og["motivo"])
            grupo_ext = frozenset(og["encontrado"])
            interno = todos - grupo_ext
            biparticoes_origem = datacao.biparticoes_do_grupo_interno(arv, grupo_ext)

            nwk = f"{prefixo}.entrada.nwk"
            rotulos = datacao.escrever_entrada(arv, datacao.ler_rotulos_do_alinhamento(aln), nwk)
            datas = datacao.escrever_arquivo_datas(crus, rotulos, f"{prefixo}.datas.txt",
                                                   interno, params["ano_pontual"])
            decimais = datas.pop("_decimais")
            rel["datas"] = datas

            minimo = datacao.minimo_de_pontas_datadas(len(interno))
            if datas["pontas_datadas_no_grupo_interno"] < minimo:
                # Sem este portão o LSD2 NÃO falha: grava timetree com date="0"
                # em toda ponta e tMRCA -266334 (VARV-49, E12 §8.1).
                return self._falha(
                    rel, timetree, rotulo.aligner,
                    f"datas de coleta insuficientes: {datas['pontas_datadas_no_grupo_interno']} "
                    f"pontas datadas no grupo interno, mínimo {minimo} "
                    f"(max(3, 10% de {len(interno)}))")

            sinal = datacao.sinal_temporal(arv, decimais, grupo_ext)
            rel["sinal_temporal"] = sinal
            inclinacao = sinal.get("inclinacao_subs_sitio_ano")
            if inclinacao is None or inclinacao <= 0:
                return self._falha(
                    rel, timetree, rotulo.aligner,
                    "sem sinal temporal no grupo interno: regressão raiz-ponta com "
                    f"inclinação {inclinacao!r} (R² {sinal.get('r2')!r}, n={sinal.get('n')}) — "
                    "datar produziria número sem lastro (E12 §4)")

            og_rotulos = [rotulos[n] for n in sorted(grupo_ext)]

            def _registrar(cmd, saida):
                tool_runs.registrar('iqtree', cmd, saida=saida)
                tool_runs.anotar('iqtree', saida,
                                 etapa="datação LSD2 (--date, topologia fixa) — E12.10",
                                 seed=params["semente"], threads=params["threads"],
                                 model=params["modelo"])

            res = datacao.rodar_lsd2(
                binario, aln, nwk, f"{prefixo}.datas.txt", prefixo, og_rotulos,
                modelo=params["modelo"], ci=params["ci"], threads=params["threads"],
                semente=params["semente"], reotimizar_ramos=params["reotimizar_ramos"],
                registrar=_registrar)
            rel["_relatorio_lsd"] = f"{prefixo}.timetree.lsd"
            rel["lsd2"] = {k: res[k] for k in (
                "linha_de_resultado_lsd2", "retorno", "erro_lsd2",
                "aviso_sem_informacao_de_data_suficiente",
                "percentual_de_ramos_internos_colapsados_lsd2",
                "sequencias_com_mais_de_50pc_gap")}
            rel["lsd2"]["comando"] = " ".join(os.path.basename(c) if os.path.isabs(c) else c
                                              for c in res["comando"])
            for k in ("taxa_subs_por_sitio_por_ano", "intervalo_da_taxa", "tmrca",
                      "tmrca_lsd2", "intervalo_do_tmrca", "solucao_unica"):
                rel[k] = res[k]

            if res["retorno"] != 0 or res["timetree"] is None:
                return self._falha(
                    rel, timetree, rotulo.aligner,
                    f"IQ-TREE/LSD2 terminou com código {res['retorno']}: "
                    f"{res['erro_lsd2'] or 'sem timetree'}")
            if not res["solucao_unica"]:
                return self._falha(
                    rel, timetree, rotulo.aligner,
                    "LSD2 sem solução única ('not enough input date information' ou ':NA' "
                    f"na linha de resultado: {res['linha_de_resultado_lsd2']!r}) — o número "
                    "que ele imprime não é estimativa (E12 §8.2)")

            por_clado, contagem = datacao.datas_por_clado(res["timetree"], todos)
            rel["datas_por_clado"] = por_clado
            # A porcentagem vem da contagem de bipartições, não do texto do
            # LSD2: ele só a imprime às vezes, e ausência não é zero (regra 5).
            # Denominador = ramos internos do grupo interno, SEM a aresta
            # grupo interno | grupo externo: ela vira a raiz da árvore de tempo
            # e não pode ser colapsada. É a mesma conta do LSD2 (conferido:
            # 210/448 = 46,875 % no ZIKV-480, igual ao que ele imprime).
            raiz = canonical_bipartition(interno, todos)
            ramos = biparticoes_origem - {raiz}
            ids_ramos = {str(canonical_item_id(b)) for b in ramos}
            datados = len(ids_ramos & set(por_clado))
            rel["percentual_de_ramos_internos_colapsados"] = (
                round(100.0 * (len(ramos) - datados) / len(ramos), 4) if ramos else None)
            rel["contagem_de_clados"] = dict(
                contagem, ramos_internos_do_grupo_interno_na_origem=len(ramos),
                ramos_internos_datados=datados)
            tool_runs.registrar_metodo(METODO_NO_MANIFESTO, "executado", saida=timetree,
                                       alinhador=rotulo.aligner)
            rel["estado"] = "executado"
            rel["timetree"] = os.path.relpath(timetree, self.destino)
            logging.info(f"Relógio molecular — {nome}: {res['linha_de_resultado_lsd2']}")
            return rel
        except Exception as e:  # noqa: BLE001 — falha declarada, nunca silenciosa
            logging.error(f"Relógio molecular — erro em {nome}: {e}", exc_info=True)
            return self._falha(rel, timetree, rotulo.aligner, f"{type(e).__name__}: {e}")

    def _escolher_principal(self, resumo: Dict) -> None:
        """
        Escolhe a árvore cuja datação vai para o topo do JSON.

        Critério declarado, não default escondido: a pedida em
        `molecular_clock_tree`; senão, a do alinhador do projeto
        (`align_method`) antes das outras, e IQ-TREE > RAxML-NG > FastTree
        (a ordem de `METODOS_DATAVEIS`). Se a pedida falhou, o topo é falha —
        trocar em silêncio pela próxima seria responder outra pergunta.
        """
        por = resumo["por_arvore"]
        pedida = getattr(self, 'molecular_clock_tree', None)
        if pedida:
            criterio = f"molecular_clock_tree = {pedida!r}"
            escolhida = por.get(pedida)
            if escolhida is None:
                resumo.update(estado="tentado_e_falhou", criterio_da_arvore_principal=criterio,
                              motivo=f"a árvore pedida {pedida!r} não foi datada "
                                     "(ausente ou fora dos métodos datáveis)")
                return
        else:
            alinhador = getattr(self, 'align_method', None)
            ordem = sorted(
                por.values(),
                key=lambda r: (r["alinhador"] != alinhador,
                               METODOS_DATAVEIS.index(r["metodo_de_inferencia"]),
                               r["arvore"]))
            ok = [r for r in ordem if r["estado"] == "executado"]
            criterio = (f"alinhador do projeto ({alinhador!r}) primeiro; depois "
                        f"{' > '.join(METODOS_DATAVEIS)}; primeira árvore datada com sucesso")
            escolhida = ok[0] if ok else None
            if escolhida is None:
                motivos = "; ".join(f"{r['arvore']}: {r['motivo']}" for r in ordem)
                resumo.update(estado="tentado_e_falhou", criterio_da_arvore_principal=criterio,
                              motivo=f"nenhuma árvore foi datada — {motivos}")
                return

        resumo["criterio_da_arvore_principal"] = criterio
        resumo["arvore_principal"] = escolhida["arvore"]
        if escolhida["estado"] != "executado":
            resumo.update(estado="tentado_e_falhou",
                          motivo=f"{escolhida['arvore']}: {escolhida['motivo']}")
            return
        for k in ("taxa_subs_por_sitio_por_ano", "intervalo_da_taxa", "tmrca", "tmrca_lsd2",
                  "intervalo_do_tmrca", "solucao_unica",
                  "percentual_de_ramos_internos_colapsados", "datas_por_clado"):
            resumo[k] = escolhida.get(k)
        origem = os.path.join(self.destino, escolhida["timetree"])
        shutil.copyfile(origem, os.path.join(self.destino, "timetree.nex"))
        resumo["timetree"] = "timetree.nex"
        resumo["estado"] = "executado"

    def _gravar(self, resumo: Dict) -> Dict:
        os.makedirs(self.destino, exist_ok=True)
        tmp = self.json_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(resumo, fh, indent=2, ensure_ascii=False)
        os.replace(tmp, self.json_path)
        logging.info(f"Relógio molecular: {resumo['estado']} — {self.json_path}")
        return resumo


# --------------------------------------------------------------------------- #
# CLI — E12.11: a mesma saída, no mesmo lugar, sem passar pelo workflow
# --------------------------------------------------------------------------- #
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Datação leve (LSD2) das árvores de um out/ já executado (E12.10/E12.11).")
    ap.add_argument("--out", required=True, help="diretório out/ do projeto")
    ap.add_argument("--outgroup", default=None, help="acessos do grupo externo, separados por vírgula")
    ap.add_argument("--outgroup-file", default=None, help="arquivo com um acesso por linha")
    ap.add_argument("--tree", default=None, help="árvore principal (nome em out/Trees)")
    ap.add_argument("--model", default=None)
    ap.add_argument("--ci", type=int, default=None)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--blfix", action="store_true", help="não reotimizar comprimentos de ramo")
    a = ap.parse_args(argv)

    cfg: Dict = {}
    backup = os.path.join(a.out, "outputs", "config_backup.json")
    if os.path.isfile(backup):
        with open(backup) as fh:
            cfg = dict(json.load(fh).get("tree_config", {}))
    cfg.update(output_path=a.out, run_molecular_clock=True, via="cli")
    if a.outgroup_file or a.outgroup:
        cfg["molecular_clock_outgroup"] = a.outgroup_file or a.outgroup
    for chave, valor in (("molecular_clock_tree", a.tree), ("molecular_clock_model", a.model),
                         ("molecular_clock_ci", a.ci), ("molecular_clock_threads", a.threads)):
        if valor is not None:
            cfg[chave] = valor
    if a.blfix:
        cfg["molecular_clock_reoptimize_branches"] = False

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    resumo = MolecularClockController(**cfg)()
    print(json.dumps({k: resumo.get(k) for k in (
        "estado", "motivo", "arvore_principal", "taxa_subs_por_sitio_por_ano",
        "intervalo_da_taxa", "tmrca", "tmrca_lsd2", "intervalo_do_tmrca", "solucao_unica",
        "percentual_de_ramos_internos_colapsados")}, indent=2, ensure_ascii=False))
    print(f"clados datados: {len(resumo.get('datas_por_clado') or {})}")
    print(f"saída: {os.path.join(a.out, 'outputs', 'molecular_clock', 'relogio_molecular.json')}")

    # Achado da fila de triagem (DEC-129): disparada por fora de workflow.py, a
    # etapa nunca escrevia no manifesto da execução original — SHA-256 e
    # tools_invoked da datação isolada não existiam em lugar nenhum. Anexa ao
    # manifesto já existente (se houver um) sem regravá-lo do zero; ver a
    # docstring de anexar_execucao_isolada sobre por que não gerar um novo.
    destino = os.path.join(a.out, "outputs", "molecular_clock")
    anexado = manifest.anexar_execucao_isolada(a.out, destino)
    if not anexado:
        logging.info(
            "Sem manifesto em %s — datação isolada não entrou na proveniência "
            "de execução (projeto anterior a M2.5, ou nunca gerado).",
            os.path.join(a.out, "outputs", manifest.MANIFEST_FILENAME))

    return 0 if resumo.get("estado") == "executado" else 1


if __name__ == "__main__":
    sys.exit(main())
