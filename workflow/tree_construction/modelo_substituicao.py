"""
Seleção e tradução do modelo de substituição — M7.3, opção C′.

Até M7.3, IQ-TREE, RAxML-NG e MrBayes rodavam **GTR+Γ4 fixo** (`-m GTR+G`,
`--model GTR+G`, `lset nst=6 rates=gamma`) e nenhum dos três selecionava
modelo (`docs/science/11-auditoria-m7-inferencia.md §3`). A opção C′ aprovada:

1. **Uma** seleção por alinhamento, numa chamada separada do IQ-TREE
   (`-m MF`: só ModelFinder, sem árvore), **antes** de qualquer método baseado
   em modelo. Não `-m MFP` dentro da busca: medido no Zika-20 que ModelFinder
   na frente da busca muda a topologia **mesmo quando escolhe o mesmo modelo**
   (RF 2, ramo de UFBoot 45 — `11-auditoria §3.3`).
2. Candidatos restritos ao que os **três** métodos conseguem rodar:
   `--mset mrbayes` (bases JC/F81, K2P/HKY, SYM/GTR — `nst` 1, 2, 6) e
   `-mrate E,I,G,I+G`. Sem o `-mrate`, o ModelFinder também testa FreeRate
   (`+R2`, `+R3`…, medido na lista de candidatos do Zika-20), que o MrBayes não
   implementa — o próprio IQ-TREE, ao escrever o bloco do MrBayes, troca `+R`
   por `+I+G` em silêncio.
3. O modelo escolhido pelo BIC é **traduzido** para cada ferramenta por uma
   tabela fechada e testada (`TABELA_TRADUCAO`). Passar a string adiante troca
   o modelo em silêncio — ver as três armadilhas abaixo.

Três armadilhas medidas (IQ-TREE 3.1.3, RAxML-NG 2.0.2, MrBayes 3.2.7)
-----------------------------------------------------------------------
- **O nome que o ModelFinder imprime não é o nome que o `-m` lê.** Na lista do
  ModelFinder com `--mset mrbayes`, `GTR`, `HKY` e `F81` *sem* `+F` são as
  variantes de **frequências iguais** (log-verossimilhança idêntica à de `SYM`,
  `K2P`, `JC`). Passados de volta a `-m`, os mesmos nomes rodam com
  frequências **empíricas** (`-m GTR` → `Model of substitution: GTR+F`).
  Medido no VARV-6: o BIC escolheu `GTR` (iguais, lnL −541 924), e o próprio
  relatório do IQ-TREE declara `GTR+F` (lnL −547 174) — e o bloco MrBayes que
  ele escreve também. Por isso a leitura aqui canoniza para `SYM`/`K2P`/`JC`.
- **RAxML-NG não usa os nomes do IQ-TREE.** `K2P` não existe (`K80`);
  `K3Pu` não existe (`K81uf`); `TIM2` existe e é **outro** modelo (frequências
  iguais). Fora do conjunto `--mset mrbayes` não há tradução aqui: recusa.
- **Frequências de base não são harmonizadas entre ferramentas.** Com `+F`, o
  IQ-TREE usa contagens empíricas; o RAxML-NG, sem sufixo, estima por ML
  (`+FO`); o MrBayes amostra sob priori Dirichlet(1,1,1,1). É o que já
  acontecia com GTR+G antes de C′, e é o ponto (1) em aberto de
  `11-auditoria §3.6` — harmonizar mudaria a árvore do RAxML-NG e exige tabela
  de diff própria. Por isso a divergência é **declarada** por ferramenta em
  `frequencias`, não escondida. Com frequências iguais, as três coincidem.
"""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional

from Bio import AlignIO

from workflow.utils import tool_runs
from workflow.utils.external_tools import require_tool

__all__ = [
    "POLITICAS_SELECAO_MODELO", "POLITICA_SELECAO_PADRAO", "TAXAS_CANDIDATAS",
    "ModeloCanonico", "ModeloIntraduzivel", "SelecaoDeModeloFalhou",
    "TABELA_TRADUCAO", "MODELO_LEGADO", "interpretar_nome_modelfinder",
    "traduzir", "comando_selecao", "ler_relatorio_modelfinder", "selecionar_modelo",
]

#: Política de modelo do `tree_config` (chave `model_selection`).
#:
#: - ``bic_mrbayes`` — C′: seleção por BIC entre os modelos que IQ-TREE,
#:   RAxML-NG e MrBayes conseguem rodar, uma vez por alinhamento.
#: - ``nenhuma`` — o comportamento anterior a M7.3, literal: `GTR+G` no IQ-TREE
#:   e no RAxML-NG, `lset nst=6 rates=gamma` no MrBayes. Existe para reproduzir
#:   execuções antigas e para medir o diff, não como alternativa metodológica.
POLITICAS_SELECAO_MODELO = ("bic_mrbayes", "nenhuma")
POLITICA_SELECAO_PADRAO = "bic_mrbayes"

#: Heterogeneidade de taxa candidata: uniforme, +I, +Γ4, +I+Γ4. Sem FreeRate.
TAXAS_CANDIDATAS = "E,I,G,I+G"

#: Teto da seleção. Medido: 1,5 s no Zika-20, 28 s no VARV-49 (1 thread).
TIMEOUT_SELECAO_PADRAO_S = 3600


class ModeloIntraduzivel(ValueError):
    """Nome de modelo fora da tabela: não se passa adiante por aproximação."""


class SelecaoDeModeloFalhou(RuntimeError):
    """A seleção não produziu um modelo utilizável (binário, tempo, relatório).

    É levantada por **cada** método baseado em modelo que dependia dela, e o
    M7.6 a registra como ``tentado_e_falhou`` com este motivo. Não há recuo
    para GTR+G: seria a mesma substituição silenciosa de D1, em outro fator.
    """


# --------------------------------------------------------------------------- #
# Forma canônica e tradução
# --------------------------------------------------------------------------- #

#: Base do ModelFinder → (nst, frequências **quando não há `+F`**).
#: `None` = a base tem frequências livres, e o que decide é a presença de `+F`
#: (sem ela, na saída do ModelFinder, são iguais — primeira armadilha).
_BASES_MODELFINDER = {
    "JC": (1, "iguais"), "F81": (1, None),
    "K2P": (2, "iguais"), "K80": (2, "iguais"), "HKY": (2, None),
    "SYM": (6, "iguais"), "GTR": (6, None),
}

#: Tokens de taxa do ModelFinder → forma canônica.
_TAXAS = {(): "", ("I",): "+I", ("G4",): "+G4", ("I", "G4"): "+I+G4"}


@dataclass(frozen=True)
class ModeloCanonico:
    """Modelo de substituição na forma que as três ferramentas compartilham.

    Attributes
    ----------
    nst : int
        Número de classes de taxa de substituição: 1 (JC/F81), 2 (K2P/HKY),
        6 (SYM/GTR) — o `nst` do MrBayes.
    frequencias : str
        ``"iguais"`` ou ``"empiricas"`` (o `+F` do IQ-TREE).
    taxas : str
        ``""``, ``"+I"``, ``"+G4"`` ou ``"+I+G4"``.
    """
    nst: int
    frequencias: str
    taxas: str


#: Nome canônico por ferramenta, por (nst, frequências). IQ-TREE com o `+F`
#: explícito (sem ele, `-m GTR` rodaria GTR+F de qualquer jeito, mas `-m SYM`
#: ≠ `-m GTR`); RAxML-NG com o nome que o `--parse` aceita.
_NOME_IQTREE = {(1, "iguais"): "JC", (1, "empiricas"): "F81+F",
                (2, "iguais"): "K2P", (2, "empiricas"): "HKY+F",
                (6, "iguais"): "SYM", (6, "empiricas"): "GTR+F"}
_NOME_RAXML = {(1, "iguais"): "JC", (1, "empiricas"): "F81",
               (2, "iguais"): "K80", (2, "empiricas"): "HKY",
               (6, "iguais"): "SYM", (6, "empiricas"): "GTR"}
_TAXA_RAXML = {"": "", "+I": "+I", "+G4": "+G", "+I+G4": "+I+G"}
_TAXA_MRBAYES = {"": "rates=equal", "+I": "rates=propinv",
                 "+G4": "rates=gamma ngammacat=4",
                 "+I+G4": "rates=invgamma ngammacat=4"}

#: O que cada ferramenta faz com as frequências de base — declarado, porque
#: com `+F` as três **não** fazem a mesma coisa (terceira armadilha).
_FREQUENCIAS = {
    "iguais": {"iqtree": "iguais (fixas)", "raxml-ng": "iguais (fixas)",
               "mrbayes": "iguais (prset statefreqpr=fixed(equal))"},
    "empiricas": {"iqtree": "empíricas (+F, contagens do alinhamento)",
                  "raxml-ng": "estimadas por ML (+FO, padrão do RAxML-NG)",
                  "mrbayes": "amostradas sob priori Dirichlet(1,1,1,1) (padrão do MrBayes)"},
}


def _traducao(m: ModeloCanonico) -> Dict[str, object]:
    chave = (m.nst, m.frequencias)
    return {
        "iqtree": _NOME_IQTREE[chave] + m.taxas,
        "raxml-ng": _NOME_RAXML[chave] + _TAXA_RAXML[m.taxas],
        "mrbayes_lset": f"lset nst={m.nst} {_TAXA_MRBAYES[m.taxas]}",
        "mrbayes_prset": ("prset statefreqpr=fixed(equal)"
                          if m.frequencias == "iguais" else None),
        "nst": m.nst,
        "frequencias": dict(_FREQUENCIAS[m.frequencias]),
        "taxas": m.taxas or "uniforme",
    }


#: A tabela inteira: 3 `nst` × 2 frequências × 4 taxas = 24 modelos, que é
#: exatamente o espaço de `--mset mrbayes -mrate E,I,G,I+G` (medido: a lista de
#: candidatos do Zika-20 só tem estas bases e estas taxas). Cada entrada tem
#: teste com os três binários reais (`test_selecao_modelo`).
TABELA_TRADUCAO: Dict[ModeloCanonico, Dict[str, object]] = {
    ModeloCanonico(nst, freq, taxa): _traducao(ModeloCanonico(nst, freq, taxa))
    for nst in (1, 2, 6) for freq in ("iguais", "empiricas")
    for taxa in ("", "+I", "+G4", "+I+G4")
}

#: Política ``nenhuma``: os literais de antes de M7.3, **byte a byte** — não a
#: tradução de GTR+F+G4 (que difere no texto: `-m GTR+F+G4`, `ngammacat=4`).
MODELO_LEGADO: Dict[str, object] = {
    "iqtree": "GTR+G",
    "raxml-ng": "GTR+G",
    "mrbayes_lset": "lset nst=6 rates=gamma",
    "mrbayes_prset": None,
}


def interpretar_nome_modelfinder(nome: str) -> ModeloCanonico:
    """
    Lê um nome **como o ModelFinder o imprime** e devolve a forma canônica.

    Atenção: a semântica é a da saída do ModelFinder, não a do `-m` — `GTR`
    sem `+F` aqui é GTR de frequências iguais (= SYM). Ver o docstring do
    módulo.

    Raises
    ------
    ModeloIntraduzivel
        Base fora de `--mset mrbayes` (TIM2, K3Pu, TN…), FreeRate (`+R`),
        frequências `+FO`/`+FQ` explícitas ou qualquer token desconhecido.
    """
    tokens = (nome or "").strip().split("+")
    base, resto = tokens[0].upper(), tokens[1:]
    if base not in _BASES_MODELFINDER:
        raise ModeloIntraduzivel(
            f"modelo {nome!r}: base {tokens[0]!r} fora do conjunto que IQ-TREE, "
            f"RAxML-NG e MrBayes compartilham ({sorted(_BASES_MODELFINDER)})")
    nst, freq_sem_f = _BASES_MODELFINDER[base]
    tem_f = "F" in resto
    if resto.count("F") > 1:
        raise ModeloIntraduzivel(f"modelo {nome!r}: '+F' repetido")
    taxas = tuple(t for t in resto if t != "F")
    if taxas not in _TAXAS:
        raise ModeloIntraduzivel(
            f"modelo {nome!r}: componentes {'+'.join(taxas)!r} fora de "
            f"{{uniforme, +I, +G4, +I+G4}} (FreeRate e frequências FO/FQ "
            f"explícitas não têm tradução para as três ferramentas)")
    frequencias = "empiricas" if tem_f else (freq_sem_f or "iguais")
    return ModeloCanonico(nst, frequencias, _TAXAS[taxas])


def traduzir(nome_modelfinder: str) -> Dict[str, object]:
    """Tradução completa de um nome do ModelFinder (ver `TABELA_TRADUCAO`)."""
    canonico = interpretar_nome_modelfinder(nome_modelfinder)
    return {"modelfinder": nome_modelfinder, **TABELA_TRADUCAO[canonico]}


# --------------------------------------------------------------------------- #
# A seleção
# --------------------------------------------------------------------------- #

def comando_selecao(binario: str, alinhamento: str, prefixo: str, semente: int) -> List[str]:
    """
    Linha de comando da seleção. `-nt 1` e `-seed` pela mesma razão do
    `iqtree_constructor` (D11, D21): o ModelFinder avalia os modelos sobre uma
    árvore inicial, e a árvore inicial é o que a paralelização e a semente
    podem mudar. `-redo`: uma reexecução no mesmo diretório (métodos que
    faltaram na anterior) não pode esbarrar no *checkpoint* da seleção antiga.
    """
    return [binario, "-s", alinhamento, "-m", "MF", "--mset", "mrbayes",
            "-mrate", TAXAS_CANDIDATAS, "-seed", str(semente), "-nt", "1",
            "-pre", prefixo, "-redo", "-quiet"]


_RE_ESCOLHIDO = re.compile(r"^Best-fit model according to BIC:\s*(\S+)", re.MULTILINE)
_RE_LINHA_CANDIDATO = re.compile(
    r"^(\S+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+[+-]\s+\S+\s+(\d+\.\d+)\s+[+-]\s+\S+"
    r"\s+(\d+\.\d+)\s+[+-]\s+(\S+)\s*$")


def ler_relatorio_modelfinder(caminho: str) -> Dict[str, object]:
    """
    Modelo escolhido e tabela de candidatos do `.iqtree` de uma seleção `-m MF`.

    Return
    ------
    dict
        ``modelo_escolhido`` (str ou ``None`` se a linha não existe — nunca um
        palpite) e ``candidatos``: ``[{modelo, lnL, AIC, AICc, BIC, w_BIC}]``
        na ordem do relatório (crescente em BIC), sem as linhas repetidas que
        o IQ-TREE imprime para `SYM`/`K2P`/`JC`.
    """
    try:
        with open(caminho, errors="ignore") as fh:
            texto = fh.read()
    except OSError:
        return {"modelo_escolhido": None, "candidatos": []}
    achado = _RE_ESCOLHIDO.search(texto)
    candidatos, vistos = [], set()
    inicio = texto.find("List of models sorted by BIC scores")
    if inicio >= 0:
        for linha in texto[inicio:].splitlines()[1:]:
            if linha.startswith("AIC, w-AIC"):
                break
            m = _RE_LINHA_CANDIDATO.match(linha)
            if not m or m.group(1) in vistos:
                continue
            vistos.add(m.group(1))
            candidatos.append({
                "modelo": m.group(1), "lnL": float(m.group(2)),
                "AIC": float(m.group(3)), "AICc": float(m.group(4)),
                "BIC": float(m.group(5)), "w_BIC": float(m.group(6))})
    return {"modelo_escolhido": achado.group(1) if achado else None,
            "candidatos": candidatos}


def selecionar_modelo(alinhamento, diretorio: str, semente: int,
                      timeout_s: Optional[int] = None) -> Dict[str, object]:
    """
    Roda a seleção C′ sobre um alinhamento e devolve o modelo traduzido.

    Parameters
    ----------
    alinhamento : Bio.Align.MultipleSeqAlignment
        O mesmo objeto que os construtores recebem; gravado em PHYLIP como o
        `iqtree_constructor` grava, para que a seleção veja o mesmo dado.
    diretorio : str
        Onde ficam o `.phylip` e os arquivos do IQ-TREE.
    semente : int
    timeout_s : int, optional

    Return
    ------
    dict
        ``comando``, ``relatorio`` (caminho do `.iqtree`), ``modelo_escolhido``,
        ``candidatos`` e ``traducao`` (saída de `traduzir`).

    Raises
    ------
    SelecaoDeModeloFalhou
        Binário ausente, código de saída ≠ 0, tempo esgotado, relatório sem a
        linha do BIC, ou modelo escolhido sem tradução.
    """
    os.makedirs(diretorio, exist_ok=True)
    caminho_aln = os.path.join(diretorio, "alinhamento.phylip")
    prefixo = os.path.join(diretorio, "modelfinder")
    timeout_s = timeout_s or TIMEOUT_SELECAO_PADRAO_S
    try:
        AlignIO.write(alinhamento, caminho_aln, "phylip")
        cmd = comando_selecao(require_tool("iqtree"), caminho_aln, prefixo, semente)
    except (OSError, ValueError) as e:
        raise SelecaoDeModeloFalhou(f"{type(e).__name__}: {e}") from e

    resultado = {"comando": cmd, "relatorio": prefixo + ".iqtree"}
    # Registrada antes de rodar, como toda chamada (`tool_runs.registrar`,
    # Notes). Sem parâmetro no nível da ferramenta: `model`/`seed` de
    # `tools_invoked.iqtree` são da inferência; o que é desta chamada vai nela.
    tool_runs.registrar("iqtree", cmd, saida=resultado["relatorio"])
    tool_runs.anotar("iqtree", resultado["relatorio"],
                     etapa="seleção de modelo (-m MF, sem árvore) — M7.3",
                     seed=semente, threads=1)
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout_s)
    except subprocess.TimeoutExpired as e:
        raise SelecaoDeModeloFalhou(
            f"ModelFinder excedeu {timeout_s} s (model_selection_timeout_s)") from e
    except subprocess.CalledProcessError as e:
        cauda = (e.stderr or e.stdout or "").strip()[-300:]
        raise SelecaoDeModeloFalhou(
            f"ModelFinder terminou com código {e.returncode}: {cauda}") from e

    lido = ler_relatorio_modelfinder(resultado["relatorio"])
    resultado.update(lido)
    if not lido["modelo_escolhido"]:
        raise SelecaoDeModeloFalhou(
            "ModelFinder terminou, mas o relatório não traz "
            "'Best-fit model according to BIC' — sem modelo, nada a propagar")
    try:
        resultado["traducao"] = traduzir(lido["modelo_escolhido"])
    except ModeloIntraduzivel as e:
        raise SelecaoDeModeloFalhou(str(e)) from e
    return resultado
