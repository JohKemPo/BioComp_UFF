"""
Biblioteca de alinhadores: o que existe, o que cada um aguenta, e o que fazer
quando não aguenta.

O pipeline oferece três alinhadores múltiplos — **MAFFT**, **Clustal Omega** e
**MUSCLE**. Eles não são intercambiáveis: têm complexidades diferentes, limites
práticos diferentes, e um deles (o Clustal Omega) já foi morto pelo *OOM killer*
neste projeto.

Este módulo existe para que essas diferenças sejam **declaradas em um lugar só**
e possam ser consultadas por quem precisa delas — a UI, ao montar o experimento;
o controlador, antes de executar; e o manifesto, ao registrar o que foi feito.

## A regra que governa a substituição

[D1](../../../docs/science/02-defeitos-que-alteram-resultado.md#d1) é o defeito
mais caro do projeto, e ele **não** é "o Clustal Omega não aguenta 250 kb".
É que o controlador, ao descobrir isso, trocava para MAFFT **mantendo o nome de
arquivo `*_clustalo_*`**. Metade dos "pipelines" dos experimentos de *Variola*
são cópias byte a byte, e o fator alinhador — que é metade do delineamento —
não existe ali.

Daí as três regras:

1. **Declarar antes.** `viability` responde, para um dado conjunto, quais
   alinhadores são viáveis e **por que** os outros não são. Isso permite avisar
   no momento da configuração, quando ainda há escolha.
2. **Nunca substituir em silêncio.** Sem autorização explícita, um alinhador
   inviável é um **erro**, não um desvio. `AlignerPolicy.on_unavailable`
   controla isso, e o padrão é falhar.
3. **Nomear pelo que executou.** Se a substituição for autorizada, o resultado
   leva o nome do alinhador que **rodou**. Um arquivo chamado `clustalo` que
   contém MAFFT é pior que uma execução que falhou.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    "Aligner",
    "ALIGNERS",
    "AlignerPolicy",
    "available_aligners",
    "viability",
    "resolve_aligner",
]


@dataclass(frozen=True)
class Measurement:
    """
    Um ponto medido de consumo, com **as condições em que foi medido**.

    Um número sem suas condições é o defeito que este projeto passou dois meses
    corrigindo: o `20000` do Clustal Omega era um limite sem procedência, e a
    substituição silenciosa que ele disparava é [D1](../../../docs/science/02-defeitos-que-alteram-resultado.md#d1).

    Attributes
    ----------
    n_sequences, max_sequence_bp : int
        Dimensões do conjunto.
    peak_rss_bytes : int or None
        Pico de memória residente. `None` quando não foi instrumentado.
    seconds : float or None
        Tempo de parede.
    outcome : str
        ``"ok"``, ``"oom"``, ``"timeout"``, ``"crash"``.
    machine : str
        Descrição curta da máquina — é isto que torna o ponto interpretável.
    machine_bytes : int or None
        Memória física **daquela** máquina. Sem ela não é possível transferir a
        observação: "morreu com 19,4 GB" só significa alguma coisa ao lado de
        "numa máquina de 31 GB".
    date : str
        Data ISO.
    """

    n_sequences: int
    max_sequence_bp: int
    outcome: str
    machine: str
    date: str
    machine_bytes: Optional[int] = None
    peak_rss_bytes: Optional[int] = None
    seconds: Optional[float] = None


@dataclass(frozen=True)
class ResourceModel:
    """
    Modelo de custo de um alinhador, separando **o que é do algoritmo** do
    **que é da máquina**.

    A distinção é o achado de [R2](../../../docs/respostasUteis/r2.md): a *lei de
    escala* é propriedade do método e transfere entre arquiteturas; o *ponto onde
    a curva cruza o orçamento* é propriedade da máquina. Declarar só o segundo,
    como um limite em pares de base, produz um número que parece universal e não é.

    Attributes
    ----------
    scaling : str
        A lei declarada, em texto — ``"n*L"``, ``"n*L^2"``, ``"L"``. **É uma
        suposição enquanto `fitted` for False.**
    fitted : bool
        Se o expoente foi **ajustado** a partir de medições. Com pontos de uma
        máquina só, expoente e deslocamento ficam confundidos: qualquer curva
        passa por dois pontos. Ajustar exige bissectar em **duas máquinas**.
    bytes_per_unit : float or None
        Constante do modelo, em bytes por unidade da lei. `None` enquanto não
        houver medição suficiente.
    measurements : tuple of Measurement
        Os pontos observados. Crescem a cada máquina nova; nunca sobrescrevem.
    """

    scaling: str
    fitted: bool = False
    bytes_per_unit: Optional[float] = None
    bytes_is_lower_bound: bool = False
    measurements: Tuple["Measurement", ...] = ()

    def unidades(self, n_sequences: int, max_sequence_bp: int) -> float:
        """Tamanho do problema na unidade da lei declarada."""
        n, L = float(n_sequences), float(max_sequence_bp)
        return {
            "n*L": n * L,
            "n*L^2": n * L * L,
            "L": L,
            "n^2*L": n * n * L,
        }.get(self.scaling, n * L)

    def estimate_bytes(self, n_sequences: int, max_sequence_bp: int) -> Optional[int]:
        """
        Memória estimada, ou `None` quando não há base para estimar.

        `None` é resposta legítima e frequente: significa *"não sabemos"*, que é
        diferente de *"cabe"*. Quem consome precisa distinguir os dois.

        Quando `bytes_is_lower_bound` é True, o valor é um **piso**: veio de um
        pico observado no momento em que o processo foi morto, e o consumo real
        seria maior se ele tivesse continuado.
        """
        if self.bytes_per_unit is None:
            return None
        return int(self.bytes_per_unit * self.unidades(n_sequences, max_sequence_bp))

    def blocking_failure(self, n_sequences: int, max_sequence_bp: int,
                         budget_bytes: Optional[int],
                         tolerance: float = 0.95) -> Optional["Measurement"]:
        """
        Falha **observada** que condena este conjunto nesta máquina.

        Uma falha vale como veto quando ocorreu num conjunto **igual ou menor**,
        numa máquina de orçamento **igual ou maior**. Observação vence
        estimativa: o pico registrado no instante em que o processo morreu é um
        piso, não o requisito, e uma estimativa construída a partir dele
        **subestima por construção** — foi assim que o MUSCLE apareceu como
        viável justamente no conjunto em que o vimos morrer.

        `tolerance` existe porque memória registrada e memória reportada nunca
        batem exatamente: o mesmo hardware informa 33.424.216.064 bytes e a
        medição foi anotada como "31 GB". Comparar com `>=` estrito faria a
        máquina não reconhecer a própria falha por 24 MB de arredondamento.
        """
        if budget_bytes is None:
            return None
        for m in self.measurements:
            if m.outcome == "ok" or m.machine_bytes is None:
                continue
            if (m.n_sequences <= n_sequences and m.max_sequence_bp <= max_sequence_bp
                    and m.machine_bytes >= budget_bytes * tolerance):
                return m
        return None

    def largest_ok(self) -> Optional[int]:
        """Maior comprimento por sequência **verificado funcionando**."""
        oks = [m.max_sequence_bp for m in self.measurements if m.outcome == "ok"]
        return max(oks) if oks else None

    def smallest_failure(self) -> Optional[int]:
        """Menor comprimento por sequência em que se **observou** falha."""
        falhas = [m.max_sequence_bp for m in self.measurements if m.outcome != "ok"]
        return min(falhas) if falhas else None


@dataclass(frozen=True)
class Aligner:
    """
    Um alinhador múltiplo e os limites que ele impõe.

    Attributes
    ----------
    key : str
        Nome usado na configuração e nos nomes de arquivo (`mafft`, `clustalo`,
        `muscle`).
    label : str
        Nome de exibição.
    binary : str
        Executável procurado no PATH.
    max_sequence_bp : int or None
        Comprimento por sequência acima do qual o alinhador é considerado
        **inviável** neste projeto. `None` significa "sem limite conhecido".
    estrategia : str or None
        Estratégia pedida ao binário. Dois alinhadores podem apontar para o
        mesmo executável e diferir só nisto.
    max_sequences : int or None
        Número de sequências acima do qual é inviável.
    note : str
        Por que o limite existe. Vai para a UI: um limite sem explicação vira
        superstição.
    """

    key: str
    label: str
    binary: str
    max_sequence_bp: Optional[int] = None
    max_sequences: Optional[int] = None
    note: str = ""
    resources: Optional["ResourceModel"] = None
    #: Estratégia pedida ao binário, quando ele aceita mais de uma. É o que
    #: permite dois braços do fator alinhador compartilharem o mesmo executável
    #: e ainda assim produzirem alinhamentos diferentes (D1).
    estrategia: Optional[str] = None

    def installed(self) -> bool:
        return shutil.which(self.binary) is not None

    def version(self) -> Optional[str]:
        """Versão detectada, ou `None` se o binário não está no PATH."""
        if not self.installed():
            return None
        # Indexado pelo BINÁRIO, não pela chave: dois alinhadores podem
        # compartilhar executável e diferir só na estratégia (D1).
        for args, padrao in _VERSAO[self.binary]:
            try:
                r = subprocess.run([self.binary, *args], capture_output=True, text=True,
                                   stdin=subprocess.DEVNULL, timeout=10, check=False)
            except (OSError, subprocess.TimeoutExpired):
                continue
            achado = re.search(padrao, (r.stdout or "") + (r.stderr or ""))
            if achado:
                return achado.group(achado.lastindex or 0)
        return None


#: Como perguntar a versão de cada um. Nenhum respeita `--version` do mesmo jeito.
_VERSAO: Dict[str, Sequence[Tuple[Sequence[str], str]]] = {
    "mafft": ((["--version"], r"v[\d.]+"),),
    "clustalo": ((["--version"], r"[\d.]+"),),
    "muscle": ((["-version"], r"[\d.]+"), (["--version"], r"[\d.]+")),
}


#: Os três alinhadores da biblioteca.
#:
#: Os limites **não são teóricos**: o de Clustal Omega vem de uma falha observada
#: neste projeto — `return code 137`, o *OOM killer* do kernel, no conjunto
#: Zika479. Os demais estão declarados como desconhecidos até que a máquina de
#: validação os meça (M7.7). **Limite não medido é `None`, nunca um palpite.**
#: Onde as medições foram feitas. Faz parte do dado: um pico de 19,4 GB só é
#: interpretável ao lado da memória total da máquina que o mediu.
MAQUINA_DEV = "i5-11400H · 12 núcleos lógicos · 31 GB"
#: Valor reportado pelo próprio sistema, não arredondado — o arredondamento
#: é o que faz uma máquina deixar de reconhecer a própria medição.
MAQUINA_DEV_BYTES = 33_424_216_064

ALIGNERS: Dict[str, Aligner] = {
    "mafft": Aligner(
        key="mafft", label="MAFFT (FFT-NS-2)", binary="mafft",
        max_sequence_bp=None, max_sequences=None,
        estrategia="retree",
        note=("Escala por estratégia: L-INS-i para conjuntos pequenos, FFT-NS até "
              "10 000 sequências, PartTree acima disso. **Nenhuma falha observada** "
              "neste projeto — alinhou genomas de 250 kb. Isso não é o mesmo que "
              "'não tem limite': é ausência de observação."),
        resources=ResourceModel(
            scaling="n*L", fitted=False, bytes_per_unit=None,
            measurements=(
                Measurement(n_sequences=20, max_sequence_bp=10_788, outcome="ok",
                            seconds=4.9, machine=MAQUINA_DEV, machine_bytes=MAQUINA_DEV_BYTES, date="2026-08-25"),
                Measurement(n_sequences=52, max_sequence_bp=228_250, outcome="ok",
                            machine=MAQUINA_DEV, machine_bytes=MAQUINA_DEV_BYTES, date="2026-08-25"),
            ),
        ),
    ),
    # Segundo braço do fator alinhador, decidido pelo usuário em 2026-08-26
    # (decisão 1 / D1 parte 2). Mesma ferramenta, mesma versão, mesmo binário:
    # o que muda é o **algoritmo**, e é esse o contraste que E4 quer medir. As
    # duas alternativas — MUSCLE e Clustal Omega — foram remedidas com o env
    # pinado e não servem em genoma de poxvírus: ver as notas de cada uma.
    "mafft_iterative": Aligner(
        key="mafft_iterative", label="MAFFT (FFT-NS-i)", binary="mafft",
        max_sequence_bp=None, max_sequences=None,
        estrategia="iterative",
        note=("MAFFT com refinamento iterativo (`--retree 2 --maxiterate 1000`), "
              "contra o `mafft` progressivo (`--maxiterate 0`). É o fator alinhador "
              "em conjuntos onde nenhuma outra ferramenta roda — e o único que "
              "existe tanto em *Variola* quanto em Zika."),
        resources=ResourceModel(
            scaling="n*L", fitted=False, bytes_per_unit=None, measurements=(),
        ),
    ),
    "clustalo": Aligner(
        key="clustalo", label="Clustal Omega", binary="clustalo",
        max_sequence_bp=20_000, max_sequences=None,
        note=("**O limite é de tempo, não de memória** — medido em 2026-08-26 sobre "
              "52 sequências de até 228 kb: não terminou em **1 h** e o pico de RSS "
              "foi de apenas **220 MB**. A afirmação anterior, de que era morto pelo "
              "OOM killer neste porte, estava errada: o código 137 observado foi em "
              "Zika479 (478 sequências curtas), que é outro regime. O limite de 20 kb "
              "por sequência continua herdado do código original."),
        resources=ResourceModel(
            scaling="n*L", fitted=False, bytes_per_unit=None,
            measurements=(
                Measurement(n_sequences=20, max_sequence_bp=10_788, outcome="ok",
                            seconds=64.0, machine=MAQUINA_DEV, machine_bytes=MAQUINA_DEV_BYTES, date="2026-08-25"),
            ),
        ),
    ),
    "muscle": Aligner(
        key="muscle", label="MUSCLE", binary="muscle",
        max_sequence_bp=10_788, max_sequences=1_000,
        note=("O limite é o **maior comprimento verificado funcionando**, não um ponto "
              "de falha conhecido: 10,8 kb roda em 34,5 s; 228 kb consumiu 19,4 GB e foi "
              "morto pelo OOM killer numa máquina de 31 GB. O limiar real está entre os "
              "dois e não foi estreitado. **Numa máquina maior o mesmo conjunto pode "
              "concluir** — ver docs/respostasUteis/r2.md."),
        resources=ResourceModel(
            scaling="n*L", fitted=False,
            # 19,4 GB / (52 x 228 250) — uma constante a partir de UM ponto de
            # falha. Serve para ordem de grandeza e nada mais; `fitted=False`
            # é o que impede que ela seja lida como previsão.
            bytes_per_unit=19_385_640 * 1024 / (52 * 228_250),
            bytes_is_lower_bound=True,
            measurements=(
                Measurement(n_sequences=20, max_sequence_bp=10_788, outcome="ok",
                            seconds=34.5, machine=MAQUINA_DEV, machine_bytes=MAQUINA_DEV_BYTES, date="2026-08-25"),
                Measurement(n_sequences=52, max_sequence_bp=228_250, outcome="oom",
                            peak_rss_bytes=19_385_640 * 1024, machine=MAQUINA_DEV, machine_bytes=MAQUINA_DEV_BYTES,
                            date="2026-08-25"),
            ),
        ),
    ),
}


def memoria_disponivel_bytes() -> Optional[int]:
    """
    Memória física da máquina **que vai executar**, não uma constante compilada.

    É o orçamento contra o qual o requisito estimado é comparado. Consultá-la em
    execução é o que faz o mesmo código dar vereditos diferentes em máquinas
    diferentes — que é o comportamento correto (R2).
    """
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


@dataclass
class AlignerPolicy:
    """
    O que fazer quando o alinhador pedido não é viável.

    Attributes
    ----------
    on_unavailable : str
        ``"fail"`` (padrão) levanta erro com o motivo. ``"fallback"`` troca pelo
        primeiro alinhador viável de `fallback_order` — e **o resultado passa a
        levar o nome do alinhador que rodou**.
    fallback_order : sequence of str
        Ordem de preferência da substituição.
    """

    on_unavailable: str = "fail"
    fallback_order: Sequence[str] = ("mafft", "muscle", "clustalo")


@dataclass
class Viability:
    """
    Veredito sobre um alinhador para um conjunto concreto **nesta máquina**.

    `estimated_bytes` e `available_bytes` são devolvidos separados de propósito:
    a mensagem útil não é "indisponível", é *"precisa de ~19 GB e a máquina tem
    31"*. A primeira é um veto sem apelação; a segunda é um requisito, e diz ao
    usuário o que mudaria o veredito.
    """

    aligner: str
    viable: bool
    installed: bool
    version: Optional[str] = None
    reasons: List[str] = field(default_factory=list)
    estimated_bytes: Optional[int] = None
    available_bytes: Optional[int] = None
    estimate_is_fitted: bool = False

    def summary(self) -> Dict[str, object]:
        return {
            "aligner": self.aligner,
            "viable": self.viable,
            "installed": self.installed,
            "version": self.version,
            "reasons": list(self.reasons),
            "estimated_bytes": self.estimated_bytes,
            "available_bytes": self.available_bytes,
            "estimate_is_fitted": self.estimate_is_fitted,
        }


def available_aligners() -> Dict[str, Optional[str]]:
    """``chave -> versão`` de cada alinhador; `None` quando não está no PATH."""
    return {chave: a.version() for chave, a in ALIGNERS.items()}


def viability(n_sequences: int,
              max_sequence_bp: int,
              available_bytes: Optional[int] = None,
              margin: float = 0.8) -> Dict[str, Viability]:
    """
    Diz, para um conjunto concreto, quais alinhadores são viáveis e por quê não.

    Serve à UI no momento de montar o experimento — que é quando ainda há
    escolha — e ao controlador, antes de executar.

    Parameters
    ----------
    n_sequences : int
        Número de sequências do conjunto.
    max_sequence_bp : int
        Comprimento da **maior** sequência, em pares de base. É o comprimento
        máximo, e não a média, porque é uma sequência só que estoura a memória.
    available_bytes : int or None, optional
        Orçamento de memória. `None` consulta **a máquina atual** — é isso que
        faz o mesmo código dar vereditos diferentes em máquinas diferentes, que é
        o comportamento correto ([R2](../../../docs/respostasUteis/r2.md)).
    margin : float, optional
        Fração do orçamento considerada utilizável. O padrão de 0,8 deixa espaço
        para o resto do sistema.

    Return
    ------
    dict of str -> Viability
    """
    resultado: Dict[str, Viability] = {}
    orcamento = available_bytes if available_bytes is not None else memoria_disponivel_bytes()

    for chave, a in ALIGNERS.items():
        motivos: List[str] = []
        instalado = a.installed()
        if not instalado:
            motivos.append(f"{a.label} não está instalado (binário `{a.binary}` ausente do PATH)")

        estimado = None
        ajustado = False
        if a.resources is not None:
            estimado = a.resources.estimate_bytes(n_sequences, max_sequence_bp)
            ajustado = a.resources.fitted

        # Falha observada tem precedência sobre qualquer estimativa.
        falha = (a.resources.blocking_failure(n_sequences, max_sequence_bp, orcamento)
                 if a.resources is not None else None)
        if falha is not None:
            pico = f", pico de {falha.peak_rss_bytes / 1e9:.1f} GB" if falha.peak_rss_bytes else ""
            motivos.append(
                f"**observado falhar** ({falha.outcome}) em {falha.n_sequences} sequências de "
                f"{falha.max_sequence_bp:,} pb numa máquina de "
                f"{falha.machine_bytes / 1e9:.0f} GB{pico}, em {falha.date}. "
                f"Esta máquina tem {(orcamento or 0) / 1e9:.0f} GB".replace(",", "."))

        # Requisito contra orçamento: é o critério que generaliza entre máquinas.
        elif estimado is not None and orcamento is not None and estimado > orcamento * margin:
            if ajustado:
                qualidade = "estimativa ajustada"
            elif a.resources is not None and a.resources.bytes_is_lower_bound:
                qualidade = ("PISO, não estimativa: veio do pico observado quando o "
                             "processo foi morto, e o consumo real seria maior")
            else:
                qualidade = "ordem de grandeza, não ajustada"
            motivos.append(
                f"precisa de ao menos ~{estimado / 1e9:.1f} GB e esta máquina tem "
                f"{orcamento / 1e9:.1f} GB ({qualidade}). Numa máquina maior pode concluir")

        # Limite absoluto: rede de segurança **apenas quando não há estimativa**.
        # Onde existe modelo de custo, ele manda — um limite em pares de base é
        # justamente o que não generaliza entre máquinas (R2), e mantê-lo aqui
        # faria a ferramenta vetar numa máquina de 128 GB pelo que mediu numa de 31.
        if (estimado is None and a.max_sequence_bp is not None
                and max_sequence_bp > a.max_sequence_bp):
            motivos.append(
                f"sequência de {max_sequence_bp:,} pb acima do maior comprimento "
                f"verificado ({a.max_sequence_bp:,} pb) — {a.note}".replace(",", "."))

        if a.max_sequences is not None and n_sequences > a.max_sequences:
            motivos.append(
                f"{n_sequences} sequências excedem o limite de {a.max_sequences} — {a.note}")

        resultado[chave] = Viability(
            aligner=chave, viable=not motivos, installed=instalado,
            version=a.version() if instalado else None, reasons=motivos,
            estimated_bytes=estimado, available_bytes=orcamento,
            estimate_is_fitted=ajustado)

    return resultado


def resolve_aligner(requested: str,
                    n_sequences: int,
                    max_sequence_bp: int,
                    policy: Optional[AlignerPolicy] = None) -> Tuple[str, Optional[str]]:
    """
    Decide qual alinhador vai efetivamente rodar.

    Parameters
    ----------
    requested : str
        Alinhador pedido pela configuração do experimento.
    n_sequences, max_sequence_bp : int
        Dimensões do conjunto.
    policy : AlignerPolicy, optional
        O padrão **falha** em vez de substituir.

    Return
    ------
    tuple
        ``(alinhador efetivo, motivo da substituição ou None)``. Quando há
        substituição, o chamador **precisa** usar o nome devolvido para nomear a
        saída — é o que D1 descreve como o defeito, não a substituição em si.

    Raises
    ------
    ValueError
        Se o alinhador pedido não existe, ou se é inviável e a política é falhar.
    """
    politica = policy or AlignerPolicy()

    if requested not in ALIGNERS:
        raise ValueError(
            f"Alinhador desconhecido: '{requested}'. Disponíveis: {sorted(ALIGNERS)}")

    vereditos = viability(n_sequences, max_sequence_bp)
    pedido = vereditos[requested]

    if pedido.viable:
        return requested, None

    motivo = "; ".join(pedido.reasons)

    if politica.on_unavailable != "fallback":
        raise ValueError(
            f"{ALIGNERS[requested].label} não é viável para este conjunto: {motivo}. "
            f"Escolha outro alinhador ou autorize a substituição com "
            f"on_unavailable='fallback' — e saiba que a saída levará o nome do "
            f"alinhador que de fato rodar.")

    for alternativa in politica.fallback_order:
        if alternativa != requested and vereditos.get(alternativa, Viability(alternativa, False, False)).viable:
            return alternativa, (
                f"{ALIGNERS[requested].label} substituído por {ALIGNERS[alternativa].label}: {motivo}")

    raise ValueError(
        f"Nenhum alinhador viável para este conjunto ({n_sequences} sequências, "
        f"maior com {max_sequence_bp} pb). Motivo do pedido original: {motivo}")
