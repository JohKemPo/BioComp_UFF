"""
Tabelas e figuras derivadas da análise de estabilidade entre pipelines.

As figuras seguem um único sistema visual: paleta categórica validada para
deficiência de visão de cores (azul `#2a78d6`, laranja `#eb6834`), rampa
sequencial de matiz único para magnitude, eixos recessivos e rótulos diretos em
vez de legendas quando há poucas séries.
"""

from __future__ import annotations

import csv
import os
from typing import Callable, Dict, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from workflow.stability.stability import StabilityAnalyzer

__all__ = ["StabilityReport", "SERIES", "INK"]

#: Slots categóricos (paleta validada; ordem fixa, nunca ciclada).
SERIES = {"canonical": "#2a78d6", "legacy": "#eb6834"}

#: Tinta de texto e elementos recessivos.
INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#8a8880", "grid": "#e4e3df"}

#: Rampa sequencial de matiz único (azul, claro -> escuro).
SEQUENTIAL = LinearSegmentedColormap.from_list(
    "phylo_blue",
    ["#f4f8fe", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
)


def _style_axes(ax) -> None:
    """Aplica o estilo comum: sem molduras supérfluas, grade recessiva."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK["grid"])
    ax.tick_params(colors=INK["secondary"], length=3, labelsize=8)
    ax.set_axisbelow(True)


class StabilityReport:
    """
    Materializa os resultados de um `StabilityAnalyzer` em CSV e PNG.

    Parameters
    ----------
    analyzer : StabilityAnalyzer
        Análise já construída.
    output_path : str
        Diretório onde tabelas e figuras serão gravadas (criado se ausente).
    label : str
        Rótulo do experimento, usado em títulos e nomes de arquivo.
    """

    def __init__(self, analyzer: StabilityAnalyzer, output_path: str, label: str) -> None:
        self.analyzer = analyzer
        self.output_path = output_path
        self.label = label
        os.makedirs(output_path, exist_ok=True)

    def _path(self, filename: str) -> str:
        return os.path.join(self.output_path, filename)

    # ------------------------------------------------------------------ #
    # Tabelas
    # ------------------------------------------------------------------ #

    def write_clade_table(self, classifier: Optional[Callable[[str], str]] = None) -> str:
        """
        Grava um CSV com um clado por linha, ordenado por suporte metodológico.

        Parameters
        ----------
        classifier : callable, optional
            Mapeia o nome de um terminal em um rótulo taxonômico; quando
            fornecido, a composição do clado é incluída.

        Return
        ------
        str
            Caminho do arquivo gravado.
        """
        path = self._path("clade_support.csv")
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            header = ["size", "n_pipelines", "support", "pipelines", "taxa"]
            if classifier:
                header += ["composition", "majority_label", "purity"]
            writer.writerow(header)

            for record in self.analyzer.clade_records():
                row = [
                    record.size,
                    len(record.pipelines),
                    round(record.support, 4),
                    ";".join(sorted(record.pipelines)),
                    ";".join(sorted(record.taxa)),
                ]
                if classifier:
                    composition: Dict[str, int] = {}
                    for taxon in record.taxa:
                        key = classifier(taxon)
                        composition[key] = composition.get(key, 0) + 1
                    majority, count = max(composition.items(), key=lambda kv: kv[1])
                    row += [
                        ";".join(f"{k}:{v}" for k, v in sorted(composition.items())),
                        majority,
                        round(count / record.size, 4),
                    ]
                writer.writerow(row)
        return path

    def write_pattern_table(self, min_support: float = 0.5) -> str:
        """
        Grava os conjuntos maximais de clados co-recuperados.

        Parameters
        ----------
        min_support : float, optional
            Suporte mínimo dos padrões exportados.

        Return
        ------
        str
            Caminho do arquivo gravado.
        """
        path = self._path("maximal_patterns.csv")
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["support", "n_pipelines", "n_clades", "clade_sizes", "pipelines"])
            for pattern in self.analyzer.maximal_patterns(min_support):
                writer.writerow(
                    [
                        round(pattern.support, 4),
                        len(pattern.pipelines),
                        pattern.size,
                        ";".join(str(len(c)) for c in pattern.clades),
                        ";".join(sorted(pattern.pipelines)),
                    ]
                )
        return path

    def write_rf_table(self) -> str:
        """
        Grava a matriz Robinson-Foulds normalizada entre pipelines.

        Return
        ------
        str
            Caminho do arquivo gravado.
        """
        path = self._path("rf_matrix.csv")
        matrix = self.analyzer.rf_matrix(normalized=True)
        names = list(matrix)
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["pipeline"] + names)
            for a in names:
                # RF indefinida vira célula vazia, nunca 0: `pandas` lê como NaN
                # e nenhuma média a confunde com "topologias idênticas" (D3).
                writer.writerow([a] + [
                    "" if matrix[a][b] is None else round(matrix[a][b], 4)
                    for b in names
                ])
        return path

    # ------------------------------------------------------------------ #
    # Figuras
    # ------------------------------------------------------------------ #

    def plot_support_profile(self) -> str:
        """
        Barras do número de clados recuperados por pelo menos k pipelines,
        contrastando a identidade canônica com a legada.

        Return
        ------
        str
            Caminho da figura gravada.
        """
        profile = self.analyzer.support_profile()
        legacy = self.analyzer.legacy_support_profile()
        levels = sorted(profile, reverse=True)

        canonical_values = [profile[k]["cumulative"] for k in levels]
        legacy_values = [legacy.get(k, 0) for k in levels]

        x = range(len(levels))
        width = 0.38
        fig, ax = plt.subplots(figsize=(6.2, 3.4), dpi=200)
        _style_axes(ax)
        ax.grid(axis="y", color=INK["grid"], linewidth=0.6)

        ax.bar([i - width / 2 for i in x], canonical_values, width * 0.94,
               color=SERIES["canonical"], label="Canonical clade identity")
        ax.bar([i + width / 2 for i in x], legacy_values, width * 0.94,
               color=SERIES["legacy"], label="Legacy 16-bit ordered identity")

        for i, (a, b) in enumerate(zip(canonical_values, legacy_values)):
            ax.text(i - width / 2, a, str(a), ha="center", va="bottom",
                    fontsize=7, color=INK["secondary"])
            ax.text(i + width / 2, b, str(b), ha="center", va="bottom",
                    fontsize=7, color=INK["secondary"])

        ax.set_xticks(list(x))
        ax.set_xticklabels([str(k) for k in levels])
        ax.set_xlabel("Recovered by at least k pipelines", fontsize=9, color=INK["secondary"])
        ax.set_ylabel("Number of clades", fontsize=9, color=INK["secondary"])
        ax.set_title(f"Cross-pipeline clade support — {self.label}",
                     fontsize=10, color=INK["primary"], loc="left")
        ax.legend(frameon=False, fontsize=8, labelcolor=INK["secondary"], loc="upper left")

        path = self._path("fig_support_profile.png")
        fig.tight_layout()
        fig.savefig(path, facecolor="white")
        plt.close(fig)
        return path

    def plot_rf_heatmap(self) -> str:
        """
        Mapa de calor das distâncias Robinson-Foulds normalizadas.

        Return
        ------
        str
            Caminho da figura gravada.
        """
        matrix = self.analyzer.rf_matrix(normalized=True)
        names = list(matrix)
        # NaN para par sem RF definida: o matplotlib deixa a célula em branco em
        # vez de pintá-la como se fosse distância zero.
        values = [[float("nan") if matrix[a][b] is None else matrix[a][b]
                   for b in names] for a in names]

        fig, ax = plt.subplots(figsize=(5.6, 4.8), dpi=200)
        image = ax.imshow(values, cmap=SEQUENTIAL, vmin=0, vmax=1)

        ax.set_xticks(range(len(names)))
        ax.set_yticks(range(len(names)))
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7, color=INK["secondary"])
        ax.set_yticklabels(names, fontsize=7, color=INK["secondary"])
        ax.tick_params(length=0)
        for side in ax.spines.values():
            side.set_visible(False)

        for i in range(len(names)):
            for j in range(len(names)):
                value = values[i][j]
                if value != value:   # NaN: par sem distância definida
                    ax.text(j, i, "—", ha="center", va="center", fontsize=6.5,
                            color=INK["secondary"])
                    continue
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6.5,
                        color="white" if value > 0.55 else INK["primary"])

        bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
        bar.set_label("Normalised Robinson–Foulds distance", fontsize=8, color=INK["secondary"])
        bar.ax.tick_params(labelsize=7, colors=INK["secondary"], length=0)
        bar.outline.set_visible(False)

        ax.set_title(f"Topological distance between pipelines — {self.label}",
                     fontsize=10, color=INK["primary"], loc="left")

        path = self._path("fig_rf_heatmap.png")
        fig.tight_layout()
        fig.savefig(path, facecolor="white")
        plt.close(fig)
        return path

    def plot_size_vs_support(self, classifier: Callable[[str], str], min_size: int = 3) -> str:
        """
        Dispersão do tamanho do clado contra seu suporte metodológico,
        separando clados taxonomicamente puros dos mistos.

        Parameters
        ----------
        classifier : callable
            Mapeia terminal em rótulo taxonômico.
        min_size : int, optional
            Tamanho mínimo do clado exibido.

        Return
        ------
        str
            Caminho da figura gravada.
        """
        rows = self.analyzer.taxonomic_purity(classifier, min_size=min_size)
        pure = [r for r in rows if r["monophyletic_group"]]
        mixed = [r for r in rows if not r["monophyletic_group"]]

        fig, ax = plt.subplots(figsize=(6.2, 3.8), dpi=200)
        _style_axes(ax)
        ax.grid(axis="both", color=INK["grid"], linewidth=0.6)

        for group, color, label in (
            (pure, SERIES["canonical"], "Single-species clade"),
            (mixed, SERIES["legacy"], "Mixed-species clade"),
        ):
            ax.scatter(
                [r["support"] for r in group],
                [r["size"] for r in group],
                s=34, color=color, alpha=0.85, edgecolor="white", linewidth=0.8, label=label,
            )

        ax.set_yscale("log")
        ax.set_xlabel("Cross-pipeline support", fontsize=9, color=INK["secondary"])
        ax.set_ylabel("Clade size (taxa, log scale)", fontsize=9, color=INK["secondary"])
        ax.set_title(f"Taxonomic coherence versus methodological support — {self.label}",
                     fontsize=10, color=INK["primary"], loc="left")
        ax.legend(frameon=False, fontsize=8, labelcolor=INK["secondary"], loc="upper left")

        path = self._path("fig_size_vs_support.png")
        fig.tight_layout()
        fig.savefig(path, facecolor="white")
        plt.close(fig)
        return path

    # ------------------------------------------------------------------ #

    def run_all(self, classifier: Optional[Callable[[str], str]] = None,
                min_support: float = 0.5) -> Dict[str, str]:
        """
        Produz todas as tabelas e figuras disponíveis.

        Parameters
        ----------
        classifier : callable, optional
            Classificador taxonômico; habilita a figura de coerência.
        min_support : float, optional
            Suporte mínimo dos padrões maximais exportados.

        Return
        ------
        dict
            Mapeamento nome do artefato -> caminho gravado.
        """
        artifacts = {
            "clade_table": self.write_clade_table(classifier),
            "pattern_table": self.write_pattern_table(min_support),
            "rf_table": self.write_rf_table(),
            "fig_support": self.plot_support_profile(),
            "fig_rf": self.plot_rf_heatmap(),
        }
        if classifier:
            artifacts["fig_purity"] = self.plot_size_vs_support(classifier)
        return artifacts
