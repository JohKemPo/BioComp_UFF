"""
Estudo de caso: estabilidade metodológica em genomas de Orthopoxvirus.

Executa a análise de estabilidade sobre um ou mais projetos do PhyloTreeMiner,
grava tabelas e figuras em `<projeto>/out/outputs/stability/` e imprime um
resumo consolidado.

Uso
---
    python -m workflow.stability.case_study \\
        --project projects/Variola_Yu_li_2007_200seq \\
        --fasta   data/workflow_dataAcquisition_li_et_al_2007_replication-RetMax200/dataset_final.fasta

Sem `--project`, os três experimentos de Variola do repositório são analisados.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Callable, Dict, List, Optional

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, "../.."))

from workflow.stability.report import StabilityReport
from workflow.stability.stability import StabilityAnalyzer, TreeSet, strip_accession_version

__all__ = ["build_classifier", "analyse_project", "main"]

#: Palavras-chave usadas para inferir a espécie a partir da descrição do FASTA.
SPECIES_KEYWORDS = (
    ("variola", "VARV"),
    ("monkeypox", "MPXV"),
    ("camelpox", "CMLV"),
    ("taterapox", "TATV"),
    ("cowpox", "CPXV"),
    ("buffalopox", "BPXV"),
    ("vaccinia", "VACV"),
    ("yoka", "YOKA"),
    ("crocodilepox", "CROC"),
)

#: Experimentos analisados quando nenhum projeto é informado explicitamente.
DEFAULT_EXPERIMENTS = (
    (
        "VARV-121",
        "projects/Variola_Yu_li_2007_200seq",
        "data/workflow_dataAcquisition_li_et_al_2007_replication-RetMax200/dataset_final.fasta",
    ),
    (
        "VARV-49",
        "projects/Variola_Yu_li_2007",
        None,
    ),
    (
        "VARV-6-noITR",
        "projects/Variola_Yu_li_2007_noITRs_6seqs",
        None,
    ),
)


def build_classifier(fasta_path: Optional[str]) -> Callable[[str], str]:
    """
    Constrói um classificador de espécie a partir dos cabeçalhos de um FASTA.

    Parameters
    ----------
    fasta_path : str or None
        Caminho do FASTA cujos cabeçalhos descrevem cada acesso. Se None, o
        classificador devolve ``"unknown"`` para todos os terminais.

    Return
    ------
    callable
        Função que mapeia o nome de um terminal em uma sigla de espécie.
    """
    descriptions: Dict[str, str] = {}
    if fasta_path and os.path.exists(fasta_path):
        with open(fasta_path, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    header = line[1:].strip()
                    # As chaves usam a mesma normalização aplicada às árvores,
                    # de modo que acessos com e sem versão se reencontrem.
                    descriptions[strip_accession_version(header.split()[0])] = header.lower()

    def classify(terminal: str) -> str:
        description = descriptions.get(strip_accession_version(terminal), terminal.lower())
        for keyword, code in SPECIES_KEYWORDS:
            if keyword in description:
                return code
        return "unknown"

    return classify


def analyse_project(label: str, project_path: str, fasta_path: Optional[str],
                    min_support: float = 0.5) -> Dict[str, object]:
    """
    Analisa um projeto e grava seus artefatos.

    Parameters
    ----------
    label : str
        Rótulo curto do experimento.
    project_path : str
        Diretório do projeto (contendo `out/Trees`).
    fasta_path : str or None
        FASTA usado para a classificação taxonômica.
    min_support : float, optional
        Suporte mínimo dos padrões maximais exportados.

    Return
    ------
    dict
        Resumo numérico do experimento.
    """
    trees_dir = os.path.join(project_path, "out", "Trees")
    tree_set = TreeSet.from_directory(trees_dir)
    analyzer = StabilityAnalyzer(tree_set)
    classifier = build_classifier(fasta_path)

    output_path = os.path.join(project_path, "out", "outputs", "stability")
    report = StabilityReport(analyzer, output_path, label)
    artifacts = report.run_all(classifier=classifier, min_support=min_support)

    audit = analyzer.legacy_audit()
    effects = analyzer.factor_effects()
    universal = analyzer.consensus_clades(min_support=1.0, min_size=2)
    purity = analyzer.taxonomic_purity(classifier, min_size=3)
    universal_purity = [row for row in purity if row["n_pipelines"] == len(tree_set)]

    summary = {
        "label": label,
        "project": project_path,
        "n_taxa": tree_set.n_taxa,
        "n_pipelines": len(tree_set),
        "pipelines": sorted(tree_set.trees),
        "mangled_labels": tree_set.mangled_labels,
        "distinct_clades": len(analyzer.clade_records()),
        "universal_clades": len(universal),
        "universal_clade_sizes": [record.size for record in universal],
        "support_profile": analyzer.support_profile(),
        "legacy_support_profile": analyzer.legacy_support_profile(),
        "legacy_audit": audit.summary(),
        "factor_effects": effects,
        "universal_clades_min3": len(universal_purity),
        "universal_clades_min3_single_species": sum(
            1 for row in universal_purity if row["monophyletic_group"]
        ),
        "artifacts": artifacts,
    }

    with open(os.path.join(output_path, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    return summary


def _print_summary(summary: Dict[str, object]) -> None:
    """Imprime o resumo de um experimento em formato legível no terminal."""
    audit = summary["legacy_audit"]
    effects = summary["factor_effects"]

    print(f"\n{'=' * 72}")
    print(f"{summary['label']}  ({summary['n_taxa']} taxa, {summary['n_pipelines']} pipelines)")
    print(f"{'=' * 72}")
    mangled = {k: v for k, v in summary["mangled_labels"].items() if v}
    if mangled:
        print(f"  rótulos adulterados pelo inferidor ... "
              f"{ {k: len(v) for k, v in mangled.items()} }")
        print(f"      táxons afetados .................. {sorted(set(sum(mangled.values(), [])))}")
    print(f"  clados distintos ..................... {summary['distinct_clades']}")
    print(f"  clados em todos os pipelines ......... {summary['universal_clades']} "
          f"(tamanhos: {summary['universal_clade_sizes'][:8]})")
    print(f"  destes, com >=3 taxa ................. {summary['universal_clades_min3']} "
          f"({summary['universal_clades_min3_single_species']} de espécie única)")
    print(f"  identidade legada: itens ............. {audit['legacy_items']} "
          f"para {audit['canonical_clades']} clados reais")
    print(f"      fragmentados por ordem ........... {audit['fragmented_clades']} "
          f"({audit['fragmentation_rate'] * 100:.1f}%)")
    print(f"      itens em colisão ................. {audit['colliding_items']} "
          f"({audit['collision_rate'] * 100:.1f}%)")
    print(f"  RF médio | trocando alinhador ........ {effects['aligner']['mean']:.3f} "
          f"(n={effects['aligner']['n']}, máx {effects['aligner']['max']:.3f})")
    print(f"  RF médio | trocando inferência ....... {effects['inference']['mean']:.3f} "
          f"(n={effects['inference']['n']}, máx {effects['inference']['max']:.3f})")


def main(argv: Optional[List[str]] = None) -> int:
    """
    Ponto de entrada de linha de comando.

    Parameters
    ----------
    argv : list of str, optional
        Argumentos; usa `sys.argv` quando omitido.

    Return
    ------
    int
        Código de saída (0 em sucesso).
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project", action="append", default=None,
                        help="Diretório do projeto a analisar (repetível).")
    parser.add_argument("--fasta", action="append", default=None,
                        help="FASTA correspondente, na mesma ordem de --project.")
    parser.add_argument("--label", action="append", default=None,
                        help="Rótulo do experimento, na mesma ordem de --project.")
    parser.add_argument("--min-support", type=float, default=0.5,
                        help="Suporte mínimo dos padrões maximais exportados.")
    args = parser.parse_args(argv)

    if args.project:
        fastas = args.fasta or [None] * len(args.project)
        labels = args.label or [os.path.basename(p.rstrip("/")) for p in args.project]
        experiments = list(zip(labels, args.project, fastas))
    else:
        experiments = list(DEFAULT_EXPERIMENTS)

    summaries = []
    for label, project, fasta in experiments:
        summaries.append(analyse_project(label, project, fasta, args.min_support))
        _print_summary(summaries[-1])

    print(f"\nArtefatos gravados em <projeto>/out/outputs/stability/.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
