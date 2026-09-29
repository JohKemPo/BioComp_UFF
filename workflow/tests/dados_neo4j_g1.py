"""
Fixture sintética para G1 (`docs/science/05-grafo-neo4j.md §G1`).

Três terminais repetidos em várias subárvores de três árvores, cobrindo os
casos que o gerador precisa tratar ao deduplicar `Metadata`:

- `T1` aparece primeiro **sem** registro GenBank (lookup falhou) e só
  depois completo, e volta a aparecer vazio numa árvore posterior;
- `T2` idem, e tem uma feature sem qualifiers e `strand` ausente;
- `T3` aparece primeiro completo, não tem nenhuma feature, e tem uma
  ocorrência vazia no fim.

A ocorrência "vazia" reproduz o que o pipeline grava quando
`fetch_local_record` devolve `{"error": ...}` — medido em
`Variola_VARV6_reexec_20260901`, onde só a primeira árvore processada tem
registro completo (o `raw_data_sequences.gb` foi sobrescrito por um
re-download com rótulo truncado, D13).

As strings trazem `;`, aspa simples e aspa dupla de propósito: são os
caracteres que já quebraram o tokenizer de blocos (C-5e, DEC-052).

`cql_legado_pre_G1.cql`, ao lado deste arquivo, é a saída do gerador
**antes** da correção de G1 sobre esta fixture (golden de caracterização).
"""

from __future__ import annotations

import copy

PROJETO = "Projeto_Fixture_G1"

_HASH = {"T1": 1111111111111111, "T2": 2222222222222222, "T3": 3333333333333333}
_NEWICK = {"T1": "MN000001", "T2": "MN000002", "T3": "MN000003"}

_REGISTRO = {
    "T1": {
        "description": "Zika virus isolate d'Ivoire; serogroup: \"Spondweni\", complete genome",
        "annotations": {
            "molecule_type": "RNA",
            "topology": "linear",
            "date": "01-JAN-2016",
            "source": "Zika virus",
            "organism": "Zika virus",
            "taxonomy": ["Viruses", "Riboviria", "Flaviviridae", "Orthoflavivirus"],
        },
        "features": [
            {
                "type": "source",
                "location": "[0:10807](+)",
                "strand": 1,
                "qualifiers": {
                    "organism": ["Zika virus"],
                    "country": ["Cote d'Ivoire: Abidjan"],
                    "note": ["passage 1; Vero cells"],
                },
            },
            {
                "type": "CDS",
                "location": "[106:10378](+)",
                "strand": 1,
                "qualifiers": {"product": ["polyprotein"]},
            },
        ],
    },
    "T2": {
        "description": "Zika virus strain ZIKV/Homo_sapiens/BRA//2016, partial genome",
        "annotations": {
            "molecule_type": "RNA",
            "topology": "linear",
            "date": "02-FEB-2017",
            "organism": "Zika virus",
            "taxonomy": ["Viruses", "Riboviria"],
        },
        "features": [
            {"type": "source", "location": "[0:10000](+)", "strand": 1,
             "qualifiers": {"country": ["Brazil"], "collection_date": ["2016"]}},
            {"type": "misc_feature", "location": "[5:50]", "qualifiers": {}},
        ],
    },
    "T3": {
        "description": "Zika virus isolate X, complete genome",
        "annotations": {"molecule_type": "RNA", "organism": "Zika virus"},
        "features": [],
    },
}


def _terminal(nome: str, completo: bool) -> dict:
    if completo:
        metadata = copy.deepcopy(_REGISTRO[nome])
    else:
        metadata = {"error": f"Acesso {_NEWICK[nome]}. não encontrado no arquivo x.gb"}
    return {"newick": _NEWICK[nome], "terminal_hash": _HASH[nome], "metadata": metadata}


def arvores_fixture() -> list:
    """Entrada no formato que `SubtreeMiner` passa a `parse_tree`.

    Returns
    -------
    list of dict
        ``[{nome_arvore: {nome_subarvore: {'supports', 'data_terminals'}}}]``
    """
    return [
        {"tree_mafft_iqtree": {
            "tree_mafft_iqtree_Inner1": {
                "supports": [0.5],
                "data_terminals": [_terminal("T1", False), _terminal("T3", True)],
            },
            "tree_mafft_iqtree_Inner2": {
                "data_terminals": [_terminal("T2", False), _terminal("T3", True)],
            },
        }},
        {"tree_mafft_fasttree": {
            "tree_mafft_fasttree_Inner1": {
                "supports": [0.5, 1.0],
                "data_terminals": [_terminal("T1", True), _terminal("T2", True),
                                   _terminal("T3", True)],
            },
            "tree_mafft_fasttree_Inner2": {
                "data_terminals": [_terminal("T1", True), _terminal("T2", True)],
            },
        }},
        {"tree_mafft_raxml": {
            "tree_mafft_raxml_Inner1": {
                "data_terminals": [_terminal("T1", False), _terminal("T3", False)],
            },
        }},
    ]


def ocorrencias_esperadas() -> list:
    """Pares (subárvore, terminal_hash) na ordem em que o gerador os visita."""
    pares = []
    for arvore in arvores_fixture():
        for _, subarvores in arvore.items():
            for nome_sub, dados in subarvores.items():
                for t in dados.get("data_terminals", []):
                    pares.append((nome_sub, t["terminal_hash"]))
    return pares
