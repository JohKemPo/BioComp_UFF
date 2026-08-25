# Proveniência — dataset_final.fasta

Conjunto derivado por remoção taxonômica. **O conjunto de origem não foi alterado.**

| Campo | Valor |
|---|---|
| Origem | `data/SMALL_li_2007_replication-RetMax200-ITRs/dataset_final.fasta` |
| Clado exigido | Orthopoxvirus (`txid10242`) |
| Sequências antes | 6 |
| Sequências depois | 5 |
| Removidas | 1 |
| Sem linhagem, mantidas | 0 |
| Gerado em (UTC) | 2026-08-25T14:12:55+00:00 |

## Removidas

| Acesso | Motivo |
|---|---|
| `NC_008030` | Nile crocodilepox virus — Crocodylidpoxvirus > Crocodylidpoxvirus nilecrocodilepox, fora de Orthopoxvirus |

## Reproduzir

```bash
cd BioComp_UFF && python ../docs/science/scripts/limpar_datasets.py
```

A conferência independente é `docs/science/scripts/auditar_taxonomia.py`.