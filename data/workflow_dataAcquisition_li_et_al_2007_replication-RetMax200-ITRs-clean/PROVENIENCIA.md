# Proveniência — dataset_final.fasta

Conjunto derivado por remoção taxonômica. **O conjunto de origem não foi alterado.**

| Campo | Valor |
|---|---|
| Origem | `data/workflow_dataAcquisition_li_et_al_2007_replication-RetMax200-ITRs/dataset_final.fasta` |
| Clado exigido | Orthopoxvirus (`txid10242`) |
| Sequências antes | 55 |
| Sequências depois | 54 |
| Removidas | 1 |
| Sem linhagem, mantidas | 2 |
| Gerado em (UTC) | 2026-08-25T14:13:19+00:00 |

## Removidas

| Acesso | Motivo |
|---|---|
| `NC_008030` | Nile crocodilepox virus — Crocodylidpoxvirus > Crocodylidpoxvirus nilecrocodilepox, fora de Orthopoxvirus |

## Mantidas sem decisão

Sem registro GenBank que permitisse decidir o clado. **Permanecem no conjunto**: retirá-las seria descartar dado por falta de metadado, que é uma decisão diferente e precisa ser tomada explicitamente.

- `DQ437594`
- `NC_003391`

## Reproduzir

```bash
cd BioComp_UFF && python ../docs/science/scripts/limpar_datasets.py
```

A conferência independente é `docs/science/scripts/auditar_taxonomia.py`.