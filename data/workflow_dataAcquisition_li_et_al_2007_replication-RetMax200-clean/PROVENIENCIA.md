# Proveniência — dataset_final.fasta

Conjunto derivado por remoção taxonômica. **O conjunto de origem não foi alterado.**

| Campo | Valor |
|---|---|
| Origem | `data/workflow_dataAcquisition_li_et_al_2007_replication-RetMax200/dataset_final.fasta` |
| Clado exigido | Orthopoxvirus (`txid10242`) |
| Sequências antes | 125 |
| Sequências depois | 121 |
| Removidas | 4 |
| Sem linhagem, mantidas | 3 |
| Gerado em (UTC) | 2026-08-25T14:13:39+00:00 |

## Removidas

| Acesso | Motivo |
|---|---|
| `MG450915` | Saltwater crocodilepox virus — Chordopoxvirinae > Crocodylidpoxvirus, fora de Orthopoxvirus |
| `MG450916` | Saltwater crocodilepox virus — Chordopoxvirinae > Crocodylidpoxvirus, fora de Orthopoxvirus |
| `NC_008030` | Nile crocodilepox virus — Crocodylidpoxvirus > Crocodylidpoxvirus nilecrocodilepox, fora de Orthopoxvirus |
| `NC_015960` | Yokapox virus — Centapoxvirus > Centapoxvirus yokapox, fora de Orthopoxvirus |

## Mantidas sem decisão

Sem registro GenBank que permitisse decidir o clado. **Permanecem no conjunto**: retirá-las seria descartar dado por falta de metadado, que é uma decisão diferente e precisa ser tomada explicitamente.

- `DQ437594`
- `HQ849551`
- `NC_003391`

## Reproduzir

```bash
cd BioComp_UFF && python ../docs/science/scripts/limpar_datasets.py
```

A conferência independente é `docs/science/scripts/auditar_taxonomia.py`.