# ZIKV-479 — o dataset 1 de Zadra *et al.* 2024

**Criado em 2026-11-17** (DEC-171, máquina `joh`). Este é **o** conjunto de Zika do BMC e do FGCS. Nenhum outro conjunto Zika entra nessas duas pesquisas.

## O que é

As **479** sequências do *dataset 1* (bruto filtrado) de:

> Zadra, N.; Rizzoli, A.; Rota-Stabelli, O. **"Comprehensive phylogenomic analysis of Zika virus: Insights into its origin, past evolutionary dynamics, and global spread."** *Virus Research*, v. 350, artigo 199490, 2024. DOI: 10.1016/j.virusres.2024.199490. Acesso aberto, CC BY 4.0.

Filtro dos autores (§3.1 do artigo): GenBank, arquivo válido, data de coleta e local presentes, sequência > 700 pb, duplicatas removidas — 1733 *hits* brutos de junho de 2021 reduzidos a 479.

## Composição e proveniência

| Fato | Valor |
|---|---|
| Registros no FASTA | **479** |
| Acessos distintos | 479 |
| `sha256` de `dataset_final.fasta` | `5e541e079b676e66f476b2babed4c344793fe04a3afaead7817f49d5268c9772` |
| Tamanho | 4 791 641 bytes |
| Comprimento das sequências | 708 a 10 812 pb |
| Ordem dos registros | **alfabética por acesso** |
| Largura de linha | 60 |
| Sequências idênticas entre si | 0 |

**A identidade com o artigo é verificada, não presumida.** A árvore publicada pelos próprios autores ([`mmc4_arvore_479_taxons.nex`](../../../docs/literatura/baseline_suplementar/mmc4_arvore_479_taxons.nex), rótulo `ACCESSION_País_Data`) tem 479 táxons, e a comparação acesso a acesso contra este conjunto dá **479 de 479, sem sobra nem falta**.

## Como foi montado

478 dos 479 vieram de `Zika479ONE/dataset_final.fasta`, com cabeçalho e sequência intactos. O que faltava foi baixado do NCBI:

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MF099651.1&rettype=fasta&retmode=text"
# >MF099651.1 Zika virus isolate Zika-GZDJ1685 polyprotein gene, partial cds
# 10 260 pb — confere com o `slen` do esummary; `Status live`, sem `ReplacedBy`
```

`MF099651` é `MF099651_China_Guizhou_08/2016/15` na árvore dos autores. Com 10 260 pb, data (08/2016) e local (China, Guizhou), satisfaz os três critérios do filtro do artigo — a ausência era da nossa aquisição, não do GenBank nem de uma exclusão dos autores.

## O conjunto anterior continua existindo

`Zika479ONE/` (478 sequências) **não foi alterado**, de propósito: é a entrada declarada no `manifest.json` de `Zika_ZIKV480_reexec_20260901` e a base da evidência do [DEC-103](../../../docs/automation/07-log-de-execucao.md). Mexer nele invalidaria o que já foi medido.

## Duas coisas que quem usar este conjunto precisa saber

1. **Cobertura de metadado, medida nos 478 registros antigos:** 478 de 478 com `/collection_date` (100%), 478 com `/geo_loc_name`, **0 com `/country`** — o NCBI renomeou o campo, e extrator que procure `/country` não acha nada. Anos de 1947 a 2018, com 341 registros de 2016: amplitude boa para sinal temporal, amostragem desequilibrada (limitação a declarar, não defeito).
2. **O conjunto é heterogêneo em comprimento** (708 a 10 812 pb): há fragmentos alinhados contra genomas quase completos. Isso é o que torna o **U4** (sobreposição mínima da p-distância) e o **P2** (limiar de *gap*) decisivos aqui — ver DEC-171 e [DM-10](../../../docs/audit/11-auditoria-ponta-a-ponta-2026-10-17.md).

## Regra de uso

Este diretório contém **apenas** `dataset_final.fasta` e este README. Qualquer FASTA extra aqui vira um segundo dataset de entrada silencioso ([AC-09](../../../docs/audit/11-auditoria-ponta-a-ponta-2026-10-17.md)).
