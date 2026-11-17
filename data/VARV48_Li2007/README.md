# VARV-48 — o conjunto de Li *et al.* 2007

**Criado em 2026-11-17** (DEC-171, máquina `joh`). Este é **o** conjunto de *Variola* do BMC e do FGCS. Nenhum outro conjunto VARV entra nessas duas pesquisas.

## O que é

As **48** sequências correspondentes aos 48 terminais de *Variola virus* da SI Fig. 4 de:

> Li, Y.; Carroll, D. S.; Gardner, S. N.; Walsh, M. C.; Vitalis, E. A.; Damon, I. K. **"On the origin of smallpox: Correlating variola phylogenics with historical smallpox records."** *PNAS*, v. 104, n. 40, p. 15787-15792, 2007. DOI: 10.1073/pnas.0609268104.

41 terminais em P-I e 7 em P-II, exatamente como a figura.

## Composição e proveniência

| Fato | Valor |
|---|---|
| Registros no FASTA | **48** |
| Acessos distintos | 48 |
| `sha256` de `dataset_final.fasta` | `0741f6107e317bdeecf8c0dd7710e3195dd0d28fc4f7dd6a36626471494860b3` |
| Tamanho | 9 087 964 bytes |
| Comprimento das sequências | 184 170 a 188 251 pb (genomas completos, **com ITRs**) |
| Ordem dos registros | **alfabética por acesso**, fixada de propósito |
| Largura de linha | 60 |

**A lista de acessos não vem do artigo.** Li *et al.* nomeiam os terminais por código de isolado (`BSH74_nur`) e não publicam tabela de acessos. A ponte código → acesso está em [`docs/science/dados_curadoria/li2007_isolados.csv`](../../../docs/science/dados_curadoria/li2007_isolados.csv) e vem da *table S1* de Esposito *et al.* 2006 (*Science* 313:807-812, DOI 10.1126/science.1125134, SOM):

- **45 de 48** com acesso confirmado por essa fonte primária;
- **3** pela ponte secundária do campo `/strain` do GenBank — `SOM77_1252` (DQ441438.1), `SOM77_1605` (DQ441439.1) e `SUD47_rum` (DQ441441.1). **Isso se declara no texto**; não é proveniência primária.

Quem posiciona cada isolado em P-I/P-II é o Li *et al.* (por nome); quem dá o acesso é o Esposito *et al.* São **duas citações**, não uma.

## Como foi montado

47 dos 48 vieram de `workflow_dataAcquisition_li_et_al_2007_replication-RetMax200/dataset_final.fasta` (a colheita local), com cabeçalho e sequência intactos. O que faltava foi baixado do NCBI:

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=DQ441421.1&rettype=fasta&retmode=text"
# >DQ441421.1 Variola virus strain Bangladesh 1974 (Shahzaman), complete genome
# 186 293 pb — confere com o `slen` do esummary
```

## Três coisas que quem usar este conjunto precisa saber

1. **`DQ441421.1` e `DQ441422.1` têm a sequência byte a byte idêntica** (ambos 186 293 pb, mesmo `sha256` de sequência). São dois isolados distintos para o Li *et al.* — `BSH74_shz` (Shahzaman) e `BSH74_sol` (Solaiman), linhas 18 e 19 da table S1 — e dois registros distintos no GenBank, mas **o mesmo genoma**.
   - **É isso que explica a ausência de `DQ441421` na colheita local:** foi removido por deduplicação por sequência, não perdido pela consulta.
   - **Consequência operacional:** `deduplicar_por_sequencia` vai descartar um dos dois, e o pipeline rodará com **47 táxons**, não 48. Isso é esperado, tem de aparecer no manifesto e precisa ser declarado no texto — nenhum método baseado em alinhamento consegue separar dois terminais com sequência idêntica.
2. **`NC_001611.1` é RefSeq**, e **fica** (decisão do usuário, DEC-171). É o acesso de `IND67_mah`. O conjunto principal, portanto, mistura um RefSeq com 47 registros INSDC — o que aciona o [D23](../../../docs/science/02-defeitos-que-alteram-resultado.md). Declarar.
3. **As ITRs estão presentes.** Não existe conjunto VARV sem ITR neste repositório, apesar do que alguns nomes de diretório antigos sugerem (`01-revisao-variola.md §1.2`). As ITRs são a origem mais provável das colunas majoritariamente *gap* — ver o limiar do **P2** (DEC-171).

## Regra de uso

Este diretório contém **apenas** `dataset_final.fasta` e este README. O pipeline lista todo arquivo do diretório de entrada antes de limpar os `*_NoPipe` ([AC-09](../../../docs/audit/11-auditoria-ponta-a-ponta-2026-10-17.md)), e qualquer FASTA extra aqui vira um segundo dataset de entrada silencioso.
