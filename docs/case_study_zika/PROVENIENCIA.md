# Proveniência de `data/Zika479ONE/`

**Por que este arquivo não está dentro de `data/Zika479ONE/`:** a etapa de leitura de input do workflow trata **todo arquivo** do diretório de dataset como uma sequência a processar (`TreeBuilderController.__init__`, `sorted(os.listdir(self.input_path))`, sem filtro de extensão) — o mesmo defeito que já forçou mover `PROVENIENCIA.md` para fora do diretório do caso Variola (DEC-078) e que reapareceu, com outro gatilho (um `*NoPipe` de cache, não um `.md`), na execução de 2026-09-22 documentada em [DEC-102](../../../docs/automation/07-log-de-execucao.md#dec-102--2026-09-23--e10--primeira-execução-real-com-achado-novo-defeito-de-pipeline--monofilia-africana-em-33-métodos-de-caráter). Até esse filtro existir (candidato de M7), nenhum arquivo de documentação deve viver dentro de `data/Zika479ONE/` — só `dataset_final.fasta`.

## O que é `data/Zika479ONE/dataset_final.fasta`

478 sequências de *Zika virus*, formato FASTA. SHA-256 (conferido em três execuções distintas, sem alteração — 2026-09-02, 2026-09-05, 2026-09-22):

```
75c3665b443891f473aa1e4522e3daa681dad6fa6ca7943cec38f38bb7220db2
```

## De onde vem — confirmado, não mais hipótese

**Documentado em 2026-09-23** ([DEC-103](../../../docs/automation/07-log-de-execucao.md#dec-103--2026-09-23--e10--proveniência-confirmada-e-oráculo-independente-rodado-rf-real-contra-a-árvore-do-próprio-artigo-2-implementações-concordam)): o usuário forneceu o material suplementar do ScienceDirect do artigo [Zadra, Rizzoli & Rota-Stabelli (2024)](../../../docs/literatura/06-zika-filogenomica-2024.md), DOI [10.1016/j.virusres.2024.199490](https://doi.org/10.1016/j.virusres.2024.199490) — em particular `mmc4.tree` (NEXUS, 479 táxons, rótulo `ACCESSION_País_Data`), a árvore/dataset publicado pelos próprios autores como seu "dataset 1 (bruto filtrado)": GenBank, junho de 2021, filtro por sequência com data e local de coleta, ≥700 pb.

Comparação accession a accession (sem versão) entre `dataset_final.fasta` e as 479 sequências de `mmc4.tree`:

- **478 de 478** accessions de `Zika479ONE` batem exatamente com o dataset do artigo.
- A única sequência do artigo **ausente** em `Zika479ONE` é `MF099651` (China, Guizhou, `08/2016/15` no rótulo do artigo) — 1 de 479, não caracterizada a fundo (pode ser diferença de data de consulta ao NCBI, deduplicação, ou atualização/retratação do registro no GenBank entre a consulta do artigo e a aquisição deste dataset).

**Conclusão:** `Zika479ONE` **é** o dataset 1 do artigo, a menos de 1 registro — não uma extração independente coincidente. Script de comparação: `comparar_proveniencia.py` (não versionado — rodado ad hoc na sessão de DEC-103; reproduzível a partir de `mmc4.tree`, disponível em `docs/literatura/baseline_suplementar/mmc4_arvore_479_taxons.nex` no repositório principal, e de `raw_data_sequences.gb` de qualquer execução deste projeto).

## O que ainda não se sabe

- **Quem, exatamente, rodou a consulta ao NCBI que produziu este `dataset_final.fasta`, e quando.** O SHA-256 é idêntico desde a primeira execução registrada no ledger (2026-09-02) — mais antigo que isso não há registro. Nenhum `manifest.json`/log de aquisição sobrevive de antes dessa data.
- **Por que `MF099651` está ausente.** Não investigado — candidato a checagem futura contra `Table S2` (`docs/literatura/baseline_suplementar/mmc3_tabela_S2_completa.pdf`), que traz uma coluna "Usable" que pode explicar a exclusão.

## Uso

Usado em [E10](../../../docs/science/04-agenda-de-pesquisa.md#e10--✅--generalização-do-resultado-central-de-m3-numa-segunda-espécie--cotejo-com-zadra-et-al-2024-zikv) para cotejar o resultado atual da ferramenta contra o baseline publicado — ver DEC-102/DEC-103/DEC-104 no ledger para o experimento completo.
