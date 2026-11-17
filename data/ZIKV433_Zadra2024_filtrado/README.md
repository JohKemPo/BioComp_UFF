# ZIKV-433 — o braço filtrado do dataset 1 de Zadra *et al.* 2024

**Criado em 2026-11-17** (DEC-171/DEC-175, protocolo [§3.3](../../../docs/science/19-protocolo-reexecucao-baselines.md)). É o **segundo braço** da decisão U5: a análise de sensibilidade do ZIKV-479.

## O que é, e de quem é o critério

As **433** sequências do `ZIKV479_Zadra2024` que têm **9 000 nt ou mais**. As 46 removidas são fragmentos de 753 a 841 nt.

**O critério não é nosso.** Remover sequências abaixo de 9 000 nt é o primeiro passo do *dataset 2* de Zadra *et al.* 2024 (§3.1 do artigo). Não inventamos filtro: usamos o dos próprios autores.

> O que **não** fizemos é o segundo passo deles, a poda por efeito de densidade de nós, que levaria a 117 sequências. Esse passo depende de decisões de amostragem que não temos como reproduzir fielmente, e o conjunto de 117 **não** é um braço desta pesquisa.

| Fato | Valor |
|---|---|
| Registros | **433** |
| `sha256` | `077de2505794e79bdba30040eff0e5196338ceb483818f2c5bd0fdaf9934380c` |
| Ordem dos registros | alfabética por acesso |
| Largura de linha | 60 |

## Como reproduzir, em qualquer máquina, sem transferir o arquivo

A derivação é **determinística** — rodada duas vezes aqui, deu o mesmo `sha256` nas duas:

```bash
python docs/science/scripts/derivar_conjunto_filtrado.py \
    BioComp_UFF/data/ZIKV479_Zadra2024/dataset_final.fasta \
    BioComp_UFF/data/ZIKV433_Zadra2024_filtrado/dataset_final.fasta
```

Se o `sha256` da saída não for o da tabela acima, **pare**: ou a entrada não é o `ZIKV479_Zadra2024` canônico, ou algo mudou no script.

## Por que este braço existe

O ZIKV-479 completo vai de 708 a 10 812 pb. Essa heterogeneidade é o que faz o **U4** recusar p-distância em 59 e mais pares e o **filtro do P2** querer apagar quase todo o alinhamento. O braço de 433 responde a pergunta "o resultado depende dos fragmentos?" sem **jamais** substituir o conjunto do baseline, que continua sendo o de 479.

## Regra de uso

Só `dataset_final.fasta` e este README. Qualquer FASTA extra aqui vira um segundo dataset de entrada silencioso ([AC-09](../../../docs/audit/11-auditoria-ponta-a-ponta-2026-10-17.md)).
