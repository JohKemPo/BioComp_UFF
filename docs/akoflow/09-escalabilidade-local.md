# 09. Escalabilidade local: o parâmetro de paralelismo não limitou a concorrência

Status: **resultado negativo**. Os dados não mostram efeito mensurável do parâmetro P sobre o tempo, e a concorrência observada não variou com P. Este documento registra o experimento, o que deu errado no desenho e o que a equipe do AkôFlow precisa esclarecer. Nenhuma execução deste lote é publicável (IN-09: manifestos com `dirty: true`).

Fontes, todas em `/home/joh/akoflow_local/`:
- `esc/grade.log`: início, fim e código de retorno de cada execução;
- `analise/<run>.json`: métricas por execução, geradas por `akoflow/analisar_run.py`;
- `esc/condicoes.log`: estado da máquina no início de cada execução (carga, memória, CPU por contêiner);
- `esc/logs/<run>.out` e `.analise.out`: saída do runtime e do analisador;
- `env/env-p<P>.yaml`, `yaml/<kind>-p<P>.yaml`: configuração de cada ponto, gerada por `esc/gerar.py`.

Os rótulos de evidência seguem a regra do lote: **medido** (lido de arquivo, com o caminho), **observado** (lido de um registro de máquina) ou **inferido** (conclusão que os dados sugerem e que não foi testada diretamente).

## 1. Pergunta e desenho

**Pergunta.** Ao variar o paralelismo disponível (P) no ambiente local, o tempo de parede do workflow granular cai até o ponto em que o caminho crítico passa a limitar o ganho?

**Desenho pretendido.**
- Dataset Zika-21, com P = 1, 2, 4, 8 e 12, três repetições por ponto (15 execuções).
- Dataset VARV-6, com P = 4 e P = 12, uma repetição por ponto (2 execuções).
- Cada ponto P é um ambiente próprio, com `cpuCores = P` e `cpuCapacity = P` no recurso `ptm-local-docker-p<P>` (`env/env-p<P>.yaml`, linhas 31–32, medido).
- Cada atividade do YAML declara `cpuLimit: '1'`. O script `esc/gerar.py` reduz esse valor para P só quando ele é maior que P. Como todas as atividades declaram 1, o ajuste nunca foi aplicado (medido: 40 ocorrências de `cpuLimit: '1'` em `yaml/zika-p1.yaml`, nenhuma alterada).

**Execuções contadas.** Foram 17 execuções com `rc=0` em `esc/grade.log`: as 15 do Zika-21 e as 2 do VARV-6. Há ainda quatro arquivos de Zika-21 (`zika21-r3` a `zika21-r6`) e um de VARV-6 (`varv6-r3`) em `analise/`. Eles são de rodadas anteriores, não estão no `grade.log` desta noite e ficam fora deste experimento.

## 2. O parâmetro variado e por que ele não fez o que se esperava

O parâmetro variado foi `cpuCapacity` do recurso. A expectativa era que ele limitasse quantas atividades rodam ao mesmo tempo.

Não limitou. Três registros mostram isso:

1. **Paralelismo máximo observado igual em todos os pontos.** O campo `paralelismo_max` de cada `analise/*.json` vai de 8 a 10 no Zika-21, inclusive com P = 1 (`analise/zika-p1-r1.json`: 9; `zika-p1-r2.json`: 9; `zika-p1-r3.json`: 9). Para P = 12 o valor é 8 ou 9. Com `cpuCapacity = 1`, o analisador registra nove atividades simultâneas.
2. **Nenhuma fila.** O campo `fila_s` é zero em todas as 680 atividades das 17 execuções (17 × 40, medido). O runtime nunca esperou por recurso, o que é incompatível com um limite de concorrência de 1 ser aplicado.
3. **Contêiner usando mais de um núcleo com limite 1.** Em `esc/condicoes.log`, no bloco de P = 1, o contêiner `akoflow-zika-p1-r2-…-arvore-mafft-iterative-raxml` aparece com `cpu=401.67%` (linha 52) e `cpu=401.45%` (linha 63). Um contêiner com `cpuLimit` 1 não deveria passar de cerca de 100%. Os contêineres `…-upgma-parsimony` e `…-nj-parsimony` aparecem com cerca de 100% cada (medido).

Consequência: os cinco pontos de P rodaram, na prática, a mesma configuração (um ambiente com a mesma concorrência efetiva), e a variação entre eles é variação de repetição e de máquina, não efeito de P. Esse é o achado que invalida o desenho.

Não se sabe ainda se o `cpuCapacity` é ignorado pelo runtime local, se o `cpuLimit` não chega ao `docker run` ou se algum outro mecanismo substitui os dois. A seção 7 lista a pergunta para a equipe do AkôFlow.

## 3. As 17 execuções

Colunas: **makespan** é o tempo do workflow medido pelo analisador (`makespan_s`). **Parede** é o tempo do `grade.log`, que inclui o polling do cliente. **Caminho crítico** é `caminho_critico_s`. **CC/makespan** é a razão entre os dois. **Par. máx.** e **Par. médio** são `paralelismo_max` e `paralelismo_medio`.

| Execução | P | Rep. | Makespan (s) | Parede (s) | Caminho crítico (s) | CC/makespan | Par. máx. | Par. médio |
|---|---|---|---|---|---|---|---|---|
| zika-p1-r1 | 1 | 1 | 434,6 | 440 | 417,2 | 96,0 % | 9 | 3,22 |
| zika-p1-r2 | 1 | 2 | 481,3 | 487 | 458,9 | 95,3 % | 9 | 2,30 |
| zika-p1-r3 | 1 | 3 | 500,4 | 504 | 482,2 | 96,4 % | 9 | 3,19 |
| zika-p2-r1 | 2 | 1 | 407,0 | 413 | 393,4 | 96,7 % | 9 | 3,22 |
| zika-p2-r2 | 2 | 2 | 516,1 | 522 | 489,5 | 94,8 % | 8 | 2,24 |
| zika-p2-r3 | 2 | 3 | 512,4 | 519 | 490,8 | 95,8 % | 10 | 3,22 |
| zika-p4-r1 | 4 | 1 | 407,1 | 412 | 391,0 | 96,0 % | 9 | 3,15 |
| zika-p4-r2 | 4 | 2 | 1380,3 | 1385 | 1358,1 | 98,4 % | 9 | 2,35 |
| zika-p4-r3 | 4 | 3 | 519,2 | 524 | 502,7 | 96,8 % | 9 | 3,24 |
| zika-p8-r1 | 8 | 1 | 405,7 | 410 | 392,7 | 96,8 % | 8 | 3,18 |
| zika-p8-r2 | 8 | 2 | 1224,4 | 1232 | 1203,3 | 98,3 % | 10 | 2,38 |
| zika-p8-r3 | 8 | 3 | 513,3 | 519 | 492,8 | 96,0 % | 9 | 3,26 |
| zika-p12-r1 | 12 | 1 | 458,3 | 465 | 440,9 | 96,2 % | 8 | 3,23 |
| zika-p12-r2 | 12 | 2 | 413,9 | 419 | 396,6 | 95,8 % | 8 | 3,14 |
| zika-p12-r3 | 12 | 3 | 505,2 | 2870 | 485,9 | 96,2 % | 9 | 3,25 |
| varv-p4-r1 | 4 | 1 | 3902,5 | 3908 | 3886,1 | 99,6 % | 7 | 2,03 |
| varv-p12-r1 | 12 | 1 | 4427,1 | 4810 | 4409,6 | 99,6 % | 7 | 2,02 |

Notas sobre a tabela:
- Todas as 17 execuções têm `status = completed`, `falha = ""` e 40 atividades, todas com estado `completed` (medido).
- Na linha **zika-p12-r3**, o parede (2870 s) excede o makespan (505 s) em cerca de 2365 s. Na **varv-p12-r1**, a diferença é de cerca de 383 s. Nas demais, a diferença é de 5 a 15 s. Os registros disponíveis não explicam essas diferenças (não verificado; ver a seção 7).
- A execução **zika-p12-r1** tem nesta tabela o registro da repetição que completou (463,2 s na saída do runtime, `esc/logs/zika-p12-r1.out`). Ver a seção 5.

## 4. Dispersão por ponto (Zika-21, makespan, três repetições)

| P | Mínimo (s) | Mediana (s) | Máximo (s) | Amplitude (s) |
|---|---|---|---|---|
| 1 | 434,6 | 481,3 | 500,4 | 65,7 |
| 2 | 407,0 | 512,4 | 516,1 | 109,1 |
| 4 | 407,1 | 519,2 | 1380,3 | 973,3 |
| 8 | 405,7 | 513,3 | 1224,4 | 818,7 |
| 12 | 413,9 | 458,3 | 505,2 | 91,3 |

As medianas dos cinco pontos vão de 458,3 s (P = 12) a 519,2 s (P = 4), uma diferença de 61 s. Dentro de cada ponto, a amplitude vai de 66 s a 973 s. Em P = 4 e P = 8 há uma execução que levou cerca de 1,2 a 1,4 mil segundos. O caminho crítico dessas execuções também é longo (1358 s e 1203 s), ou seja, as atividades ficaram mais lentas, e não só o fim do workflow foi atrasado. A dispersão dentro dos pontos é maior que a diferença entre eles, então os dados não separam os pontos.

Para o VARV-6, com uma repetição por ponto: P = 4 levou 3902,5 s e P = 12 levou 4427,1 s (medido). Nessa comparação, P = 12 foi mais lento que P = 4, o oposto do que se esperaria de mais paralelismo, mas há uma repetição só por ponto, sem medida de dispersão, e a execução de P = 12 coincidiu com outras execuções na mesma máquina (seção 6). Não é possível interpretá-la como efeito de P.

## 5. Resultado

1. **Não há efeito mensurável do parâmetro P sobre o tempo de parede** nos dados coletados. A variação dentro de um mesmo ponto é maior que a diferença entre os pontos.
2. **A concorrência observada não variou com P** (paralelismo máximo de 8 a 10 no Zika-21 em todos os pontos, inclusive P = 1), e não houve fila em nenhuma atividade. O parâmetro variado não limitou a concorrência (seção 2).
3. **O caminho crítico ocupa 94,8 % a 98,4 % do makespan no Zika-21 e 99,6 % no VARV-6**, em todas as 17 execuções (medido). O workflow, nas execuções registradas, termina quando termina a sequência de atividades longas, e não quando falta concorrência. Isso é consistente com o caminho crítico como limite, mas o experimento não o isola: o mesmo padrão aparece nas execuções lentas, nas quais a interferência externa é provável (inferido).

**Execução truncada, que não é resultado.** A tentativa inicial de `zika-p12-r1` foi identificada como truncada na análise do lote: um registro com tempo de cerca de 162 s, *speedup* de 4,13 e caminho crítico de quatro atividades terminando em `subarvores`, sem chegar a `empacotar-saida`. Essa execução não completou o workflow e o número não entra como resultado. **Não foi possível reconferir esse registro a partir dos arquivos atuais**: o `analise/zika-p12-r1.json` que existe hoje é o da repetição que completou (458,3 s, caminho de dez atividades até `empacotar-saida`), e a linha do `grade.log` de 00:03:55 mostra apenas `rc=1` com 0 s para a primeira tentativa, sem registro de métricas. O registro truncado foi sobrescrito. A conclusão de que a tentativa estava truncada depende da descrição do lote, e deve ser conferida no log de origem antes de qualquer citação.

Com isso, o "speedup" de 4,13 não é um resultado. Os *speedups* registrados em `analise/*.json` (campo `speedup_vs_mono`) dependem de um tempo de monólito de referência que não faz parte deste experimento, e não são usados aqui.

## 6. Ameaças à validade

- **Máquina compartilhada.** Durante o experimento rodaram o Neo4j (`phylotree_neo4j`, com amostras de até 220 % de CPU em `esc/condicoes.log`, linha 38, medido), um Postgres (`mbmscraping-postgres-mbm-1`) e o servidor do AkôFlow. A carga média de 1 minuto no início das execuções ficou entre cerca de 10 e 16 nas janelas de 00:06 a 00:35 e de 02:04 (medido em `esc/condicoes.log`), em uma máquina que o manuscrito declara com 12 núcleos lógicos (não reconferido aqui).
- **Execuções sobrepostas.** Dois loops de execução foram iniciados na mesma noite, `esc/grade.sh` (INÍCIO às 23:30:20) e `esc/grade2.sh` (INÍCIO às 00:03:55), gravando no mesmo `grade.log` e nos mesmos identificadores de execução. Há entradas `rc=1` com 0 a 1 s de duração, seguidas de execuções com `rc=0` do mesmo identificador (medido em `grade.log`). Os blocos de `condicoes.log` mostram contêineres de uma execução (por exemplo, `zika-p4-r2`) listados no início de outra (`zika-p8-r2`, 00:15:02), e de `zika-p12-r2` durante o início de `varv-p12-r1` (00:35:43) (medido). A conclusão é que execuções rodaram ao mesmo tempo na mesma máquina, e que o registro de `analise/` pode ter sobrescrito uma execução por outra com o mesmo identificador (inferido, não testado diretamente). Essa é a principal ameaça: ela afeta os pontos com execuções lentas (P = 4 e P = 8, com 1385 s e 1232 s) e a comparação do VARV-6.
- **Três repetições por ponto**, com dispersão de até 973 s dentro de um ponto. Três repetições não bastam para estimar uma mediana com intervalo útil.
- **Uma máquina só**, sem medida em outro hardware.
- **Conjunto pequeno**: Zika-21 (vinte sequências, conforme o manuscrito) e VARV-6 (seis genomas de cerca de 250 kb). Uma atividade de alinhamento maior poderia mudar o caminho crítico.
- **Não há controle do ambiente**: nenhuma execução foi feita com a máquina ociosa, e a versão do pacote do AkôFlow usada é a do worktree `AKOFLOW_Phylotree_fgcs` no momento da execução (não reconferida neste lote).

## 7. Perguntas para a equipe do AkôFlow

1. **Qual é o parâmetro correto para limitar a concorrência local?** Entre `cpuCapacity` do recurso, `cpuLimit` por atividade e um eventual limite de execuções simultâneas do runtime, qual deve ser usado para fixar P, e como cada um é convertido em limite de contêiner no runtime local (`docker run`)?
2. **Por que `cpuCapacity = 1` não limitou?** Todas as atividades declaram `cpuLimit: '1'`, e um contêiner chegou a 401 % de CPU (`esc/condicoes.log`, linhas 52 e 63). O `cpuLimit` é aplicado como `--cpus` no `docker run`? Há um caminho em que ele é ignorado no runtime local?
3. **Como obter o tempo por atividade separado do tempo de fila?** Todas as 680 atividades registram `fila_s = 0`. O runtime local tem fila, ou a métrica não é medida nesse runtime? Há como registrar, por atividade, o instante de início do contêiner, o de fim e o de liberação do recurso, para separar a espera do tempo de execução?
4. **De onde vêm as diferenças entre o tempo de parede e o makespan** em zika-p12-r3 (cerca de 2365 s) e varv-p12-r1 (cerca de 383 s)? O polling do cliente explica, ou algum estágio antes da primeira atividade (por exemplo, a criação do workspace) entra na conta?

## 8. Próximos passos deste lote (fora do escopo desta tarefa)

- Refazer o experimento com um único executor por vez, com identificadores que não se repetem entre loops, e com um bloqueio de execução na máquina.
- Antes de refazer, confirmar a resposta à pergunta 2 da seção 7, porque sem ela o parâmetro continua sem efeito.
- Repetir pelo menos cinco vezes por ponto e reportar mediana e intervalo, como o manuscrito já prevê.
- Não citar nenhum número deste documento como resultado de escalabilidade enquanto as execuções não forem refeitas.
