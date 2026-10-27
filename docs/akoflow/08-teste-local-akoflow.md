# 08 - Teste local do AkôFlow com o workflow granular do PhyloTreeMiner

Máquina: `joh` (Ubuntu, kernel 6.8.0-138, 12 núcleos lógicos, 31 GB de RAM, Docker local). Data: 2026-10-27. AkôFlow `v1.0.11`. Workflow: `akoflow/workflow.yaml` da G4b (39 atividades de desenho, 40 com `empacotar-saida`), ajustado por `akoflow/gerar_yaml.py`. Etiquetas: **medido** (comando e saída vistos), **observado**, **inferido**, **não testado**.

> **Ressalva geral.** Máquina local única, conjuntos pequenos (Zika-21: 20 registros no `.gb`; VARV-6: 6 genomas de ~250 kb), um só run completo por conjunto para a conclusão de tempo. Nenhum número daqui é resultado científico publicável: os artefatos não são de execução limpa (`dirty: true`, IN-09) e o VARV-6 não é conjunto pós-D28 catalogado. Os resultados reais virão de nuvem ou HPC.

## 1. Instalação (medido, ordem executada)

Tudo em `/home/joh/akoflow_local/` (fora do repositório). Release oficial, opção B de 07 §3.2 (servidor Docker sem Desktop).

```bash
# 1. artefatos da release v1.0.11 (daemon e buildkit), verificados por SHA-256
cd /home/joh/akoflow_local/release && sha256sum -c akoflow-runtime-v1.0.11-linux-amd64.sha256
docker load -i akoflow-daemon-v1.0.11-linux-amd64.tar
docker load -i akoflow-buildkit-v1.0.11-linux-amd64.tar
# 2. servidor: compose oficial + um ajuste local (compose.override.yaml)
cd /home/joh/akoflow_local/server        # .env: AKOFLOW_RELEASE_TAG=v1.0.11, AKOFLOW_API_TOKEN=<aleatório>, AKOFLOW_PORT=8080
docker compose -f compose.yaml -f compose.override.yaml up -d      # porta 127.0.0.1:8080 (loopback)
# 3. imagem do workflow (3,64 GB; micromamba + as 7 ferramentas pinadas; build.log = 104 KB)
cd <worktree> && docker build -f akoflow/Dockerfile -t phylotreeminer-akoflow:local .
# 4. ambiente local (runtime driver=local, 12 cores, 24 GiB), conjunto, execução
export AKOFLOW_API_URL=http://127.0.0.1:8080/akoflow-api AKOFLOW_API_TOKEN=...
python akoflow/gerar_yaml.py --nome ptm-zika21 --serializar-saidas --empacotar --saida zika21.yaml
python akoflow/rodar_akoflow.py --yaml zika21.yaml --ambiente akoflow/environment.local.yaml --id zika21-r6 --coleta coleta/zika21-r6
python akoflow/analisar_run.py zika21-r6 --mono 672 --json analise/zika21-r6.json
```

`compose.override.yaml` (único desvio do compose oficial): o daemon inicia cada atividade como `docker run` no Docker do host via `docker.sock`, então o `AKOFLOW_LOCAL_WORKSPACE_ROOT` precisa ser o **mesmo caminho no host e no daemon** (`/home/joh/akoflow_local/workspace`, montado nos dois).

Fluxo da API (necessário; o corpo só com `run` falha com `FOREIGN KEY constraint failed`, medido): ambiente, escopo, topologia, definição do workflow, sessão de planejamento HEFT, candidato selecionado, `POST /execution-runs/` com o *bundle* completo. `akoflow/rodar_akoflow.py` faz tudo e, num *rerun* com o mesmo nome, reaproveita a definição existente (409/422/500 tolerados nas criações idempotentes).

## 2. Ajustes ao YAML e por quê

As divergências de 07 §4 foram confirmadas ou refinadas na execução. O que `gerar_yaml.py` faz, e o que o teste mostrou:

| Ajuste | Motivo (07 §4 / medido) |
|---|---|
| `command.executable` (OCI, `phylotreeminer-akoflow:local`, `delivery: auto`) no lugar de `spec.image` | `spec.image` é forma legada (07 §4.1) |
| `runtime: ptm-local-runtime` e `environment.local.yaml` (runtime, recurso, binding) | o `runtime:` da atividade não escolhe nada; sem ambiente, escopo e plano não há execução (07 §4.3 itens 2-3). |
| `expectedOutputs` em toda atividade: fragmento `out/outputs/prov/<atividade>.json` + saídas de contrato | afirma existência e checksum das saídas. **Medido:** o AkôFlow **cancela** a atividade se um `expectedOutput` não for observado (r2: `consolidar-manifesto` cancelada porque não grava fragmento próprio; por isso essa atividade ficou sem o fragmento) |
| `dataDependencies` conferidas contra `dependsOn`; `sizeBytes` do alinhamento reescalado | a doc exige `dependsOn` correspondente (07 §4.1); conferido por assert no gerador |
| **`rsync` e conflito de caminho (07 §4.3 item 1):** todas as atividades gravam em `out/` e cada uma só escreve arquivos de **nome próprio** (fragmento de proveniência por atividade, árvores e subárvores por nome de método). Nenhuma atividade paralela reescreve o mesmo caminho com conteúdo distinto | o sucessor recebe o workspace inteiro de cada predecessor e o `rsync` rejeita conteúdo diferente no mesmo caminho. Com a convenção de nomes, **os 86 `rsync` de cada run passaram sem conflito** (medido, `dataTransfers`: 86, todos concluídos). Os conflitos que apareceram foram outros, abaixo |
| `--serializar-saidas`: `graficos-suporte` passa a depender também de `gerar-cql` | **Medido (r1 do Zika e do VARV-6, falha):** duas atividades consumindo `minerar-fpmax` ao mesmo tempo disputam a liberação do *lease* do workspace; o daemon falha com `release workspace leases ... lstat .../out/Subtrees/...: no such file or directory`, e a seguinte com `context canceled`. Serializar elimina a corrida (inferido: é corrida do daemon, não do nosso código) |
| `--empacotar`: atividade final `empacotar-saida` (`tar -cf resultados/out.tar out`) | **Medido:** o daemon poda do workspace final as entradas herdadas; só sobrevive o que a **última atividade grava**. Sem ela, não há como recolher `out/` completo |
| `--sem-gb` (VARV-6): remove `--gb-congelado` | não há `.gb` congelado na imagem para o VARV-6; `baixar-metadados` roda em modo online (NCBI) |
| `USER root` no `akoflow/Dockerfile` | **Medido:** o daemon cria o workspace de cada atividade como root no host e o monta no contêiner; com o usuário padrão do micromamba (uid 57439) a atividade falha com `PermissionError: '.'` |
| `memoryLimit`/`cpuLimit` explícitos (2 Gi nas árvores; 1 Gi no modelo; `cpuLimit 4` nos alinhamentos) | padrões da doc são 0,1 núcleo e 16 Mi (07 §4.3 item 5). **Medido (VARV-6 r2):** `alinhar-mafft` e `alinhar-mafft-iterative` saíram com código 137 (OOM), em 122 s; o `varv6.yaml` usado no r3 traz `memoryLimit: 2Gi` e `cpuLimit: 4` nos alinhamentos, e o r3 completou (o valor exato do r2 não foi preservado nos arquivos) |

## 3. Execuções realizadas

Todos os runs na mesma máquina. O run de referência para cada conjunto é o **último completo** (Zika-21 r6, VARV-6 r3); os anteriores são o caminho de depuração.

| Run | Resultado | Motivo / observação (medido) |
|---|---|---|
| zika21-r1 | falhou (36 ok, 2 falhas) | corrida de `rsync` dos dois consumidores de `minerar-fpmax` (§2) |
| zika21-r2 | 38 ok, 1 cancelada | `expectedOutput` `prov/consolidar-manifesto.json` não observado |
| zika21-r3 | sem coleta | interrompido (log vazio); sem dado |
| zika21-r4 | completo, 39 atividades, 455,1 s | ainda sem `empacotar-saida` |
| zika21-r5 | completo, 40, 346,8 s | com `empacotar-saida` (equivalência em `equiv/zika_r5_*`) |
| **zika21-r6** | **completo, 40, 409,5 s** | run de referência do Zika-21 |
| varv6-r1 | falhou (36 ok, 2 falhas), 3 197 s | a mesma corrida de `rsync` |
| varv6-r2 | falhou (2 ok, 2 falhas), 122 s | OOM (código 137) nos dois alinhamentos |
| **varv6-r3** | **completo, 40, 3 905,2 s** (3 908 s de parede) | run de referência do VARV-6 |

A variação entre r4, r5 e r6 do Zika-21 (455, 347, 410 s) vem sobretudo da concorrência entre os RAxML/NJ-parcimônia/UPGMA-parcimônia, que disputam 12 núcleos (inferido; um único run de cada).

## 4. Tabela atividade x tempo x estado (medido, `POST /provenance/sql/` do próprio daemon)

Duração = `runtime_seconds` (tempo dentro do contêiner; exclui preparação do workspace e *rsync*). Ordenado pelo início no Zika-21 r6.

| Atividade | Zika-21 duração (s) | Zika estado | VARV-6 duração (s) | VARV-6 estado |
|---|---:|---|---:|---|
| `preparar-entrada` | 1.43 | completed | 2.07 | completed |
| `alinhar-mafft` | 8.03 | completed | 3548.53 | completed |
| `baixar-metadados` | 0.35 | completed | 6.86 | completed |
| `alinhar-mafft-iterative` | 13.73 | completed | 3781.16 | completed |
| `modelo-mafft` | 4.64 | completed | 7.33 | completed |
| `arvore-mafft-upgma-parsimony` | 3.42 | completed | 74.73 | completed |
| `arvore-mafft-upgma-distance` | 3.30 | completed | 3.99 | completed |
| `arvore-mafft-nj-distance` | 3.03 | completed | 4.12 | completed |
| `arvore-mafft-nj-parsimony` | 3.11 | completed | 126.79 | completed |
| `arvore-mafft-fasttree` | 7.02 | completed | 39.55 | completed |
| `subarvores-mafft-nj-parsimony` | 0.53 | completed | 2.06 | completed |
| `subarvores-mafft-upgma-parsimony` | 0.53 | completed | 2.56 | completed |
| `subarvores-mafft-upgma-distance` | 0.56 | completed | 3.34 | completed |
| `subarvores-mafft-nj-distance` | 0.46 | completed | 2.19 | completed |
| `arvore-mafft-iterative-upgma-distance` | 3.13 | completed | 3.11 | completed |
| `arvore-mafft-iterative-upgma-parsimony` | 136.48 | completed | 54.05 | completed |
| `modelo-mafft-iterative` | 4.20 | completed | 5.53 | completed |
| `arvore-mafft-iterative-nj-parsimony` | 374.99 | completed | 91.40 | completed |
| `arvore-mafft-iterative-nj-distance` | 3.24 | completed | 2.86 | completed |
| `arvore-mafft-iterative-fasttree` | 8.44 | completed | 31.07 | completed |
| `arvore-mafft-raxml` | 334.37 | completed | 41.79 | completed |
| `arvore-mafft-iqtree` | 13.38 | completed | 5.85 | completed |
| `subarvores-mafft-fasttree` | 1.07 | completed | 1.94 | completed |
| `subarvores-mafft-iterative-upgma-distance` | 1.09 | completed | 2.43 | completed |
| `subarvores-mafft-iterative-nj-distance` | 1.40 | completed | 1.42 | completed |
| `arvore-mafft-iterative-raxml` | 329.00 | completed | 30.07 | completed |
| `arvore-mafft-iterative-iqtree` | 12.24 | completed | 5.15 | completed |
| `subarvores-mafft-iterative-fasttree` | 1.34 | completed | 1.34 | completed |
| `subarvores-mafft-iqtree` | 1.24 | completed | 2.87 | completed |
| `subarvores-mafft-iterative-iqtree` | 1.29 | completed | 3.21 | completed |
| `subarvores-mafft-iterative-upgma-parsimony` | 1.00 | completed | 1.91 | completed |
| `subarvores-mafft-iterative-raxml` | 0.74 | completed | 1.48 | completed |
| `subarvores-mafft-raxml` | 0.62 | completed | 2.10 | completed |
| `subarvores-mafft-iterative-nj-parsimony` | 0.63 | completed | 1.74 | completed |
| `heatmap-rf` | 2.99 | completed | 3.24 | completed |
| `minerar-fpmax` | 1.63 | completed | 3.92 | completed |
| `gerar-cql` | 0.72 | completed | 2.06 | completed |
| `graficos-suporte` | 1.51 | completed | 2.38 | completed |
| `consolidar-manifesto` | 1.48 | completed | 1.46 | completed |
| `empacotar-saida` | 0.03 | completed | 0.08 | completed |

Resumo por run (de `analisar_run.py`, `analise/*.json`):

| Run | Atividades | Makespan | Trabalho total | Caminho crítico medido | Paralelismo máx. / médio | Fila |
|---|---:|---:|---:|---:|---:|---:|
| zika21-r6 | 40 completas | 409,5 s | 1 288,4 s | 396,2 s | 8 / 3,15 | 0 s |
| varv6-r3 | 40 completas | 3 905,2 s | 7 909,8 s | 3 886,3 s | 7 / 2,03 | 0 s |

Caminho crítico do Zika-21 r6: `preparar-entrada` > `alinhar-mafft-iterative` > `arvore-mafft-iterative-nj-parsimony` (375 s) > `subarvores-...-nj-parsimony` > `minerar-fpmax` > `gerar-cql` > `graficos-suporte` > `consolidar-manifesto` > `empacotar-saida`. No VARV-6 r3 o tempo é dominado por `alinhar-mafft-iterative` (3 781 s) e `alinhar-mafft` (3 549 s), **em paralelo**.

**Por que o alinhamento do VARV-6 leva ~1 h (e não era travamento).** Observado: ao inspecionar os contêineres, os dois alinhamentos ainda avançavam; terminaram em 3 549 s e 3 781 s. O monolito gastou 8 303 s só na fase "Building trees" (2 h 18 min). Genomas de ~250 kb (alinhamento de 250 452 colunas) com MAFFT, e `--maxiterate 2` acima de 50 kb (`alignmentSeq.py:408-410`), tornam o alinhamento o gargalo (inferido: a causa fina, memória ou iterações, não foi perfilada).

## 5. Speedup medido contra o monolito e contra o teórico

Monolito de referência: o mesmo `workflow.py` dentro da mesma imagem, no Docker local (`ref/rodar_ref.sh`), código de saída 0: **Zika-21 671,7 s; VARV-6 8 314,7 s** (medido, `ref/*/inicio` e `fim`).

| | Monolito | Makespan AkôFlow | Speedup medido | Caminho crítico | Speedup teórico "monolito / caminho crítico" desta máquina | Teórico do 02 §4.3 (`geomesh`) |
|---|---:|---:|---:|---:|---:|---:|
| Zika-21 | 671,7 s | 409,5 s | **1,64x** | 396,2 s | 1,70x | 3,12x |
| VARV-6 | 8 314,7 s | 3 905,2 s | **2,13x** | 3 886,3 s | 2,14x | n/a (02 tem VARV-49: 1,78x) |

Leitura (medido para os números, inferido para a explicação):

- O speedup medido fica a **1-3 %** do limite do caminho crítico medido: o *overhead* do AkôFlow é pequeno frente às durações. Não é o limite do 02, porque este usa durações de outra máquina (o caminho crítico de 299 s do 02 vs 396 s aqui: aqui o NJ-parcimônia/RAxML disputam núcleos e vão a 330-375 s).
- Overhead do AkôFlow por run (soma sobre as 40 atividades, `run.json`): transferência de workspace 68 s (Zika) e 80 s (VARV-6); espera de partida do contêiner 19 e 26 s; preparação 18 e 21 s; fila 0. Volume transferido pelos 86 `rsync`: 109 MB (Zika) e **734 MB (VARV-6)**. Ao longo do caminho crítico (9 atividades), custa em torno de 15-25 s: coerente com o modelo "~5 ativações x *o* s" do 02 (inferido).
- O paralelismo médio (3,15 e 2,03) é menor que o máximo (8 e 7): o trabalho se concentra em poucos ramos longos (RAxML e NJ-parcimônia). Com mais núcleos, o teto é o caminho crítico, não o número de máquinas.
- O 02 §4.5 propõe quebrar o `--all` do RAxML-NG; não foi testado aqui.

## 6. Equivalência com o monolito (`python -m workflow.equivalencia`)

Comando (medido): `python -m workflow.equivalencia <ref/out> <saida/out> --json equiv/<x>.json`. Classes: H1 alinhamento, H2 árvores de distância/parcimônia, H3 árvores ML, H4 subárvores (SHA-256); H5 `.gb`; P padrões FPMax; S1-S5 semânticas.

**Zika-21 r6 contra o monolito local (`ref/zika`): equivalente nos artefatos científicos.**

| H1 | H2 | H3 | H4 | P | S1-S5 | H5 (`.gb`) |
|---|---|---|---|---|---|---|
| 4/4 | 4/4 | 16/16 | 182/182 | 13/13 | 5/5 | 0/1 |

A única diferença é a ordem dos registros do `.gb` (mesmo tamanho, 128 975 bytes; mesmo conjunto de 20 linhas `LOCUS` ordenadas: `cmp` sem diferença). Não afeta nada a jusante (alinhamento, árvores e padrões idênticos por hash).

**VARV-6 r3 contra o monolito local (`ref/varv6`): equivalente nas árvores de recorte `mafft`, divergente no braço `mafft_iterative`.**

| H1 | H2 | H3 | H4 | P | S2 (FPMax) | S4 (CQL) | S1, S3, S5 |
|---|---|---|---|---|---|---|---|
| 2/4 | 6/8 | 8/16 | 42/60 | 17/17 | igual | igual | diferem |

- O alinhamento `dataset_final_mafft_iterative.aln` (e a cópia em `tmp/`) difere: mesmos 6 táxons, mesmo tamanho (250 452 colunas), com posições distintas em **5 das 6 sequências** (8 colunas em `NC_008030.1` e 6 em cada uma de `NC_001611.1`, `L22579.1`, `DQ437592.1`, `DQ437591.1`; `NC_008291.1` idêntica; medido por comparação coluna a coluna). Todas as árvores derivadas desse alinhamento diferem (NJ-distância, NJ-parcimônia, FastTree, IQ-TREE, RAxML do braço iterativo e suas subárvores); as do braço `mafft` (retree) são idênticas. Os **17 padrões FPMax e o CQL, porém, coincidem**.
- Causa **inferida, não isolada**: o MAFFT iterativo com `--thread` é não determinístico entre execuções ou entre contextos de CPU (a G4b/02 já anota "sem semente" para o alinhador). Não rodei o MAFFT duas vezes no monolito para separar "ruído do MAFFT" de "efeito do AkôFlow" [não testado]. É o passo de maior prioridade para a equipe.
- **A referência VARV-6 está degradada** (medido): o monolito local escreveu em `raw_data_sequences.gb` o texto de erro `Failed to understand id NC_008030. ...` (accessions truncados, 210 bytes, 0 registros). O run AkôFlow (modo online, 6 registros, 2,6 MB) obteve os metadados. Logo S3 (metadados) e S1 (seleção de modelo, que depende do alinhamento) divergem **por defeito da referência e pelo MAFFT**, não por falha do AkôFlow; S5 (`outputs_sha256`) diverge por consequência. [Isso reproduz o D13 (rótulos truncados) no caminho do monolito.]
- O Zika r5 foi equivalente (H1-H4, P, S) a duas referências (`equiv/zika_r5_vs_R1.json`, `zika_r5_vs_Rlocal.json`).

## 7. O que a proveniência do AkôFlow captura (observado, `run.json` de 3,4-4,2 MB por run)

Seções do `GET /execution-runs/<id>/`: `activities` (por atividade: recurso planejado e alocado, `attempt`, `status`, `readyAt`, `queuedAt`, `startedAt`, `finishedAt`, `runtimeSeconds`, `queueSeconds`, `transferSeconds`, `transferBytes`, `overheadSeconds`, `cost`, metadados de partida e preparação), `activityAttempts`, `dataObjects` e `dataLocations` (294 por run VARV-6), `dataTransfers` (86; com estratégia, número de arquivos e segundos do `rsync`), `workspaces` e `workspaceLeases`, `events`, `timeline`, `handles`, `artifactMaterializations`, `cloudBilling` (zerado, local). O SQL `POST /provenance/sql/` expõe `execution_runs`, `task_executions` e `activity_dependencies`; foi a base de `analisar_run.py`.

Não capturou (medido): `auditoria.json`, `entidades.json` e `linhagem.json` voltaram vazios (2-19 bytes) nesta versão local; **entradas, versões das ferramentas, sementes e a proveniência *científica*** (qual ferramenta, qual comando, hash do alinhamento) não estão no AkôFlow: vêm dos fragmentos `out/outputs/prov/<atividade>.json` do próprio PhyloTreeMiner (comando, `entradas_sha256`, ambiente, versão). As duas camadas são complementares: o AkôFlow dá **tempo, recurso, transferência e dependência**; o fragmento dá **o que foi calculado**. O `rsync` registra bytes e arquivos, não hash do conteúdo [observado].

## 8. Problemas encontrados

1. Corrida do daemon na liberação de *leases* com dois consumidores simultâneos de uma mesma atividade (r1 dos dois conjuntos). Contorno: serializar. **Reportar à equipe.**
2. `expectedOutputs` cancelam a atividade quando falta um arquivo; mensagem clara, mas a granularidade é tudo ou nada.
3. O daemon poda as entradas herdadas do workspace final; exige uma atividade de empacotamento.
4. Workspace criado como root no host: a imagem precisa rodar como root (ou o daemon respeitar o `uid`).
5. Limites de memória precisam ser reais: o padrão de 16 Mi é inutilizável, e 2 Gi estoura no alinhamento de genomas de 250 kb.
6. O corpo de `POST /execution-runs/` exige o *bundle* completo; a doc não diz isso explicitamente (07 §5). Erro obtido: `FOREIGN KEY constraint failed`.
7. A coleta (`rodar_akoflow.py`) falhou uma vez ao ler um JSON vazio (`JSONDecodeError`, run zika21-r1); corrigida.
8. `rsync` de 734 MB no VARV-6 (alinhamentos de 1,5 MB por árvore, mas o workspace herdado inclui os de todos os predecessores): em nuvem ou HPC com rede lenta, isso vira custo relevante (inferido).
9. Reprodutibilidade científica: MAFFT iterativo não reprodutível bit a bit (§6).

## 9. O que levar para a equipe do Wesley

- Corrida no *lease* do workspace com dois consumidores (§8.1): vale um *issue* com o log do daemon do zika21-r1/varv6-r1.
- Como expor o token e a API a scripts sem o Desktop; o *bundle* de `execution-runs` documentado; `expectedOutputs` parciais (aviso em vez de cancelamento).
- Opção de preservar o workspace herdado no último nó (evitar a atividade de empacotamento), e de não transferir o workspace inteiro do predecessor (`dataDependencies` selecionando arquivos): pesa 734 MB no VARV-6.
- Registrar hash do conteúdo dos arquivos transferidos na proveniência.
- Perguntas abertas do 04: *overhead* por atividade em Kubernetes e em SLURM (aqui local: ~0,5 s de partida e ~1-2 s de *rsync* por atividade).

## 10. Ressalvas

- Uma máquina, 12 núcleos, runtime `local` (Docker). Não houve rede entre nós, nem Kubernetes, nem SLURM [não testado]. O *overhead* real de nuvem ou HPC será maior.
- Conjuntos pequenos; o único conjunto grande com alinhamento pesado é o VARV-6, com 6 genomas. Um run completo por conjunto para a conclusão de tempo; sem repetições, sem intervalo de confiança. A diferença entre os tempos de r4/r5/r6 (§3) mostra variação de até ~30 %.
- A concorrência com outros serviços da máquina (Neo4j e Postgres em execução) não foi isolada.
- O monolito de comparação rodou **na mesma imagem** e no mesmo Docker, mas com o `num_threads` do config (8); ver §6 sobre o MAFFT.
- Nenhuma execução aqui é publicável (IN-09); o VARV-6 não está entre os projetos pós-D28 catalogados.

## 11. Disco e encerramento

Medido ao final: dados em `/home/joh/akoflow_local/` ~1,2 GB (arquivos 685 MB + workspace 467 MB, parte sobreposta com `saida/` 304 MB); imagens `phylotreeminer-akoflow:local` 3,64 GB, `akoflow/daemon` 1,05 GB, `akoflow/buildkit` 242 MB; volumes do servidor ~30 MB. **Total do teste ~5-6,5 GB, abaixo de 15 GB.** Disco livre: 38 GB de 457 GB (92 % usado).

Containers `server-akoflow-server-1` e `server-buildkitd-1` parados (`docker compose stop`) e **mantidos** (as imagens ficam). Para retomar: `cd /home/joh/akoflow_local/server && docker compose -f compose.yaml -f compose.override.yaml up -d`. Para liberar tudo: `docker compose down -v` e `docker rmi phylotreeminer-akoflow:local`.
