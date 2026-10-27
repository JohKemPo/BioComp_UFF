# 05 · Insumos para a refatoração em biblioteca + workflow agnóstico

Meta G4a. Este documento é **só mapeamento**: a refatoração fica para depois dos dois manuscritos (painel da frente, seção "Fora desta frente"). Ele junta o que o desenho granular revelou e liga cada achado às portas do *blueprint* de infraestrutura agnóstica: o artefato `SzU2c9obZKcYhpfL9UdtT6` e `docs/audit/99-futuro-infra-agnostica.md`, com as portas `ComputePort`, `StoragePort` e afins, e a migração em fases *strangler-fig*.

## 1. A tese, em uma linha

**O paralelismo de hoje é pior do que o modelável** porque o domínio e a orquestração estão no mesmo objeto.

- **Hoje.** O único paralelismo está dentro das ferramentas: `mafft --thread`, `raxml-ng --threads`. O laço alinhadores × métodos é sequencial (`treeBuilderController.py:682-714`), `parallelWorkflow.py:17-18` é `NotImplemented`, e o `parsl` está nas dependências sem uso (DEC-127).
- **No modelo granular.** São 15 ativações simultâneas possíveis no Zika-21, e um speedup teórico de 3,1× (Zika-21) e 1,8× (VARV-52), sem mudar nenhum cálculo ([02 §4](02-desenho-granular.md#4-análise-de-paralelismo)).

## 2. O que o desenho obrigou a fazer, e o que isso diz sobre a biblioteca

| # | Achado | Hoje (arquivo:linha) | O que a biblioteca precisa ter | Porta / fase do *blueprint* |
|---|---|---|---|---|
| R1 | Controladores misturam orquestração, cálculo e I/O | `TreeBuilderController` (1 359 linhas) faz laço, cache, alinhamento, inferência, PNG e heatmap | Funções puras por passo, com contrato de arquivo, que são o núcleo das CLIs da G4b: `alinhar(fasta, alinhador, threads) → aln`; `selecionar_modelo(aln, semente) → selecao`; `inferir(aln, metodo, selecao, semente, threads) → arvore`; `subarvores(arvore, gb) → json`; `minerar(jsons, cfg) → resultado`; `gerar_cql(resultado, projeto)` | núcleo de domínio (fase 1) |
| R2 | Proveniência num registro **global do processo** | `tool_runs.py:38-81`; drenado em `workflow.py:166` | Cada passo devolve ou grava o próprio registro (fragmento). O manifesto vira uma **função de agregação** de fragmentos. `anexar_execucao_isolada` (`manifest.py:144-229`) já é esse caso em pequena escala | `JobStorePort` / proveniência |
| R3 | Estado em memória entre passos | `_selecoes_modelo` (`TBC:1066`), `multi_trees` (`TBC:387`), `downloaded_sequences_cache` (`subtreeBuilderController.py:89`) | Artefatos explícitos: `selecao.json`, o conjunto de árvores lido do disco, o `.gb` completo | `StoragePort` |
| R4 | Caminhos derivados por manipulação de texto | `split('/Trees')` (`builder.py:452`, `:527`, `:573`); `replace('Trees','outputs/Plots')` (`TBC:589`, `:775`, `:839`); nome do projeto deduzido do caminho (`neo4jProcessing.py:6-15`) | Um objeto de layout (`LayoutDoProjeto`) que **dá** os caminhos. O nome do projeto vai como parâmetro | `StoragePort` (fase 2: `fsspec`) |
| R5 | Reuso por "o arquivo existe" | `TBC:747-754`, `:809-816`, `:911-917`; `subtree_construction/builder.py:135-140` (AC-04) | Cache por **conteúdo**: hash da entrada + parâmetros + versões. Ou nenhum cache, e o orquestrador decide | `StoragePort` + cache por hash (*trade-off* 03 do *blueprint*) |
| R6 | Escrita no diretório de entrada | `dataValidation.py:133`; `clean_NoPipe` (`TBC:188`, `:474`) (AC-09) | Entrada somente leitura. A deduplicação é um passo com saída própria | `StoragePort` |
| R7 | Rede no meio do pipeline | `download_sequences` dentro da extração de subárvores (`subtree_construction/builder.py:94-105`); `Entrez.email` fixo (`treeUtils.py:480`) | Aquisição de metadados como passo separado, no começo, com modo *offline* (`.gb` congelado) para reprodução | `SecretsPort` (chave NCBI) + `ComputePort` |
| R8 | Ordem dependente do sistema de arquivos | `os.listdir` sem `sorted` (`subtreeBuilderController.py:81`) | Ordem canônica declarada em toda agregação | domínio |
| R9 | Releitura do `.gb` por terminal e por subárvore | `fetch_local_record` faz `SeqIO.parse` inteiro a cada chamada (`treeUtils.py:407-424`), chamado em `subtree_construction/builder.py:148`. É daí que vêm 200–335 s por árvore no VARV-52 [inferido] | Índice do `.gb` (acesso → registro) carregado uma vez por passo. **Zona sagrada** (extração de metadado): exige o protocolo, junto com F2.5 | domínio |
| R10 | Configuração sem esquema | `setattr` de todo `kwarg` (`TBC:117-120`) (AC-02) | Esquema tipado, que é a fonte única das CLIs e do gerador do DAG | domínio |
| R11 | Log e stdout globais por processo | `run_logging.py:88-123`; `sys.stdout` redirecionado (`workflow.py:132-133`) | Um *logger* por passo. A concatenação é responsabilidade de quem agrega | `EventBusPort` (progresso) |
| R12 | O ambiente registrado vê o host, não a cota | `psutil.cpu_count` em `_environment()` [inferido] | Ler a cota do cgroup. O orquestrador informa o recurso alocado | `ComputePort` |
| R13 | Heatmap de RF com o RF errado e com árvores duplicadas | `metrics.py:13-58`; `TBC:809-816` ([01 §5.6](01-workflow-atual.md#56-heatmap-acumulado-nunca-sai-e-o-rf-do-heatmap-tem-defeito-achado-lateral)) | Usar a identidade canônica de `stability/` também para a figura. **É mudança de número numa figura**: vai para a triagem | domínio (zona sagrada) |

## 3. O workflow agnóstico que sai daqui

A G4b vai produzir as três peças que viram, na refatoração, o "workflow agnóstico de infraestrutura":

1. **CLIs por passo**, que são a unidade de trabalho de qualquer executor;
2. **um gerador de DAG a partir da configuração** (hoje o modelo `dag.py` do [03, anexo A](03-dataflow.md#anexo-a--script-do-modelo-dag));
3. **um consolidador de proveniência**.

Os executores viram adaptadores da `ComputePort`, todos consumindo o mesmo DAG:

| Adaptador | Estado | Observação |
|---|---|---|
| sequencial local | existe hoje: o `workflow.py` | continua como referência de equivalência |
| *pool* local (processos) | não existe | ganho imediato numa máquina só, sem infraestrutura nova [inferido] |
| Parsl | dependência instalada, sem uso | `python_app`/`bash_app` por passo; os *providers* cobrem SLURM e Kubernetes |
| AkôFlow | esboço YAML ([03](03-dataflow.md)) | o DAG vira `workflow.yaml` |

O critério de aceite de todo adaptador novo é o **mesmo teste de equivalência** de [02 §5](02-desenho-granular.md#5-equivalência-com-o-monolito). Ele vira o *golden test* da refatoração, exigido pela regra 4 do `CLAUDE.md`.

## 4. Ordem sugerida (para o planejamento futuro, não agora)

1. **Fechar antes F2.2 (AC-09), F2.3 (AC-04), F2.4 (D13) e F2.5 (`.gb`) no monolito.** Esses quatro são os que impedem a equivalência de metadados no VARV ([02 §5.2](02-desenho-granular.md#52-de-onde-pode-vir-não-determinismo-e-o-que-o-desenho-faz-com-cada-fonte)).
2. **Extrair R1 e R2** (passos puros + fragmentos). O `workflow.py` vira um adaptador sequencial sobre eles, e a saída tem de ser bit a bit igual à de hoje.
3. **R3 a R8**, que são os contratos de arquivo, e o gerador de DAG.
4. **Os adaptadores** de *pool* local, Parsl e AkôFlow.
5. **R9 e R13**, que mudam número: só com o protocolo de zona sagrada e parecer no ledger.
