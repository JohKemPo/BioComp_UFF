# 02 · Desenho granular do workflow para o AkôFlow

Meta G4a. Este documento propõe o DAG: as atividades atômicas, o contrato de cada uma, a análise de paralelismo, o plano de equivalência com o `workflow.py` monolítico e a proveniência. É **só desenho**. A implementação é a meta G4b e depende da aprovação deste texto. As etiquetas de evidência ([medido], [observado], [inferido], [não testado]) seguem a convenção de [01](01-workflow-atual.md#escopo-e-convenções).

## 1. Princípios

1. **Mesmo cálculo, outra orquestração.** Cada atividade chama as **mesmas funções** de biblioteca que os controladores chamam hoje. Nenhuma fórmula, parâmetro de ferramenta, semente ou política muda:
   - `AlignmentSeqs.align_sequences_mafft` e `resolve_aligner`;
   - `selecionar_modelo`;
   - os `TreeBuilder.*_constructor`;
   - `SubtreeBuilder.subtree_constructor`;
   - `SubtreeMiner.process_group` e `consolidate_fpmax_results`;
   - `parse_tree`, `process_histogram_frequence` e `ExecutionManifest`.

   O protocolo de zona sagrada (`04-rigor-cientifico.md §3`) só seria exigido se algum número mudasse, e o plano de equivalência (§5) é o oráculo que mostra que nenhum mudou. O monolito `workflow.py` **continua existindo** e é a referência.
2. **Uma ativação por chave de dados.** Há uma ativação por alinhador; uma por par (alinhamento × método de inferência); uma por árvore na extração de subárvores; e uma por árvore na datação, quando ela está ligada.
3. **Tudo que passa de uma atividade a outra é arquivo.** O estado em memória de [01 §5.1](01-workflow-atual.md#51-estado-em-memória-entre-etapas) vira arquivo:
   - a seleção de modelo vira `selecao.json`;
   - `raw_data` vira um JSON por árvore;
   - o cache de acessos vira o `.gb` completo, baixado uma vez;
   - o registro `tool_runs` vira um fragmento de proveniência por atividade.
4. **Caminhos únicos por ativação.** No AkôFlow, antes de uma atividade começar, o *workspace* inteiro de cada predecessora é copiado por `rsync`; arquivos idênticos no mesmo caminho são pulados e **conteúdo conflitante no mesmo caminho é rejeitado**. A fonte é `docs/docs/internal/workflow-spec.md` do repositório `UFFeScience/akoflow`, commit `aeea81c` (2026-10-23), seção *Top-level fields*. Daí a regra: duas ativações nunca escrevem o mesmo caminho relativo. O layout `out/` do monolito já é quase todo único por chave; os acréscimos são `out/outputs/prov/<atividade>.json`, `out/outputs/subarvores/<árvore>.json` e `out/outputs/stdout/<atividade>.txt`.
5. **Workspace limpo e entrada intocada.** Cada execução começa num *workspace* novo, o que fecha o AC-04 por construção. A deduplicação grava em `out/entrada/`, não em `data/`, o que fecha o AC-09 neste caminho. O monolito não muda, e o AC-09 continua aberto nele (F2.2).
6. **D17 e D21 não são negociáveis.** RAxML-NG sempre com `--threads N --workers 1`, nunca `auto`, e com o **mesmo N** do monolito no teste de equivalência. IQ-TREE e ModelFinder com `-nt 1`. A semente vem do `config_backup.json`, igual para todas as ativações (12345 nos dois projetos).
7. **Falha declarada, não engolida.** Com `on_method_failure: continue`, uma ativação de inferência que falha **sai com código 0** e grava o fragmento com estado `tentado_e_falhou` e o motivo, sem árvore. A extração de subárvores dessa chave vira *no-op* declarado. Com `fail` (o padrão), a ativação sai com código diferente de zero. O AkôFlow não inicia as sucessoras de uma atividade que falhou (`workflow-spec.md`, "A failed transfer leaves the successor unstarted"; para código de saída não-zero, ver a pergunta Q3 do [04](04-handoff-equipe-akoflow.md)).

## 2. As atividades atômicas

Ponto de entrada proposto: `python -m workflow.akoflow <subcomando> --config out/outputs/config_backup.json --out <projeto>/out [chave]`. É um módulo novo, entregue na G4b. Os **parâmetros científicos vêm do `config_backup.json`**, que é o mesmo arquivo do monolito; a linha de comando só diz a chave da ativação. Assim o YAML não tem como divergir da configuração.

O `--out` precisa terminar em `<nome_do_projeto>/out`, porque o `.cql` deriva o nome do projeto da estrutura de diretórios ([01 §5.2](01-workflow-atual.md#52-caminhos-fixos-e-derivados-por-texto)).

### 2.1 Tabela de atividades

As colunas CPU, memória e tempo trazem a fonte. "Zika" e "VARV" são os dois projetos de referência ([01 §3](01-workflow-atual.md#3-tempos-medidos-por-etapa)).

| Atividade | Subcomando e chave | Função de biblioteca reaproveitada | CPU (`cpuLimit`) | Tempo Zika / VARV | Memória |
|---|---|---|---|---|---|
| `preparar-entrada` | `preparar --entrada <dir>` | `eh_fasta` (`TBC:54-73`), `duplicate_*`, `validate_sequences`, `deduplicar_por_sequencia` (`dataValidation.py:8-155`) | 1 [inferido: Python puro] | 0,003 s / 0,085 s [medido, log] | [DADO]; proposta 512Mi [SUPOSIÇÃO] |
| `baixar-metadados` | `baixar-metadados --acessos …` (ou `--gb-congelado <arq>`) | `download_sequences` (`treeUtils.py:467-534`), chamado **uma vez** e em modo de acréscimo | 1 | 4,51 s / 35,05 s [medido, log; o tempo é dominado pelos `sleep(3)` entre lotes e pela latência do NCBI] | [DADO]; 512Mi [SUPOSIÇÃO] |
| `alinhar-<alg>` | `alinhar --alinhador <alg>` | `resolve_aligner`, `AlignmentSeqs.align_sequences_mafft(estrategia=…)` (`alignmentSeq.py:439-521`) | `num_threads` (4 / 16) [config] | ≤ 6,09 e ≈ 10,06 s / ≈ 602,7 e ≈ 2 674,2 s [medido/inferido, 01 §3] | **não medida** (`aligners.py:285`, `bytes_per_unit=None`); 2Gi [SUPOSIÇÃO] |
| `modelo-<alg>` | `selecionar-modelo --alinhador <alg>` | `selecionar_modelo` (`modelo_substituicao.py:289`), `tool_runs.registrar_selecao_modelo` | 1 (`-nt 1`) | 1,26 e 1,15 s / 9,95 e 9,79 s [medido, `modelfinder.log`] | ≤ 14 MB, declarado pelo próprio ModelFinder [observado]; 1Gi de folga [SUPOSIÇÃO] |
| `arvore-<alg>-{nj,upgma}-distance` | `inferir --alinhador <alg> --metodo nj_distance` | `TreeBuilder.distance_matrix` + `distance_constructor` (`builder.py:350-405`), `save_tree`, `save_tree_image` | 1 | ~0,1 s / ~25,7 s [medido, log] | [DADO] |
| `arvore-<alg>-{nj,upgma}-parsimony` | idem, `--metodo nj_parsimony` | `parsimony_constructor` (`builder.py:407-443`) | 1 | 246,4 e 78,8 s (braço iterativo do Zika) / não roda no VARV (`ignore_mode`) [medido, log] | [DADO] |
| `arvore-<alg>-fasttree` | `--metodo fasttree` | `fasttree_constructor` (`builder.py:520-564`) | 1 [observado: binário sem OpenMP] | 2,8–3,0 s / 233,8–237,0 s [medido, log do FastTree] | [DADO] |
| `arvore-<alg>-iqtree` | `--metodo iqtree` | `iqtree_constructor` (`builder.py:445-518`) | 1 (`-nt 1`, D21) | 5,1–6,2 s / 46,7 s [medido, log do IQ-TREE] | ≤ 22 MB, declarado pelo IQ-TREE [observado] |
| `arvore-<alg>-raxml` | `--metodo raxml` | `raxml_ng_constructor` (`builder.py:566-664`) | `raxml_threads` (4 / 8) [config] | 271,3–285,5 s / 2 039,7–2 129,0 s [medido, `.raxml.log`] | [DADO]; o RAxML-NG 2.0.2 não imprime estimativa |
| `arvore-<alg>-mrbayes` | `--metodo mrbayes` | `mrbayes_constructor` (`builder.py:777-940`) | 1 [inferido; MPI não verificado] | sem medição nestes projetos; timeout de 3 600 s (`builder.py:78`) | [DADO] |
| `subarvores-<alg>-<metodo>` | `subarvores --arvore out/Trees/<árvore>.nexus` | `SubtreeBuilder.subtree_constructor` (`subtree_construction/builder.py:67-177`), com o cache **pré-preenchido** com os rótulos da árvore, para não baixar nada | 1 | 0,28–0,42 s / 204,4–334,5 s com o `.gb` completo [medido, 2 árvores; 01 §3.2] | [DADO]. No VARV, o `metadata.json` agregado tem 244 MB; no VARV-49, que é pré-D28 e tem o `.gb` completo, chega a 990 MB [medido] |
| `heatmap-rf` | `heatmap-rf --stem <stem>` | `flatten_multi_trees`, `process_rf_distance`, `plot_heatmap_distances` (`metrics.py:108-232`), sobre `multi_trees` **reconstruído como no laço do monolito** (inclusive a duplicação, [01 §5.1](01-workflow-atual.md#51-estado-em-memória-entre-etapas)) | 1 | 1,6 s / ~2,4 s [medido, log] | [DADO] |
| `minerar-fpmax` | `minerar` | `SubtreeMiner.group_data_by_tree_base` + `process_group` (`miner.py:81-177`) sobre os JSON por árvore, em **ordem canônica** (nome do arquivo) | 1 | ≤ 0,56 s / ≤ 20 s [medido: FPMax + `.cql`, sem separação no log] | [DADO]; VARV na casa de GB [inferido do tamanho do JSON] |
| `gerar-cql` | `gerar-cql` | `parse_tree(trees=grupo, path, mode=base_name)` (`neo4jProcessing.py:17-46`), lendo o `metadata.json` minerado | 1 | contido no tempo de `minerar` | [DADO] |
| `graficos-suporte` | `graficos-suporte` | `process_histogram_frequence` (`metrics.py:420-459`) | 1 | 0,16 s / — [medido, log Zika 11:30:20,278→,439] | [DADO] |
| `datar-<alg>-<metodo>` (se `run_molecular_clock`) | `datar --arvore …` | `MolecularClockController._datar` (`molecularClockController.py:309-430`) | `molecular_clock_threads` (1) | **desligado** nos dois projetos | [DADO] |
| `datar-consolidar` (se ligado) | `datar-consolidar` | `_escolher_principal` + `_gravar` (`molecularClockController.py:433-495`) | 1 | — | — |
| `estabilidade` (extensão, **fora do monolito**) | `python -m workflow.stability.case_study --project … --fasta …` (CLI que já existe, `case_study.py:194-237`) | `StabilityAnalyzer`, `StabilityReport.run_all` | 1 | não medido | [DADO] |
| `consolidar-manifesto` | `consolidar` | `ExecutionManifest` (`manifest.py:273-650`): `register_*`, `register_outputs`, `finish` | 1 | não medido; no monolito acontece dentro do `finally` | [DADO] |

### 2.2 Contrato de entradas e saídas

Todos os caminhos são relativos ao *workspace* da ativação. `<p>` é o diretório do projeto, por exemplo `Zika_21seq_akoflow_<data>`; `<s>` é o stem do FASTA (`dataset_final`).

| Atividade | Lê | Escreve |
|---|---|---|
| `preparar-entrada` | `/app/data/<conjunto>/*.fasta` (da imagem) ou o volume de entrada; `/app/config/<cfg>.json` | `<p>/out/outputs/config_backup.json` (com `#` → nome do projeto, como em `workflow.py:52-75`); `<p>/out/entrada/<s>.fasta`, ou `<s>_NoPipe` se houve deduplicação; `<p>/out/entrada/acessos.txt`; `prov/preparar-entrada.json` (com `inputs_sha256` e o registro de `deduplicacao`) |
| `baixar-metadados` | `acessos.txt`; opcionalmente `--gb-congelado` | `<p>/out/outputs/raw_data_sequences.gb`; `prov/baixar-metadados.json` |
| `alinhar-<alg>` | `out/entrada/<s>*` | `out/tmp/<s>_<alg>.aln`; `out/Align/<s>_<alg>.aln` (cópia, como em `TBC:475`); `out/tmp/tmpFilesAlignment/`; `prov/alinhar-<alg>.json` |
| `modelo-<alg>` | `out/tmp/<s>_<alg>.aln` | `out/tmp/modelfinder_<s>_<alg>/{alinhamento.phylip, modelfinder.*}` (7 arquivos); `out/tmp/modelfinder_<s>_<alg>/selecao.json` (**novo**: a tradução usada pelos métodos); `prov/modelo-<alg>.json` |
| `arvore-<alg>-<m>` | `.aln`, e `selecao.json` para `iqtree`, `raxml` e `mrbayes` | `out/Trees/tree_<s>_<alg>_<m>.nexus`; `out/outputs/Plots/tree_<s>_<alg>_<m>.png`; `out/tmp/<ferramenta>_tree_<s>_<alg>_<m>/*`; `prov/arvore-<alg>-<m>.json` (estado M7.6, `tools_invoked`) |
| `subarvores-<alg>-<m>` | a árvore e `raw_data_sequences.gb` | `out/Subtrees/tree_<s>_<alg>_<m>_InnerN.nexus`; `out/outputs/subarvores/tree_<s>_<alg>_<m>.json` (o `result_dict` dessa árvore); `prov/subarvores-….json` |
| `heatmap-rf` | todas as `out/Trees/tree_<s>_*.nexus` e os fragmentos `arvore-*` (para saber quais foram executadas) | `out/outputs/Plots/tree_<s>_heatmap_of_distances.png` |
| `minerar-fpmax` | `out/outputs/subarvores/*.json` | `out/outputs/all_results_fpmax.csv` (ou `resul_of_fpmax_<x>.csv`); `out/outputs/metadata.json`; `out/outputs/metadata.csv` |
| `gerar-cql` | `metadata.json` | `out/outputs/neo4j_commands_tree_dataset.cql` (`neo4j_commands.cql` fora do modo OFST) |
| `graficos-suporte` | `metadata.json` | `out/outputs/Plots/support_frequence_{bar_plot,frequency_polygon}.png` |
| `consolidar-manifesto` | tudo (staging transitivo; ver a pergunta Q1 do 04) | `out/outputs/manifest.json`; `out/outputs/log_setup_<data>_<run_id>.log` (concatenação dos logs por atividade, na ordem canônica); `out/outputs/output_log.txt` (concatenação de `stdout/*`); opcionalmente `manifest.prov.json` (depende do F2.6 ou de `gerar_prov.py`) |

### 2.3 Imagem

Proposta: **uma imagem só**, como no exemplo `infrastructure-selection` (`spec.image`, um valor para todas as atividades).

- **Base:** a receita pinada da raiz do repositório (`environment.yml` com MAFFT 7.526, IQ-TREE 3.1.3, RAxML-NG 2.0.2, FastTree 2.2.0, MrBayes 3.2.7 e Clustal Omega 1.2.4, mais o `conda-lock.yml`), materializada com `mambaorg/micromamba:1.5.10-jammy`, como já faz `Backend/Dockerfile`.
- **Código:** `BioComp_UFF` no commit da G4b.
- **Dados de teste:** a entrada do Zika-21 e o `.gb` congelado, para o modo *offline*.
- **Sem `ENTRYPOINT`.** O exemplo diz "The image intentionally has no competing `ENTRYPOINT`". O `BioComp_UFF/Dockerfile` atual tem `ENTRYPOINT ["conda","run",…,"workflow.py"]` e um `environment.yml` **sem pinos**. Não serve.
- A versão de cada ferramenta faz parte do resultado (DEC-044). A imagem tem de passar no equivalente a `check_dependencies.sh --strict` durante o *build*.

Esta máquina (`joh`) **não** tem o env pinado: o RAxML-NG resolvido é o 1.2.2, o IQ-TREE o 3.0.1 e o MAFFT o v7.525 [medido: `raxml-ng --version`, `iqtree3 --version` e `mafft --version` no env `Phylotreeminer`]. Ela não serve como referência de equivalência.

Uma imagem por ferramenta (biocontainers) reduziria o tamanho, mas trocaria a pilha e o *wrapper* Python. Fica para depois da equivalência.

## 3. O DAG

```
preparar-entrada ─┬─► baixar-metadados ───────────────────────────────┐
                  └─► alinhar-<alg> ─┬─► modelo-<alg> ─┬─► arvore-<alg>-iqtree ──┐  │
                     (×A)            │                 ├─► arvore-<alg>-raxml ───┤  │
                                     │                 └─► arvore-<alg>-mrbayes ─┤  │
                                     ├─► arvore-<alg>-fasttree ──────────────────┤  │
                                     ├─► arvore-<alg>-{nj,upgma}-distance ───────┤  │
                                     └─► arvore-<alg>-{nj,upgma}-parsimony ──────┤  │
                                                                                 ▼  ▼
                                               subarvores-<alg>-<m> (uma por árvore; também depende de baixar-metadados)
                                               heatmap-rf (todas as árvores) │      │
                                               datar-<alg>-<m> (opcional) ───┘      ▼
                                                                         minerar-fpmax ─┬─► gerar-cql ────────┐
                                                                                         └─► graficos-suporte ─┤
                                                         heatmap-rf ─────────────────────────────────────────►┤
                                                                                          consolidar-manifesto ◄┘
```

Os diagramas em Mermaid, geral e expandido para o Zika-21, estão em [03](03-dataflow.md).

## 4. Análise de paralelismo

### 4.1 Método

O DAG foi montado por um script de modelo, `dag.py` (íntegra no [03, anexo A](03-dataflow.md#anexo-a--script-do-modelo-dag)), com `networkx` 3.2.1:

- a **largura** é a maior anticadeia, pelo teorema de Dilworth: *n* menos o emparelhamento máximo no fecho transitivo;
- o **caminho crítico** é o caminho mais longo, ponderado pelas durações medidas de [01 §3](01-workflow-atual.md#3-tempos-medidos-por-etapa).

Atividades sem medição, como `consolidar-manifesto`, `gerar-cql` e `graficos-suporte` no VARV, entram com 0 s, e o resultado é um **limite inferior** do caminho crítico. O modelo **não** inclui o *overhead* do AkôFlow: subir o contêiner, o `rsync` entre *workspaces* e a espera na fila. Esse é o [DADO] que a equipe mede (§4.5).

### 4.2 Contagem de ativações

| | Zika-21 | VARV-52 |
|---|---:|---:|
| `preparar-entrada` | 1 | 1 |
| `baixar-metadados` | 1 | 1 |
| `alinhar-*` | 2 | 2 |
| `modelo-*` | 2 | 2 |
| `arvore-*` | 14 (2 × [4 básicos + 3 ML]; MrBayes ignorado) | 10 (2 × [2 distância + 3 ML]; parcimônia e MrBayes ignorados) |
| `subarvores-*` | 14 (4 *no-op*: as árvores básicas do braço `mafft` falham por D28) | 10 |
| `heatmap-rf`, `minerar-fpmax`, `gerar-cql`, `graficos-suporte`, `consolidar-manifesto` | 5 | 5 |
| **Total de ativações** | **39** (84 arestas) | **31** (64 arestas) |
| **Largura máxima (anticadeia)** | **15** | **11** |
| Pico de ativações simultâneas no escalonamento ASAP, com recurso ilimitado e as durações medidas | 8 | 6 |
| Pico de CPU (soma dos `cpuLimit`) no ASAP | 13 | 33 |

Saída literal do modelo:

```
$ python dag.py zika
 "ativacoes": 39, "arestas": 84, "largura_max_anticadeia": 15,
 "caminho_critico_s": 298.6,
 "caminho": [["preparar-entrada",0.003],["alinhar-mafft-iterative",10.06],
             ["modelo-mafft-iterative",1.152],["arvore-mafft-iterative-raxml",285.761],["heatmap-rf",1.6]],
 "trabalho_total_s": 931.7, "T_mono_s": 931.566, "speedup_teorico_vs_mono": 3.12,
 "pico_ativacoes_simultaneas_ASAP": 8, "pico_cpu_ASAP": 13
$ python dag.py varv
 "ativacoes": 31, "arestas": 64, "largura_max_anticadeia": 11,
 "caminho_critico_s": 5078.4,
 "caminho": [["preparar-entrada",0.085],["alinhar-mafft-iterative",2674.2],["modelo-mafft-iterative",9.787],
             ["arvore-mafft-iterative-raxml",2039.65],["subarvores-mafft-iterative-raxml",334.457],["minerar-fpmax",20.234]],
 "trabalho_total_s": 11534.6, "T_mono_s": 9039.43, "speedup_teorico_vs_mono": 1.78,
 "pico_ativacoes_simultaneas_ASAP": 6, "pico_cpu_ASAP": 33
```

(A saída foi reformatada em menos linhas; os valores são os impressos.)

Com a datação ligada, o Zika-21 ganharia 6 `datar-*` (iqtree, raxml e fasttree × 2 braços) e 1 `datar-consolidar`. A `estabilidade` acrescenta 1 ativação.

### 4.3 Caminho crítico e speedup teórico

| | Zika-21 | VARV-52 |
|---|---|---|
| Caminho crítico | `alinhar-mafft-iterative` (10,1 s) → `modelo` (1,2 s) → `arvore-…-raxml` (285,8 s) → `heatmap-rf` (1,6 s) = **298,6 s** | `alinhar-mafft-iterative` (2 674,2 s) → `modelo` (9,8 s) → `arvore-…-raxml` (2 039,7 s) → `subarvores-…-raxml` (334,5 s) → `minerar` (20,2 s) = **5 078,4 s** |
| Monolito, medido | 931,6 s | 9 039,4 s |
| **Speedup teórico** (monolito ÷ caminho crítico) | **3,12×** | **1,78×** |
| O que limita | o RAxML-NG: 95 % do caminho crítico | o MAFFT iterativo (53 %) e o RAxML-NG (40 %) |

Ressalvas:

- **No VARV, o trabalho total sobe de 9 039 s para 11 535 s.** No monolito, 8 das 10 extrações de subárvores foram baratas (30–48 s) **porque o `.gb` tinha sido truncado para 3 registros** (F2.5, [01 §5.4](01-workflow-atual.md#54-o-genbank-sobrescrito-no-meio-da-extração-f25-já-conhecido)). No desenho granular, todas leem o `.gb` completo, e o modelo usa para todas os 334,5 s medidos na árvore processada com o `.gb` completo. É uma suposição conservadora [inferido]. O caminho crítico não depende dela: ele usa a medição da própria árvore `mafft_iterative_raxml`.
- O `.gb` completo é relido **por terminal, por subárvore** (`treeUtils.py:407-424`, chamado em `subtree_construction/builder.py:148`). É daí que vêm os 200–335 s por árvore no VARV contra 0,3 s no Zika [inferido]. Indexar o `.gb` uma vez por ativação baixaria esse tempo em ordens de grandeza, mas é mudança de código na zona de metadados: fica como insumo ([05](05-insumos-refatoracao.md)), não como parte deste desenho.
- **Os números são do `geomesh`** (Threadripper 2970WX). Em outro hardware, as durações mudam, e talvez a ordem do caminho crítico [não testado].

### 4.4 Granularidade de nível 2: o que daria mais paralelismo, e por que fica fora

| Ideia | Ganho potencial | Por que não entra na equivalência |
|---|---|---|
| Quebrar o `--all` do RAxML-NG em busca de ML + *k* blocos de `--bootstrap` + `--support` | Tira os ~2 000 s do bootstrap do caminho crítico [inferido] | As réplicas mudam com a semente de cada bloco, e o suporte FBP muda. Exige decisão científica (seed por bloco, e se o ML de `--search` é idêntico ao de `--all` com a mesma semente [não testado]) |
| `--threads` maior no RAxML-NG | Reduz o tempo do braço crítico | Que `N` não muda a árvore com `--workers 1` é afirmado no D26 (`02-defeitos…md:1073`), mas **não encontrei a medição**. A equivalência fica no **mesmo N**; variar `N` é um experimento à parte |
| `-nt N` no IQ-TREE | — | Proibido pelo D21 (3 topologias em 3 repetições com `-nt 4`) |
| Mais threads no MAFFT | Reduz o braço crítico do VARV | É seguro para a equivalência só se o alinhamento continuar idêntico. No Zika, `--thread 4` e `--thread 8` deram o mesmo SHA-256 [observado, §5.2], mas não testamos em genoma longo |

### 4.5 O que a equipe do AkôFlow mede, e que este documento não tem

O *overhead* por ativação (*pull* da imagem, partida do *pod*, `rsync`), multiplicado pelas ~5 ou 6 ativações do caminho crítico, entra direto no *makespan*. Com 39 ativações e um caminho crítico de 299 s no Zika, um *overhead* de *o* segundos por ativação custa cerca de 5*o* no *makespan* [inferido]. É um [DADO] do FGCS.

## 5. Equivalência com o monolito

### 5.1 Pergunta e hipótese

*O workflow granular, executado pelo AkôFlow, produz os mesmos artefatos científicos que o `workflow.py` monolítico, a partir da mesma entrada, da mesma configuração, da mesma imagem e do mesmo tipo de nó?*

H0 (aceite): em toda classe declarada determinística, a diferença é zero. A tolerância está fixada **aqui, antes de qualquer execução** (§5.4).

### 5.2 De onde pode vir não determinismo, e o que o desenho faz com cada fonte

| Fonte | Existe? Evidência | Tratamento |
|---|---|---|
| **Ordem e concorrência entre atividades** | Não muda resultado. Cada ativação é um processo com entradas explícitas e semente própria (`random_seed` repassada por chamada: `TBC:1016-1019`, `builder.py:25-51`). Nenhum RNG global é consumido em sequência entre métodos [inferido da leitura; o único `random` sem semente está no subamostrador do Clustal Omega, `alignmentSeq.py:79`, que não roda com MAFFT] | Nenhum. A concorrência muda o tempo, não o resultado |
| **Threads dentro de uma ferramenta** | O IQ-TREE com `-nt 4` deu 3 topologias em 3 repetições; o RAxML-NG com `--threads 4 --workers 1` deu 1 em 3 (D21, tabela de medição). O `auto` do RAxML-NG deu RF = 8 (D17) | `-nt 1` e `--workers 1` com N igual ao do monolito. MAFFT: `--thread 4` (execução de 2026-10-17) e `--thread 8` (2026-09-01) deram SHA-256 idênticos nos dois alinhamentos do Zika [observado: `325cc07c…` e `7c0751c3…`] |
| **Ordem de agregação** | `os.listdir` sem `sorted` ([01 §5.5](01-workflow-atual.md#55-ordem-dependente-do-sistema-de-arquivos)) define a ordem do `metadata.json`, do `metadata.csv`, do `.cql` e das linhas do CSV do FPMax | O granular usa a ordem canônica (nome do arquivo). A comparação **canonicaliza** antes de comparar (§5.4, classe S) |
| **Hardware (SIMD)** | O IQ-TREE escolhe o *kernel* sozinho ("Kernel: AVX+FMA", no log). O RAxML-NG tem `--simd none\|sse3\|avx\|avx2`, padrão automático [observado no `--help` do 1.2.2 desta máquina; conferir no 2.0.2]. No `--help` do IQ-TREE 3.0.1 não achei opção de SIMD (`grep simd\|kernel\|avx`) [observado]. Ótimos quase empatados (D17, D21: ΔlnL de 0,001) podem trocar de lado com outro caminho de ponto flutuante [não testado] | No teste E1, **mesmo tipo de nó** para o monolito e para o granular. Em nós heterogêneos (E2), declarar `--simd` no RAxML-NG e fixar o tipo de nó do IQ-TREE |
| **Entrada externa que muda com o tempo** | Registros GenBank baixados ao vivo | Modo *offline*: `baixar-metadados --gb-congelado` usa o `.gb` da execução de referência. O modo *online* é diagnóstico (§5.5, passo 6) |
| **F2.5 e D13** | O metadado depende da ordem e dos rótulos truncados; no VARV-52, o `.gb` final tem 3 registros [medido] | O Zika-21 não manifesta o defeito: 1 download, 20 registros. Por isso **a equivalência de metadados se faz no Zika-21**. No VARV, só as classes H valem até F2.4 e F2.5 fecharem |
| **Artefatos com carimbo de tempo** | `.iqtree`, `.log`, `.raxml.log`, `.ckp.gz` e o manifesto têm datas, durações ou caminhos | Excluídos da comparação por hash (classe X) |

Evidência adicional de estabilidade: as árvores de **FastTree e IQ-TREE** do Zika-21 saíram com **SHA-256 idêntico** em duas execuções do monolito com um mês de intervalo (2026-09-01 e 2026-10-17, mesma máquina) [observado]:

```
6ba5936a37c8 6ba5936a37c8 tree_dataset_final_mafft_fasttree.nexus
c3b6369f451b c3b6369f451b tree_dataset_final_mafft_iqtree.nexus
08d0fda8a91f 08d0fda8a91f tree_dataset_final_mafft_iterative_fasttree.nexus
3f38f36385c1 3f38f36385c1 tree_dataset_final_mafft_iterative_iqtree.nexus
```

A execução de 2026-09-01 é **anterior ao D28**: as árvores NJ, UPGMA e de parcimônia dela diferem por causa da correção, e as de RAxML-NG porque ela não usava `--all`. Ela só serve de controle para os métodos que o D28 não toca.

### 5.3 Execuções do teste (Zika-21, ~16 min por execução monolítica no `geomesh`)

| Rótulo | O quê | Para quê |
|---|---|---|
| R1, R2 | `python workflow.py -p /app/config/zika21.json`, **dentro da imagem**, no mesmo tipo de nó das ativações, cada uma num diretório novo | Referência e **controle nulo**: R1 × R2 mede o não determinismo intrínseco do monolito. Se ele não for zero numa classe, a classe muda de categoria antes de olhar o granular |
| G1, G2 | Workflow granular no AkôFlow, mesma imagem, mesmo tipo de nó, `--gb-congelado` = `.gb` de R1 | Teste e réplica |
| X | o projeto existente `Zika_21seq_validacao_20261017` (`geomesh`, fora do contêiner) | Diagnóstico de hardware e ambiente (E2), sem valor de aceite |

### 5.4 Classes, comparação e tolerância (fixadas antes de rodar)

| Classe | Artefatos (Zika-21) | Como comparar | Tolerância E1 (mesmo tipo de nó) | E2 (hardware diferente; diagnóstico) |
|---|---|---|---|---|
| **H1** | `out/Align/*.aln` (2), `out/tmp/*.aln` (2) | SHA-256 | **0 diferenças** | 0 |
| **H2** | `out/Trees/*_{nj,upgma}_{distance,parsimony}.nexus` (4) | SHA-256 | **0** | RF = 0 |
| **H3** | `out/Trees/*_{fasttree,iqtree,raxml}.nexus` (6), mais `.raxml.support`, `.raxml.bestTree`, `.treefile`, `.contree` em `tmp/` | SHA-256 | **0** | RF = 0, não enraizado, com rótulos normalizados por `strip_accession_version`. Suporte e comprimento de ramo: relatar Δ, sem portão |
| **H4** | `out/Subtrees/*.nexus` (182) | SHA-256; nomes iguais | **0** | consequência de H2 e H3 |
| **H5** | `out/outputs/raw_data_sequences.gb` | SHA-256 (no modo *offline* é cópia) | **0** | — |
| **S1** | seleção de modelo: `modelo_escolhido`, a lista `candidatos` com `lnL`, `AIC`, `AICc`, `BIC` e `w_BIC`, e `traducao` | igualdade campo a campo dos valores impressos pelo ModelFinder | **0** | modelo escolhido igual |
| **S2** | `all_results_fpmax.csv` | sem a coluna de índice; `itemsets` convertido em tupla ordenada; conjunto de linhas `(itemset, support, min_support_threshold, max_support_threshold, n_trees)` | **igualdade de conjuntos** | idem, se RF = 0 |
| **S3** | `metadata.json`, `metadata.csv` | JSON canônico: lista ordenada por nome de árvore, chaves ordenadas, números exatos; CSV sem índice, linhas ordenadas | **igualdade** | idem |
| **S4** | `neo4j_commands_tree_dataset.cql` | multiconjunto de comandos (separados por `;`), com espaço em branco normalizado e `project: '<nome>'` trocado por um marcador | **igualdade** | idem |
| **S5** | `manifest.json` | conjunto de `inference_methods` `(metodo, alinhador, estado, saida)`, mais `motivo` normalizado; `model_selection`; `tools_invoked` com caminhos normalizados (esperada uma diferença: o caminho de entrada do MAFFT, `data/…` contra `out/entrada/…`); chaves e valores de `outputs_sha256`; `params`, salvo `input_path`; `reproducibility`; `execution_mode`; `outcome.status` | **igualdade**, com a exceção nomeada | — |
| **P** | PNG (13) | mesmo conjunto de nomes; arquivo não vazio | nomes iguais | — |
| **X** (excluídos) | logs, `output_log.txt`, `.iqtree`, `.log` das ferramentas, `.ckp.gz`, `.rba`, `run_id`, carimbos de tempo, `environment`, `git` | — | — | — |

O `outputs_sha256` do monolito já cobre H1, H2, H3 (só as `.nexus`), H4 e os dois CSV ([01 §5.3](01-workflow-atual.md#53-um-processo-só-escreve-o-manifesto)). A comparação de base é, portanto, `manifest_R1["outputs_sha256"]` contra `manifest_G1["outputs_sha256"]`, chave a chave. Os CSV entram canonicalizados (S2, S3), porque a ordem das linhas depende da ordem de agregação.

**Regra de decisão (pré-registrada).**

1. Se R1 × R2 tiver diferença numa classe H ou S, essa classe é reclassificada e o achado vai para o ledger **antes** de olhar G.
2. G1 é aceito se, em toda classe H e S que passou no controle nulo, a diferença for **zero**, salvo a exceção nomeada em S5.
3. G2 repete G1.
4. Qualquer diferença em H2 ou H3 dispara o diagnóstico com o oráculo independente: `docs/science/scripts/oraculo_rf_dendropy.py` para RF, e ΔlnL a partir dos logs. O diagnóstico **não** muda o veredito.

### 5.5 Procedimento

1. Construir a imagem e registrar o *digest* e as versões das ferramentas dentro dela (equivalente ao `check_dependencies.sh --strict`).
2. R1 e R2 no tipo de nó *N*, com `/app/config/zika21.json` igual ao `config_backup.json` do Zika-21 de referência, trocando só os caminhos.
3. Copiar `R1/out/outputs/raw_data_sequences.gb` para a entrada do granular como `--gb-congelado`.
4. G1 e G2 no AkôFlow, com todas as ativações no tipo de nó *N*.
5. Rodar `comparar_equivalencia.py R1 G1` (entrega da G4b), que imprime uma tabela por classe com as contagens de iguais e diferentes. Fazer o mesmo com R1 × R2 e R1 × G2.
6. Diagnóstico opcional: G3 com `baixar-metadados` *online*, contra R1, para ver a deriva do GenBank; e X × R1 (E2).
7. Registrar no ledger: comando e saída literal de cada comparação.

## 6. Cobertura das saídas (C3)

Cada um dos **274 arquivos** do `Zika_21seq_validacao_20261017/out`, com a mesma contagem de [01 §6](01-workflow-atual.md#6-saídas-de-um-projeto-real), aponta para exatamente uma atividade.

| Padrão | Zika-21 | Atividade que produz |
|---|---:|---|
| `out/Align/<aln>.aln` | 2 | `alinhar-<alg>` |
| `out/tmp/<aln>.aln` | 2 | `alinhar-<alg>` |
| `out/tmp/modelfinder_<aln>/*` (7 por braço) | 14 | `modelo-<alg>` |
| `out/Trees/<aln>_{nj,upgma}_{distance,parsimony}.nexus` | 4 | `arvore-<alg>-<m>` |
| `out/Trees/<aln>_{fasttree,iqtree,raxml}.nexus` | 6 | `arvore-<alg>-<m>` |
| `out/tmp/iqtree_*/*` (9 por árvore) | 18 | `arvore-<alg>-iqtree` |
| `out/tmp/fasttree_*/*` (2 por árvore) | 4 | `arvore-<alg>-fasttree` |
| `out/tmp/raxml_*/*` (10 por árvore) | 20 | `arvore-<alg>-raxml` |
| `out/outputs/Plots/<árvore>.png` | 10 | `arvore-<alg>-<m>` |
| `out/outputs/Plots/tree_dataset_final_heatmap_of_distances.png` | 1 | `heatmap-rf` |
| `out/Subtrees/*_InnerN.nexus` | 182 | `subarvores-<alg>-<m>` |
| `out/outputs/raw_data_sequences.gb` | 1 | `baixar-metadados` |
| `out/outputs/all_results_fpmax.csv` | 1 | `minerar-fpmax` |
| `out/outputs/metadata.json` | 1 | `minerar-fpmax` |
| `out/outputs/metadata.csv` | 1 | `minerar-fpmax` |
| `out/outputs/neo4j_commands_tree_dataset.cql` | 1 | `gerar-cql` |
| `out/outputs/Plots/support_frequence_*.png` | 2 | `graficos-suporte` |
| `out/outputs/config_backup.json` | 1 | `preparar-entrada` |
| `out/outputs/manifest.json` | 1 | `consolidar-manifesto` |
| `out/outputs/log_setup_<data>_<run_id>.log` | 1 | `consolidar-manifesto`, por concatenação dos logs de atividade |
| `out/outputs/output_log.txt` | 1 | `consolidar-manifesto`, por concatenação de `stdout/*` |
| **Soma** | **274** | sem resto |

Conferência da soma: 2 + 2 + 14 + 4 + 6 + 18 + 4 + 20 + 10 + 1 + 182 + 1 + 1 + 1 + 1 + 1 + 2 + 1 + 1 + 1 + 1 = 274.

Além dos 274:

- o diretório `out/tmp/tmpFilesAlignment/` vem de `alinhar-<alg>`;
- o transitório `data/<conjunto>/dataset_final_NoPipe` é substituído por `out/entrada/<s>_NoPipe`, gravado por `preparar-entrada`;
- as saídas condicionais de [01 §6](01-workflow-atual.md#6-saídas-de-um-projeto-real) também têm dono:
  - `mrbayes` → `arvore-<alg>-mrbayes`;
  - `molecular_clock/**` → `datar-*` e `datar-consolidar`;
  - `resul_of_fpmax_<s>.csv` → `minerar-fpmax`;
  - `neo4j_commands.cql` → `gerar-cql`;
  - `tree_Acumulate_heatmap…` → `heatmap-rf`;
  - Clustal Omega e MUSCLE → `alinhar-<alg>`.

Os arquivos **novos** que só o granular produz são `out/outputs/prov/*.json`, `out/outputs/subarvores/*.json`, `out/outputs/stdout/*`, `out/entrada/*` e `selecao.json`. Nenhum deles tem sufixo que `register_outputs` hasheia (`.nexus`, `.nwk`, `.aln`, `.csv`), e por isso não alteram as chaves de `outputs_sha256` [inferido de `manifest.py:413-430`].

## 7. Proveniência

### 7.1 Fragmento por atividade

Cada ativação roda num processo novo, então o registro `tool_runs` já começa vazio. Ao terminar, mesmo com falha, o *wrapper* grava `out/outputs/prov/<atividade>.json` com:

```json
{
  "fragmento_versao": 1,
  "atividade": "arvore-mafft-raxml",
  "chave": {"stem": "dataset_final", "alinhador": "mafft", "metodo": "raxml"},
  "inicio_utc": "…", "fim_utc": "…",
  "desfecho": {"status": "concluido|falhou", "motivo": null},
  "entradas_sha256": {"out/tmp/dataset_final_mafft.aln": "…", "out/tmp/modelfinder_dataset_final_mafft/selecao.json": "…"},
  "saidas_sha256": {"out/Trees/tree_dataset_final_mafft_raxml.nexus": "…"},
  "tools_invoked": "<tool_runs.execucoes()>",
  "inference_methods": "<tool_runs.metodos()>",
  "model_selection": "<tool_runs.selecoes_modelo()>",
  "reproducibility": "<reproducibility_settings(tree_config)>",
  "tools_available": "<manifest.tool_versions()>",
  "environment": "<manifest._environment()> + cota de CPU do cgroup, quando legível",
  "akoflow": {"run_id": "[se exposto ao contêiner — Q2]", "atividade": "arvore-mafft-raxml"}
}
```

Os caminhos são gravados **relativos à raiz do projeto**, já no fragmento, pela mesma regra de `ExecutionManifest._relativo` e `_fora_da_raiz`. Assim, nenhum caminho absoluto de *workspace* ou de *home* entra no manifesto (D15).

`psutil.cpu_count()` dentro de um contêiner devolve as CPUs **do host**, não a cota [inferido]. O ambiente registrado precisa ler a cota do cgroup, ou o manifesto vai declarar 48 CPUs onde havia 4.

### 7.2 Recomposição do `manifest.json` (v2, mesma forma)

`consolidar-manifesto` cria um `ExecutionManifest(project_root=<p>, params=config_backup)`. Depois:

1. `register_input` com o SHA-256 da entrada, tirado do fragmento de `preparar-entrada`;
2. `register_reproducibility` e `register_execution_mode`, pelo mesmo código de `workflow.py:100-128`, sobre a configuração e as ferramentas da imagem;
3. para cada fragmento, **na ordem canônica do monolito** (por alinhador: `nj_distance`, `nj_parsimony`, `iqtree`, `fasttree`, `raxml`, `mrbayes`, `upgma_distance`, `upgma_parsimony`; por último `relogio_molecular`, como em `TBC:684-712` e no manifesto observado), chama `register_tool_run`, `register_method_state` e `register_model_selection`;
4. `register_outputs(out)` e `register_outputs(out/outputs/molecular_clock, suffixes=('.nex','.json'))`, iguais a `workflow.py:167-171`;
5. `register_outcome`. `started_at` é o menor `inicio_utc` e `finished_at` o maior `fim_utc`;
6. acrescenta um campo **aditivo**, `execucao_granular`, com a lista de atividades (tempos, desfecho, ambiente e ids do AkôFlow). `manifest_version` continua 2, pelo mesmo critério dos campos de M7.6: "acrescentado sem mudar a forma de nenhum outro".

O mecanismo já existe em escala menor: `anexar_execucao_isolada` (`manifest.py:144-229`) anexa a etapa de datação rodada em outro processo, marcando `"via": "isolado"`. A consolidação generaliza esse caminho.

### 7.3 Mapeamento para a proveniência do AkôFlow

Fontes:

- o artigo do SBBD 2024, §3.2: o modelo tem as classes *Workflow*, *Activity*, *Metrics* e *Logs*, guarda tudo em SQLite e coleta métricas pela *Metrics API* com janela mínima de 15 s;
- a documentação no commit `aeea81c`: `guides/data/provenance-and-audit.md` (runs, plans, task attempts, transfers, artifact manifests, linhagem, SQL somente leitura) e `internal/activity-telemetry.md` (CPU, memória e I/O de bloco por cgroup v2 a cada 10 s, com amostra inicial e final);
- o README do exemplo `infrastructure-selection`: a tabela *Generated files* mostra nome, tamanho e SHA-256.

**A versão do AkôFlow a citar no FGCS é a que rodar o experimento**, e o artigo de 2024 descreve uma arquitetura anterior.

| Manifesto PTM | AkôFlow | Ligação proposta |
|---|---|---|
| `run_id` | id da *execution run* | `consolidar` grava os dois. Depende de o AkôFlow expor o id dentro do contêiner (Q2) |
| `execucao_granular.atividades[].inicio/fim` | *task attempt*: tempo de execução, fila, *startup* | Os do AkôFlow incluem fila e *staging*. A diferença entre os dois **é** o *overhead* (§4.5) |
| `outputs_sha256` | *Generated files* (SHA-256) e *artifact manifest* | Checagem cruzada: o mesmo caminho tem de ter o mesmo hash nos dois registros. É uma segunda prova de integridade, independente |
| `inputs_sha256` e as entradas dos fragmentos | *transfers* (bytes, origem, destino) e `dataDependencies` | A linhagem arquivo → atividade fica nos dois |
| `tools_invoked[].command` (comando **da ferramenta**) | `command` da atividade (o *wrapper*) | São complementares: o AkôFlow vê o *wrapper*, o PTM vê o `raxml-ng …` real |
| `tools_available`, `environment` | *runtime* e *resource* (tipo de máquina) | O PTM registra o que o processo vê; o AkôFlow, onde ele foi alocado |
| — (o PTM não mede) | *activity metrics*: CPU, pico de memória, I/O | **É a fonte das memórias de §2.1**, hoje [DADO] |
| `log_file`, `output_log.txt` | *Logs* (stdout) da atividade | Iguais em conteúdo; o PTM concatena |
| `manifest.prov.json` (DEC-147: `gerar_prov.py`; nativo no F2.6) | grafo de linhagem (exportável em JSON) | O PROV do PTM passa a ter `startTime` e `endTime` por passo, que hoje faltam ("O manifesto não registra horário por chamada", `gerar_prov.py:28-30`) |

## 8. Decisões que este desenho pede

| # | Decisão | Opções | Recomendação |
|---|---|---|---|
| D-a | Separar `gerar-cql` de `minerar-fpmax` | (1) separado, como pedido; (2) junto, como em `miner.py:64-79` | (1). Mas a G4b prova com S4 que `parse_tree` sobre o `metadata.json` relido dá o mesmo `.cql` [inferido; não testado] |
| D-b | Como `baixar-metadados` obtém os acessos | (1) cabeçalhos do FASTA; (2) terminais das árvores, como o monolito | (1). No Zika os conjuntos coincidem [observado: 20 acessos, 1 download]. No VARV o comportamento muda junto com F2.4 e F2.5, que estão no plano M8 |
| D-c | Onde fica a entrada | (1) dentro da imagem, como no exemplo (`COPY data/`); (2) volume ou *storage* do AkôFlow | (1) para o Zika (55 820 B) e o VARV (10,3 MB); (2) para o ZIKV-480 e maiores |
| D-d | Nível 2 do RAxML-NG (§4.4) | — | Fora da G4b. Exige parecer científico |
| D-e | `estabilidade` no DAG | (1) atividade opcional; (2) fora | (1), marcada como extensão e fora da equivalência |
