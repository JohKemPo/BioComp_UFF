# 01 · Como o `workflow.py` executa hoje

Meta G4a da frente de manuscritos de 2026-10. Este documento descreve o pipeline **como ele é**, para que o desenho granular ([02](02-desenho-granular.md)) possa ser conferido contra ele, linha por linha.

## Escopo e convenções

- **Código lido:** o submódulo `BioComp_UFF` no commit `29833ff`, no worktree `AKOFLOW_Phylotree_fgcs`. Todo `arquivo:linha` abaixo se refere a esse commit. Os caminhos são relativos à raiz do submódulo.
- **Execuções de referência.** As duas únicas posteriores ao D28:
  - `projects/Zika_21seq_validacao_20261017`, `run_id fee226ea…`;
  - `projects/Variola_VARV52_reexec_20261016`, `run_id 11938c6a…`.

  As duas rodaram na máquina `geomesh` (48 CPUs lógicas, 50,5 GB, segundo o `manifest.json["environment"]`). As duas foram feitas a partir do commit `6dafd42` com `dirty: true`, e por isso **nenhuma é publicável** (IN-09). Os tempos que aparecem aqui servem para dimensionar o desenho, não para o artigo.
- **Etiquetas de evidência:**
  - **[medido]**: um número saído de um log, de um manifesto ou de um comando rodado nesta meta;
  - **[observado]**: um fato visto num artefato ou na saída de um comando;
  - **[inferido]**: uma conclusão tirada da leitura do código, sem execução;
  - **[não testado]**: uma hipótese.
- **Fonte dos tempos.** Os tempos por etapa saem da diferença entre os carimbos do `log_setup_<data>_<run_id>.log` de cada projeto, ou do próprio log da ferramenta (`Total wall-clock`, `Elapsed time`, `Total time`). A fonte vem indicada em cada linha.

## 1. Visão geral

`workflow.py` é um único processo Python. Ele lê a configuração, abre o manifesto e chama três controladores **em série**:

1. `TreeBuilderController`, que faz alinhamento e inferência;
2. `SubtreeBuilderController`, que extrai as subárvores, coleta os metadados GenBank e roda a mineração FPMax;
3. `MolecularClockController`, a datação.

```
workflow.py ──► TreeBuilderController() ──► SubtreeBuilderController() ──► MolecularClockController() ──► manifest.finish()
                 (alinhar, modelo, 14       (subárvores, NCBI, FPMax,        (LSD2; desligado nos
                  árvores, PNG, heatmap)     metadata, .cql, gráficos)         dois projetos)
```

O DEC-127 registra isso como decisão: "o pipeline é inteiramente sequencial". O `parsl` está em `requirements.txt` e não é chamado em lugar nenhum. `workflow/optimization/parallelWorkflow.py` tem só um esboço, `parallel_tree_construction`, cujo corpo é `NotImplemented` (`parallelWorkflow.py:17-18`) [observado]. O laço alinhadores × métodos do modo `advanced` é sequencial (`TBC:682-714`).

## 2. Etapas, uma a uma (C1)

`TBC` abrevia `workflow/controller/treeBuilderController.py`.

| # | Etapa | Onde (arquivo:linha) | Entradas | Saídas (caminhos reais, Zika-21) | Parâmetros | Ferramenta externa | Semente / threads |
|---|---|---|---|---|---|---|---|
| E0 | Carregar a configuração; trocar `#` pelo nome do projeto | `workflow.py:13-67` | `-p config.json`, ou `-cw` com o JSON inline, mais `-iData` e `-pName` | — | `project_name`, `output_log`, `tree_config`, `subtree_config` | — | — |
| E1 | Criar `out/`, gravar `config_backup.json` | `workflow.py:71-75` | params | `out/outputs/config_backup.json` | — | — | — |
| E2 | Abrir o manifesto: commit, ambiente, versões das ferramentas, SHA-256 da entrada, reprodutibilidade, modo de execução; primeira gravação | `workflow.py:80-130`, `workflow/utils/manifest.py:273-650` | params, `input_path` | `out/outputs/manifest.json` (parcial) | `mode`, `ignore_mode` | 7 ferramentas chamadas com `--version` (`manifest.py:54-106`) | registra `random_seed`, `raxml_threads` e `iqtree_threads` (`builder.py:25-51`) |
| E3 | Log da execução e redirecionamento de stdout | `workflow.py:93-94`, `:132-133`; `workflow/utils/run_logging.py:88-123` | — | `out/outputs/log_setup_2026-10-17_fee226ea3bbf.log`, `out/outputs/output_log.txt` | `log_file: true` | — | — |
| E4 | Iniciar o controlador de árvores: diretórios, lista de FASTA, `AlignmentSeqs`, políticas; `clean_NoPipe(input_path)` | `TBC:101-188` (diretórios em `136-140`, FASTA em `142-148`, `AlignmentSeqs` em `152-159`, limpeza em `188`) | `input_path` | `out/{outputs,tmp,Align,Trees,outputs/Plots}/`, `out/tmp/tmpFilesAlignment/` (`alignmentSeq.py:49-51`) | `num_threads`, `on_method_failure`, `model_selection`, `ignore_mode` | — | — |
| E5 | Validar e deduplicar o FASTA | `TBC:405`, `:507-545`; `workflow/utils/dataValidation.py:58-155` | `data/<conjunto>/dataset_final.fasta` | se houver duplicata: `data/<conjunto>/dataset_final_NoPipe`, **gravado no diretório de entrada** (`dataValidation.py:133`, AC-09). No VARV-52: 54 → 51 registros | — | — | — |
| E6 | Laço por arquivo × alinhador × {NJ, UPGMA}, no modo `advanced` | `TBC:390-466`, `:662-714` | FASTA | — | `aligners` (padrão `('mafft','mafft_iterative')`, `TBC:1224`), `ignore_mode` | — | — |
| E7 | Alinhar, na primeira árvore de cada braço que precisa de alinhamento | `TBC:911-921` (distância), `:978-987` (parcimônia), `:1265-1299` (métodos avançados); `_resolver_alinhador` `:1182-1209`; `_alinhar` `:1244-1263`; `alignmentSeq.py:439-521` | FASTA (ou `_NoPipe`) | `out/tmp/dataset_final_mafft.aln`, `out/tmp/dataset_final_mafft_iterative.aln` | estratégia `retree` (`--retree 2 --maxiterate 0`) ou `iterative` (`--maxiterate 1000`, ou `2` acima de 50 kb, `alignmentSeq.py:408-410`, `:472-483`) | `mafft --thread N …` (`alignmentSeq.py:504`) | `num_threads` (4 no Zika-21, 16 no VARV-52); sem semente |
| E8 | Árvore por distância: p-distância D28, depois NJ ou UPGMA | `TBC:716-779`, `:876-940`; `builder.py:350-405`; `distancia_p.py:191-220` | `.aln` | `out/Trees/tree_dataset_final_<alg>_{nj,upgma}_distance.nexus` | `construct_tree_method` interno (`nj`/`upgma`) | — (Biopython) | determinístico, 1 thread [inferido] |
| E9 | Árvore por parcimônia: Fitch com ambiguidade IUPAC e NNI, partindo da árvore NJ ou UPGMA | `TBC:942-1005`; `builder.py:407-443` | `.aln` | `out/Trees/tree_…_{nj,upgma}_parsimony.nexus` | — | — (Biopython) | 1 thread [inferido] |
| E10 | Seleção de modelo (ModelFinder, BIC), uma vez por alinhamento, sob demanda e guardada em memória | `TBC:1035-1113`; `modelo_substituicao.py:234-245`, `:289-340` | `.aln` | `out/tmp/modelfinder_dataset_final_<alg>/{alinhamento.phylip, modelfinder.*}` | `model_selection` (padrão `C′`), timeout de 3600 s (`modelo_substituicao.py:79`) | `iqtree3 -m MF --mset mrbayes -mrate E,I,G,I+G -nt 1 -seed S -redo` | semente `random_seed`, `-nt 1` |
| E11 | IQ-TREE | `TBC:1115-1127`; `builder.py:445-518` (comando em `474-480`) | `.aln` e o modelo | `out/tmp/iqtree_tree_…_iqtree/*` (9 arquivos), `out/Trees/tree_…_iqtree.nexus` | `-m <BIC> -bb 1000` | `iqtree3` | `-seed 12345 -nt 1` (D21); `iqtree_threads` só é declarado |
| E12 | FastTree | `TBC:1129-1140`; `builder.py:520-564` (comando em `533`) | `.aln` | `out/tmp/fasttree_tree_…_fasttree/{.fasta,.log}`, `out/Trees/tree_…_fasttree.nexus` | `-nt -gtr` (GTR+CAT20) | `FastTree` | sem semente; o log mostra a versão `Double precision`, sem OpenMP [observado] |
| E13 | RAxML-NG | `TBC:1142-1154`; `builder.py:566-664` (comando em `604-611`) | `.aln` e o modelo | `out/tmp/raxml_tree_…_raxml/*` (10 arquivos), `out/Trees/tree_…_raxml.nexus` (lido de `.raxml.support`) | `--all --tree rand{10} --bs-trees 1000` | `raxml-ng` | `--seed 12345 --threads raxml_threads --workers 1` (D17) |
| E14 | MrBayes | `TBC:1156-1168`; `builder.py:777-940` | `.aln` e o modelo | `out/tmp/mrbayes_tree_…/*`, `out/Trees/tree_…_mrbayes.nexus` | `mrbayes_*` (`builder.py:70-79`), timeout de 3600 s | `mb` | `seed`/`swapseed` = `random_seed` (`builder.py:766`); **ignorado** nos dois projetos (`ignore_mode`) |
| E15 | PNG de cada árvore | `TBC:587-590`, `:773-776`, `:837-840`, `:1301-1359` | árvore em memória | `out/outputs/Plots/tree_dataset_final_<alg>_<metodo>.png` | — | matplotlib | — |
| E16 | Heatmap de RF por arquivo; heatmap acumulado | `TBC:449-454`, `:467-471`, `:847-874`; `workflow/utils/metrics.py:13-232` | `multi_trees`, **em memória** | `out/outputs/Plots/tree_dataset_final_heatmap_of_distances.png`. O acumulado não sai (ver §5.6) | — | seaborn | — |
| E17 | Fechar o controlador de árvores: `clean_NoPipe`, copiar `tmp/` para `Align/` | `TBC:474-475`; `dataCleaning.py:84-113` | `out/tmp/*.aln` | `out/Align/dataset_final_<alg>.aln` | — | — | — |
| E18 | Iniciar o controlador de subárvores: lista de árvores **sem ordenar** e cache de acessos em memória | `subtreeBuilderController.py:49-90` (`os.listdir` em `81`, cache em `89`) | `out/Trees/` | `out/Subtrees/` | `subtree_config` | — | — |
| E19 | Baixar o GenBank do NCBI, na primeira árvore em que aparecem acessos novos | `subtree_construction/builder.py:94-105`; `workflow/utils/treeUtils.py:467-530` | rótulos dos terminais | `out/outputs/raw_data_sequences.gb`, **aberto com `"w"`** (`treeUtils.py:486`) | lotes de 10, 3 tentativas | Entrez `efetch` (rede) | — |
| E20 | Extrair as subárvores e o metadado de cada terminal | `subtree_construction/builder.py:67-177` (escrita só se o arquivo não existir, em `134-140`; leitura do `.gb` por terminal em `148`); `treeUtils.py:407-465` | árvore, `.gb` | `out/Subtrees/tree_…_InnerN.nexus` (182 no Zika-21, 492 no VARV-52) | `resume_infos` | — | — |
| E21 | Gravar o `metadata.json` antes da mineração | `subtreeBuilderController.py:112-122` | `raw_data` em memória | `out/outputs/metadata.json`, sobrescrito em E23 | `save_metadata` | — | — |
| E22 | FPMax no modo OFST, com o suporte variando de 0,1 a 1,0; consolidar; gerar o `.cql` | `subtreeBuilderController.py:131-137`; `subtree_mining/miner.py:64-226`; `workflow/utils/neo4jProcessing.py:17-46`, `:183` | `raw_data` em memória | `out/outputs/all_results_fpmax.csv`, `out/outputs/neo4j_commands_tree_dataset.cql` | `mode: OFST`, `support_fpmax: auto` | mlxtend `fpmax` | — |
| E23 | Gravar os metadados minerados | `subtreeBuilderController.py:138-146` | resultado do FPMax | `out/outputs/metadata.csv`, `out/outputs/metadata.json` | — | — | — |
| E24 | Histogramas de suporte | `subtreeBuilderController.py:151-157`; `metrics.py:420-459` | resultado do FPMax | `out/outputs/Plots/support_frequence_{bar_plot,frequency_polygon}.png` | — | matplotlib | — |
| E25 | Datação LSD2, uma por árvore datável (`iqtree`, `raxml`, `fasttree`), mais a escolha da árvore principal | `molecularClockController.py:177-241`, `:309-430`, `:433-495`; `molecular_clock/datacao.py` | `Trees/`, `Align/`, `.gb` | `out/outputs/molecular_clock/{relogio_molecular.json, timetree.nex, <arvore>/…}` | `run_molecular_clock` e `molecular_clock_*` | `iqtree3 --date` (LSD2) | `-seed`, `molecular_clock_threads=1`. **Desligada nos dois projetos**: no manifesto, `ignorado_por_configuracao` |
| E26 | Fechar o manifesto: desfecho, drenar `tool_runs`, SHA-256 das saídas, `finish` | `workflow.py:157-172`; `manifest.py:413-430`, `:513-537`, `:644-650` | registro global `tool_runs` | `out/outputs/manifest.json` (final) | — | — | — |

## 3. Tempos medidos por etapa

### 3.1 Zika-21

20 sequências de cerca de 10,8 kb. Configuração `advanced`, MrBayes ignorado, `on_method_failure: continue`. Tempo total pelo manifesto: 14:14:48,885 → 14:30:20,451 UTC = **931,6 s** [medido].

| Etapa | mafft | mafft_iterative | Fonte |
|---|---:|---:|---|
| Alinhamento | ≤ 6,09 s (inclui a p-distância que falhou) | ≈ 10,06 s (10,17 s menos a NJ, de ~0,11 s) | log, linhas 14→15 e 94→99 |
| NJ / UPGMA por distância | falha em < 0,01 s (`SemSitioComparavel`, D28) | ~0,11 s cada | log |
| Parcimônia NJ / UPGMA | falha em < 0,01 s | 246,35 s / 78,79 s | log, 11:19:48,452→11:23:54,805 e 11:28:51,164→11:30:09,953 |
| ModelFinder | 1,256 s | 1,152 s | `modelfinder.log`, "Total wall-clock" |
| IQ-TREE | 6,216 s | 5,141 s | `.log` do IQ-TREE |
| FastTree | 2,77 s | 2,97 s | `.log` do FastTree, "Total time" |
| RAxML-NG | 271,348 s | 285,481 s | `.raxml.log`, "Elapsed time" |
| PNG por árvore | 0,2–0,43 s | | log |
| Heatmap por arquivo | 1,60 s | | log, 11:30:10,170→11:30:11,770 |
| Download NCBI (20 acessos, 2 lotes) | 4,51 s | | log, linhas 160→165 |
| Subárvores (por árvore, sem o download) | 0,28–0,42 s | | log, linhas 318–1873 |
| FPMax + `.cql` + CSV/JSON + gráficos | 0,73 s | | log, 11:30:19,713→11:30:20,439 |

As quatro árvores básicas do braço `mafft` falham no Zika-21 porque 4 pares do alinhamento progressivo não têm sítio comparável (D28). Com `on_method_failure: continue`, a execução segue e produz **10 árvores de 14 pipelines** [observado em `manifest.json["inference_methods"]`].

### 3.2 VARV-52

54 registros, 51 depois da deduplicação, genomas de cerca de 200 kb. `ignore_mode: [mrbayes, parsimony]`, `num_threads` 16, `raxml_threads` 8. Tempo total: 19:54:53,573 → 22:25:33,003 UTC = **9 039,4 s** [medido].

| Etapa | mafft | mafft_iterative | Fonte |
|---|---:|---:|---|
| Alinhamento | ≈ 602,7 s [inferido: 628,4 s com a NJ, menos 25,7 s da UPGMA medida] | ≈ 2 674,2 s [inferido: 2 699,9 − 25,7] | log, linhas 16→18 e 62→64 |
| UPGMA por distância (com releitura do `.aln`) | 25,69 s | 25,65 s | log |
| ModelFinder | 9,947 s | 9,787 s | `modelfinder.log` |
| IQ-TREE | 46,673 s | 46,747 s | `.log` do IQ-TREE |
| FastTree | 233,79 s | 236,97 s | `.log` do FastTree |
| RAxML-NG | 2 129,041 s | 2 039,650 s | `.raxml.log` |
| Download NCBI (51 acessos, 6 lotes) | 35,05 s | | log, linhas 113→130 |
| Subárvores, 10 árvores | 872,3 s ao todo; 334,5 s na árvore `mafft_iterative_raxml`, que foi processada com o `.gb` completo; de 30 a 48 s nas que vieram depois de o `.gb` ser truncado (§5.4) | | log |
| FPMax + `.cql` + CSV/JSON + gráficos | 20,2 s | | log, 19:25:12,649→19:25:32,883 |

### 3.3 Memória

Nenhum manifesto nem log registra o **pico de RSS** de nenhuma etapa [observado: `grep -i "memory\|RAM"` nos logs das ferramentas]. O que existe:

- o IQ-TREE e o ModelFinder declaram a própria estimativa: "NOTE: 2/14/22 MB RAM is required" e "ModelFinder requires 14 MB RAM" [observado];
- o RAxML-NG 2.0.2 não imprime estimativa no `.raxml.log`;
- o MAFFT não tem medição de memória no `aligners.py`: `bytes_per_unit=None`, na linha 285.

## 4. Paralelismo de hoje

| Onde | Como | Evidência |
|---|---|---|
| Entre etapas, métodos, alinhadores e árvores | **nenhum**: laços Python em série | `TBC:390`, `:682-714`; `subtreeBuilderController.py:106`; DEC-127 |
| Dentro do MAFFT | `--thread num_threads` | `alignmentSeq.py:504` |
| Dentro do RAxML-NG | `--threads raxml_threads --workers 1`: pthreads, uma busca por vez | `builder.py:607`; `.raxml.log`: "1 worker(s) x 4 thread(s)" |
| Dentro do IQ-TREE e do ModelFinder | **desligado de propósito**: `-nt 1` (D21) | `builder.py:479`; `modelo_substituicao.py:244`; log: "Kernel: AVX+FMA - 1 threads (48 CPU cores detected)" |
| FastTree, NJ, UPGMA, parcimônia, subárvores, FPMax | uma thread [inferido: Biopython e Python puros; o FastTree sem OpenMP] | — |
| MrBayes | um processo `mb`; não verificamos se o pacote do bioconda é MPI [não testado] | `builder.py:855-860` |

Etapas que rodam em série **sem precisar**:

- os dois braços de alinhamento;
- os até 7 métodos sobre cada alinhamento: FastTree, distância e parcimônia não dependem do ModelFinder;
- a extração de subárvores de cada árvore;
- o download do NCBI, que só depende da lista de acessos e espera todas as árvores;
- a datação de cada árvore.

No VARV-52, durante os 8 146,6 s da fase de árvores (log, linha 102), no máximo 16 das 48 CPUs lógicas estiveram ocupadas: 16 no MAFFT, 8 no RAxML-NG, 1 nas demais etapas. Na parcimônia do Zika-21, a máquina ficou ~5,4 min com 1 thread ocupada [inferido dos parâmetros e dos tempos acima].

## 5. Acoplamentos que impedem quebrar o workflow

### 5.1 Estado em memória entre etapas

Cada item abaixo vira um arquivo no desenho granular.

- **Seleção de modelo.** Fica em `self._selecoes_modelo`, por alinhamento (`TBC:1066-1072`, `:1112`). IQ-TREE, RAxML-NG e MrBayes recebem a tradução em memória, e não de um arquivo. O relatório `.iqtree` fica em `tmp/`, mas o resultado que se usa não sai da memória.
- **Árvores do heatmap.** `multi_trees` acumula as árvores em memória para o heatmap de RF (`TBC:387`, `:441-454`).
- **Árvores avançadas visitadas duas vezes.** No modo `advanced`, cada método avançado é visitado **duas vezes por braço**: uma no laço da NJ e outra no da UPGMA (`TBC:684-712`). Na segunda visita a árvore é lida do disco (`TBC:809-816`) e **entra de novo** em `multi_trees`, de modo que o heatmap compara cópias da mesma árvore [inferido da leitura; o registro `reaproveitado` dessa segunda visita é descartado porque `tool_runs.registrar_metodo` mantém o primeiro desfecho, `tool_runs.py:198-199`].
- **Cache de acessos baixados.** `downloaded_sequences_cache` é um `set` compartilhado entre as árvores (`subtreeBuilderController.py:89`, `:182`).
- **Dados para a mineração.** `raw_data`, a lista que a mineração recebe, também fica em memória (`subtreeBuilderController.py:109`, `:136`).
- **Registro de proveniência do processo.** O manifesto depende de um **registro global do processo**: `tool_runs._EXECUCOES`, `_METODOS` e `_SELECOES_MODELO` (`tool_runs.py:38-81`), drenados uma vez no fim (`workflow.py:166`). Esse módulo diz de si mesmo: "O registro é **de processo**: uma execução do pipeline é um processo". Uma etapa que roda em outro processo não aparece no manifesto, a menos que use `anexar_execucao_isolada` (`manifest.py:144-229`), o caminho já usado pela datação isolada (DEC-129).
- **Log único por processo.** `run_logging.configurar` e `garantir` abrem um log só por processo (`run_logging.py:88-123`), com o `run_id` no nome (D22).

### 5.2 Caminhos fixos e derivados por texto

- `tmp/`, `Align/`, `Trees/` e `outputs/Plots` são fixos em relação a `output_path` (`TBC:136-140`).
- O builder **deriva** o diretório temporário cortando o caminho da árvore em `/Trees` (`builder.py:452`, `:527`, `:573`).
- O PNG é gravado em `path.replace('Trees', 'outputs/Plots')` (`TBC:589`, `:775`, `:839`). Isso troca **toda** ocorrência de `Trees` no caminho absoluto: um diretório de projeto com `Trees` no nome desviaria os PNG [inferido].
- O nome do projeto gravado no `.cql` vem da **estrutura de diretórios** (`neo4jProcessing.py:6-15`: `basename(dirname(dirname(path)))`). Num *workspace* que não siga `<projeto>/out/outputs`, o `.cql` sai com o nome errado [inferido].
- A datação deduz o `.aln` a partir do nome do arquivo da árvore (`molecularClockController.py:322-325`).

### 5.3 Um processo só escreve o manifesto

`ExecutionManifest` é criado no começo (`workflow.py:81-88`) e fechado no `finally` (`:157-172`). Ele hasheia as saídas percorrendo `out/` e excluindo `tmp/`, só para os sufixos `.nexus`, `.nwk`, `.aln` e `.csv` (`manifest.py:413-430`). No Zika-21 isso dá 196 entradas: 2 alinhamentos, 10 árvores, 182 subárvores e 2 CSV [observado]. **`metadata.json` e `raw_data_sequences.gb` ficam fora do hash**, por decisão (`manifest.py:417-419`).

### 5.4 O GenBank sobrescrito no meio da extração (F2.5, já conhecido)

`download_sequences` abre `raw_data_sequences.gb` com `"w"` (`treeUtils.py:486`). O cache em memória compara rótulos sem normalizar a versão. Quando aparecem acessos "novos" numa árvore posterior, o arquivo é **truncado** para esses acessos.

- **VARV-52.** O `.gb` final tem **3 registros** (`grep -c ^LOCUS` → 3, para 51 folhas) [medido]. O segundo download, de 3 sequências, acontece na 3.ª árvore, às 19:20:14 (log, linha 1398).
- **Zika-21.** Houve um download só, e o `.gb` tem os 20 registros [medido]. O defeito não se manifesta nesse conjunto.

Isso já está na triagem de 2026-09-25 e no plano F2.5 (`docs/automation/21-plano-de-fechamento.md`). Como consequência, o metadado de cada árvore depende da **ordem** em que as árvores são lidas, e essa ordem é a do `os.listdir` (§5.5).

### 5.5 Ordem dependente do sistema de arquivos

`SubtreeBuilderController` lista `out/Trees` com `os.listdir`, **sem `sorted`** (`subtreeBuilderController.py:81`). A primeira árvore processada foi `…_mafft_iterative_upgma_parsimony` no Zika-21 e `…_mafft_iqtree` no VARV-52 [observado nos logs]. Essa ordem define:

- a ordem de `metadata.json`, de `metadata.csv` e do `.cql`;
- com F2.5, **quais árvores recebem metadado**;
- no `.cql`, qual ocorrência de cada `Metadata` prevalece (`DeduplicadorMetadata`: "mantido o primeiro", `neo4jProcessing.py:37-41`).

O `TreeBuilderController`, por outro lado, ordena a entrada (`TBC:143`).

### 5.6 Heatmap acumulado nunca sai, e o RF do heatmap tem defeito (achado lateral)

O acumulado começa como matriz 8×8 (`TBC:386`). `somarMatrizes` devolve `None` quando as dimensões diferem (`TBC:866-867`), e com `distance_matrix=None` a função de plotagem retorna sem desenhar (`metrics.py:213-214`). Com 10 ou mais árvores por arquivo, o heatmap acumulado nunca é gerado [inferido; observado que não há `tree_Acumulate_*` nos dois projetos].

Mais grave: `metrics.rf_distance` (`metrics.py:13-58`) indexa os terminais pela **posição em cada árvore** (`terminals.index(t)`). Duas árvores com a mesma topologia e outra ordem de folhas dão RF > 0 [medido nesta meta]:

```
$ python -c "…; t=Phylo.read('…/Trees/tree_dataset_final_mafft_iqtree.nexus','nexus'); u=deepcopy(t); u.ladderize(); …"
mesma topologia, ordem das folhas trocada por ladderize:
  ordem igual? False
  metrics.rf_distance(t,u) = 30
  metrics.rf_distance(t,t) = 0
```

Isso afeta só o PNG `tree_dataset_final_heatmap_of_distances.png`. O RF da zona sagrada é outro, e é conferido pelo oráculo dendropy (`stability/`, `oraculo_rf_dendropy.py`). **Não está no catálogo D1..D28** que conferimos: vai para a fila de triagem do orquestrador. Esta meta não corrige nada.

### 5.7 Escrita no diretório de entrada (AC-09)

A deduplicação grava `data/<conjunto>/dataset_final_NoPipe` (`TBC:535-536`, `dataValidation.py:133`), e `clean_NoPipe` apaga o arquivo no começo e no fim (`TBC:188`, `:474`). O MAFFT do VARV-52 recebeu `dataset_final_NoPipe` [observado em `tools_invoked.mafft.runs[].command`]. Se duas execuções sobre o mesmo conjunto rodarem em paralelo, uma apaga a entrada da outra [inferido].

### 5.8 Reaproveitamento silencioso de artefatos (AC-04)

Antes de construir, cada pipeline procura a árvore no disco e a reaproveita se ela existir: `TBC:747-754` e `:809-816`, desfecho `reaproveitado`. O mesmo vale para o alinhamento:

- nos métodos de distância e parcimônia, basta existir (`TBC:911-917`, `:978-981`);
- nos avançados, o arquivo precisa existir **e** não estar vazio (`TBC:1279-1291`).

As subárvores só são escritas se ainda não existirem (`subtree_construction/builder.py:135-140`). Reexecutar no mesmo `out/` mistura, portanto, artefatos de execuções diferentes.

### 5.9 Rede no meio do pipeline

A extração de subárvores (E19) faz chamadas Entrez com e-mail fixo `email@email.com` (`treeUtils.py:480`) e espera de 3 s entre lotes. Sem saída para a internet, o pipeline não falha: cada lote é tentado 3 vezes, a falha só vira `log.error("Falha definitiva no lote …")`, e o `.gb` fica sem aqueles registros. Os terminais saem então com `{"error": …}` no lugar do metadado [inferido de `treeUtils.py:518-534` e `:407-424`].

### 5.10 Configuração sem esquema

Os controladores fazem `setattr` de todo `kwarg` (`TBC:117-120`, `subtreeBuilderController.py:60-61`). Uma chave com erro de digitação é ignorada em silêncio (AC-02).

## 6. Saídas de um projeto real

Contagem por padrão em `projects/Zika_21seq_validacao_20261017/out`, com `find out -type f`, nomes normalizados por `sed`. Total de **274 arquivos** em 16 diretórios, 15 961 195 bytes [medido]. O VARV-52 tem os mesmos padrões e **584 arquivos**, 518 639 835 bytes: lá são 492 subárvores, e o `metadata.json` tem 243 843 555 bytes [medido].

| Padrão | Zika-21 | Etapa |
|---|---:|---|
| `out/Align/<aln>.aln` | 2 | E17 |
| `out/Trees/<aln>_{nj,upgma}_{distance,parsimony}.nexus` | 4 | E8, E9 |
| `out/Trees/<aln>_{fasttree,iqtree,raxml}.nexus` | 6 | E11–E13 |
| `out/Subtrees/*_InnerN.nexus` | 182 | E20 |
| `out/tmp/<aln>.aln` | 2 | E7 |
| `out/tmp/modelfinder_<aln>/*` (7 por braço) | 14 | E10 |
| `out/tmp/iqtree_*/*` (9 por árvore) | 18 | E11 |
| `out/tmp/fasttree_*/*` (2 por árvore) | 4 | E12 |
| `out/tmp/raxml_*/*` (10 por árvore) | 20 | E13 |
| `out/outputs/Plots/<árvore>.png` | 10 | E15 |
| `out/outputs/Plots/tree_dataset_final_heatmap_of_distances.png` | 1 | E16 |
| `out/outputs/Plots/support_frequence_*.png` | 2 | E24 |
| `out/outputs/config_backup.json` | 1 | E1 |
| `out/outputs/log_setup_*.log`, `out/outputs/output_log.txt` | 2 | E3 |
| `out/outputs/manifest.json` | 1 | E2, E26 |
| `out/outputs/raw_data_sequences.gb` | 1 | E19 |
| `out/outputs/metadata.json` | 1 | E21, E23 |
| `out/outputs/metadata.csv` | 1 | E23 |
| `out/outputs/all_results_fpmax.csv` | 1 | E22 |
| `out/outputs/neo4j_commands_tree_dataset.cql` | 1 | E22 |
| **Total** | **274** | |

Existe ainda um diretório vazio, `out/tmp/tmpFilesAlignment/` (E4), e, no VARV-52, um arquivo transitório fora de `out/`, `data/<conjunto>/dataset_final_NoPipe` (E5, apagado em E17).

Saídas **condicionais**, que não aparecem nesses dois projetos:

- `out/tmp/mrbayes_*` e `tree_*_mrbayes.nexus` (E14);
- `out/outputs/molecular_clock/**` (E25);
- `resul_of_fpmax_<s>.csv`, com `support_fpmax` fixo (`miner.py:175`);
- `neo4j_commands.cql`, fora do modo OFST (`neo4jProcessing.py:27`);
- `tree_Acumulate_heatmap_of_distances.png`, só com exatamente 8 árvores (§5.6);
- alinhamentos de `clustalo` e `muscle`, se declarados em `aligners`.

O mapeamento de cada saída para uma atividade granular está em [02 §6](02-desenho-granular.md#6-cobertura-das-saídas-c3).
