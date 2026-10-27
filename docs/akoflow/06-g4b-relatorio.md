# 06 · Relatório da G4b (implementação do workflow granular)

Máquina `joh`, 2026-10-26/27. Nada foi commitado. Tudo [medido] salvo indicação. **Nenhuma execução é publicável** (manifesto `dirty: true`, IN-09), e esta máquina está **fora da receita**.

## 1. O que foi entregue

| Caminho (em `BioComp_UFF/`) | O quê |
|---|---|
| `workflow/akoflow.py` | uma CLI por atividade (`preparar`, `baixar-metadados`, `alinhar`, `selecionar-modelo`, `inferir`, `subarvores`, `heatmap-rf`, `minerar`, `gerar-cql`, `graficos-suporte`, `estabilidade`, `consolidar`); chamam as funções do monolito (`_process_*_tree_method`, `SubtreeBuilder`, `SubtreeMiner.process_group`, `parse_tree`, `ExecutionManifest`) |
| `workflow/equivalencia.py` | comparador das classes H1-H5, S1-S5, P do 02 §5.4 |
| `workflow/utils/manifest.py` | `absorver_registros` (aditivo; `drain_tool_runs` passou a delegar a ele) |
| `akoflow/workflow.yaml`, `config/zika21.json`, `data/`, `Dockerfile`, `conferir_ferramentas.py`, `rodar_local.py` | YAML final (39 ativações), imagem, checagem de versões da receita, executor local (serial ou paralelo, `--retomar`) |
| D29 | patch de 4 arquivos + `remine.py` + `test_item_biparticao.py` trazidos do checkout principal |

`git diff --stat` (com os novos em intent-to-add): 17 arquivos, 6 325 inserções, 18 remoções (D29: 4 arquivos, 199 +/18 -, dos quais `manifest.py` 89 inclui a minha refatoração).

## 2. Testes

- `test_item_biparticao`: 26 testes, OK (1 xfail esperado).
- Suíte do submódulo por lista explícita, com `projects/` ligado por link simbólico só nas cópias em /tmp (os testes só leem): **antes 380 testes, depois 401 (+21 novos). As mesmas 4 falhas nas duas** (3 de `TestMrBayesDeVerdade` e `test_recupera_as_datas_de_no`, do LSD2): dependem das ferramentas, e o env desta máquina difere da receita. Sem regressão; não verificado se passam no env pinado.
- Novos: `test_akoflow` (preparar/dedup sem tocar a entrada, baixar, subarvores -> minerar -> gerar-cql == `SubtreeBuilderController` nas mesmas árvores, no-op declarado, falha com fragmento, ordem de visita, `absorver_registros` == `drain`, comparador) e `test_akoflow_yaml` (campos do exemplo real, DAG acíclico, `dataDependencies` ⊂ `dependsOn`, cada comando aceito pelo parser, `cpuLimit` = config).

## 3. Equivalência no Zika-21 (02 §5)

Ambiente: env `Phylotreeminer` desta máquina: MAFFT v7.525, IQ-TREE 3.0.1, RAxML-NG 1.2.2, FastTree 2.2.0, MrBayes 3.2.7, Clustal Omega 1.2.4 (receita: 7.526, 3.1.3, 2.0.2). Mesmo env e mesmo PATH para R e G (E1). R1, R2 = `workflow.py` em diretórios novos (1 315 s e 1 316 s, rodando juntos); G1 = paralelo (12 CPU), G2 = serial, ambos com o `.gb` de R1.

| Classe | R1×R2 (controle nulo) | R1×G1 | R1×G2 |
|---|---|---|---|
| H1 alinhamentos (4) | 0 dif. | 0 | 0 |
| H2 árvores básicas (4) | 0 | 0 | 0 |
| H3 árvores ML + `.support/.bestTree/.treefile/.contree` (16) | 0 | 0 | 0 |
| H4 subárvores (182) | 0 | 0 | 0 |
| H5 `.gb` | 0 | 0 | 0 |
| S1 modelo / S2 FPMax / S3 metadados / S4 `.cql` / S5 manifesto | 0 | 0 | 0 |
| P PNG (13) | 0 | 0 | 0 |

R2×G1 também deu 0. Em S5 entram 194 hashes de `outputs_sha256` e 17 estados de método. **Veredito: dentro da tolerância (zero) em todas as classes**, no E1. Ressalvas:
1. Teste sem o AkôFlow: o executor local compartilha **um** diretório e dispensa `rsync` entre *workspaces* e contêiner. A conferência de conflito de caminho do `rsync` não foi exercida [não testado].
2. Versões fora da receita; a imagem não foi construída [não testado, `docker build` não rodado].
3. Piloto com o `.gb` da `geomesh` (G0): igual em tudo salvo H5 (0 iguais, 1 dif.): o NCBI devolve os registros em **ordem diferente** entre downloads, com o mesmo conteúdo (S3 e S4 iguais). Daí o `.gb` congelado.
4. 4 árvores básicas do braço `mafft` falham por D28, como no monolito; viram `tentado_e_falhou` e *no-op* nas extrações.

## 4. Tempos (s) [medido; 12 CPU, sem overhead de contêiner]

| | monolito R1 | granular serial G2 | granular paralelo G1 |
|---|---:|---:|---:|
| makespan | 1 315 (com R2 ao lado) | 712,4 | **408,9** |
| trabalho total | | 712,4 | 1 328,3 (inflado por contenção) |

Por atividade em G2 (soma das 2 ativações): `nj-parsimony` 203,8; `raxml` 374,5 (máx. 187,4); `upgma-parsimony` 67,6; `alinhar` 18,8; `iqtree` 12,2; `fasttree` 7,6; `modelo` 4,6; `heatmap-rf` 2,8; `minerar` 1,5; `consolidar` 1,4; `gerar-cql` 0,7; `graficos` 1,4; cada `subarvores` <= 0,7; cada `nj/upgma-distance` <= 1,7. Em G1 o caminho crítico foi a parcimônia iterativa do `nj` (387,7 s), com `raxml` em 337 e 331 s. O 02 previa o RAxML-NG no caminho crítico (298,6 s teóricos, `geomesh`); aqui a parcimônia ganhou porque o RAxML-NG 1.2.2 e a CPU diferem [observado]. Speedup medido sobre o monolito de 1 315 s: 3,2x, mas o monolito rodou com R2 ao lado; contra G2 (mesma carga) o paralelo dá 1,74x. Não é número para o artigo.

## 5. Em aberto

- `memoryLimit`, `runtime`, `image`, `namespace`: SUPOSIÇÃO. Memória não medida.
- `estabilidade` implementada como invólucro, **não testada**. `datar-*` não implementado (`consolidar` recusa `run_molecular_clock`).
- Rodar R1/R2/G1/G2 na imagem e no nó da equipe, na receita pinada (E1 de verdade).
- A mensagem do 1º método que depende de uma seleção de modelo que falhou difere da do monolito (o monolito repete a mensagem original só no 1º). Sem efeito nos resultados [inferido; não exercitado].
- Achados fora de escopo, não corrigidos: o erro "acesso não encontrado" no `metadata.json` grava o caminho absoluto do `.gb` (`treeUtils.py`, D15); `process_histogram_frequence` quebra se `outputs/Plots` não existe.
