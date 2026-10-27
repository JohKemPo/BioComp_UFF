# 04 · Handoff para a equipe do AkôFlow

**Para:** Wesley Ferreira e a equipe do AkôFlow (e-Science, IC/UFF)
**De:** equipe do PhyloTreeMiner
**Assunto:** workflow filogenético granular para rodar no AkôFlow e medir para o artigo do FGCS

> **Estado:** implementado (meta G4b), testado localmente sem o AkôFlow e sem construir a imagem; ver [06](06-g4b-relatorio.md). Em `BioComp_UFF/`: `workflow/akoflow.py` (CLIs), `akoflow/workflow.yaml`, `akoflow/Dockerfile`, `workflow/equivalencia.py`. Nenhuma execução do PhyloTreeMiner é publicável ainda: todo manifesto atual tem `dirty: true`. Os tempos abaixo servem para planejar, não para o artigo.

## 1. Em duas frases

O PhyloTreeMiner alinha um conjunto de genomas virais com mais de um alinhador, infere uma árvore por par (alinhamento × método) e minera com FPMax os clados que se repetem entre as árvores. Hoje isso roda num único processo Python, **em série**. Quebramos esse processo em atividades independentes, com contrato de arquivo, para que o AkôFlow rode em paralelo tudo o que não depende de outra coisa.

## 2. As atividades

No Zika-21, que é o conjunto do teste, são **39 ativações**. A lista completa, com os contratos, está em [02 §2](02-desenho-granular.md#2-as-atividades-atômicas). O grafo está em [03](03-dataflow.md).

| Atividade | O que faz | Ativações (Zika-21) | CPU | Tempo de referência* |
|---|---|---:|---:|---:|
| `preparar-entrada` | valida o FASTA, remove sequências idênticas, grava a configuração | 1 | 1 | < 0,1 s |
| `baixar-metadados` | baixa os registros GenBank dos acessos (NCBI), ou copia um `.gb` congelado | 1 | 1 | 4,5 s (online) |
| `alinhar-<alg>` | MAFFT progressivo (`mafft`) e iterativo (`mafft_iterative`) | 2 | 4 | 6 s e 10 s |
| `modelo-<alg>` | escolhe o modelo de substituição (ModelFinder, BIC) | 2 | 1 | 1,2 s |
| `arvore-<alg>-<método>` | uma árvore: NJ e UPGMA por distância, parcimônia, FastTree, IQ-TREE, RAxML-NG | 14 | 1 (RAxML-NG: 4) | 0,1 a 285 s |
| `subarvores-<alg>-<método>` | extrai as subárvores de uma árvore e anexa os metadados GenBank | 14 | 1 | 0,3–0,4 s |
| `heatmap-rf` | figura de distâncias RF entre as árvores | 1 | 1 | 1,6 s |
| `minerar-fpmax` | FPMax sobre os clados de todas as árvores | 1 | 1 | < 0,6 s |
| `gerar-cql` | comandos Cypher para carga no Neo4j | 1 | 1 | — |
| `graficos-suporte` | histogramas de suporte | 1 | 1 | 0,2 s |
| `consolidar-manifesto` | junta a proveniência de todas as atividades num `manifest.json` | 1 | 1 | — |

\* Medido na execução monolítica de 2026-10-17, máquina `geomesh` (AMD Threadripper 2970WX, 48 CPUs lógicas). Num nó de vocês será diferente.

Quatro das 14 árvores **falham de propósito** no Zika-21, porque o alinhamento progressivo deixa pares de sequências sem sítio comparável (defeito D28, tratado). Essas ativações terminam com código 0 e um registro `tentado_e_falhou`, e as subárvores delas viram *no-op*. Isso é esperado e faz parte do teste.

Os números do paralelismo, pelo nosso modelo e sem o *overhead* de vocês:

- **Zika-21:** largura máxima 15; caminho crítico de 298,6 s contra 931,6 s do monolito, ou **3,1× teórico**. O RAxML-NG é 95 % do caminho.
- **VARV-52** (51 genomas de varíola, ~200 kb): 31 ativações; largura 11; caminho crítico de 5 078 s contra 9 039 s, ou **1,8×**. Ali mandam o MAFFT iterativo e o RAxML-NG.

## 3. Construir a imagem

**Uma imagem** para todas as atividades, como no exemplo `infrastructure-selection` (`spec.image`). A receita fica com a G4b. O esboço:

```dockerfile
# Esboço; o arquivo real é BioComp_UFF/akoflow/Dockerfile (contexto de build = raiz do PhyloTreeMiner; não foi construído).
FROM mambaorg/micromamba:1.5.10-jammy
# compilador para a seção pip, como em Backend/Dockerfile do PhyloTreeMiner
USER root
RUN apt-get update && apt-get install -y --no-install-recommends build-essential && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
# receita PINADA da raiz do PhyloTreeMiner (MAFFT 7.526, IQ-TREE 3.1.3, RAxML-NG 2.0.2,
# FastTree 2.2.0, MrBayes 3.2.7, Clustal Omega 1.2.4; Python 3.10)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml requirements.txt requirements-dev.txt /tmp/receita/
RUN micromamba install -y -n base -f /tmp/receita/environment.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER BioComp_UFF/ /app/BioComp_UFF/
ENV PYTHONPATH=/app/BioComp_UFF
# dados do teste de equivalência (pequenos): FASTA, configuração e o .gb congelado
COPY --chown=$MAMBA_USER:$MAMBA_USER akoflow/data/  /app/data/
COPY --chown=$MAMBA_USER:$MAMBA_USER akoflow/config/ /app/config/
# sem ENTRYPOINT: o comando de cada atividade vem do workflow.yaml
```

Pontos que importam para o resultado:

1. **As versões das ferramentas fazem parte do resultado.** Não troquem por `latest` nem por outro canal. A G4b inclui no *build* uma checagem estrita de versões, e a imagem tem de passar nela.
2. **Não usem o `BioComp_UFF/Dockerfile` atual.** Ele tem um `environment.yml` sem pinos e um `ENTRYPOINT` que roda o monolito.
3. **Registrem o *digest* da imagem** e passem para nós. Ele entra no manifesto e no artigo.

## 3.1 Como rodar (G4b)

```bash
# no env do projeto (ou na imagem)
cd BioComp_UFF
python akoflow/conferir_ferramentas.py          # versões contra a receita; sai 1 se divergir
python akoflow/rodar_local.py --workdir /tmp/G1/Zika_G1 --modo paralelo --jobs 12 --tempos /tmp/G1/tempos.json
#   --modo serial | --retomar (pula atividades com prov/ concluído) | --app-dir <dir com config/ e data/>
# uma atividade só, como o AkôFlow a chama (cwd = diretório do projeto):
python -m workflow.akoflow alinhar --config out/outputs/config_backup.json --alinhador mafft --out out
# equivalência (classes H/S/P do 02 §5.4); sai 0 só se H e S forem iguais:
python -m workflow.equivalencia R1/Zika_R1/out G1/Zika_G1/out
```

- `workflow.yaml` usa `/app/config` e `/app/data` (dentro da imagem). O `.gb` da imagem é o da `geomesh`; para o teste de equivalência troque-o pelo `raw_data_sequences.gb` de R1 (o NCBI muda a ordem dos registros entre downloads).
- Cada atividade grava `out/outputs/prov/<atividade>.json` e `stdout/<atividade>.txt`; `consolidar` junta tudo em `manifest.json` (+ campo aditivo `execucao_granular`). Ele precisa ver as saídas de **todas** as atividades (pergunta Q1).
- Resultado local no Zika-21 (executor local, 12 CPU): zero diferença em todas as classes contra o monolito; ver [06](06-g4b-relatorio.md).
- Dados e config: `akoflow/data/Zika479_Test_large/` (FASTA e `.gb`) e `akoflow/config/zika21.json`.

## 4. Dados de entrada

| Conjunto | Arquivo | Tamanho | Registros | SHA-256 |
|---|---|---:|---:|---|
| **Zika-21** (teste) | `BioComp_UFF/data/Zika479_Test_large/dataset_final.fasta` | 55 820 B | 20 sequências, ~10,8 kb | `656aa379314fb6aabe09875b8c5d0a0711ef5de67d2e0f39d77447ed66901a64` |
| VARV-52 (escala) | `BioComp_UFF/data/workflow_dataAcquisition_li_et_al_2007_replication-RetMax200-ITRs-clean/dataset_final.fasta` | 10 316 575 B | 54 registros (51 após deduplicação), ~200 kb | `8435653302b86fd8b7a31dc3ac3ffe4b905efd8fabcd68208a7d4efbf8a077fc` |

Os hashes foram medidos com `sha256sum` e conferem com o `inputs_sha256` dos manifestos de referência.

Volume gerado por projeto, medido no monolito (`du -sb out/`):

- **Zika-21:** 16 MB em 274 arquivos;
- **VARV-52:** 519 MB em 584 arquivos, quase tudo o `metadata.json`, de 244 MB.

Num projeto de varíola com os metadados completos, o `metadata.json` chega a **~1 GB**: medido no VARV-49 de 2026-09-01, que é uma execução anterior à correção D28. Com *staging* por `rsync` entre *workspaces*, esse arquivo atravessa as atividades finais. Ver Q1.

**Rede.** `baixar-metadados` online chama o NCBI Entrez. Para o teste de equivalência e para o padrão *offline* dos *showcases* (`SHOWCASE_STANDARD.md`, item 2), usamos `--gb-congelado`: o `.gb` da execução de referência vai dentro da imagem e nenhuma atividade sai para a internet.

## 5. O teste de equivalência

**Objetivo.** Mostrar que o workflow granular no AkôFlow dá **o mesmo resultado** que o monolito `workflow.py`. Fazemos isso antes de qualquer medição de tempo, porque speedup de um resultado diferente não vale nada. O plano completo e as tolerâncias, fixadas antes de rodar, estão em [02 §5](02-desenho-granular.md#5-equivalência-com-o-monolito).

O que precisamos que vocês rodem, todos **no mesmo tipo de nó** e com a **mesma imagem**:

1. **R1 e R2.** O monolito, duas vezes, dentro da imagem e em diretórios novos:
   ```
   python /app/BioComp_UFF/workflow.py -p /app/config/zika21.json
   ```
   Cerca de 16 min cada, pela medição no `geomesh`. Pode ser uma atividade única de AkôFlow ou um `docker run` no mesmo tipo de nó.
2. **G1 e G2.** O `akoflow/workflow.yaml` (a versão final da G4b, derivada do esboço [`workflow-zika21.esboco.yaml`](workflow-zika21.esboco.yaml)), duas vezes, com o `.gb` congelado de R1.
3. **Comparar.** Rodar o comparador da G4b, `python -m workflow.equivalencia R1/out G1/out`, o mesmo para R1 × R2 e para R1 × G2, e nos mandar a saída inteira e os dois `out/`.

**Critério.** Zero diferença em alinhamentos, árvores e subárvores (hash SHA-256), na escolha de modelo, no resultado do FPMax, nos metadados e no `.cql` (comparados depois de ordenar). As figuras PNG e os logs ficam fora da comparação. Se R1 × R2 já mostrar diferença numa classe, essa classe sai do critério, e nós registramos o achado.

**Por que o mesmo tipo de nó.** IQ-TREE e RAxML-NG escolhem as instruções vetoriais (AVX2, AVX-512) conforme a CPU. Duas árvores quase empatadas em verossimilhança podem trocar de lado com outro caminho de ponto flutuante. Isso nunca foi testado entre CPUs diferentes neste projeto. Um experimento separado com nós heterogêneos é bem-vindo, mas como diagnóstico.

## 6. O que esperamos medir para o FGCS

Proposta, para vocês ajustarem ao que o AkôFlow já coleta:

| Medida | Como | Observação |
|---|---|---|
| Tempo por ativação | *task attempt* do AkôFlow (execução, fila, *startup*, *staging*), mais os carimbos dos nossos fragmentos `prov/*.json` | A diferença entre os dois é o *overhead* por ativação |
| *Makespan* e speedup | granular no AkôFlow ÷ monolito no **mesmo** tipo de nó | Mostrar contra o teórico (3,1× no Zika, 1,8× no VARV) |
| Eficiência e cenários de alocação | como no artigo de 2024 (cenários C1–C4): por exemplo, 1 vCPU por ativação de uma thread, 4 vCPU no alinhamento e no RAxML-NG, e variantes | `cpuLimit` precisa ser **maior ou igual** às threads da ferramenta (Q9) |
| Custo | `pricePerSecond` × recursos × tempo, pelo modelo de custo do *environment* de vocês; separar processamento e armazenamento, como na Tabela 2 de 2024 | O Zika é pequeno; o custo interessa no VARV |
| Dados movidos | *transfers* do AkôFlow (bytes, arquivos) | O `metadata.json` do VARV pesa |
| Pico de memória por atividade | telemetria cgroup v2 | **Nós não temos essa medida.** Ela vai para os `memoryLimit`, hoje SUPOSIÇÃO |
| Proveniência | Para cada arquivo de `outputs_sha256` do nosso manifesto, o SHA-256 do *Generated files* do AkôFlow tem de ser igual; linhagem exportada em JSON | É a segunda prova de integridade, independente da nossa |
| Repetições | 5 execuções por cenário, como em 2024; mediana e dispersão | Pré-registrar os cenários antes de rodar |

Um plano de escalonamento (HEFT ou PRISM) usaria as nossas durações medidas como perfil de recurso. Pela especificação, `simulation.durationSeconds` torna a atividade **só de simulação**, então não pode ir no YAML real. Ver Q8.

## 7. Perguntas em aberto para a equipe

| # | Pergunta | Por que importa |
|---|---|---|
| Q1 | O *staging* copia o *workspace* da predecessora **com o que ela recebeu das predecessoras dela**, ou só o que ela produziu? | O `consolidar-manifesto` precisa ver todas as saídas. Se a cópia não for transitiva, ele tem de depender diretamente das 39 atividades |
| Q2 | O contêiner recebe, por variável de ambiente, o id da *run*, o nome da atividade e o número da tentativa? | Para ligar o nosso `manifest.json` à proveniência do AkôFlow sem casar por horário |
| Q3 | Código de saída ≠ 0 numa atividade: a *run* inteira falha, ou só os descendentes ficam parados? Existe "permitir falha"? | Hoje, a política `on_method_failure: continue` sai com 0 e um estado declarado. Queremos saber se há uma forma nativa |
| Q4 | Há saída para a internet nos nós? Como se passa um segredo (a chave `NCBI_API_KEY`) sem pô-lo no YAML? `credentialRef`? | Para o `baixar-metadados` online (fora do teste de equivalência) |
| Q5 | Dá para fixar o **tipo de nó** de atividades específicas (`resourceSelector`? perfis no `environment.yaml`?) | Para o determinismo do IQ-TREE e do RAxML-NG em ambiente heterogêneo |
| Q6 | O que recomendam, imagem única ou por atividade? Há limite de tamanho, ou *cache* de imagem por nó? | A imagem com conda e as 7 ferramentas deve ter alguns GB [não medido] |
| Q7 | Onde ficam as saídas finais e por quanto tempo? Devemos usar *promote data* para o registro científico? | Para o artigo apontar onde está cada artefato |
| Q8 | Como passar as nossas durações medidas para o planejador numa execução real (`activityResourceProfiles`?) | `simulation` desliga a execução real |
| Q9 | O `cpuLimit` é cota rígida de cgroup? Com MAFFT `--thread 16` e `cpuLimit: "4"`, o que acontece? | Superalocação mexe no tempo, e o tempo é a medida do artigo |
| Q10 | Uma atividade de 0,1 s gera métricas, com amostra inicial e final? | Metade das nossas ativações dura menos de 1 s |
| Q11 | Qual versão do AkôFlow citar (*release* ou commit)? Pretendem publicar o `environment.yaml`, o `scope.yaml` e o `topology.yaml` usados? | Reprodutibilidade do lado de vocês |
| Q12 | Querem incluir o VARV-52 na primeira rodada, ou só o Zika-21 até a equivalência fechar? | O VARV leva cerca de 2 h 30 min no monolito |

## 8. Contato e próximos passos

1. O PhyloTreeMiner aprova este desenho (G4a) e implementa a G4b: CLIs, YAML final, imagem e comparador.
2. Entregamos a vocês a imagem (ou a receita), o YAML, o `environment.yaml` sugerido e o `.gb` congelado.
3. Vocês rodam R1, R2, G1 e G2 (§5) e nos devolvem a saída do comparador.
4. Com a equivalência aceita, vêm as medições do §6.
