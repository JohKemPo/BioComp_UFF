# 07 - Documentação do AkôFlow: Tutorials e How-to guides

Consulta para o usuário e para o orquestrador do PhyloTreeMiner. Escopo: as seções "Tutorials" e "How-to guides" da barra lateral de https://www.akoflow.com/docs/.

## 0. Como esta leitura foi feita (leia antes de confiar)

- O `WebFetch` das páginas do site devolveu **só resumos parafraseados** e se recusou a copiar blocos literais. Por isso, **o texto-fonte foi lido no repositório** `UFFeScience/akoflow`, pasta `docs/` (clone raso, commit `aeea81c`, o mesmo que o `workflow.yaml` do nosso worktree cita). A estrutura da barra lateral vem de `docs/sidebars.ts`.
- Os comandos abaixo foram copiados desse fonte. O site publicado pode estar à frente de `aeea81c`; não comparei com o site.
- Níveis de leitura usados abaixo: **[L]** lida inteira; **[P]** lida em parte (só o início, cerca de 2,8 KB); **[R]** só resumo do `WebFetch`; **[N]** não lida (só o título na barra lateral).
- Nada aqui foi executado. Nenhum comando foi testado.

## 1. Índice das páginas

URL base: `https://www.akoflow.com/docs`. O caminho vem do `id`/arquivo no fonte (`docs/docs/<caminho>.md`); o slug no site pode diferir. Só `getting-started` tem slug declarado (`/getting-started`).

### Tutorials

| Página | URL | Leitura |
|---|---|---|
| Getting started | /getting-started | L (fonte) e R |
| Installation | /installation | L |
| Downloads | /downloads | L |
| Interface tour | /guides/interface-tour | L |
| First workflow > Run through the API (SimGrid) | /guides/workflows/first-run | L |
| Showcases > Overview | /showcase | L (resumo da página) |
| Showcases > Scientific Workflows > Montage 050 locally | /showcase/scientific-workflows/montage-050 | L (a parte textual; os componentes React não foram lidos) |
| Showcases > Machine Learning > Infrastructure selection | /showcase/ai/machine-learning/infrastructure-selection | P (só o cabeçalho) |
| Showcases: as outras 26 páginas (ai-assisted-analysis, surrogate-model, experiment-feedback-loop; training-with-mlflow, hyperparameter-search, batch-inference, model-comparison, anomaly-detection, computer-vision, model-ensemble, model-distillation, inference-benchmark, lora-fine-tuning, model-conversion; rag-indexing, synthetic-dataset, audio-transcription, document-extraction, ai-video-studio, local-rag-factory; document-intelligence; video-intelligence; blender-render-farm; bounded-agent, supervised-tool-loop, multi-agent-review) | /showcase/... | **N** |

### How-to guides

| Subseção | Página | URL | Leitura |
|---|---|---|---|
| Workflows | Define a workflow | /guides/workflows/definitions | L |
| Workflows | Plan a workflow | /guides/workflows/planning | L |
| Workflows | Execute and monitor a workflow | /guides/workflows/executions | L |
| Infrastructure | Environments | /guides/infrastructure/environments | L |
| Infrastructure | Model a SimGrid environment | /guides/infrastructure/simgrid | L |
| Infrastructure | Connect a Kubernetes environment | /guides/infrastructure/kubernetes | L |
| Infrastructure | Define execution scopes and network links | /guides/infrastructure/execution-scopes | L |
| Infrastructure | Browse and manage storage | /guides/infrastructure/storage | L |
| Infrastructure | Cloud capacity and machine configuration | /guides/infrastructure/cloud-capacity | P |
| Infrastructure | Configure Google Cloud | /guides/infrastructure/gcp | P |
| Infrastructure | AWS and S3 support | /guides/infrastructure/aws | L |
| Infrastructure | Connect an HPC and SLURM cluster | /guides/infrastructure/hpc-slurm | L |
| Data and evidence | Artifacts, storage, and builds | /guides/data/artifacts | L |
| Data and evidence | Provenance and audit | /guides/data/provenance-and-audit | L |
| Operations | Instance management | /guides/operations/instance-management | L |
| Operations | Run the AkôFlow server on a Linux instance | /guides/operations/server-instance | L |
| Operations | Search and notifications | /guides/operations/search-and-notifications | P |
| Operations | Interactive console and commands | /guides/operations/interactive-console | P |
| Operations | Credentials and SSH service keys | /guides/operations/credentials-and-ssh | P |
| Operations | Troubleshooting | /guides/operations/troubleshooting | L |

### Páginas fora do escopo, lidas só porque as de dentro mandam usá-las

| Página | Arquivo no fonte | Por que foi lida |
|---|---|---|
| Set up API access for the tutorials | `docs/docs/tutorials/api-access.md` (URL: /tutorials/api-access) | linkada por definitions, planning, executions, execution-scopes. Não está na barra lateral. |
| Portable workflow specification | `docs/docs/internal/workflow-spec.md` (URL: /internal/workflow-spec) | é a referência dos campos do `workflow.yaml`; fica em "Developing AkôFlow". |
| `docker-compose.yml` e `.env.example` da raiz do repositório | raiz de `UFFeScience/akoflow` | descrevem a pilha de desenvolvimento a que a first-run se refere. |

Existem no fonte, **não lidas** e fora da barra lateral: `tutorials/register-hpc.md`, `tutorials/connect-cloud.md`, `tutorials/compare-plans.md`, `guides/workflows/first-local-run.md`, `guides/data/{artifact-locations,audit-events,build-executable,provenance}.md`, `guides/infrastructure/{cloud-support,machine-configurations}.md`, `guides/operations/{find-records,follow-notifications,personal-preferences}.md`. Também **não lidas**: as seções Explanations, Reference (incl. `reference/environment-yaml`, `reference/api-overview`) e API endpoints.

## 2. Resumo por página

Convenção da API: base `http://127.0.0.1:8080/akoflow-api` (com `/akoflow-api`, sem barra final), cabeçalho `Authorization: Bearer $AKOFLOW_API_TOKEN`. `202 Accepted` significa "enfileirado", não "concluído": sempre faça polling.

### 2.1 Getting started [L]

- Ensina: o mapa. AkôFlow define um workflow (atividades e dependências de dados), decide onde cada atividade roda, compara planos, executa o plano escolhido e guarda resultados e proveniência.
- O que dá para fazer hoje: simular localmente com SimGrid; rodar em Kubernetes ou HPC/SLURM "após configurar o ambiente"; integração com nuvem com limites (a descoberta e o provisionamento de EC2 na AWS **não** estão implementados).
- Caminho do primeiro uso: (1) instalar o AkôFlow Desktop e confirmar que conecta; (2) fazer o tour da interface; (3) "se você tem a pilha de desenvolvimento e credenciais de API", rodar o exemplo SimGrid que verifica 3 atividades e 2 transferências.
- Armadilha explícita da página: o tutorial SimGrid **usa arquivos do repositório e a API; não é um tutorial só de Desktop**. Um passo a passo completo só com o Desktop "ainda precisa ser validado e documentado".
- Conceito: um **plano** descreve a execução proposta; uma **execução (run)** registra a observada. O Desktop e a API operam sobre os mesmos registros persistidos. A API é preferível para experimentos repetíveis.

### 2.2 Installation [L]

- Instala o **AkôFlow Desktop**, que carrega o runtime no Docker local no primeiro uso.
- Pré-requisitos:

| Plataforma | Requisito |
|---|---|
| macOS | Docker Desktop e o instalador Desktop da mesma Release |
| Windows | Docker Desktop com contêineres Linux e o instalador Desktop |
| Linux | Docker Engine, plugin Docker Compose v2 e o instalador Desktop |

- Passos: abrir https://github.com/UFFeScience/akoflow/releases/latest; baixar o instalador do seu SO e arquitetura; instalar e abrir; esperar a tela **Overview** e o indicador de conexão com o daemon.
- O primeiro uso baixa os arquivos de runtime da mesma Release e os carrega no Docker. Não é preciso copiar token do daemon para o Desktop.
- Armadilhas: o Docker precisa estar rodando antes de abrir o AkôFlow; não misturar ativos de Releases diferentes; atualizar/reverter sempre por Release nomeada, exportando a instância antes.
- Para servidor numa instância Linux, ver 2.20.

### 2.3 Downloads [L]

| Plataforma | Ativo | Ação |
|---|---|---|
| macOS (Intel e Apple silicon) | `.dmg` universal | abrir e arrastar para Applications |
| Windows x64 | `.exe` | seguir as notas da Release |
| Debian/Ubuntu x64 | `.deb` | `sudo apt install ./downloaded-file.deb` (trocar o nome do arquivo) |
| Linux x64 | `.AppImage` | dar permissão de execução e abrir |

- Ignorar: `.blockmap` e `latest*.yml` (metadados do atualizador); ZIP/TAR de código-fonte; os `.tar` de daemon/BuildKit e os `.sha256` (são do guia de servidor).
- Download automatizado pela API pública do GitHub (precisa de `curl` e `jq`):

```bash
AKOFLOW_RELEASE_API=https://api.github.com/repos/UFFeScience/akoflow/releases/latest
curl --fail-with-body --fail --location --silent --show-error \
  "$AKOFLOW_RELEASE_API" \
  -o release.json

# Select the exact asset name shown on the release page.
AKOFLOW_ASSET='<desktop-asset-filename>'

AKOFLOW_DOWNLOAD_URL=$(jq -er --arg name "$AKOFLOW_ASSET" \
  '.assets[] | select(.name == $name) | .browser_download_url' release.json) || exit 1
curl --fail-with-body --fail --location --show-error \
  "$AKOFLOW_DOWNLOAD_URL" -o "$AKOFLOW_ASSET"
```

- Verificação no Linux (esperar `OK`):

```bash
AKOFLOW_DIGEST=$(jq -er --arg name "$AKOFLOW_ASSET" \
  '.assets[] | select(.name == $name) | .digest | select(startswith("sha256:"))' \
  release.json) || exit 1
printf '%s  %s\n' "${AKOFLOW_DIGEST#sha256:}" "$AKOFLOW_ASSET" | sha256sum --check -
```

- Armadilha: a disponibilidade de ativos e o formato do instalador mudam entre Releases; ler as notas, sobretudo para Windows e ARM64. A página diz que a validação de download é separada da validação de instalação no SO.

### 2.4 Interface tour [L]

- Mapa Desktop para API:

| Área | Para que serve | Domínio da API |
|---|---|---|
| Overview | ver ambientes, workflows e runs concluídos | vários endpoints de leitura |
| Workflows | definições versionadas e DAG | `/workflow-definitions/` |
| Planning sessions | comparar algoritmos e promover candidato a plano | `/planning-sessions/` |
| Infrastructure | ambientes, recursos, escopos | `/environments/`, `/resources/`, `/execution-scopes/` |
| Runs | iniciar e inspecionar execuções | `/execution-runs/` |
| Artifacts | registrar, construir, localizar, materializar executáveis | `/artifacts/`, `/artifact-builds/` |
| Provenance | evidência e linhagem | `/provenance/` |
| Audit | histórico de eventos operacionais | `/audit-events/` |

- Hierarquia de propriedade: Workflow -> Plan -> Execution -> Activity; Environment -> Inventory -> Resource -> Provisioning -> Operation.
- Busca global: `Ctrl+K` (ou `⌘K`). Terminal global fixo no rodapé (sessões pertencem ao daemon). Algumas instâncias são somente leitura (bloqueiam criação, edição, execução e ações destrutivas).

### 2.5 First workflow: "Your first workflow with SimGrid" [L]

- Ensina: registrar o exemplo versionado, executá-lo e verificar computação e rede. **Nada vai para Kubernetes, SLURM ou nuvem.** Não ensina escalonamento automático (importa um plano fixo).
- Pré-requisitos: clone do repositório do AkôFlow; daemon rodando com o runner SimGrid; `curl` e `jq`; o token da API, se a autenticação estiver ligada. A pilha de desenvolvimento escuta a porta 8080.

```bash
export AKOFLOW_API_URL="http://127.0.0.1:8080/akoflow-api"
export AKOFLOW_API_TOKEN="<token>"
```

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/preflight/" | jq
```

Só continuar com `server.available` igual a `true`. O runner SimGrid precisa estar no contêiner do daemon ou em `AKOFLOW_SIMGRID_BINARY`.

- Ordem fixa dos seis `POST`, cada um referenciando IDs do anterior. Cabeçalho de todos: `Content-Type: application/yaml`, `--data-binary @arquivo`.

| # | Endpoint | Arquivo (em `examples/simulation/`) |
|---|---|---|
| 1 | `POST /environments/` | `environment.yaml` |
| 2 | `POST /execution-scopes/` | `scope.yaml` |
| 3 | `POST /network-topologies/` | `topology.yaml` |
| 4 | `POST /workflow-definitions/` | `workflow.yaml` |
| 5 | `POST /schedule-plans/` | `plan-request.yaml` |
| 6 | `POST /execution-runs/` | `execution-request.yaml` |

Exemplo do passo 1 e sua verificação:

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/yaml' \
  --data-binary @examples/simulation/environment.yaml \
  "$AKOFLOW_API_URL/environments/"
```

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/environments/simulation-example/" |
  jq '{environment: .environment.id,
       runtime: .runtimes[0].driver,
       resources: [.resources[].id]}'
```

IDs esperados: `simulation-example`, `simgrid`, `simulated-edge`, `simulated-cloud`.

- Polling do run (criação é assíncrona):

```bash
while :; do
  projection=$(curl --fail-with-body --silent --show-error \
    -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
    "$AKOFLOW_API_URL/execution-runs/simulation-example-run-v1/") || exit 1
  run_state=$(printf '%s' "$projection" | jq -r '.run.status')
  printf '%s\n' "$projection" | jq \
    '{status: .run.status,
      completed: .run.completedActivityCount,
      activities: .run.activityCount}'
  case "$run_state" in completed|failed|cancelled) break ;; esac
  sleep 1
done
```

- Critério de aprovação (invariantes): `status` = `completed`; 3 de 3 atividades; 2 registros de transferência; `transferredBytes` = `120000000`; `breakdown.computeSeconds` e `breakdown.transferSeconds` maiores que zero. Referência de 2026-09-11: makespan de cerca de 21,593 s.
- Atalho de um comando, e artefatos do run em `storage/simgrid/<run-id>-<instance>/` (`platform.xml`, `simulation.json`, `result.json`, `runner.log`):

```bash
AKOFLOW_API_URL="$AKOFLOW_API_URL" \
AKOFLOW_API_TOKEN="$AKOFLOW_API_TOKEN" \
sh examples/simulation/run.sh
```

- Armadilhas: IDs fixos (`simulation-example`, `simulation-example-run-v1`); repetir só parte da sequência dá `422` ou erro de chave estrangeira, então use instância nova ou troque os IDs em todos os seis arquivos; `401` = faltou token; `FOREIGN KEY constraint failed` no plano = ID de atividade do plano não bate com o workflow registrado; `akoflow-simgrid-runner: executable file not found` = falta o runner; run em `pending` = ver eventos e log do daemon; zero bytes transferidos = faltam `dataDependencies` ou produtor e consumidor ficaram no mesmo recurso. Unidades: banda em **bits/s**, tamanhos em **bytes**.
- Equivalente no Desktop (passos 1-6 da página): Environments, Execution scopes, Workflows > Definitions > **Import YAML**, aba **Plans**, iniciar o plano, abrir o run. "A run is complete only when the header says `completed` and the activity summary says `3/3 settled`."
- **Inconsistência da doc**: este tutorial manda usar `examples/simulation/...`; as páginas simgrid/executions usam `akoflow-examples/simulation/...`. No clone `aeea81c` não existe pasta `examples/` na raiz; os arquivos estão no repositório separado `UFFeScience/akoflow-examples` (a página kubernetes e o showcase apontam para ele).

### 2.6 Showcases [L só o overview e a Montage; o resto N]

- 28 workflows executáveis. Cada página explica DAG, atividades, entradas, saídas e proveniência, **sem código**; as imagens Docker, programas Python e validadores ficam em `UFFeScience/akoflow-examples`. Os 27 exemplos pequenos rodam CPU real em perfis reduzidos; a Montage 050 (58 atividades, 138 transferências de arquivo nomeadas) é o fluxo científico completo e longo.
- Comando local de cada showcase (gerado pelo componente `WorkflowLocalRun.tsx`; precisa de Go e, para a maioria, Python 3.12+; Docker para Montage e Blender):

```bash
git clone https://github.com/UFFeScience/akoflow-examples.git &&
cd akoflow-examples &&
cd <examplePath> &&
./prepare.sh &&
./run.sh
```

(há também um `git checkout --detach <sourceCommit>` entre o `cd akoflow-examples` e o `cd <examplePath>`, quando o contrato do exemplo define `sourceCommit`). Para o Desktop, a partir do diretório do exemplo, com o daemon local e `AKOFLOW_API_TOKEN` definidos:

```bash
go run ../../tools/daemon-local -root ../.. -example <examplePath>
```

(com `-wait=false` para os `scientific-workflows/*`). Segundo o componente, esse comando valida a conexão com o daemon, constrói a imagem declarada, registra o ambiente local e o workflow e submete um run real.
- Montage 050: dois "slots de CPU" locais, sem GCP/Kubernetes/HPC; saída em `outputs/manifest.json`, `outputs/dag-validation.json`, `outputs/stage-<activity>.json`; `./validate.sh` confere as 58 atividades. Os perfis de duração herdados foram medidos em outra CPU. Não foi lido o artigo de blog sobre a Montage.
- Infrastructure selection (única página de ML relevante para nós, pois nosso `workflow.yaml` cita esse exemplo): só li o cabeçalho (título, resumo, estágios e saídas `decision.json` e `cost-comparison.csv`). Não li o YAML do exemplo.
- Aviso da página: execução em HPC, Kubernetes ou nuvem "ainda exige o ambiente descrito no guia de infraestrutura correspondente".

### 2.7 Define a workflow [L]

- Ensina: criar a definição (Desktop) ou importar YAML portátil (API). A definição é **versionada e imutável**; planos e runs apontam para a versão usada.
- Desktop: Workflows > Workflow definitions > **Create workflow** > nome e namespace > atividades (artefato executável, comando, limites de CPU e memória, predecessoras separadas por vírgula) > prévia do DAG > criar. Há **Import YAML**, exportar e duplicar. A exportação remove IDs gerados e caminhos de runtime resolvidos.
- API: exige o `AKOFLOW_API_URL`/`TOKEN` do 2.23. Para execução real, `command.entrypoint` e `command.executable` são obrigatórios. `spec.image`, `image` por atividade e `run` são formas legadas "ainda aceitas"; as novas definições devem usar `command` e `command.executable`. YAML mínimo da página (só simulação, `command: {}`):

```yaml
name: portable-simulation-demo
spec:
  namespace: examples
  activities:
    - name: prepare
      command: {}
      cpuLimit: "1"
      memoryLimit: 128Mi
      simulation:
        model: deterministic
        durationSeconds: 4
    - name: analyze
      command: {}
      cpuLimit: "2"
      memoryLimit: 1Gi
      dependsOn: [prepare]
      simulation:
        model: deterministic
        durationSeconds: 12
```

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/yaml' \
  --data-binary @workflow.yaml \
  "$AKOFLOW_API_URL/workflow-definitions/"

# List definitions
curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/workflow-definitions/"

# Read one definition
curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/workflow-definitions/portable-simulation-demo/"

# Export portable YAML
curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -o exported-workflow.yaml \
  "$AKOFLOW_API_URL/workflow-definitions/portable-simulation-demo/export/"
```

- A resposta da criação traz `version.id`, que se usa para abrir a sessão de planejamento. Armadilha: o exemplo SimGrid versionado usa a forma legada.
- Os campos exatos estão na especificação portátil (2.24).

### 2.8 Plan a workflow [L]

- Ensina: plano manual ou candidatos automáticos (PRISM Time, PRISM Cost, HEFT) comparados antes de executar.
- Desktop: abrir a definição > **Generate plan** > **Automatic planning** > escolher escopo (execução ou simulação) e topologia > algoritmos (PRISM: contagem de opções e beam width) > opcionais deadline, budget, matriz de interferência > **Generate candidate plans** > inspecionar Gantt > selecionar. PRISM Cost e PRISM Time são objetivos exclusivos; HEFT é baseline e **não usa** a matriz de interferência. Os modelos de predição de HEFT e PRISM são diferentes: compare execuções observadas para saber qual plano rende mais.
- Plano manual: escolher escopo e topologia, atribuir **cada atividade** a um runtime e recurso compatíveis e informar a duração esperada.
- API:

```bash
curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/planning-algorithms/"
```

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/json' \
  -d '{
    "id": "planning-simulation-example",
    "workflowVersionId": "simulation-example-workflow-v1",
    "executionScopeId": "simulation-example-v1-scope",
    "networkTopologyId": "simulation-network-v1",
    "algorithms": [
      {"id": "prism-time", "configuration": {"optionCount": 25, "beamWidth": 120}},
      {"id": "prism-cost", "configuration": {"optionCount": 25, "beamWidth": 120}},
      {"id": "heft", "configuration": {}}
    ],
    "deadlineSeconds": 0,
    "budget": 0
  }' \
  "$AKOFLOW_API_URL/planning-sessions/"
```

Retorna `202`. Depois: `GET /planning-sessions/<id>/` até `completed`, `GET .../candidates/`, e selecionar:

```bash
curl --fail-with-body -X POST \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/planning-sessions/planning-simulation-example/candidates/$AKOFLOW_CANDIDATE_ID/select/"
```

`201` com o plano salvo. Plano manual: `POST /schedule-plans/` com o envelope de validação (`plan-request.yaml`); plano já montado: `POST /schedule-plans/import/` com `{ "plan": ... }` (a origem passa a `imported`).
- Armadilhas: `feasible` indica que o candidato passou na validação do plano, **não** que os vínculos de runtime estão prontos; a validação de plano importado **não verifica vínculos de runtime** (confirme um binding habilitado e compatível com o modo para cada recurso atribuído antes de executar); `POST /schedule-plans/` e `/import/` salvam as métricas previstas que você manda, sem recalcular (já a página first-run diz que o servidor reavalia planos importados com o modelo atual; as duas frases coexistem na doc); uma sessão congela o escopo, e mudar o escopo depois exige nova sessão.

### 2.9 Execute and monitor a workflow [L]

- Ensina: iniciar um run a partir de um plano salvo e inspecionar observações.
- Modo **derivado do escopo do plano**: escopo de simulação gera run de simulação; escopo de execução gera run real. "Real" envia atividades para a máquina local, Kubernetes ou SLURM.
- Estados do run: `created` -> `running` -> `completed`/`failed` (cancelamento só antes da conclusão). Atividades: `blocked` -> `ready` -> `preparing` -> `running` -> `completed`/`failed`.
- Em runs reais há três instantes: submitted (entregue ao runtime), started (runtime alocou) e container-started (código do usuário pode começar). Fila, transferência, interferência e overhead aparecem na decomposição; "makespan" é relógio de parede e o "tempo de estágio acumulado" pode passar dele quando há sobreposição.
- API:

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/yaml' \
  --data-binary @akoflow-examples/simulation/execution-request.yaml \
  "$AKOFLOW_API_URL/execution-runs/"

curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/execution-runs/simulation-example-run-v1/"

curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/execution-runs/"
```

- Resposta de detalhe: `run`, `activities`, `dataTransfers`, `handles`, `events`.
- Armadilhas: um `404` logo após submeter pode ser só a fila ainda não processada; se o run nunca aparecer, a validação do worker falhou antes de salvar (ver log do servidor). Para falhas: começar por `run.failureReason`, depois atividade, `handle.failure`, código de saída, log e eventos em ordem; se a falha foi na preparação, olhar preparação de executável/workspace e transferências (pode nem haver log de runtime).

### 2.10 Environments [L]

- Um ambiente descreve onde o AkôFlow planeja ou roda trabalho. **Real** (runtime local, SSH, Kubernetes, SLURM ou nuvem) ou **simulação** (SimGrid). Definições são versionadas; escopos e planos apontam para a **versão** do ambiente.
- Ciclo de vida (`environment.status`): `defined` -> `connecting` -> `connected` -> `discovering` -> `ready` (-> `degraded`/`unreachable`). Planejamento e execução esperam `ready`.
- Corpo de criação (é um `EnvironmentDefinition`, não só um nome). Exemplo local mínimo da página:

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/json' \
  -X POST "$AKOFLOW_API_URL/environments/" \
  -d '{
    "environment":{"id":"local-lab","name":"Local lab","status":"defined"},
    "version":{"id":"local-lab-v1","environmentId":"local-lab","version":1,"status":"published"},
    "runtimes":[{"environmentVersionId":"local-lab-v1","id":"local-lab-local","name":"Local execution","driver":"local","mode":"execution","capabilities":{"container":true}}],
    "resources":[{"id":"local-lab-machine","environmentVersionId":"local-lab-v1","executionTarget":"direct","type":"local_machine","name":"Local machine","cpuCores":4,"memoryBytes":8589934592}],
    "resourceRuntimeBindings":[{"resourceId":"local-lab-machine","runtimeId":"local-lab-local"}]
  }'
```

- Teste de conexão (não cria ambiente), saúde e descoberta (operações diferentes: saúde verifica acesso, descoberta coleta inventário):

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/json' \
  -X POST "$AKOFLOW_API_URL/connection-tests/" \
  -d '{"id":"hpc-ssh","name":"HPC login","type":"ssh","endpoint":"login.example.org:22","username":"researcher","credentialRef":"<ssh-key-reference>"}'

curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -X POST "$AKOFLOW_API_URL/environment-connections/hpc-ssh/health/"
curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -X POST "$AKOFLOW_API_URL/environment-connections/hpc-ssh/discover/"
```

- Armadilha: segredos ficam no daemon; o ambiente guarda só `credentialRef`. Um ambiente de simulação não valida endpoint SSH/Kubernetes.

### 2.11 Model a SimGrid environment [L]

- Ensina: modelar recursos, duração e rede para simulação. Não serve para validar conexão SSH, Kubernetes, nuvem ou SLURM.
- Campos de recurso (`environment.yaml`): `cpuCores` (núcleos do host SimGrid), `cpuCapacity` (manter coerente com `cpuCores`), `computeSpeedup` (duração = `baseRuntimeSeconds / computeSpeedup`), `pricePerSecond`, `bootOverheadSeconds`, `containerOverheadSeconds`, `schedulable`.

```yaml
resources:
  - id: simulated-edge
    cpuCores: 2
    cpuCapacity: 2
    computeSpeedup: 1
    pricePerSecond: 0
    containerOverheadSeconds: 0.1
    schedulable: true
  - id: simulated-cloud
    cpuCores: 8
    cpuCapacity: 8
    computeSpeedup: 4
    pricePerSecond: 0.001
    bootOverheadSeconds: 0.5
    containerOverheadSeconds: 0.2
    schedulable: true
```

- Atividade: `runtime: simgrid` e `simulation: {model: deterministic, durationSeconds: N}` (ou `flops`). Dê um perfil **por atividade**; um padrão único esconde diferenças.
- Rede (`topology.yaml`): `links` com `sourceResourceId`, `targetResourceId`, `bandwidthBitsPerSecond`, `latencySeconds`, `pricePerByte`, `bidirectional`, `sharingPolicy` (`shared`, `independent`, `fatpipe`). Transferência = `latencySeconds + sizeBytes / (bandwidthBitsPerSecond / 8)`. Transferência só existe se produtor e consumidor estão em recursos diferentes. PRISM rejeita rota ausente; o baseline HEFT pode estimar transferência zero sem link.
- O formulário de escopo do Desktop cria topologia **vazia** e não expõe criação de links: use a API.
- Bundle completo: `sh akoflow-examples/simulation/run.sh` depois de `git clone --depth 1 https://github.com/UFFeScience/akoflow.git akoflow-simgrid` (ver a inconsistência de caminhos em 2.5).
- Tabela de erros comuns da página: todas as atividades com duração minúscula igual (faltou perfil por atividade); nenhuma transferência (`sizeBytes` > 0, recursos diferentes, rota na topologia); transferência rápida demais (unidades); atividades paralelas em uma raia (`cpuCores`, `coreId`); recurso ausente nos candidatos (`schedulable: true` e versão do ambiente no escopo).

### 2.12 Connect a Kubernetes environment [L]

- Ensina: ligar um cluster existente para o AkôFlow rodar atividades de contêiner como **Jobs**. Recomenda começar por um cluster Kind descartável.
- Pré-requisitos: API do Kubernetes alcançável **a partir da máquina do daemon**; `kubectl` de administrador só para criar namespace e conta de serviço; namespace dedicado; imagens que todos os nós consigam baixar (o runtime passa a imagem direto ao Kubernetes e **não configura credenciais de registro nem `imagePullSecrets`**); storage class ou armazenamento compartilhado se houver troca de arquivos.
- Passo 1: RBAC mínimo (Namespace, ServiceAccount, Role, RoleBinding). A Role dá `create/get/list/delete` em `jobs` (grupo `batch`) e em `pods`, `services`, `persistentvolumeclaims`; `get` em `pods/log`; `create` em `pods/exec`. Para descoberta de nós, um `ClusterRole` à parte com `get,list` em `nodes`.

```bash
kubectl apply -f akoflow-runtime-rbac.yaml
kubectl auth can-i --as=system:serviceaccount:akoflow:akoflow-runtime create jobs -n akoflow
kubectl auth can-i --as=system:serviceaccount:akoflow:akoflow-runtime get pods/log -n akoflow
```

- Passo 2: guardar o token fora da definição do ambiente:

```bash
set -o pipefail
kubectl -n akoflow create token akoflow-runtime --duration=1h \
  | jq -R '{id:"research-kubernetes",token:.}' \
  | curl --fail-with-body \
      -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
      -H 'Content-Type: application/json' --data-binary @- \
      "$AKOFLOW_API_URL/kubernetes-tokens/" -o kubernetes-reference.json || exit 1
jq -e '.credentialRef' kubernetes-reference.json
```

Use o `credentialRef` retornado. O cliente aceita `file:<path>`, `env:<variable>` ou `bearerToken` na configuração. Em produção, `configuration.caFile`; `insecureSkipTlsVerify: true` só para Kind.
- Passo 3: conexão, runtime e recursos (trecho do `akoflow-examples/kind/environment.yaml`):

```yaml
connections:
  - id: kind-akoflow-connection
    type: kubernetes
    endpoint: https://host.docker.internal:63395
    credentialRef: file:storage/credentials/kubernetes/kind-akoflow.token
    configuration:
      namespace: akoflow
      insecureSkipTlsVerify: true
runtimes:
  - id: kind-kubernetes
    driver: kubernetes
    mode: execution
    configuration:
      namespace: akoflow
      connectionId: kind-akoflow-connection
    capabilities:
      batch: true
      container: true
      cancellation: true
resources:
  - id: kind-worker
    type: kubernetes_machine
    providerId: akoflow-worker
    cpuCores: 4
    cpuCapacity: 4
    schedulable: true
    metadata:
      kubernetesNode: akoflow-worker
resourceRuntimeBindings:
  - resourceId: kind-worker
    runtimeId: kind-kubernetes
    enabled: true
```

O `namespace` do runtime prevalece sobre o da conexão. `providerId`/`metadata.kubernetesNode` deve ser o valor exato de `kubernetes.io/hostname` (vira `nodeSelector` do Pod).
- Passo 4: armazenamento (PVC):

```yaml
storages:
  - id: kind-akoflow-data
    type: pvc
    endpoint: pvc://akoflow/akoflow-data
    shared: true
    readOnly: false
    configuration:
      namespace: akoflow
      claimName: akoflow-data
      mountPath: /akoflow/data
    runtimeBindings:
      - runtimeId: kind-kubernetes
        default: true
        containerPath: /akoflow/data
```

PVC `ReadWriteOnce` monta em um nó só; planos que fixam atividades em nós diferentes precisam de outro modo de acesso ou co-localização.
- Passo 5: validar no Desktop (testar conexão, descobrir, importar workflow pequeno, escopo, plano, **Real execution**) ou seguir o README do Kind em `UFFeScience/akoflow-examples/tree/main/kind`. Limpeza: `kubectl -n akoflow get jobs,pods,services,pvc -l app.kubernetes.io/managed-by=akoflow` e, no Kind descartável, `kind delete cluster --name akoflow`.
- Armadilhas: o console interativo exige `kubectl` no host do daemon; Job `Pending` = ver `kubectl describe pod`; `403` = faltou permissão específica.

### 2.13 Define execution scopes and network links [L]

- Um escopo diz ao planejador quais versões de ambiente ele pode usar e é um "universo congelado" por sessão. Com uma máquina só, topologia vazia basta.
- Desktop: Infrastructure > Execution scopes > **Create scope** (ID, nome, versões publicadas de ambiente, opção de criar topologia vazia).
- API:

```bash
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/json' \
  -X POST "$AKOFLOW_API_URL/execution-scopes/" \
  -d '{
    "id":"hybrid-research",
    "name":"Hybrid research",
    "networkTopologyId":"hybrid-network-v1",
    "environmentVersionIds":["hpc-v1","cloud-v1"]
  }'

curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H 'Content-Type: application/json' \
  -X POST "$AKOFLOW_API_URL/network-topologies/" \
  -d '{
    "id":"hybrid-network-v1",
    "name":"HPC to cloud",
    "version":1,
    "executionScopeId":"hybrid-research",
    "links":[{
      "id":"hpc-cloud",
      "sourceResourceId":"hpc-cluster",
      "targetResourceId":"cloud-capacity-small",
      "bandwidthBitsPerSecond":1000000000,
      "latencySeconds":0.03,
      "pricePerByte":0,
      "bidirectional":true
    }]
  }'
```

- Armadilhas: a API **não valida** que os recursos do link pertencem ao escopo nem mede a rede real; os valores de topologia só afetam estimativas.

### 2.14 Browse and manage storage [L]

- Navegação limitada a **raízes aprovadas** (descoberta ou conectores configurados); ações dependem de capacidade (somente leitura, não saudável ou sem acesso verificado a partir dos nós de computação desabilitam escrita).
- Endpoints: `GET /environments/<env>/storages/`, `GET /storages/<id>/roots/`, `GET /storages/<id>/entries/?path=...&limit=...` (paginação por `nextCursor`), `POST /storages/<id>/downloads/` (ou `/archives/` para diretório), `/copies/`, `/checksum/`, `GET /storage-downloads/<id>/` e `.../content/`. Promoção: `POST /storages/<id>/promote-data/` (campos `path`, `id`, `workflowVersionId`, `runId`, `activityId`) e `/promote-artifact/` (`path`, `id`, `name`, `version`, `scope`, `scopeId`). Promover **registra** o caminho, não envia nem move o arquivo.
- Armadilha: não construa caminhos fora das raízes devolvidas.

### 2.15 Cloud capacity [P], Google Cloud [P], AWS and S3 [L]

- Cloud capacity: separa catálogo do provedor, alvos de capacidade, configurações de máquina (playbooks Ansible) e instâncias provisionadas. Matriz v1.0: GCP descobre e provisiona (Terraform); **AWS não descobre nem provisiona EC2**; transferência por GCS ou S3 e endpoints compatíveis. Li só o início (catálogo e alvo de capacidade).
- GCP: conta de serviço com JSON, Compute Engine API (e Billing API para preços); o token pede o escopo `cloud-platform`; o Terraform cria uma instância, disco de boot e regra de firewall SSH. Li só o início.
- AWS: S3 por `s3://bucket/prefix` com `AWS_ACCESS_KEY_ID` e `AWS_SECRET_ACCESS_KEY` (e `AWS_SESSION_TOKEN` se temporárias) **no processo do servidor**. Salvar credencial AWS no AkôFlow não a disponibiliza para S3. O navegador de storage do Desktop envia requisições S3 sem assinatura. "Um run contra conta AWS ou bucket de teste ainda não foi verificado." Para computar em máquinas AWS existentes, conecte-as por SSH ou Kubernetes.

### 2.16 Connect an HPC and SLURM cluster [L]

- Conecta ao **nó de login** por SSH e submete atividades via `sbatch`; os testes locais não emulam o SLURM e não validam SSH, alocação, contabilização nem política do site. Não trate o nó de login como recurso de computação.
- Pré-requisitos: conta com SSH não interativo e `sbatch`, `squeue`, `sacct`, `scancel`; chave pública autorizada em cada salto; fingerprint do host, porta, partição, conta/QoS, caminho de workspace visível pelos nós de computação; runtime de contêiner aprovado (o exemplo usa Apptainer e não o instala); alocação de teste aprovada.
- Conexão (a porta vai em `configuration.port`, `endpoint` é só o host):

```yaml
connections:
  - id: research-hpc-ssh
    environmentId: research-hpc
    name: Research HPC login node
    type: ssh
    endpoint: login.example.org
    username: researcher
    credentialRef: file:storage/credentials/ssh/research-hpc
    configuration:
      port: 22
      hostKeyAlias: research-hpc-login
      knownHostsFile: storage/credentials/ssh/known_hosts
      proxyCommand: ssh gateway.example.org -W login.example.org:22
      scriptDirectory: /scratch/researcher/akoflow/scripts
```

- Runtime e recursos (`examples/slurm/environment.yaml`): runtime `driver: slurm`, `mode: execution`, `configuration.partition`, `containerRuntime: apptainer`; recursos `cluster` (não escalonável), `hpc_partition` (`executionTarget: batch`, `providerId` = nome da partição) e `hpc_machine` (nó de login, `executionTarget: direct`). Storages `lustre` e `nfs` com `runtimeBindings` (`hostPath`, `containerPath`).
- Validação em ordem: saúde da conexão, descoberta (usa `sinfo`), `sbatch` de sonda, escrita/leitura no armazenamento a partir de um nó alocado, job Apptainer mínimo, um run de uma atividade.
- Armadilhas: testar SSH **a partir do host do daemon**; job que some do `squeue` não é falha automática (o AkôFlow usa `sacct`, com `squeue`/`scontrol` como reserva); cancelar chama `scancel`; fechar a sessão do console para liberar a alocação `srun`.

### 2.17 Artifacts, storage, and builds [L]

- Separa **dados científicos** de **artefatos executáveis** (imutáveis, versionados). Arquivos produzidos podem ser promovidos ao registro científico; o **manifesto de artefatos** é o que torna um run citável.
- Construir executável a partir de imagem Docker (gera SIF; o pull e a conversão rodam no serviço de build):

```bash
REGISTERED=$(curl -sS -X POST \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "artifactId":"busybox",
    "version":"1.36",
    "image":"docker.io/library/busybox:1.36",
    "architecture":"amd64"
  }' \
  "$AKOFLOW_API_URL/artifacts/docker/")
# Read .build.id from REGISTERED, then start it:
curl --fail-with-body -X POST -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/artifact-builds/$BUILD_ID/runs/"
```

- Contexto de build próprio: `POST /build-contexts/` multipart (`-F "context=@context.tar.gz"`), depois `/artifact-builds/` com `id`, `artifactVersionId`, `contextDigest`, `recipeDigest`, `cacheKey`. Materialização: `planned`, `reconciling`, `transferring`, `verifying`, `committed`, `failed`; só é `committed` quando `verifiedDigest` = `digest`.
- Armadilha: nunca coloque credenciais de registro, SSH ou nuvem nos payloads; caminhos locais do navegador não são contexto de build do servidor.

### 2.18 Provenance and audit [L]

- **Proveniência** responde "como este resultado foi produzido" (workflows, planos, runs, atividades, transferências, dados); **auditoria** responde "que operação ocorreu, quando, em qual alvo, com que resultado". São armazenadas separadamente.
- Consulta de proveniência: `GET /provenance/entities/` (catálogo; não invente nomes de entidade ou campo), `GET /provenance/entities/runs/?q=...&filterField=...&filterValue=...&page=...&pageSize=...&sortField=...&sortOrder=...`, linhagem `GET /provenance/lineage/runs/$RUN_ID/?direction=both&depth=2&maxNodes=300` (resposta com `root`, `nodes`, `edges`, `truncated`; teto de 300 nós).
- SQL somente leitura (`SELECT`/`WITH`, timeout de 10 s, 200 linhas por página, esquema seguro):

```bash
curl --fail-with-body -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/provenance/sql/schema/"
curl --fail-with-body -X POST -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "sql":"SELECT id, status, created_at FROM execution_runs WHERE status = :status ORDER BY created_at DESC",
    "parameters":{"status":"completed"},
    "page":1,
    "pageSize":200
  }' \
  "$AKOFLOW_API_URL/provenance/sql/"
```

- Auditoria: `GET /audit-events/` com filtros `eventType`, `environmentId`, `resourceId`, `connectionId`, `sessionId`, `executionId`, `outcome` (`started`, `succeeded`, `failed`), `limit`.
- Armadilha: favoritos e histórico de SQL do Desktop ficam no perfil do navegador; o texto de resumo de auditoria pode conter erros de sistemas externos (revisar antes de publicar).

### 2.19 Instance management [L]

- Uma **instância** é uma instalação do plano de controle com seu catálogo (SQLite). Operações: `GET/PUT /instance/` (`id` e `name` obrigatórios; `transferBufferBytes` entre 5 e 64 MiB, padrão 8 MiB), exportar ZIP saneado (`GET /instances/default/export/?includeArtifacts=false`; credenciais **não** vão), importar como snapshot somente leitura (`POST /instances/import/`, até 8 GiB compactado, 10 000 entradas, 64 GiB expandido), ativar (`POST /instance-activations/<id>/`) e reset de fábrica (`POST /factory-reset/`, irreversível).
- Em snapshot, toda requisição não-`GET` retorna `423 Locked` com `{"error":"the selected instance is a read-only snapshot"}`; a única exceção é voltar para `default`.

### 2.20 Run the AkôFlow server on a Linux instance [L]

Passo a passo completo está na seção 3.2 (opção B). Pontos-chave: Linux amd64/arm64 com Docker Engine e Compose v2; usa os arquivos `akoflow-daemon-<tag>-linux-<arch>.tar`, `akoflow-buildkit-<tag>-linux-<arch>.tar` e `akoflow-runtime-<tag>-linux-<arch>.sha256` de uma Release; sem registro de contêiner; o Compose monta o `/var/run/docker.sock` do host e roda o BuildKit privilegiado (**não usar em host não confiável ou multiusuário**); API presa a `127.0.0.1`, acesso remoto por túnel SSH ou proxy reverso; CORS só com origens HTTPS exatas em `AKOFLOW_API_ALLOWED_ORIGINS` (sem curinga).

### 2.21 Search and notifications [P], Interactive console [P], Credentials and SSH keys [P]

- Busca: `GET /search/?q=...&types=workflow,execution,artifact&limit=20` (limite padrão 30, máximo 100); tipos: `workflow`, `execution`, `artifact`, `environment`, `resource`, `plan`, `scope`, `materialization`. Notificações são do Desktop, sem endpoint próprio. Li só o início.
- Console: comando único (`POST /console-sessions/` com `{"resourceId":...,"actorId":...}` cria uma sessão interativa, WebSocket); `503` se o motor foi iniciado sem suporte a terminal. Li só o início.
- Credenciais: o segredo é enviado uma vez, guardado com modo `0600`, e só `credentialRef` circula. Chave de serviço Ed25519 gerada por `POST /ssh-keys/` com `{"id":"...","comment":"..."}`; resposta com `id`, `credentialRef`, `publicKey`, `fingerprint`. Li só o início (importação de chave e atribuição não lidas).

### 2.22 Troubleshooting [L]

- Árvore de decisão por fronteira: motor/pré-requisitos, autenticação, modo somente leitura, conexão/saúde/descoberta, busca, terminal, armazenamento, planejamento/execução, importação/troca de instância, coleta de evidência.
- Diagnóstico base: `curl --fail-with-body "${AKOFLOW_API_URL%/akoflow-api}/"` (saúde do raiz) e `GET /preflight/` (daemon, Docker do host e BuildKit). `401` = token ausente ou inválido; `403` pode ser fora do limite de loopback; `423` = snapshot somente leitura (`POST /instance-activations/default/`); `503 interactive console is unavailable`.
- Fluxo de falha de run: primeira atividade falha ou bloqueada, motivo, recurso/runtime resolvidos e logs, comparar previsto e observado, materialização de artefatos, correlacionar IDs em Audit e Provenance. Reset de fábrica é último recurso: exporte antes. Redigir tokens, chaves e segredos antes de compartilhar evidência.

### 2.23 API access for the tutorials (fora do escopo, lida) [L]

```bash
set -o pipefail
export AKOFLOW_API_URL='http://127.0.0.1:8080/akoflow-api'
read -rsp 'Akoflow API token: ' AKOFLOW_API_TOKEN; printf '\n'
export AKOFLOW_API_TOKEN
```

```bash
curl --fail-with-body "$AKOFLOW_API_URL/preflight/" | jq
curl --fail-with-body \
  -H "Authorization: Bearer $AKOFLOW_API_TOKEN" \
  "$AKOFLOW_API_URL/environments/" | jq
```

O primeiro deve trazer `server.available: true`; o segundo, um catálogo (vazio é válido). Para o Desktop empacotado, a conexão é gerenciada por um proxy dele; a doc diz que os tutoriais gráficos não exigem extrair a credencial interna.

### 2.24 Portable workflow specification (fora do escopo, lida) [L]

Campos do `workflow.yaml` (fonte em `docs/docs/internal/workflow-spec.md`; um arquivo mais antigo e **desatualizado** com o mesmo nome existe em `docs/internal/workflow-spec.md`, descreve `spec.storageClassName`/`mountPath` e afirma que atividades rodam em sequência; ignorei-o).

| Campo | Tipo | Obrigatório | Significado |
|---|---|---:|---|
| `name` | string | sim | nome de exibição; origem do ID estável (minúsculas com hífen) |
| `spec.namespace` | string | sim | namespace lógico |
| `spec.image` | string | não | imagem OCI padrão das atividades reais sem executável |
| `spec.activities` | array | sim | atividades |
| `spec.dataDependencies` | array | não | arestas de dados e tamanhos lógicos |
| `spec.storageClassName`, `storageSize`, `storagePolicy`, `mountPath` | | não | **aceitos mas não persistidos nem aplicados**; configure armazenamento no ambiente e no runtime |

Campos de atividade: `name` (obrigatório, único), `run` (atalho: vira `sh -c <run>`), `command` (objeto: `executable`, `entrypoint`, `arguments`, `environment`, `workingDirectory`), `image`, `runtime` (**dica retida como metadado**), `cpuLimit` (padrão `0.1`; aceita `0.5` ou `250m`), `memoryLimit` (padrão `16Mi`; `Ki`, `Mi`, `Gi` ou bytes), `dependsOn`, `resourceSelector` (dica), `keepDisk` (dica), `mountPath` (dica), `simulation` (presença torna a atividade **somente simulação**), `expectedOutputs` (afirma que os caminhos relativos existem no workspace após a atividade; ausência faz o run falhar).

Fontes de executável (`command.executable.source.type`): `catalog`, `oci`, `local-container-image`, `build`, `local-file`, `remote-file`, `object-storage`, `http`; estratégias de entrega: `auto`, `managed`, `use-in-place`, `destination-pull`, `gateway-transfer`, `build-and-transfer`, `prefer-in-place`.

Dados: `dataDependencies` com `producerActivity`, `consumerActivity`, `logicalName` (rótulo de planejamento, não precisa ser caminho), `sizeBytes` (positivo); a tupla (produtor, consumidor, nome) é única. **`dataDependencies` não cria aresta de execução**: precisa do `dependsOn` correspondente, senão o consumidor não espera e não há transferência. Antes de iniciar um sucessor, o AkôFlow faz `rsync` do **workspace completo** de cada predecessor, ignora arquivos idênticos e **rejeita conteúdo conflitante no mesmo caminho relativo** em vez de sobrescrever. Falha de transferência deixa o sucessor sem iniciar. Num run com falha, instâncias de nuvem com saídas efêmeras são mantidas e podem gerar cobrança.

## 3. Como instalar e rodar localmente o AkôFlow para testar um workflow

Montado só com o que as páginas lidas dizem.

### 3.1 Resposta direta

- **Roda localmente com Docker**, não com kind/minikube. O Desktop "carrega o runtime no Docker local" e o servidor é um par de contêineres Docker (daemon e BuildKit). Os requisitos são Docker Desktop (macOS, Windows com contêineres Linux) ou Docker Engine com Compose v2 (Linux).
- **Kind** só aparece como **alvo** opcional (o exemplo `kind` do `akoflow-examples`, para testar o runtime Kubernetes). **Minikube não é citado em nenhuma página lida.** Para o nosso teste local, Kubernetes não é necessário: o runtime `local` (`driver: local`, `type: local_machine`) executa contêineres no Docker da máquina.
- Serviços que sobem:
  - Pela Release (opção A, Desktop, ou opção B, servidor): daemon (`akoflow-server`, API na 8080, `/akoflow-api`, SQLite) e BuildKit (`buildkitd`). Volumes nomeados guardam SQLite, credenciais gerenciadas, workspaces de simulação, artefatos e cache do BuildKit. O Compose monta `/var/run/docker.sock`.
  - O Desktop é a interface por cima do daemon.
  - O SimGrid runner vem na imagem do servidor.
  - Pilha de desenvolvimento (opção C): o `docker-compose.yml` da raiz do repositório sobe `akoflow-server` (compilado com `go run ./cmd/server`) e `buildkitd`, porta 8080, token padrão `akoflow-development-token` se `AKOFLOW_API_TOKEN` não estiver definido. Isso veio do arquivo, não de uma página de guia.

### 3.2 Passo a passo

**Opção A: Desktop (estação de trabalho).**
1. Instalar e iniciar o Docker (Docker Desktop, ou Docker Engine + plugin Compose v2).
2. Abrir https://github.com/UFFeScience/akoflow/releases/latest e baixar o ativo do seu SO (Linux x64: `.deb` ou `.AppImage`; ver 2.3). Opcional: verificar o digest SHA-256 como em 2.3.
3. Instalar: `sudo apt install ./arquivo.deb` (Debian/Ubuntu) ou `chmod +x` no `.AppImage`.
4. Abrir o Desktop e esperar a tela **Overview** com o indicador "daemon conectado". O primeiro uso baixa os `.tar` de runtime da mesma Release e os carrega no Docker.
5. Se não chegar ao Overview, abrir os detalhes da falha e consultar o Troubleshooting (2.22).
6. Para usar a API (scripts), é preciso o URL e o token. A doc diz que o Desktop empacotado os gerencia via proxy e **não explica onde ler o token**; ver seção 5.

**Opção B: servidor numa instância Linux (API sem Desktop).** Comandos exatos da página:

```bash
export AKOFLOW_RELEASE_TAG="<release-tag>" # replace with the tag on Releases
case "$(uname -m)" in
  x86_64) export AKOFLOW_ARCH="amd64" ;;
  aarch64|arm64) export AKOFLOW_ARCH="arm64" ;;
  *) echo "Unsupported architecture: $(uname -m)" >&2; exit 1 ;;
esac
```

```bash
mkdir -p ~/akoflow-server/releases
cd ~/akoflow-server/releases
export AKOFLOW_RELEASE_URL="https://github.com/UFFeScience/akoflow/releases/download/${AKOFLOW_RELEASE_TAG}"
curl --fail-with-body --location --remote-name \
  "${AKOFLOW_RELEASE_URL}/akoflow-daemon-${AKOFLOW_RELEASE_TAG}-linux-${AKOFLOW_ARCH}.tar"
curl --fail-with-body --location --remote-name \
  "${AKOFLOW_RELEASE_URL}/akoflow-buildkit-${AKOFLOW_RELEASE_TAG}-linux-${AKOFLOW_ARCH}.tar"
curl --fail-with-body --location --remote-name \
  "${AKOFLOW_RELEASE_URL}/akoflow-runtime-${AKOFLOW_RELEASE_TAG}-linux-${AKOFLOW_ARCH}.sha256"
sha256sum --check "akoflow-runtime-${AKOFLOW_RELEASE_TAG}-linux-${AKOFLOW_ARCH}.sha256"
```

```bash
docker image load --input "akoflow-daemon-${AKOFLOW_RELEASE_TAG}-linux-${AKOFLOW_ARCH}.tar"
docker image load --input "akoflow-buildkit-${AKOFLOW_RELEASE_TAG}-linux-${AKOFLOW_ARCH}.tar"
docker image inspect \
  "akoflow/daemon:${AKOFLOW_RELEASE_TAG}" \
  "akoflow/buildkit:${AKOFLOW_RELEASE_TAG}" \
  --format '{{.RepoTags}}'
```

```bash
cd ~/akoflow-server
curl --fail-with-body --location --remote-name \
  "https://raw.githubusercontent.com/UFFeScience/akoflow-examples/main/server-instance/compose.yaml"

umask 077
read -r -s -p "AkôFlow API token: " AKOFLOW_API_TOKEN
printf '\n'
printf 'AKOFLOW_RELEASE_TAG=%s\nAKOFLOW_API_TOKEN=%s\nAKOFLOW_PORT=8080\n' \
  "$AKOFLOW_RELEASE_TAG" "$AKOFLOW_API_TOKEN" > .env
unset AKOFLOW_API_TOKEN
```

```bash
docker compose -f compose.yaml up -d
docker compose -f compose.yaml ps
docker compose -f compose.yaml logs --tail=100 akoflow-server

set -o pipefail
curl --fail-with-body --silent http://127.0.0.1:8080/akoflow-api/preflight/ | jq .
```

Esperar `server.available: true`. Provar o token com `GET /environments/` (lista vazia é válida; `401` = token diferente do `.env`). Parar sem perder estado: `docker compose -f compose.yaml down` (**nunca** `--volumes`, que apaga banco, artefatos, credenciais e cache). Atualizar: baixar/verificar/carregar a nova Release, mudar só `AKOFLOW_RELEASE_TAG` no `.env` e repetir `up -d`.

**Opção C: pilha de desenvolvimento a partir do código-fonte** (a first-run a pressupõe). A doc não traz um guia para isso; o que há é o `docker-compose.yml` do repositório. Clonar `https://github.com/UFFeScience/akoflow` e subir esse Compose (não testei; o comando não está em página de guia).

**Rodar um workflow de teste (a partir de qualquer opção, via API):**
1. Definir `AKOFLOW_API_URL` e `AKOFLOW_API_TOKEN` (2.23) e conferir `GET /preflight/`.
2. Conferir o SimGrid com o tutorial 2.5 (clonar `UFFeScience/akoflow-examples` para os YAML de `simulation/`, ver a inconsistência de caminho), até `status: completed`, 3/3 atividades, `transferredBytes` = `120000000`.
3. Para um workflow real local: registrar um ambiente com runtime `local` (corpo de 2.10, trocando IDs e recursos), criar o escopo (2.13), importar o workflow (2.7), gerar ou criar o plano manual (2.8), iniciar o run com escopo de execução (2.9) e acompanhar. Esse caminho **não está demonstrado de ponta a ponta** em nenhuma página lida; o showcase usa `tools/daemon-local` (2.6) para isso, mas essa ferramenta é do `akoflow-examples` e não li o código dela.

## 4. Diferenças entre a doc e o que já desenhamos

Comparado com `akoflow/workflow.yaml` (40 atividades; cabeçalho cita `akoflow-examples @ d890c3171e36` e `akoflow docs/docs/internal/workflow-spec.md @aeea81c`) e `docs/akoflow/02-desenho-granular.md`. Li o `workflow.yaml` por amostragem (cabeçalho, duas atividades, o fim) e as seções 2.3 e 8 do 02; as demais seções do 02 não foram relidas.

### 4.1 Campos

| Item | Nosso `workflow.yaml` | Doc | Veredito |
|---|---|---|---|
| Raiz `name`, `spec.namespace`, `spec.activities` | presentes (`namespace: phylotreeminer`, suposição) | iguais | compatível |
| `command.entrypoint` + `command.arguments` | `entrypoint: python`, `arguments: [-m, workflow.akoflow, ...]` | mesmos campos | compatível |
| `command.executable` | **ausente**; usamos `spec.image` | `spec.image` vira executável OCI com `delivery.strategy: auto`, mas é forma "legada"; as novas definições "devem usar `command.executable`" | funciona, porém legado. Considerar `command.executable.source.type: oci` + `reference` |
| `cpuLimit: "4"`, `memoryLimit: 2Gi` | presentes | `cpuLimit` aceita `0.5`/`250m`; padrão `0.1`; `memoryLimit` padrão `16Mi` | compatível |
| `dependsOn` | presente | presente | compatível |
| `runtime: phylotreeminer-local` por atividade | presente, marcado SUPOSIÇÃO | `runtime` é **dica retida como metadado**, não aloca recurso nem seleciona modo | **diverge**: o runtime real vem do ambiente (`runtimes[].driver`), dos `resourceRuntimeBindings`, do escopo e do plano; o campo da atividade não o escolhe |
| `spec.dataDependencies` | presente (producer/consumer/logicalName/sizeBytes) | mesmo formato; `logicalName` é rótulo de planejamento, `sizeBytes` > 0, e **exige `dependsOn` correspondente** | verificar que cada par tem o `dependsOn` (não conferi) |
| Nomes de ID | `phylotreeminer-zika21-granular` | ID derivado de `name`, minúsculas com hífen | compatível |
| `expectedOutputs` | **ausente** | opcional; afirma existência de arquivos relativos; falha o run se faltar | **falta**; útil para verificar saídas |
| `memoryLimit`/`cpuLimit` como "limite" | usamos para dimensionar | também alimentam o planejamento (`cpuLimit` é demanda) | semântica a confirmar com a equipe |

Nomes de campo que **mudaram**: não encontrei nenhum campo do nosso YAML que tenha sido renomeado; o que mudou é o **status**: `spec.image`, `image` por atividade e `run` passaram a legados. Nossos comentários dizem "Formato do workflow.yaml REAL de akoflow-examples"; a doc usa um formato "portátil" e o importador faz a normalização.

### 4.2 Como se declaram runtime, volumes e dados

- **Runtime**: no **ambiente** (`environment.yaml`, uma `EnvironmentDefinition`): `runtimes[]` (`id`, `driver`: `local`, `kubernetes`, `slurm`, `simgrid`...; `mode`: `execution` ou simulação; `configuration`; `capabilities`), `resources[]` (`type`: `local_machine`, `kubernetes_machine`, `hpc_partition`...; `cpuCores`, `memoryBytes`, `schedulable`) e `resourceRuntimeBindings[]` (`resourceId`, `runtimeId`, `enabled`). Não é o `workflow.yaml` que o define. Nosso `runtime: phylotreeminer-local` precisa de um ambiente da equipe que declare esse ID (nossos comentários já dizem isso), mas, mesmo assim, só como dica.
- **Volumes/armazenamento**: **não** no workflow. `storageClassName`, `storageSize`, `storagePolicy` e `mountPath` do workflow são aceitos e **ignorados**. No ambiente: `storages[]` (`type`: `pvc`, `nfs`, `lustre`..., `endpoint`, `shared`, `readOnly`, `configuration`, `runtimeBindings[]` com `hostPath`/`containerPath`). Nosso 02 escolheu (D-c) a entrada **dentro da imagem** (`COPY data/`) para Zika e VARV e "volume ou storage do AkôFlow" para ZIKV-480 e maiores: o segundo caminho passa pelo ambiente, não pelo workflow.
- **Dados entre atividades**: não há volume compartilhado declarado no workflow; o AkôFlow faz `rsync` do **workspace completo** do predecessor para o sucessor antes de iniciar. `dataDependencies` serve ao planejamento (tamanhos e transferências previstas), não define quais arquivos movem. Executáveis e imagens: `command.executable.source`.

### 4.3 Divergências críticas para o nosso desenho

1. **Conflito de caminho no `rsync`.** O sucessor recebe o workspace inteiro de cada predecessor e a doc diz que conteúdo diferente no mesmo caminho relativo é **rejeitado**. O nosso desenho usa `--out out` em todas as atividades e fragmentos como `out/outputs/config_backup.json`; atividades paralelas que escrevem o mesmo arquivo com conteúdo diferente (ex.: 14 pipelines com o mesmo nome de arquivo) travariam o sucessor que depende de ambas. Verificar nos nossos caminhos de saída (`03-dataflow.md`).
2. **`runtime:` na atividade não escolhe nada.** Sem um `environment.yaml` com runtime, recurso e binding (e escopo e plano), o workflow não executa. O que o `workflow.yaml` sozinho entrega é a definição versionada.
3. **Execução real exige plano.** O fluxo é ambiente -> escopo -> topologia -> workflow -> plano -> run. Não há "submit workflow.yaml" direto. O modo (real ou simulação) vem do escopo do plano, e uma atividade com `simulation` é só simulação; uma sem `simulation` é só real.
4. **Validação de plano importado não checa runtime binding**; um plano pode parecer viável e falhar na execução.
5. **`cpuLimit`/`memoryLimit` têm padrões baixos** (`0.1` núcleo, `16Mi`): omitir é perigoso para MAFFT/IQ-TREE/RAxML. Nossos valores são todos explícitos, mas marcados SUPOSIÇÃO.
6. **Suposições em aberto no nosso YAML**: `image` (`ghcr.io/SUPOSICAO/...`), `namespace`, `runtime`, memória. A doc de Kubernetes avisa que não configura `imagePullSecrets`; imagem privada exige configuração no cluster.
7. **Imagem sem `ENTRYPOINT`**: nosso 02 já diz isso (citando o exemplo). A doc lida não traz essa frase; não pude reconfirmá-la.
8. **Mismatch de caminhos da doc** (`examples/` x `akoflow-examples/`): qualquer instrução nossa que cite o exemplo SimGrid deve usar o repositório `akoflow-examples`.

## 5. O que a doc não responde

1. Onde ler o token e a URL da API quando se usa **apenas o Desktop empacotado**, para scripts. Ele "gerencia a conexão por proxy" e a doc nega que se precise extraí-la, mas não diz como chamar a API a partir do host.
2. Qual é o procedimento de pilha de desenvolvimento (clonar, `docker compose up`, token padrão). A first-run pressupõe "a pilha de desenvolvimento" sem guia; só o `docker-compose.yml` do repositório a mostra.
3. Se existe um tutorial só com o Desktop que rode um workflow **real local** (runtime `local`). A doc diz que ele "ainda precisa ser validado".
4. Contrato exato do runtime `local`: como a imagem é puxada ou construída, onde ficam os workspaces, como se monta entrada de dados do host, como o `command` é executado no contêiner (shell, `WORKDIR`, `ENTRYPOINT`).
5. Onde fica o diretório de trabalho da atividade e como se referenciam as saídas do predecessor (o caminho exato do workspace herdado via `rsync`) e se ele é o `workingDirectory` padrão.
6. Se o `rsync` pode ser restringido a arquivos específicos ou se o workspace herdado sempre vai inteiro (impacto para fan-out com FASTA grande, ZIKV-480).
7. Limites de escala do `dependsOn`/`dataDependencies` (quantas atividades, quantas arestas) e o custo de overhead por atividade (o nosso 02 estima 5 vezes o overhead no caminho crítico).
8. Como passar variáveis de ambiente por atividade fora de `command.environment`, e como tratar segredos (ex.: `NCBI_API_KEY`) nas atividades.
9. Cancelar, repetir ou retomar um run após falha de uma atividade (reaproveitar sucessos): a doc só diz "rerun after correcting the issue".
10. Exportar os artefatos de um run para fora (formato do manifesto de artefatos) e como obter a hash/versão da imagem como proveniência (relevante para a nossa reprodutibilidade).
11. Se `spec.dataDependencies` sem `dependsOn` realmente existe só como rótulo (a doc diz que é armazenado e ignorado na execução).
12. Compatibilidade entre versões: a doc não diz qual `contract`/versão do daemon aceita qual versão do formato de `workflow.yaml`; nosso arquivo foi gerado contra `aeea81c`.
13. Se o showcase `infrastructure-selection` (nosso modelo) usa os campos portáteis atuais ou os legados; não li o YAML dele.
14. Minikube: nenhuma página lida o cita; não sei se é suportado ou testado.
15. As 26 páginas de showcase, as páginas parciais ([P]) e as seções Reference/Explanations não foram lidas; o que elas dizem é desconhecido para este relatório.
