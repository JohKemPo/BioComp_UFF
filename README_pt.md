# PhyloTreeMiner - Um workflow para análise de árvores filogenéticas frequentes

Este projeto fornece um workflow completo e padronizado para análises filogenéticas, encapsulado em um ambiente Docker. A solução permite a execução de pipelines de alinhamento e construção de árvores filogenéticas de forma reprodutível e eficiente em qualquer máquina que possua Docker.

O script principal `workflow.py` é configurado através de um arquivo `templates/config.json`, permitindo total flexibilidade na escolha dos métodos e parâmetros.

## Por que Docker?

* **Reprodutibilidade:** O ambiente é 100% definido pelo arquivo `environment.yml`. Qualquer pessoa que construir esta imagem terá as mesmas versões de todos os softwares (Python, MAFFT, IQ-TREE, etc.), garantindo que os resultados sejam consistentes.
* **Portabilidade:** Funciona em qualquer sistema operacional (Linux, macOS, Windows) que tenha o Docker instalado.
* **Isolamento de Dependências:** Você não precisa instalar MAFFT, RAxML-NG ou qualquer pacote Python em sua máquina local. A imagem Docker contém tudo o que é necessário, evitando conflitos com outros projetos.
* **Facilidade de Teste:** Perfeito para testes de performance, permitindo que você execute o mesmo workflow em diferentes máquinas (com mais ou menos CPUs) simplesmente rodando o mesmo comando.

---

## Ferramentas Incluídas

Este ambiente Docker vem pré-instalado com todas as ferramentas de bioinformática e pacotes Python necessários, definidos no arquivo `environment.yml`. Os principais componentes incluem:

**Ferramentas de Linha de Comando:**
* MAFFT
* Clustal Omega (clustalo)
* IQ-TREE
* RAxML-NG

**Pacotes Python:**
* Python 3.10
* Biopython
* DendroPy
* Pandas, NumPy, Matplotlib

---

## Guia de Uso

Siga estes passos para construir a imagem e executar o workflow.

### Pré-requisitos

O único pré-requisito para rodar este projeto é ter o **[Docker](https://www.docker.com/products/docker-desktop/)** instalado e em execução no seu sistema.

## 1. Usando Docker (Recomendável)

### Passo 1: Construir a Imagem Docker

Antes de executar o workflow pela primeira vez, você precisa construir a imagem. Este processo lê o `Dockerfile`, baixa todas as dependências do Conda e copia seu código-fonte para dentro da imagem.

**Atenção:** Este passo pode levar vários minutos na primeira vez, pois o Conda precisa resolver e baixar todas as dependências.

Abra um terminal na raiz do projeto e execute:

```bash
docker build -t meu-workflow:latest .
```

- `-t meu-workflow:latest` : Define um nome (meu-workflow) e uma tag (latest) para a sua imagem.

- `.` : Indica ao Docker para procurar o Dockerfile no diretório atual.

### Passo 2: Configurar sua Análise

O workflow é controlado pelo arquivo `templates/config.json`. Antes de rodar, edite este arquivo para definir os parâmetros da sua análise.

### Diretório de Configuração

Abaixo está um guia detalhado do arquivo de configuração JSON para o workflow. Ele descreve cada uma das opções e o propósito de cada configuração. O arquivo `.json` está disponível em `templates/config.json`.

```
{
    "log_file": true,
    "project_name": "ABC_3",
    "output_log": "./projects/#/out",
    "tree_config": {
        "mode": "auto",
        "ignore_mode": "",
        "construct_tree_method": "distance",
        "align_method": "mafft",
        "num_threads": 1,
        "input_path": "./data/Zika479_Test",
        "output_path": "./projects/#/out",
        "output_format": "nexus"
    },
    "subtree_config": {
        "construct_tree_method": "distance",
        "input_path": "./projects/#/out/Trees",
        "output_path": "./projects/#/out",
        "input_format": "nexus",
        "output_format": "nexus",
        "resume_infos": true,
        "save_metadata": true,
        "subtree_miner": true,
        "subtree_miner_configs": {
            "mode": "OFST",
            "save_fpmax": true,
            "output_path": "./projects/#/out",
            "support_fpmax": "auto"
        }
    }
}
```

### Explicação das Configurações

**1. log_file:**

- **Descrição:** Define se o workflow irá gerar um arquivo de log durante a execução.
- **Tipo:** booleano
- **Valores:** true significa que o log será gerado.

**2. output_log:**

- **Descrição:** Caminho onde o arquivo de log será salvo.
- **Tipo:** string
- **Exemplo:** "./projects/test_artigo_fulldataset2/out"

**3. project_name:**

- **Descrição:** Nome do projeto.
- **Tipo:** string
- **Exemplo:** "test_artigo_fulldataset2"

**4. tree_config:**

- **Descrição:** Configuração para a construção da árvore principal.

    **4.1. mode:**

    - **Descrição:** Modo de operação para construção de árvores.
    - **Tipo:** string
    - **Valores:** 
        - "auto" (automático), 
        - "OFST" (Apenas conjuntos da mesma árvore), 
        - "parsimony" (Usar parcimônia), 
        - "distance" (Usar matriz de distância).

    **4.2. construct_tree_method:**

    - **Descrição:** Método usado para construir a árvore.
    - **Tipo:** string
    - **Valores:** 
        - "upgma",
        - "nj".

    **4.3. input_path:**

     - **Descrição:** Caminho para os dados de entrada usados para construção da árvore.
     - **Tipo:** string
     - **Exemplo:** "./data/testset"

    **4.4. output_path:**

    - **Descrição:** Caminho para salvar a árvore gerada.
    - **Tipo:** string
    - **Exemplo:** "./projects/test_artigo_fulldataset2/out"
    
    **4.5. output_format:**

    - **Descrição:** Formato de saída da árvore gerada.
    - **Tipo:** string
    - **Valores:** 
        - "nexus",
        - "nwk".

    **4.6. align_method:**

    - **Descrição:** Método de alinhamento de sequências para construção de árvores.
    - **Tipo:** string
    - **Valores:** 
        - "mafft",
        - "clustalw".

**5. subtree_config:**

- **Descrição:** Configuração para construção e mineração de subárvores.

    **5.1. construct_tree_method:**

    - **Descrição:** Método usado para construção de árvores durante a mineração de subárvores.
    - **Tipo:** string
    - **Valor:** "nj" (usa o método Neighbor-Joining)

    **5.2. input_path:**

    - **Descrição:** Caminho para os arquivos de árvore de entrada a serem minerados.
    - **Tipo:** string
    - **Exemplo:** "./projects/test_artigo_fulldataset2/out/Trees"

    **5.3. output_path:**

    - **Descrição:** Caminho para salvar as subárvores geradas.
    - **Tipo:** string
    - **Exemplo:** "./projects/test_artigo_fulldataset2/out"

    **5.4. input_format:**

    - **Descrição:** Formato de entrada das árvores.
    - **Tipo:** string
    - **Valor:** "nexus"

    **5.5. output_format:**

    - **Descrição:** Formato de saída das subárvores geradas.
    - **Tipo:** string
    - **Valor:** "nexus"

    **5.6. resume_infos:**

    - **Descrição:** Define se informações anteriores serão reutilizadas ou se o processo começará do zero.
    - **Tipo:** booleano
    - **Valor:** false (o processo começará do zero)

    **5.7. save_metadata:**

    - **Descrição:** Define se os metadados para subárvores mineradas serão salvos.
    - **Tipo:** booleano
    - **Valor:** true (os metadados serão salvos)

    **5.8. subtree_miner:**

    - **Descrição:** Ativa ou desativa o processo de mineração de subárvores.
    - **Tipo:** booleano
    - **Valor:** true (o minerador de subárvores está ativo)

**6. subtree_miner_configs:**

- **Descrição:** Configurações específicas para o minerador de subárvores.

    **6.1. mode:**

    - **Descrição:** Modo de mineração de subárvores.
    - **Tipo:** string
    - **Valor:** "OFST" (apenas da mesma árvore base)

    **6.2. save_fpmax:**

    - **Descrição:** Define se os resultados do FPMax serão salvos.
    - **Tipo:** booleano
    - **Valores:** 
        - false,
        - true

    **6.3. output_path:**

    - **Descrição:** Caminho para salvar os resultados da mineração de subárvores.
    - **Tipo:** string
    - **Exemplo:** "./projects/test_artigo_fulldataset2/out"

    **6.4. support_fpmax:**

    - **Descrição:** Configuração de suporte do FPMax.
    - **Tipo:** string
    - **Valores:** 
        - "auto" (todos os valores), 
        - Qualquer número entre 0.1 e 0.9

### Parâmetros Chave:

`input_path`: O caminho dentro do container para seu arquivo FASTA. Mantenha o prefixo `data/`.

`output_path`: O caminho dentro do container onde os resultados serão salvos. Mantenha o prefixo `projects/`.

`num_threads`: O mais importante para performance. Define o número de threads que as ferramentas (MAFFT, IQ-TREE, etc.) devem usar.

`align_method`: O método de alinhamento a ser usado (ex: "mafft", "clustalw").

`construct_tree_method`: O método de construção de árvore (ex: "upgma", "nj").

### Passo 3: Executar o Workflow

Com a imagem construída e o `config.json` pronto, execute o container. O comando abaixo usa volumes `(-v)` para criar um "espelho" das suas pastas locais dentro do container.

**Isso garante que:**

1. O container possa ler seus arquivos em `data/`.

2. Os resultados escritos em `projects/` sejam salvos diretamente no seu computador.

3. O container leia o `config.json` mais recente do seu computador.

**Para Linux e macOS:**

```bash
docker run --rm \
  -v "$(pwd)/data":/app/data \
  -v "$(pwd)/projects":/app/projects \
  -v "$(pwd)/templates/config.json":/app/templates/config.json \
  meu-workflow:latest
```

**Para Windows (usando PowerShell):**

```powershell
docker run --rm `
  -v "${pwd}/data":/app/data `
  -v "${pwd}/projects":/app/projects `
  -v "${pwd}/templates/config.json":/app/templates/config.json `
  meu-workflow:latest
```

### Resultados

- `--rm`: Remove o container automaticamente após a execução, mantendo seu sistema limpo.

- `-v "$(pwd)/...":/app/...`: Mapeia as pastas locais (`$(pwd)/...`) para as pastas dentro do container (`/app/...`).

- `meu-workflow:latest`: O nome da imagem que você quer usar.

Após a execução, seus arquivos de resultado (árvores, alinhamentos, logs) estarão disponíveis na sua pasta `projects/` local.

## 2. Executando localmente

### Execução

Para executar o workflow, simplesmente informe o arquivo com as configurações do projeto através da CLI:

```bash
python3 workflow.py --path "templates/config.json"
```

ou

```bash
python workflow.py --path "templates/config.json"
```

## Documentação

Toda a documentação do projeto pode ser acessada executando o comando, se o diretório `html` não foi gerado após executar `setupWorkflow.py`:

```bash
pdoc3 --force --html workflow/
```

Após executar, simplesmente abra `index.html` no navegador, disponível em `html/workflow/index.html`.

## Licença

Este projeto está licenciado sob a [Licença MIT](https://opensource.org/licenses/MIT).