
# PhyloTreeMiner - A workflow for analysis of frequent phylogenetics trees

If you prefer, access the documentation in Portuguese, [Documentation in PT-BR](README_pt.md).


This project provides a complete and standardized workflow for phylogenetic analyses, encapsulated in a Docker environment. The solution allows execution of alignment and phylogenetic tree construction pipelines in a reproducible and efficient manner on any machine with Docker.

The main script `workflow.py` is configured through a `templates/config.json` file, allowing full flexibility in choosing methods and parameters.

## Why Docker?

* **Reproducibility:** The environment is 100% defined by the `environment.yml` file. Anyone building this image will have the same versions of all software (Python, MAFFT, IQ-TREE, etc.), ensuring consistent results.
* **Portability:** Works on any operating system (Linux, macOS, Windows) with Docker installed.
* **Dependency Isolation:** You don't need to install MAFFT, RAxML-NG, or any Python packages on your local machine. The Docker image contains everything needed, avoiding conflicts with other projects.
* **Easy Testing:** Perfect for performance testing, allowing you to run the same workflow on different machines (with more or fewer CPUs) by simply running the same command.

---

## Included Tools

This Docker environment comes pre-installed with all necessary bioinformatics tools and Python packages, defined in the `environment.yml` file. The main components include:

**Command Line Tools:**
* MAFFT
* Clustal Omega (clustalo)
* IQ-TREE
* RAxML-NG

**Python Packages:**
* Python 3.10
* Biopython
* DendroPy
* Pandas, NumPy, Matplotlib

---

## Usage Guide

Follow these steps to build the image and execute the workflow.

### Prerequisites

The only prerequisite to run this project is having **[Docker](https://www.docker.com/products/docker-desktop/)** installed and running on your system.

## 1. Using Docker (Recommended)

### Step 1: Build the Docker Image

Before running the workflow for the first time, you need to build the image. This process reads the `Dockerfile`, downloads all Conda dependencies, and copies your source code into the image.

**Note:** This step may take several minutes the first time, as Conda needs to resolve and download all dependencies.

Open a terminal at the project root and execute:

```bash
docker build -t meu-workflow:latest .
```

- `-t meu-workflow:latest` : Defines a name (meu-workflow) and tag (latest) for your image.

- `.` : Tells Docker to look for the Dockerfile in the current directory.

### Step 2: Configure Your Analysis

The workflow is controlled by the `templates/config.json` file. Before running, edit this file to define your analysis parameters.

### Configuration Directory

Below is a detailed guide to the JSON configuration file for the workflow. It describes each of the options and the purpose of each setting. The `.json` file is available at `templates/config.json`.

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

### Configuration Explanation

**1. log_file:**

- **Description:** Defines whether the workflow will generate a log file during execution.
- **Type:** boolean
- **Values:** true means the log will be generated.

**2. output_log:**

- **Description:** Path where the log file will be saved.
- **Type:** string
- **Example:** "./projects/test_artigo_fulldataset2/out"

**3. project_name:**

- **Description:** Name of project.
- **Type:** string
- **Example:** "test_artigo_fulldataset2"

**4. tree_config:**

- **Description:** Configuration for the main tree construction.

    **4.1. mode:**

    - **Description:** Operation mode for tree construction.
    - **Type:** string
    - **Values:** 
        - "auto" (automatic), 
        - "OFST" (Only sets of the same tree), 
        - "parsimony" (Use parsimony), 
        - "distance" (Use distance matrix).

    **4.2. construct_tree_method:**

    - **Description:** Method used to construct the tree.
    - **Type:** string
    - **Values:** 
        - "upgma",
        - "nj".

    **4.3. input_path:**

     - **Description:** Path to the input data used for tree construction.
     - **Type:** string
     - **Example:** "./data/testset"

    **4.4. output_path:**

    - **Description:** Path to save the generated tree.
    - **Type:** string
    - **Example:** "./projects/test_artigo_fulldataset2/out"
    
    **4.5. output_format:**

    - **Description:** Output format of the generated tree.
    - **Type:** string
    - **Values:** 
        - "nexus",
        - "nwk".

    **4.6. align_method:**

    - **Description:** Sequence alignment method for tree construction.
    - **Type:** string
    - **Values:** 
        - "mafft",
        - "clustalw".

**5. subtree_config:**

- **Description:** Configuration for subtree construction and mining.

    **5.1. construct_tree_method:**

    - **Description:** Method used for tree construction during subtree mining.
    - **Type:** string
    - **Value:** "nj" (uses Neighbor-Joining method)

    **5.2. input_path:**

    - **Description:** Path to the input tree files to be mined.
    - **Type:** string
    - **Example:** "./projects/test_artigo_fulldataset2/out/Trees"

    **5.3. output_path:**

    - **Description:** Path to save the generated subtrees.
    - **Type:** string
    - **Example:** "./projects/test_artigo_fulldataset2/out"

    **5.4. input_format:**

    - **Description:** Input format of the trees.
    - **Type:** string
    - **Value:** "nexus"

    **5.5. output_format:**

    - **Description:** Output format of the generated subtrees.
    - **Type:** string
    - **Value:** "nexus"

    **5.6. resume_infos:**

    - **Description:** Defines whether previous information will be reused or the process will start from scratch.
    - **Type:** boolean
    - **Value:** false (the process will start from scratch)

    **5.7. save_metadata:**

    - **Description:** Defines whether metadata for mined subtrees will be saved.
    - **Type:** boolean
    - **Value:** true (metadata will be saved)

    **5.8. subtree_miner:**

    - **Description:** Enables or disables the subtree mining process.
    - **Type:** boolean
    - **Value:** true (subtree miner is active)

**6. subtree_miner_configs:**

- **Description:** Specific configurations for subtree miner.

    **6.1. mode:**

    - **Description:** Subtree mining mode.
    - **Type:** string
    - **Value:** "OFST" (only from the same base tree)

    **6.2. save_fpmax:**

    - **Description:** Defines whether FPMax results will be saved.
    - **Type:** boolean
    - **Values:** 
        - false,
        - true

    **6.3. output_path:**

    - **Description:** Path to save the subtree mining results.
    - **Type:** string
    - **Example:** "./projects/test_artigo_fulldataset2/out"

    **6.4. support_fpmax:**

    - **Description:** FPMax support setting.
    - **Type:** string
    - **Values:** 
        - "auto" (all values), 
        - Any number between 0.1 and 0.9

### Key Parameters:

`input_path`: The path inside the container to your FASTA file. Keep the `data/` prefix.

`output_path`: The path inside the container where results will be saved. Keep the `projects/` prefix.

`num_threads`: Most important for performance. Sets the number of threads that tools (MAFFT, IQ-TREE, etc.) should use.

`align_method`: The alignment method to use (ex: "mafft", "clustalw").

`construct_tree_method`: The tree construction method (ex: "upgma", "nj").

### Step 3: Execute the Workflow

With the image built and `config.json` ready, run the container. The command below uses volumes `(-v)` to create a "mirror" of your local folders inside the container.

**This ensures that:**

1. The container can read your files in `data/`.

2. Results written to `projects/` are saved directly to your computer.

3. The container reads the latest `config.json` from your computer.

**For Linux and macOS:**

```bash
docker run --rm \
  -v "$(pwd)/data":/app/data \
  -v "$(pwd)/projects":/app/projects \
  -v "$(pwd)/templates/config.json":/app/templates/config.json \
  meu-workflow:latest
```

**For Windows (using PowerShell):**

```powershell
docker run --rm `
  -v "${pwd}/data":/app/data `
  -v "${pwd}/projects":/app/projects `
  -v "${pwd}/templates/config.json":/app/templates/config.json `
  meu-workflow:latest
```

### Results

- `--rm`: Automatically removes the container after execution, keeping your system clean.

- `-v "$(pwd)/...":/app/...`: Maps local folders (`$(pwd)/...`) to folders inside the container (`/app/...`).

- `meu-workflow:latest`: The name of the image you want to use.

After execution, your result files (trees, alignments, logs) will be available in your local `projects/` folder.

## 2. Running Locally

### Execution

To execute the workflow, simply enter the file with the project settings through the CLI:

```bash
python3 workflow.py --path "templates/config.json"
```

or

```bash
python workflow.py --path "templates/config.json"
```

## 3. Cross-Pipeline Stability Analysis

The `workflow/stability/` package analyses a set of trees that were built from the same data by *different* pipelines, and measures how much of the topology survives the change of method. It complements `workflow/subtree_mining/`: instead of mining subtrees within a tree, it treats each pipeline as a transaction and each clade as an item, so the support of a clade is the fraction of aligner × inference-method combinations that recover it.

It provides:

* **Canonical clade identity** — an order-invariant, 128-bit identity that replaces the ordered 16-bit hash, plus an audit that quantifies how much the legacy scheme fragments or collides clades.
* **Label reconciliation** — IQ-TREE and RAxML-NG truncate the version digit of RefSeq accessions when writing trees; affected labels are reconciled and reported per pipeline.
* **Exact pattern mining** — with M pipelines the closed-pattern lattice is enumerated exactly in O(2^M·|C|), so maximal patterns are exact rather than heuristic. No `mlxtend` dependency.
* **Factor decomposition** — mean Robinson–Foulds distance across pipeline pairs that differ only in the aligner versus only in the inference method.

Run the bundled Orthopoxvirus case study over all three Variola projects:

```bash
python -m workflow.stability.case_study
```

Or a single project:

```bash
python -m workflow.stability.case_study \
    --project projects/Variola_Yu_li_2007_200seq \
    --fasta   data/workflow_dataAcquisition_li_et_al_2007_replication-RetMax200/dataset_final.fasta \
    --label   VARV-121
```

Tables (`clade_support.csv`, `maximal_patterns.csv`, `rf_matrix.csv`, `summary.json`) and figures are written to `projects/<project>/out/outputs/stability/`. The case study, its findings and the accompanying manuscript draft are in [`docs/case_study_variola/`](docs/case_study_variola/).

Tests:

```bash
python -m unittest workflow.tests.test_stability
```

## Documentation

All project documentation can be accessed by running the command, if the `html` directory has not been generated after running `setupWorkflow.py`:

```bash
pdoc3 --force --html workflow/
```

After running, simply open `index.html` in the browser, available at `html/workflow/index.html`.

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).