# A workflow for analysis of frequent phylogenetics trees

## Prerequisites

- Python 3.10.12
- Clustalw (version 2.1)
- FASTA file with protein sequences

## Installing Dependencies

Before running the project, it is necessary to install the Python dependencies specified in the `requirements.txt` file. To do this, run the following command in the terminal:

1. **Execute the script `workflow/setupWorkflow.py` or follow the steps below**

```
python3 workflow/setupWorkflow.py
```

2. **Install requeirements `manually`**:

:exclamation:  Not necessary if you have performed the previous step.

```
pip install --ignore-installed -r requirements.txt
```

To install the dependencies on Linux (Ubuntu):

```
sudo apt update
sudo apt-get install clustalw
```

### Installation of other necessary tools:

**Phyml**
```
sudo apt install phyml
```

<!-- **Muscle**
```
wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar -xzvf muscle3.8.31_i86linux64.tar.gz
``` -->
<!-- 
```
conda install -c etetoolkit ete3 ete_toolchain
ete3 build check
``` -->
**Clustalw**
```
sudo apt install clustalw 
``` 

**Mafft** 

``` 
sudo apt install mafft 
```

## Workflow directory

```
workflow/
├── alignment
│   └── alignmentSeq.py
│
├── controller
│   ├── subtreeBuilderController.py
│   ├── subtreeMinerController.py
│   └── treeBuilderController.py
│
├── data_handler
│
├── optimization
│
├── subtree_analysis
│
├── subtree_construction
│   └── builder.py
│
├── subtree_mining
│   └── miner.py
│
├── tests
│
├── tree_construction
│   └── builder.py
│
├── utils
│   ├── dataCleaning.py
│   ├── dataValidation.py
│   ├── messages.py
│   ├── metrics.py
│   ├── neo4jProcessing.py
│   └── treeUtils.py
│
├── setupWorkflow.py
└── visualization
```

### Configuration directory

Below is a detailed guide to the JSON configuration file for the workflow. It describes each of the options and the purpose of each setting. The `.json` file is available at `templates/config.json`.



```
{
    "log_file": true,
    "output_log": "./projects/test_artigo_fulldataset2/out",
    "tree_config": {
        "mode": "auto",
        "construct_tree_method": "upgma",
        "input_path": "./data/testset",
        "output_path": "./projects/test_artigo_fulldataset2/out",
        "output_format": "nexus",
        "align_method": "mafft"
    },
    "subtree_config": {
        "construct_tree_method": "nj",
        "input_path": "./projects/test_artigo_fulldataset2/out/Trees",
        "output_path": "./projects/test_artigo_fulldataset2/out",
        "input_format": "nexus",
        "output_format": "nexus",
        "resume_infos": false,
        "save_metadata": true,
        "subtree_miner": true,
        "subtree_miner_configs": {
            "mode": "OFST",
            "save_fpmax": false,
            "output_path": "./projects/test_artigo_fulldataset2/out",
            "support_fpmax": "auto"
        }
    }
}
```

### Configuration Explanation

**1. log_file:**

- **Description:** Defines whether the workflow will generate a log file during execution.
- **Type:** boolean.
- **Values:** true means the log will be generated.

**2. output_log:**

- **Description:** Path where the log file will be saved.
- **Type:** string
- **Example:** "./projects/test_artigo_fulldataset2/out"

**3. tree_config:**

- **Description:** Configuration for the main tree construction.

    **3.1. mode:**

    - **Description:** Operation mode for tree construction.
    - **Type:** string
    - **Value:** 
        - "auto" (automatic), 
        - "OFST" (Only sets of the same tree), 
        - "parsimony" (Use parsimony), 
        - "distance" (Use distance matrix).

    **3.2. construct_tree_method:**

    - **Description:** Method used to construct the tree.
    - **Type:** string
    - **Value:** 
        - "upgma",
        -  "nj".

    **3.3. input_path:**

     - **Description:** Path to the input data used for tree construction.
     - **Type:** string
     - **Example:** "./data/testset"

    **3.4. output_path:**

    - **Description:** Path to save the generated tree.
    - **Type:** string
    - **Example:** "./projects/test_artigo_fulldataset2/out"
    
    **3.5. output_format:**

    - **Description:** Output format of the generated tree.
    - **Type:** string
    - **Value:** 
        - "nexus",
        - "nwk".

    **3.6. align_method:**

    - **Description:** Sequence alignment method for tree construction.
    - **Type:** string
    - **Value:** 
        - "mafft",
        - "clustalw".

**4. subtree_config:**

- **Description:** Configuration for subtree construction and mining.

    **4.1. construct_tree_method:**

    - **Description:** Method used for tree construction during subtree - **mining.**
    - **Type: st**ring
    Value: "nj" (uses Neighbor-Joining method)

    **4.2. input_path:**

    - **Description:** Path to the input tree files to be mined.
    - **Type:** string
    - **Exampl**e: "./projects/test_artigo_fulldataset2/out/Trees"

    **4.3. output_path:**

    - **Description:** Path to save the generated subtrees.
    - **Type:** string
    - **Exampl**e: "./projects/test_artigo_fulldataset2/out"

    **4.4. input_format:**

    - **Description:** Input format of the trees.
    - **Type:** string
    - **Value:** "nexus"

    **4.5. output_format:**

    - **Description:** Output format of the generated subtrees.
    - **Type:** string
    - **Value:** "nexus"

    **4.6. resume_infos:**

    - **Description:** Defines whether previous information will be reused - **or the** process will start from scratch.
    - **Type: b**oolean
    Value: false (the process will start from scratch)

    **4.7. save_metadata:**

    - **Description:** Defines whether metadata for mined subtrees will be saved.
    - **Type:** boolean
    - **Value:** true (metadata will be saved)

    **4.8. subtree_miner:**

    - **Description:** Enables or disables the subtree mining process.
    - **Type:** boolean
    - **Value:** true (subtree miner is active)

**5. subtree_miner_configs:**

- **Description:** Specific configurations for subtree miner.

    **5.1. mode:**

    - **Description:** Subtree mining mode.
    - **Type:** string
    - **Value:** "OFST" (only from the same base tree)

    **5.2. save_fpmax:**

    - **Description:** Defines whether FPMax results will be saved.
    - **Type:** boolean
    - **Value:** 
        - false,
        - true

    **5.3. output_path:**

    - **Description:** Path to save the subtree mining results.
    - **Type:** string
    - **Example**: "./projects/test_artigo_fulldataset2/out"

    **5.4. support_fpmax:**

    - **Description:** FPMax support setting.
    - **Type:** string
    - **Value:** 
        - "auto" (all values), 
        - Any number between 0.1 at 0.9

### Data directory

```
data/
├── full_dataset
│   
├── full_dataset_plasmodium
│   
├── plasm
│   
├── plasmodium
│   
└── testset
```

## Execution

To execute the workflow, simply enter the file with the project settings through the CLI:

```
python3 workflow.py  --path "templates/config.json"
```

or

```
python workflow.py  --path "templates/config.json"
```

## Documentation

All project documentation can be accessed by running the command, if the directory `html` has not been generated after running `setupWorkflow.py`:

```
pdoc3 --force --html workflow/
```

After running, simply open `index.html` in the browser, available at `html/workflow/index.html`.


## Licença

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
