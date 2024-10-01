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

```
templates/
└── config.json
```

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
