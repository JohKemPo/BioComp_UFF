import json, sys, os, argparse

from workflow.controller.treeBuilderController import TreeBuilderController
from workflow.controller.subtreeBuilderController import SubtreeBuilderController

BASE_PATH = os.path.dirname(os.path.abspath(__file__))

# Setup argument parser
parser = argparse.ArgumentParser(description="Execute the workflow with specified config file.")
parser.add_argument("--path", "-p", required=False, help="Path to the JSON configuration file.")
parser.add_argument("--config_workflow", "-cw", required=False, help="Params of workflow.")
parser.add_argument("--inputData", "-iData", required=False, help="Path to the input sequences.")
parser.add_argument("--projectName", "-pName", required=False, help="Name of project.")
args = parser.parse_args()

try:
    if args.inputData:
        os.path.exists(args.inputData)
except FileExistsError:
    print(f"{args.inputData} not exists.")
    sys.exit(1)

# Load the config file from the provided path
path = args.path
configs = args.config_workflow
params = None
try:
    if args.path:
        with open(args.path, 'r') as config_file:
            params = json.load(config_file)
    elif args.config_workflow:
        params = json.loads(args.config_workflow)
    else:
        print("Erro: Nenhuma configuração fornecida. Use --path ou --config_workflow.", file=sys.stderr)
        sys.exit(1)

    if args.inputData:
        params.setdefault('tree_config', {})['input_path'] = args.inputData
    if args.projectName:
        params['project_name'] = args.projectName
        
except FileNotFoundError:
    print(f"Config file at {path} not found.")
    sys.exit(1)


def replace_project_name(config):
    project_name = config['project_name']
    
    def recursive_replace(value):
        if isinstance(value, str):
            return value.replace('#', project_name)
        elif isinstance(value, dict):
            return {k: recursive_replace(v) for k, v in value.items()}
        elif isinstance(value, list):
            return [recursive_replace(item) for item in value]
        return value
    
    return recursive_replace(config)

if not args.config_workflow:
    params = replace_project_name(params)



# Create output directories if they do not exist
os.makedirs(params['output_log'], exist_ok=True)
os.makedirs(os.path.join(params['output_log'], 'outputs'), exist_ok=True)

# Save backup params
with open(os.path.join(params['output_log'], 'outputs','config_backup.json'),'w') as f:
    json.dump(params, f,  indent=4)
    
# Redirect output to log file if enabled
if params.get('log_file'):
    sys.stdout = open(os.path.join(params['output_log'], 'outputs', 'output_log.txt'), "w")

# Initialize and execute the TreeBuilderController
tree_builder_controller = TreeBuilderController(**params["tree_config"])
tree_builder_controller()

# # Initialize and execute the SubtreeBuilderController
subtree_builder_controller = SubtreeBuilderController(**params["subtree_config"])
subtree_builder_controller()