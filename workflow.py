import json, sys, os, argparse

from workflow.controller.treeBuilderController import TreeBuilderController
from workflow.controller.subtreeBuilderController import SubtreeBuilderController

# Setup argument parser
parser = argparse.ArgumentParser(description="Execute the workflow with specified config file.")
parser.add_argument("--path", "-p", required=True, help="Path to the JSON configuration file.")
args = parser.parse_args()

# Load the config file from the provided path
path = args.path
try:
    with open(path, 'rb') as configs:
        params = json.load(configs)
except FileNotFoundError:
    print(f"Config file at {path} not found.")
    sys.exit(1)

# Create output directories if they do not exist
os.makedirs(params['output_log'], exist_ok=True)
os.makedirs(os.path.join(params['output_log'], 'outputs'), exist_ok=True)

# Redirect output to log file if enabled
if params.get('log_file'):
    sys.stdout = open(os.path.join(params['output_log'], 'outputs', 'output_log.txt'), "w")

# Initialize and execute the TreeBuilderController
tree_builder_controller = TreeBuilderController(**params["tree_config"])
tree_builder_controller()

# Initialize and execute the SubtreeBuilderController
subtree_builder_controller = SubtreeBuilderController(**params["subtree_config"])
subtree_builder_controller()