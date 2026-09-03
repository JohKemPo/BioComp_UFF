import json, sys, os, argparse

from workflow.controller.treeBuilderController import TreeBuilderController
from workflow.controller.subtreeBuilderController import SubtreeBuilderController
from workflow.utils.manifest import ExecutionManifest
from workflow.utils import run_logging
from workflow.utils import external_tools
from workflow.tree_construction.builder import reproducibility_settings

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



os.makedirs(params['output_log'], exist_ok=True)
os.makedirs(os.path.join(params['output_log'], 'outputs'), exist_ok=True)

with open(os.path.join(params['output_log'], 'outputs','config_backup.json'),'w') as f:
    json.dump(params, f,  indent=4)

# Manifesto de execução (M2.5 / D11): commit, versões de ferramenta, ambiente,
# sementes e SHA-256 de entradas e saídas. Gravado ANTES de rodar, para que uma
# execução que morra no meio ainda diga em que ambiente morreu.
project_root = os.path.dirname(os.path.abspath(params['output_log'].rstrip(os.sep)))
manifest = ExecutionManifest(
    project_root=project_root,
    params=params,
    repos={
        'BioComp_UFF': BASE_PATH,
        'PhyloTreeMiner': os.path.dirname(BASE_PATH),
    },
)
# D22 — o log desta execução carrega o `run_id` no nome. Antes ele se chamava
# `log_setup_{ano}_{mês}_{dia}.log` e era aberto em append, de modo que duas
# execuções no mesmo dia caíam no mesmo arquivo e ninguém conseguia separá-las
# depois. Configurado ANTES dos controladores, que é quem escrevia primeiro.
manifest.register_log(
    run_logging.configurar(os.path.join(params['output_log'], 'outputs'), manifest.run_id))

entrada = params.get('tree_config', {}).get('input_path')
if entrada and os.path.exists(entrada):
    manifest.register_input(entrada)

# Semente e paralelização efetivas, resolvidas pela MESMA função que o builder
# usa — o manifesto não pode declarar um valor e o pipeline usar outro.
manifest.register_reproducibility(reproducibility_settings(params.get('tree_config', {})))

# D18 — "disponível" não é "executado". Modo básico/auto roda só distância e
# parcimônia; nada além do `mode` cru em config_backup.json dizia isso, e o
# denominador de todo suporte metodológico (`M`) ficava sem proveniência.
# Calculado ANTES de rodar: disponibilidade e ignore_mode já são conhecidos, e
# não depende do resultado da execução.
_tree_cfg = params.get('tree_config', {})
_modo_solicitado = _tree_cfg.get('mode', 'advanced')
_ignore_bruto = _tree_cfg.get('ignore_mode') or []
if isinstance(_ignore_bruto, str):
    _ignore = [m.strip().lower() for m in _ignore_bruto.split(',') if m.strip() and m.strip() != 'none']
else:
    _ignore = [m.lower() for m in _ignore_bruto if m and m != 'none']
# Métodos avançados do controlador (`treeBuilderController._process_advanced_mode`)
# -> chave de `external_tools.CANDIDATOS`, que usa "raxml-ng", não "raxml".
_MAPA_METODO_FERRAMENTA = {'iqtree': 'iqtree', 'fasttree': 'fasttree',
                           'raxml': 'raxml-ng', 'mrbayes': 'mrbayes'}
_avancados_disponiveis = [
    metodo for metodo, ferramenta in _MAPA_METODO_FERRAMENTA.items()
    if external_tools.resolve_tool(ferramenta)
]
_avancados_executados = (
    [] if _modo_solicitado in ('auto', 'basic')
    else [m for m in _avancados_disponiveis if m not in _ignore]
)
manifest.register_execution_mode(_modo_solicitado, _avancados_disponiveis, _avancados_executados)

manifest.write()

if params.get('log_file'):
    sys.stdout = open(os.path.join(params['output_log'], 'outputs', 'output_log.txt'), "w")

try:
    tree_builder_controller = TreeBuilderController(**params["tree_config"])
    tree_builder_controller()

    subtree_builder_controller = SubtreeBuilderController(**params["subtree_config"])
    subtree_builder_controller()
finally:
    # Também no caminho de erro: o manifesto de uma execução que falhou é o que
    # permite diagnosticar a falha depois, e D17 mostrou que elas acontecem.
    #
    # `drain_tool_runs` recolhe o que alinhadores e inferidores registraram
    # durante a execução. Vem antes de `finish` e dentro do `finally` pelo mesmo
    # motivo: numa execução interrompida, saber qual comando estava rodando é
    # metade do diagnóstico.
    manifest.drain_tool_runs()
    manifest.register_outputs(os.path.join(params['output_log']))
    manifest.finish()