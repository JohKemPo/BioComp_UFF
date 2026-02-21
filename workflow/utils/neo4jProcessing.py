import os, json
from neo4j import GraphDatabase

USER_PLACEHOLDER = "<<USER_UID>>"

def parse_tree(trees: list, path: str = None, mode: str = "completo"):
    if all(isinstance(i, list) for i in trees):
        cypher_commands = []
        for tree_list in trees:
            for tree in tree_list:
                for tree_name, subtrees in tree.items():
                    cypher_commands.extend(generate_cypher(tree_name, subtrees))
        
        # Salvar em um arquivo
        with open(os.path.join(path, 'neo4j_commands.cql'), "w") as f:
            f.write("\n".join(cypher_commands))
    else:
        cypher_commands = []
        for tree in trees:
            for tree_name, subtrees in tree.items():
                cypher_commands.extend(generate_cypher(tree_name, subtrees))
        
        # Salvar em um arquivo
        with open(os.path.join(path, f'neo4j_commands_{mode}.cql'), "w") as f:
            f.write("\n".join(cypher_commands))

def sanitize_json_string(json_obj):
    """Transforma o dicionário JSON em string corretamente escapada para Cypher."""
    json_str = json.dumps(json_obj, separators=(",", ":"))  
    return json_str.replace("'", "\\'")

def generate_cypher(tree_name, subtrees):
    cypher_statements = []
    
    cypher_statements.append(f"""
    MERGE (u:User {{uid: '{USER_PLACEHOLDER}'}})
    CREATE (t:Tree {{name: '{tree_name}', uid: '{USER_PLACEHOLDER}'}})
    MERGE (u)-[:OWNS]->(t);
    """)
    
    for subtree_name, subtree_data in subtrees.items():
        cypher_statements.extend(create_subtree(tree_name, subtree_name, subtree_data))
    return cypher_statements

def create_subtree(parent_name, subtree_name, subtree_data):
    cypher_statements = []
    
    # Criar nó da Subtree (filho) ligado ao pai
    cypher_statements.append(f"""
    MATCH (parent {{name:'{parent_name}', uid: '{USER_PLACEHOLDER}'}})
    WHERE parent:Tree OR parent:Subtree
    CREATE (child:Subtree {{name: '{subtree_name}', uid: '{USER_PLACEHOLDER}'}})
    CREATE (parent)-[:HAS_SUBTREE]->(child);
    """)
    
    # Adicionar relacionamentos para cada suporte
    supports = subtree_data.get('supports', [])
    for support in supports:
        cypher_statements.append(f"""
        MATCH (child:Subtree {{name: '{subtree_name}', uid: '{USER_PLACEHOLDER}'}})
        MERGE (s:Support {{value: {support}}})
        CREATE (child)-[:HAS_SUPPORT]->(s);
        """)
        
    # Processamento de Metadata e Features
    metadatas = subtree_data.get('data_terminals', [])
    for metadata_entry in metadatas:
        # Acessar os dados internos
        inner_meta = metadata_entry.get('metadata', {})
        annotations = inner_meta.get('annotations', {})
        
        # Extração dos campos solicitados
        props = {
            'molecule_type': annotations.get('molecule_type'),
            'topology': annotations.get('topology'),
            'date': annotations.get('date'),
            'source': annotations.get('source'),
            'organism': annotations.get('organism'),
            'taxonomy': annotations.get('taxonomy'), # Lista
            'description': inner_meta.get('description'),
            'newick_id': metadata_entry.get('newick'),
            'terminal_hash': metadata_entry.get('terminal_hash')
        }
        
        # Construir string de propriedades para o nó Metadata
        prop_str_parts = []
        for k, v in props.items():
            if v is not None:
                # json.dumps garante a formatação correta de strings e listas para o Cypher
                prop_str_parts.append(f"{k}: {json.dumps(v)}")
        
        prop_body = ", ".join(prop_str_parts)
        
        # Início do bloco de comando para este metadado
        query_parts = []
        query_parts.append(f"MATCH (child:Subtree {{name: '{subtree_name}', uid: '{USER_PLACEHOLDER}'}})")
        # Cria o nó Metadata com os campos extraídos
        query_parts.append(f"CREATE (m:Metadata {{{prop_body}}})")
        query_parts.append(f"CREATE (child)-[:HAS_METADATA]->(m)")
        
        # Processamento das Features
        features = inner_meta.get('features', [])
        for i, feature in enumerate(features):
            f_var = f"f_{i}" # Variável única para a feature dentro da query
            
            f_type = feature.get('type')
            f_loc = feature.get('location')
            f_strand = feature.get('strand')
            
            f_props = {}
            if f_type: f_props['type'] = f_type
            if f_loc: f_props['location'] = f_loc
            if f_strand is not None: f_props['strand'] = f_strand
            
            f_prop_str = ", ".join([f"{k}: {json.dumps(v)}" for k, v in f_props.items()])
            
            # Criar nó Feature ligado ao Metadata
            query_parts.append(f"CREATE (m)-[:HAS_FEATURE]->({f_var}:Feature {{{f_prop_str}}})")
            
            # Processamento dos Qualifiers (Sub-atributos da Feature)
            qualifiers = feature.get('qualifiers', {})
            for j, (q_key, q_val) in enumerate(qualifiers.items()):
                q_var = f"q_{i}_{j}"
                # q_val é geralmente uma lista de strings no GenBank/Biopython
                q_val_json = json.dumps(q_val)
                
                # Criar nó Qualifier ligado à Feature
                query_parts.append(f"CREATE ({f_var})-[:HAS_QUALIFIER]->({q_var}:Qualifier {{key: '{q_key}', value: {q_val_json}}})")
        
        query_parts.append(";") # Finaliza o comando Cypher deste bloco
        
        # Adiciona o bloco completo à lista
        cypher_statements.append("\n".join(query_parts))
        
    # Recursão para subárvores filhas
    for key, value in subtree_data.items():
        if isinstance(value, dict):
            cypher_statements.extend(create_subtree(subtree_name, key, value))
            
    return cypher_statements

class Neo4jUploader:

    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        if self.driver:
            self.driver.close()

    def execute_cql_file(self, cql_file_path):
        with self.driver.session() as session:
            with open(cql_file_path, 'r') as file:
                cql_commands = file.read()
                session.run(cql_commands)
                print(f"Executed commands from {cql_file_path}")

