import os, json
from neo4j import GraphDatabase

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
    cypher_statements.append(f"CREATE (t:Tree {{name: '{tree_name}'}});")
    for subtree_name, subtree_data in subtrees.items():
        cypher_statements.extend(create_subtree(tree_name, subtree_name, subtree_data))
    return cypher_statements

def create_subtree(parent_name, subtree_name, subtree_data):
    cypher_statements = []
    
    # Extrair os suportes, se existirem
    supports = subtree_data.get('supports', [])
    metadatas = subtree_data.get('data_terminals', [])
    
    # Criar o nó da subárvore
    cypher_statements.append(f"""
    MATCH (parent:Tree {{name: '{parent_name}'}})
    CREATE (child:Subtree {{name: '{subtree_name}'}})
    CREATE (parent)-[:HAS_SUBTREE]->(child);
    """)
    
    # Adicionar relacionamentos para cada suporte
    for support in supports:
        cypher_statements.append(f"""
        MATCH (child:Subtree {{name: '{subtree_name}'}})
        MERGE (s:Support {{value: {support}}})
        CREATE (child)-[:HAS_SUPPORT]->(s);
        """)
        
    for metadata in metadatas:
        json_escaped = sanitize_json_string(metadata)
        cypher_statements.append(f"""
        MATCH (child:Subtree {{name: '{subtree_name}'}})
        MERGE (m:Metadata {{value: '{json_escaped}'}})
        CREATE (child)-[:HAS_METADATA]->(m);
        """)
        
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

