import os, json, logging
from neo4j import GraphDatabase

USER_PLACEHOLDER = "<<USER_UID>>"

def project_name_from_output_path(path: str) -> str:
    """Deriva o nome do projeto a partir do caminho de saída.

    `path` é sempre `.../projects/<project_name>/out/outputs` (é o que
    `SubtreeMiner.output_path` + '/outputs' produz) — não existe hoje nenhum
    campo `project_name` propagado até aqui via config, então a identidade
    do projeto vem da própria estrutura de diretórios que o pipeline já usa.
    """
    return os.path.basename(os.path.dirname(os.path.dirname(path)))

def parse_tree(trees: list, path: str = None, mode: str = "completo"):
    project_name = project_name_from_output_path(path)
    # G1 — um deduplicador por arquivo gerado: cada `.cql` é carregável sozinho
    # e cria o `Metadata` de cada terminal na primeira vez que o vê.
    dedup = DeduplicadorMetadata(project_name)
    if all(isinstance(i, list) for i in trees):
        cypher_commands = []
        for tree_list in trees:
            for tree in tree_list:
                for tree_name, subtrees in tree.items():
                    cypher_commands.extend(generate_cypher(tree_name, subtrees, project_name, dedup))
        nome_arquivo = 'neo4j_commands.cql'
    else:
        cypher_commands = []
        for tree in trees:
            for tree_name, subtrees in tree.items():
                cypher_commands.extend(generate_cypher(tree_name, subtrees, project_name, dedup))
        nome_arquivo = f'neo4j_commands_{mode}.cql'

    if dedup.conflitos:
        logging.warning(
            "neo4j: %d terminal(is) com registro GenBank divergente entre ocorrências; "
            "mantido o primeiro. terminal_hash: %s",
            len(dedup.conflitos), sorted(set(dedup.conflitos))[:20])

    # Salvar em um arquivo
    with open(os.path.join(path, nome_arquivo), "w") as f:
        f.write("\n".join(cypher_commands))

def sanitize_json_string(json_obj):
    """Transforma o dicionário JSON em string corretamente escapada para Cypher."""
    json_str = json.dumps(json_obj, separators=(",", ":"))  
    return json_str.replace("'", "\\'")


# --------------------------------------------------------------------------
# G1 — `Metadata` é entidade de terminal, não de ocorrência
# (docs/science/05-grafo-neo4j.md §G1).
#
# Antes, cada ocorrência de um terminal em cada subárvore de cada árvore
# recriava com CREATE o nó `Metadata` e todo o registro GenBank abaixo dele
# (`Feature`, `Qualifier`). Medido nos `.cql` das reexecuções de 2026-09-01:
# 5,29 M de `Qualifier` para 49 táxons em VARV-49, 4,63 M para 121 em
# VARV-121. O fator é o número de subárvores que contêm o terminal, isto é,
# depende da topologia — não é corrigível por uma constante numa consulta.
#
# Agora cada terminal tem UM nó `Metadata`, chaveado por
# (uid, project, terminal_hash), e cada ocorrência vira só uma aresta
# `HAS_METADATA` da subárvore para ele. `uid` entra na chave para não fundir
# nós de usuários diferentes (a partição por `uid` é o único isolamento do
# grafo); `project` entra porque projetos diferentes podem ter baixado versões
# diferentes do mesmo registro (reexecuções pós-M1).
# --------------------------------------------------------------------------

_CHAVES_SO_IDENTIDADE = ("newick_id: ", "terminal_hash: ")


def chave_metadata(project_name: str, terminal_hash) -> str:
    """Corpo do mapa de propriedades que identifica o `Metadata` de um terminal.

    Parameters
    ----------
    project_name : str
        Projeto de origem (mesmo valor gravado em `Tree.project`).
    terminal_hash : int or str
        `terminal_hash` do terminal (`canonical_item_id([newick])`). Aceita o
        texto já serializado, que é o que o reparo de `.cql` legado tem.

    Returns
    -------
    str
        Ex.: ``uid: '<<USER_UID>>', project: "P", terminal_hash: 123``.
    """
    return (f"uid: '{USER_PLACEHOLDER}', project: {json.dumps(project_name)}, "
            f"terminal_hash: {terminal_hash}")


def instrucao_metadata(chave: str, prop_body: str, linhas_feature: list) -> str:
    """Instrução que cria (ou completa) o `Metadata` e seu registro GenBank.

    `MERGE` + `SET +=` em vez de `CREATE`: o executor em lote do Backend roda
    blocos em paralelo (`max_workers`), então a ligação de uma ocorrência
    pode chegar antes desta instrução — ela cria o nó só com a chave, e este
    `SET` o completa. As `Feature` só são criadas se o nó ainda não tem
    nenhuma, o que torna a instrução reexecutável sem duplicar o registro.
    """
    partes = [f"MERGE (m:Metadata {{{chave}}})", f"SET m += {{{prop_body}}}"]
    if linhas_feature:
        partes.append("WITH m WHERE NOT EXISTS { (m)-[:HAS_FEATURE]->() }")
        partes.extend(linhas_feature)
    partes.append(";")
    return "\n".join(partes)


def instrucao_ligacao(subtree_name: str, chave: str) -> str:
    """Uma ocorrência do terminal na subárvore: só a aresta `HAS_METADATA`.

    `MERGE` no `Metadata` (e não `MATCH`) para que a ligação não se perca se
    for executada antes da instrução que cria o nó — ver `instrucao_metadata`.
    """
    return "\n".join([
        linha_match_subtree(subtree_name),
        f"MERGE (m:Metadata {{{chave}}})",
        "CREATE (child)-[:HAS_METADATA]->(m)",
        ";",
    ])


def linha_match_subtree(subtree_name: str) -> str:
    return f"MATCH (child:Subtree {{name: '{subtree_name}', uid: '{USER_PLACEHOLDER}'}})"


class DeduplicadorMetadata:
    """Decide, ocorrência a ocorrência, quais instruções emitir para um terminal.

    Compartilhado pelo gerador (`create_subtree`) e pelo reparo de `.cql` já
    gerados (`workflow/utils/deduplicar_metadata_cql.py`), para que os dois
    produzam o mesmo texto.

    Regra: a instrução de criação sai na primeira ocorrência do terminal e,
    se essa primeira ocorrência veio sem registro GenBank (lookup falhou —
    `{"error": ...}`), sai de novo na primeira ocorrência completa. Toda
    ocorrência emite a ligação. Duas ocorrências completas com conteúdo
    diferente são registradas em `conflitos` (fica o primeiro registro).

    Parameters
    ----------
    project_name : str
        Projeto de origem, parte da chave do nó.
    """

    def __init__(self, project_name: str):
        self.project_name = project_name
        # terminal_hash (str) -> None (só a versão vazia foi emitida) ou o
        # conteúdo completo emitido (prop_body, linhas de feature).
        self.vistos: dict = {}
        self.conflitos: list = []

    @staticmethod
    def eh_completo(prop_body: str, linhas_feature: list) -> bool:
        """Se a ocorrência traz registro GenBank (algo além de newick_id/hash)."""
        return bool(linhas_feature) or not prop_body.startswith(_CHAVES_SO_IDENTIDADE)

    def instrucoes(self, subtree_name: str, terminal_hash, prop_body: str,
                   linhas_feature: list) -> list:
        """Instruções Cypher para uma ocorrência do terminal na subárvore.

        Returns
        -------
        list of str
            ``[criação, ligação]`` ou só ``[ligação]``.
        """
        h = str(terminal_hash)
        chave = chave_metadata(self.project_name, h)
        completo = self.eh_completo(prop_body, linhas_feature)
        conteudo = (prop_body, tuple(linhas_feature))

        saida = []
        if h not in self.vistos or (completo and self.vistos[h] is None):
            saida.append(instrucao_metadata(chave, prop_body, linhas_feature))
            self.vistos[h] = conteudo if completo else None
        elif completo and self.vistos[h] != conteudo:
            self.conflitos.append(h)
        saida.append(instrucao_ligacao(subtree_name, chave))
        return saida


def generate_cypher(tree_name, subtrees, project_name, dedup=None):
    cypher_statements = []
    if dedup is None:
        dedup = DeduplicadorMetadata(project_name)

    cypher_statements.append(f"""
    MERGE (u:User {{uid: '{USER_PLACEHOLDER}'}})
    CREATE (t:Tree {{name: '{tree_name}', uid: '{USER_PLACEHOLDER}', project: '{project_name}'}})
    MERGE (u)-[:OWNS]->(t);
    """)

    for subtree_name, subtree_data in subtrees.items():
        cypher_statements.extend(create_subtree(tree_name, subtree_name, subtree_data, dedup))
    return cypher_statements

def create_subtree(parent_name, subtree_name, subtree_data, dedup):
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

        terminal_hash = props['terminal_hash']
        if terminal_hash is None:
            # Sem chave não há como deduplicar; recriar o nó por ocorrência é
            # justamente o defeito G1.
            raise ValueError(
                f"data_terminal sem terminal_hash na subárvore {subtree_name!r} "
                f"(newick={props['newick_id']!r})")

        linhas_feature = []
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
            linhas_feature.append(f"CREATE (m)-[:HAS_FEATURE]->({f_var}:Feature {{{f_prop_str}}})")
            
            # Processamento dos Qualifiers (Sub-atributos da Feature)
            qualifiers = feature.get('qualifiers', {})
            for j, (q_key, q_val) in enumerate(qualifiers.items()):
                q_var = f"q_{i}_{j}"
                # q_val é geralmente uma lista de strings no GenBank/Biopython
                q_val_json = json.dumps(q_val)
                
                # Criar nó Qualifier ligado à Feature
                linhas_feature.append(f"CREATE ({f_var})-[:HAS_QUALIFIER]->({q_var}:Qualifier {{key: '{q_key}', value: {q_val_json}}})")

        cypher_statements.extend(
            dedup.instrucoes(subtree_name, terminal_hash, prop_body, linhas_feature))

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

