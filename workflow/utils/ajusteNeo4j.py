import re
import json
from pathlib import Path

# Caminho de entrada/saída
INPUT_CQL = "/home/hilai360/Documents/Joao - IC/Zika/projects/Zika 479 Test/out/outputs/neo4j_commands_tree_dataset.cql"
OUTPUT_CQL = "adjusted.cql"

def sanitize_json_string(json_obj):
    """Transforma o dicionário JSON em string corretamente escapada para Cypher."""
    json_str = json.dumps(json_obj, separators=(",", ":"))  # remove espaços extras
    return json_str.replace("'", "\\'")

def ajustar_metadata_com_contexto(cql_text):
    linhas = cql_text.splitlines()
    resultado = []
    i = 0

    while i < len(linhas):
        linha = linhas[i].strip()

        # Captura o contexto do MATCH
        if linha.startswith("MATCH (child:Subtree"):
            resultado.append(linhas[i])  # mantém o indentado original
            i += 1
            continue

        # Processa bloco com MERGE de Metadata
        if linha.startswith("MERGE (m:Metadata"):
            json_linha = linha
            while ";" not in linhas[i]:
                i += 1
                json_linha += " " + linhas[i].strip()
            i += 1  # pula o ponto e vírgula

            # Extrai o JSON da string Cypher
            match = re.search(r"value:\s*'(.*)'", json_linha, re.DOTALL)
            if match:
                try:
                    raw_json_str = match.group(1).replace("\\'", "'")  # desscapa aspas simples
                    parsed = json.loads(raw_json_str)
                    json_escaped = sanitize_json_string(parsed)

                    resultado.append(f"WITH '{json_escaped}' AS json")
                    resultado.append("CALL apoc.convert.fromJsonMap(json) YIELD value AS data")
                    resultado.append("MERGE (m:Metadata {value: data})")
                except Exception as e:
                    resultado.append("// ERRO DE PARSE JSON")
                    print(f"⚠️ Erro ao converter JSON na linha {i}: {e}")
            else:
                resultado.append("// ERRO: JSON não encontrado")
                print(f"⚠️ JSON não encontrado na linha {i}")
            continue

        # Mantém a linha de relacionamento
        if "CREATE (child)-[:HAS_METADATA]->(m)" in linha:
            resultado.append(linhas[i])
            i += 1
            continue

        # Qualquer outra linha é copiada normalmente
        resultado.append(linhas[i])
        i += 1

    return "\n".join(resultado)

# Execução principal
if __name__ == "__main__":
    cql_path = Path(INPUT_CQL)
    if not cql_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {INPUT_CQL}")

    texto_original = cql_path.read_text(encoding="utf-8")
    texto_ajustado = ajustar_metadata_com_contexto(texto_original)
    Path(OUTPUT_CQL).write_text(texto_ajustado, encoding="utf-8")

    print(f"✅ Arquivo corrigido salvo como: {OUTPUT_CQL}")
