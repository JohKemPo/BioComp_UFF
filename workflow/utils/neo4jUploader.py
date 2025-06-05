from neo4j import GraphDatabase

# Configurações da instância
URI = "neo4j+s://795fbce6.databases.neo4j.io"
USERNAME = "neo4j"
PASSWORD = "W5oxKVUAezR78m2WIsblPwAwcNlb-L1amRIlGH6xbd8"
DATABASE = "neo4j" 

# Caminho para o arquivo .cql
CQL_FILE = "/home/hilai360/Documents/Joao - IC/Zika/adjusted.cql"

# Conectar ao banco
driver = GraphDatabase.driver(URI, auth=(USERNAME, PASSWORD))

def executar_cql(caminho_arquivo):
    with driver.session(database=DATABASE) as session:
        with open(caminho_arquivo, "r", encoding="utf-8") as file:
            comandos = file.read()
            for comando in comandos.strip().split(";"):
                if comando.strip():
                    print("Executando:", comando.strip()[:80], "...")
                    session.run(comando.strip())

try:
    executar_cql(CQL_FILE)
    print("✅ Upload executado com sucesso.")
except Exception as e:
    print("❌ Erro durante o upload:", e)
finally:
    driver.close()