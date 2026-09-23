import os
import sys

from neo4j import GraphDatabase

# Configurações da instância — lidas do ambiente (DEC-052: credencial Neo4j
# Aura em texto puro achada em jun/2025; corrigido para nunca mais gravar
# segredo em código). Sem `NEO4J_URI`/`NEO4J_PASSWORD`, o script recusa
# rodar em vez de cair para um valor hardcoded — mesmo padrão de
# `Backend/src/seguranca.py` (ADMIN_TOKEN) e de `.env.example` na raiz do
# projeto principal.
URI = os.environ.get("NEO4J_URI")
USERNAME = os.environ.get("NEO4J_USERNAME", "neo4j")
PASSWORD = os.environ.get("NEO4J_PASSWORD")
DATABASE = os.environ.get("NEO4J_DATABASE", "neo4j")

if not URI or not PASSWORD:
    raise SystemExit(
        "NEO4J_URI e NEO4J_PASSWORD precisam estar no ambiente — nenhuma "
        "credencial fica hardcoded neste arquivo (DEC-052)."
    )

# Caminho para o arquivo .cql — era um caminho absoluto de outra máquina
# (`/home/hilai360/...`), nunca portável; agora é argumento de linha de
# comando, sem valor padrão hardcoded.
if len(sys.argv) < 2:
    raise SystemExit("uso: python neo4jUploader.py <caminho para o .cql>")
CQL_FILE = sys.argv[1]

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