# Use uma imagem base com Python
FROM python:3.10.12-slim
ENV PYTHONUNBUFFERED 1

# Cria o ambiente virtual e define o PATH
RUN python3 -m venv BioComp
ENV VIRTUAL_ENV=/BioComp
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Instale compiladores e outras bibliotecas necessárias para pacotes como o pysqlite3
RUN apt-get update && apt-get install -y gcc libsqlite3-dev


# Instala as dependências do sistema e ferramentas bioinformáticas
RUN apt-get update && apt-get install -y \
    clustalw \
    phyml \
    mafft \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Define o diretório de trabalho
WORKDIR /app

# Copia todos os arquivos da aplicação para o contêiner
COPY . /app


# Instala as dependências do Python no ambiente virtual
RUN pip install -r requirements.txt

# # Gera a documentação usando pdoc3 se ainda não foi gerada
# RUN if ! [ -d "html" ]; then \
#         pip install --no-cache-dir pdoc3 && \
#         pdoc3 --force --html workflow/ && \
#         echo "Documentação gerada em html/workflow/index.html"; \
#     fi

# Define o comando de entrada padrão para executar o workflow com o ambiente virtual ativado
ENTRYPOINT ["python3", "workflow.py"]

# Define o comando padrão que pode ser sobrescrito
CMD ["--path", "templates/config.json"]
