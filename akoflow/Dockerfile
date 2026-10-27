# syntax=docker/dockerfile:1
#
# Imagem única do workflow granular no AkôFlow (como o exemplo `infrastructure-selection`).
#
# CONTEXTO DE BUILD = raiz do PhyloTreeMiner (onde estão environment.yml, requirements*.txt
# e o submódulo BioComp_UFF/). A receita PINADA do ambiente é a de lá, e não é duplicada aqui:
#
#   docker build -f BioComp_UFF/akoflow/Dockerfile -t phylotreeminer-akoflow:$(git -C BioComp_UFF rev-parse --short HEAD) .
#
# Sem ENTRYPOINT: o AkôFlow define o comando de cada atividade (`python -m workflow.akoflow ...`),
# e um ENTRYPOINT concorrente a quebraria. (BioComp_UFF/Dockerfile tem um e um environment.yml
# sem pinos: não serve para isto.)
FROM mambaorg/micromamba:1.5.10-jammy

# compilador para a seção pip, como em Backend/Dockerfile
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# MAFFT 7.526, IQ-TREE 3.1.3, RAxML-NG 2.0.2, FastTree 2.2.0, MrBayes 3.2.7, Clustal Omega 1.2.4
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml requirements.txt requirements-dev.txt /tmp/receita/
RUN micromamba install -y -n base -f /tmp/receita/environment.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /app/BioComp_UFF
COPY --chown=$MAMBA_USER:$MAMBA_USER BioComp_UFF/workflow/ /app/BioComp_UFF/workflow/
COPY --chown=$MAMBA_USER:$MAMBA_USER BioComp_UFF/akoflow/conferir_ferramentas.py /app/BioComp_UFF/akoflow/
# Entrada do teste (Zika-21: FASTA de 54,5 KB e o .gb congelado do modo offline) e configuração.
COPY --chown=$MAMBA_USER:$MAMBA_USER BioComp_UFF/akoflow/config/ /app/config/
COPY --chown=$MAMBA_USER:$MAMBA_USER BioComp_UFF/akoflow/data/ /app/data/
ENV PYTHONPATH=/app/BioComp_UFF MPLBACKEND=Agg

# A versão da ferramenta faz parte do resultado (DEC-044): o build falha se divergir da receita.
RUN python /app/BioComp_UFF/akoflow/conferir_ferramentas.py

# O daemon do AkôFlow cria o workspace de cada atividade como root no host (runtime local) e o
# monta no contêiner; com o usuário padrão do micromamba (uid 57439) a atividade falha com
# `PermissionError: '.'` (medido, docs/akoflow/08). Por isso a imagem roda como root.
USER root
WORKDIR /work
