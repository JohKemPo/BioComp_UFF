import subprocess
import sys
import os
import logging
import datetime

date = datetime.datetime.now()
logging.basicConfig(level=logging.INFO, 
                    filename=f"log_setup_{date.year}_{date.month}_{date.day}.log",
                    format='%(asctime)s - %(levelname)s - %(message)s')
logging.info("Iniciando a processo de criação do ambiente e instalação de dependências...")
def create_virtualenv_and_install_deps():
    """
    Cria um ambiente virtual e instala dependências listadas em requirements.txt.

    Esta função verifica se um ambiente virtual com o nome `BioComp` já existe,
    e se não existir, cria um novo ambiente virtual. Ela também atualiza o pip,
    instala as dependências listadas no arquivo `requirements.txt` e exibe as
    bibliotecas instaladas/atualizadas. Um aviso final é impresso com a instrução
    para ativar o ambiente manualmente.

    Parameters
    ----------
    None

    Return
    ------
    `None`: 
        Esta função não retorna nada, mas executa os comandos de criação
        e configuração do ambiente virtual, incluindo a instalação de
        dependências e exibição de bibliotecas instaladas.

    Raises
    ------
    subprocess.CalledProcessError
        Se algum comando de subprocesso falhar, como ao tentar atualizar o pip ou
        instalar dependências.
    
    FileNotFoundError
        Se o arquivo `requirements.txt` não for encontrado no diretório atual.

    Examples
    --------
    >>> create_virtualenv_and_install_deps()
    Ambiente virtual criado em: BioComp
    Dependências instaladas.
    ------------------------------
    Bibliotecas instaladas/atualizadas:
    <lista de bibliotecas>
    Para ativar o ambiente:
    $ source BioComp/bin/activate
    """
    # Nome do ambiente virtual
    venv_dir = "BioComp"

    try:
        logging.info("     Atualizando pip...")
        subprocess.check_call(['python', '-m', 'pip', 'install', '--upgrade', 'pip'])
    except Exception as error:
        logging.error(error)

    # Verificar se o ambiente virtual já existe
    if not os.path.exists(venv_dir):
        # Criar o ambiente virtual
        subprocess.check_call([sys.executable, "-m", "venv", venv_dir])
        print(f"Ambiente virtual criado em: {venv_dir}")
        logging.info(f"     Ambiente virtual criado em: {venv_dir}")
    else:
        print(f"Ambiente virtual já existe em: {venv_dir}")
        logging.info(f"     Ambiente virtual já existe em: {venv_dir}")

    # Ativar o ambiente virtual e instalar dependências
    activate_venv = os.path.join(venv_dir, 'Scripts', 'activate') if os.name == 'nt' else os.path.join(venv_dir, 'bin', 'activate')

    # Instalar dependências do requirements.txt
    if os.path.exists('requirements.txt'):
        subprocess.check_call([os.path.join(venv_dir, 'bin', 'python'), '-m', 'pip', 'install', '--upgrade', '-r', 'requirements.txt'])
        print("Dependências instaladas.")
        logging.info("     Dependências instaladas.")
    else:
        print("Arquivo requirements.txt não encontrado. Nenhuma dependência foi instalada.")
        logging.error("     Arquivo requirements.txt não encontrado. Nenhuma dependência foi instalada.")

    print(f"\n{'-'*30}\n\nBibliotecas instaladas/atualizadas:\n")
    subprocess.check_call(['pip','freeze'])

    print(f"\n{'-'*30}\n\nPara ativar o ambiente:\n\n   $ source {activate_venv}\n")


def install_bioinformatics_tools():
    """
    Instala as ferramentas ClustalW, MAFFT e PhyML.

    Esta função verifica se as ferramentas bioinformáticas `ClustalW`, `MAFFT` 
    e `PhyML` estão disponíveis no sistema e tenta instalá-las usando o 
    gerenciador de pacotes `apt`. Um log é mantido para registrar o sucesso 
    ou falha das instalações.

    Parameters
    ----------
    None

    Return
    ------
    None
        Esta função não retorna nada, mas executa comandos de instalação
        das ferramentas bioinformáticas necessárias.

    Raises
    ------
    subprocess.CalledProcessError
        Se algum comando de subprocesso falhar durante a instalação das ferramentas.

    Examples
    --------
    >>> install_bioinformatics_tools()
    ClustalW instalado com sucesso.
    MAFFT instalado com sucesso.
    PhyML instalado com sucesso.
    """

    tools = ['clustalw', 'mafft', 'phyml']

    for tool in tools:
        try:
            logging.info(f"     Instalando {tool}...")
            print(f"     Instalando {tool}...")
            subprocess.check_call(['sudo', 'apt', 'install', '-y', tool])
            logging.info(f"          {tool.capitalize()} instalado com sucesso.")
            print(f"     {tool.capitalize()} instalado com sucesso.")
        except subprocess.CalledProcessError as error:
            logging.error(     f"Erro ao instalar {tool}: {error}")
            print(f"          Erro ao instalar {tool}: {error}")

def gen_doc():
    """
    Generates the project's documentation using pdoc3.

    This function runs the pdoc3 command to generate HTML documentation 
    for the project located in the 'workflow/' directory.

    Logging is used to track the progress and any errors that occur during 
    the generation process.
    """
    print("Gerando documentação atualizada...")
    logging.info("Gerando documentação atualizada...")

    try:
        subprocess.check_call(['pdoc3', '--force', '--html', 'workflow/'])
        logging.info("Documentação gerada com sucesso.")
    except subprocess.CalledProcessError as error:
        logging.error(f"Erro ao gerar documentação: {error}")
    except Exception as error:
        logging.error(f"Erro inesperado: {error}")

if __name__ == "__main__":
    install_bioinformatics_tools()
    gen_doc()
    create_virtualenv_and_install_deps()