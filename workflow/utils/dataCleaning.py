import os
import shutil

def clean_Trees(output_path: str) -> None:
    """
    Remove todos os arquivos no diretório 'Trees', exceto o arquivo 'file.gitkeep'.

    Parameters
    ----------
    output_path : str
        O caminho para o diretório onde o diretório 'Trees' está localizado.

    Returns
    -------
    None
    """
    dir_Trees = os.path.join(output_path, 'Trees')
    file_trees = os.listdir(dir_Trees)

    for name_file_trees in file_trees:
        path_trees = os.path.join(dir_Trees, name_file_trees)
        if name_file_trees != "file.gitkeep":
            os.remove(path_trees)

def clean_tmp(output_path: str) -> None:
    """
    Remove todos os arquivos no diretório 'tmp', exceto o arquivo 'file.gitkeep'.

    Parameters
    ----------
    output_path : str
        O caminho para o diretório onde o diretório 'tmp' está localizado.

    Returns
    -------
    None
    """
    dir_tmp = os.path.join(output_path, 'tmp')
    files_tmp = os.listdir(dir_tmp)

    for name_file_tmp in files_tmp:
        path_tmp = os.path.join(dir_tmp, name_file_tmp)
        if name_file_tmp != "file.gitkeep":
            os.remove(path_tmp)

def clean_dir(path_clean: str) -> None:
    """
    Remove todos os arquivos e o diretório especificado, exceto o arquivo 'file.gitkeep'.

    Parameters
    ----------
    path_clean : str
        O caminho para o diretório que será limpo e removido.

    Returns
    -------
    None
    """
    if os.path.exists(path_clean):
        for name in os.listdir(path_clean):
            if name != "file.gitkeep":
                file_path = os.path.join(path_clean, name)
                os.remove(file_path)
        shutil.rmtree(path_clean)

def clean_files(path: str) -> None:
    """
    Remove todos os arquivos no diretório especificado, exceto o arquivo 'file.gitkeep'.

    Parameters
    ----------
    path : str
        O caminho para o diretório cujos arquivos serão removidos.

    Returns
    -------
    None
    """
    dir_tmp = os.path.join(path)
    arquivos_tmp = os.listdir(dir_tmp)
    for name_file in arquivos_tmp:
        if name_file != "file.gitkeep":
            os.remove(os.path.join(dir_tmp, name_file))

def clean_NoPipe(path: str) -> None:
    """
    Remove todos os arquivos no diretório especificado que contêm 'NoPipe' no nome.

    Parameters
    ----------
    path : str
        O caminho para o diretório onde os arquivos serão limpos.

    Returns
    -------
    None
    """
    dir_NoPipe = os.path.join(path)
    for file_name in os.listdir(dir_NoPipe):
        if 'NoPipe' in file_name:
            file_path = os.path.join(dir_NoPipe, file_name)
            os.remove(file_path)

def copiar_arquivos(origem, destino):
    """
    Copia todos os arquivos da pasta 'origem' para a pasta 'destino'.
    Cria a pasta destino se não existir.
    """
    os.makedirs(destino, exist_ok=True)
    if not os.path.exists: return

    for arquivo in os.listdir(origem):
        caminho_origem = os.path.join(origem, arquivo)
        caminho_destino = os.path.join(destino, arquivo)

        if os.path.isfile(caminho_origem):
            shutil.copy2(caminho_origem, caminho_destino)  