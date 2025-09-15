import os, shutil, fnmatch

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

def copiar_arquivos(origem, destino, pattern='*'):
    """
    Copia arquivos com padrão específico, mantendo estrutura de diretórios.
    """
    os.makedirs(destino, exist_ok=True)
    for item in os.listdir(origem):
        if fnmatch.fnmatch(item, pattern):
            src_path = os.path.join(origem, item)
            dst_path = os.path.join(destino, item)
            if os.path.isfile(src_path):
                shutil.copy2(src_path, dst_path)