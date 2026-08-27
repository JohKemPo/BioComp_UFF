from Bio import SeqIO
import logging
import os

# Função que verifica se todas as sequências são proteínas válidas no formato FASTA
def validate_sequences(file_path):
    """
    """
    # Define os caracteres válidos para uma sequência de proteína.
    valid_characters = set('ACDEFGHIKLMNPQRSTVWY')
    try:
        for record in SeqIO.parse(file_path, 'fasta'):
            sequence = str(record.seq).upper()
            if not sequence:  
                return False
            if not set(sequence).issubset(valid_characters):
                return False
    except FileNotFoundError:
        print(f"O arquivo '{file_path}' não foi encontrado.")
        return False

    return True

def duplicate_names(file_path):
    """
    """
    name_count = {}
    try:
        for record in SeqIO.parse(file_path, 'fasta'):
            name = record.id
            name_count[name] = name_count.get(name, 0) + 1
            if name_count[name] > 1:
                return True
    except FileNotFoundError:
        print(f"O arquivo '{file_path}' não foi encontrado.")
        return False

    return False

def duplicate_seq(file_path):
    """
    """
    list_seq = list()
    try:
        for record in SeqIO.parse(file_path, 'fasta'):
            seq = record.seq
            list_seq.append(seq)
            if(list_seq.count(seq) > 1):
                return (True,record.id)
    except FileNotFoundError:
        print(f"O arquivo '{file_path}' não foi encontrado.")
        return (False,"")

    return (False,"")

def deduplicar_por_sequencia(name, path, outputpath):
    """
    Grava um FASTA com uma entrada por **sequência distinta**, e devolve o que
    foi descartado.

    Esta função chamava-se `remove_pipe` e não removia pipe nenhum: o nome dizia
    uma coisa e o corpo fazia outra — deduplicação por conteúdo, guardando a
    primeira ocorrência. O chamador a invocava sob a mensagem "iniciando remoção
    de pipes", e **nada registrava quais registros saíam**.

    Isso não é detalhe de nomenclatura. Nos conjuntos de *Variola* a aquisição
    baixa, para o mesmo genoma, o registro do GenBank **e** a cópia curada do
    RefSeq (`NC_008291` = `DQ437594`, `NC_003391` = `AF438165`), de modo que o
    projeto chamado VARV-49 recebe **52 registros** e produz árvores com **49
    folhas** — sem que a diferença apareça em log, em `summary.json` ou em
    *Methods*. `n` é um número publicado. Ver
    [D23](../../../docs/science/02-defeitos-que-alteram-resultado.md#d23).

    Pior: como se guarda a **primeira** ocorrência e a ordem do arquivo difere
    entre conjuntos, o acesso que representa o grupo externo **muda de
    experimento para experimento** — `DQ437594` em VARV-52, `NC_008291` em
    VARV-121.

    A composição não é alterada aqui: a decisão do usuário em 2026-08-26 foi
    **declarar agora e corrigir na aquisição depois**. O que muda é que o
    descarte deixa de ser silencioso.

    Return
    ------
    tuple of (str, list)
        Caminho do FASTA deduplicado e lista de ``(descartado, mantido)``.
    """
    sequences = list(SeqIO.parse(path, "fasta"))
    unique_sequences = {}
    descartados = []
    for sequence in sequences:
        chave = str(sequence.seq)
        if chave not in unique_sequences:
            unique_sequences[chave] = sequence
        else:
            descartados.append((sequence.id, unique_sequences[chave].id))

    output_file_tmp = os.path.join(outputpath, f'{name}_NoPipe')
    SeqIO.write(list(unique_sequences.values()), output_file_tmp, "fasta")

    if descartados:
        logging.warning(
            f"{len(sequences)} registros → {len(unique_sequences)} sequências distintas. "
            f"Descartados por sequência idêntica: "
            + "; ".join(f"{perdido} (idêntico a {mantido})" for perdido, mantido in descartados)
        )
    return output_file_tmp, descartados


#: Nome antigo, mantido para não quebrar chamador de fora do repositório. Ele
#: descreve o que a função **não** faz; use `deduplicar_por_sequencia`.
def remove_pipe(name, path, outputpath):
    caminho, _ = deduplicar_por_sequencia(name, path, outputpath)
    return caminho
