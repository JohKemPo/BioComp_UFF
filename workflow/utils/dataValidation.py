from Bio import SeqIO
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

def remove_pipe(name, path, outputpath):
    sequences = list(SeqIO.parse(path, "fasta"))
    # Criar um dicionário para armazenar as sequências únicas
    unique_sequences = {}
    # Iterar pelas sequências do arquivo de entrada
    for sequence in sequences:
        # Verificar se a sequência já existe no dicionário de sequências únicas
        if str(sequence.seq) not in unique_sequences:
            # Se a sequência é única, armazená-la no dicionário
            unique_sequences[str(sequence.seq)] = sequence
    # Criar uma lista de sequências únicas
    unique_sequences_list = list(unique_sequences.values())
    # Salvar as sequências únicas em um arquivo de saída
    output_file_tmp = os.path.join(outputpath,f'{name}_NoPipe')
    SeqIO.write(unique_sequences_list, output_file_tmp, "fasta")
    return output_file_tmp
