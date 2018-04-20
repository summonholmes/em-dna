from Bio import SeqIO
from check_if_fasta import check_if_fasta


def input_start_fasta():
    fasta_file_seq = []
    user_fasta_path = input("Please specify the path of the fasta file: ")
    fasta_file_seq = check_if_fasta(fasta_file_seq, user_fasta_path)
    return fasta_file_seq