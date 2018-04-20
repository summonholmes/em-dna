from Bio import SeqIO


def check_if_fasta(fasta_file_seq, user_fasta_path):
    try:
        for record in SeqIO.parse(user_fasta_path, "fasta"):
            fasta_file_seq.append(str(record.seq))
        return fasta_file_seq
    except IOError:
        print("Unable to open file.  Try again.")
        exit(0)