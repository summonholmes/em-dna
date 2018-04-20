from Bio import SeqIO


def check_if_int(input_start_align):
    try:
        input_start_align = int(input_start_align)
        return input_start_align
    except:
        print("You didn't provide an integer ya turkey!")
        exit(0)


def check_if_fasta(fasta_file_seq, user_fasta_path):
    try:
        for record in SeqIO.parse(user_fasta_path, "fasta"):
            fasta_file_seq.append(str(record.seq))
        return fasta_file_seq
    except IOError:
        print("Unable to open file.  Try again.")
        exit(0)