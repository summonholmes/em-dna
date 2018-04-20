from new_zeros_matrix import new_zeros_matrix_em


def em_motif(motif_width, fasta_file_seq, len_list, len_seq):
    em_seq_motif = new_zeros_matrix_em(len_list, len_seq, motif_width)
    for i in range(len_list):
        for j in range(len_seq[i] - motif_width):
            x = fasta_file_seq[i]
            z = j + motif_width
            em_seq_motif[i][j] = x[j:z]
    return em_seq_motif