def new_zeros_matrix(amount, motif_width):
    zeros_matrix = []
    for i in range(amount):
        zeros_matrix.append([])
        for j in range(motif_width):
            zeros_matrix[i].append(0)
    return zeros_matrix

def new_zeros_matrix_em(len_list, len_seq, motif_width):
    em_seq_motif = []
    for i in range(len_list):
        em_seq_motif.append([])
        for j in range(len_seq[i] - motif_width):
            em_seq_motif[i].append(0)
    return em_seq_motif