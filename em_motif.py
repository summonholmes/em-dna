def em_motif(motif_width, fasta_file_seq, len_list, len_seq):
    em_seq_motif = []
    for i in range(len_list):
        em_seq_motif.append([])
        for j in range(len_seq[i] - motif_width):
            em_seq_motif[i].append(0)
    for i in range(len_list):
        for j in range(len_seq[i] - motif_width):
            x = fasta_file_seq[i]
            z = j + motif_width
            em_seq_motif[i][j] = x[j:z]
    return em_seq_motif