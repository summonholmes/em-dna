from random import randint, seed
from em_matrices import blank_matrix, new_em_zeros_matrix


def em_motif(motif_width, fasta_file_seq, len_list, len_seq):
    em_seq_motif = new_em_zeros_matrix(len_list, len_seq, motif_width)
    for i in range(len_list):
        for j in range(len_seq[i] - motif_width):
            x = fasta_file_seq[i]
            z = j + motif_width
            em_seq_motif[i][j] = x[j:z]
    return em_seq_motif


def init_background_motif_counts(len_list, count_bases, motif):
    background_actg = {
        "background_a": count_bases[0],
        "background_c": count_bases[1],
        "background_t": count_bases[2],
        "background_g": count_bases[3]
    }
    for i in range(len_list):
        background_actg["background_a"] -= motif[i].count('A')
        background_actg["background_c"] -= motif[i].count('C')
        background_actg["background_t"] -= motif[i].count('T')
        background_actg["background_g"] -= motif[i].count('G')
    return background_actg


def init_motifs(motif_width, fasta_file_seq, len_list, motif_start_pos):
    motif = []
    for i in range(len_list):
        x = fasta_file_seq[i]
        y = motif_start_pos[i]
        z = motif_start_pos[i] + motif_width
        motif.append(x[y:z])
    return motif


def init_scores_pos_motifs(user_iter):
    return {
        "max_scores": blank_matrix(user_iter),
        "max_pos": blank_matrix(user_iter),
        "max_motifs": blank_matrix(user_iter)
    }


def start_rand(len_list, len_seq, motif_width):
    seed(a=None)
    return [randint(0, (i - motif_width)) for i in len_seq]