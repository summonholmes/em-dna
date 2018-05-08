from em_matrices import motif_actg_matrix, new_zeros_matrix


def add_pseudocounts(motif_width, score_matrix):
    score_matrix_pseudo = new_zeros_matrix(4, motif_width + 1)
    for i in range(4):
        for j in range(motif_width + 1):
            score_matrix_pseudo[i][j] = score_matrix[i][j] + 1
    return score_matrix_pseudo


def count_all_bases(fasta_file_seq, len_list):
    count_actg = {"count_a": 0, "count_c": 0, "count_t": 0, "count_g": 0}
    for i in range(len_list):
        count_actg["count_a"] += fasta_file_seq[i].count('A')
        count_actg["count_c"] += fasta_file_seq[i].count('C')
        count_actg["count_t"] += fasta_file_seq[i].count('T')
        count_actg["count_g"] += fasta_file_seq[i].count('G')
    return count_actg


def count_motif_bases(motif_width, len_list, motif):
    motif_pos_actg = motif_actg_matrix(len_list)
    for i in range(len_list):
        for j in range(motif_width):
            motif_pos_actg["motif_pos_a"][i].append(motif[i].find(
                'A', j, j + 1))
            motif_pos_actg["motif_pos_c"][i].append(motif[i].find(
                'C', j, j + 1))
            motif_pos_actg["motif_pos_t"][i].append(motif[i].find(
                'T', j, j + 1))
            motif_pos_actg["motif_pos_g"][i].append(motif[i].find(
                'G', j, j + 1))
    return motif_pos_actg


def init_background_motif_counts(len_list, count_bases, motif):
    background_actg = {
        "background_a": count_bases["count_a"],
        "background_c": count_bases["count_c"],
        "background_t": count_bases["count_t"],
        "background_g": count_bases["count_g"]
    }
    for i in range(len_list):
        background_actg["background_a"] -= motif[i].count('A')
        background_actg["background_c"] -= motif[i].count('C')
        background_actg["background_t"] -= motif[i].count('T')
        background_actg["background_g"] -= motif[i].count('G')
    return background_actg


def normalize_counts(motif_width, len_list, count_all_motif_bases):
    for i in count_all_motif_bases.keys():
        for j in range(len_list):
            for k in range(motif_width):
                count_all_motif_bases[i][j][k] += 1
    return count_all_motif_bases