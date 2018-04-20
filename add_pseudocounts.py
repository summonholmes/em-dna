from new_zeros_matrix import new_zeros_matrix


def add_pseudocounts(motif_width, score_matrix):
    score_matrix_pseudo = new_zeros_matrix(4, motif_width + 1)
    for i in range(4):
        for j in range(motif_width + 1):
            score_matrix_pseudo[i][j] = score_matrix[i][j]
    for i in range(4):
        for j in range(1, motif_width+1):
            score_matrix_pseudo[i][j] += 1
    return score_matrix_pseudo