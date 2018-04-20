from math import log
from new_zeros_matrix import new_zeros_matrix


def entropy_matrix(motif_width, score_matrix_freq):
    score_matrix_entropy = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_entropy[i][j] = round(
                score_matrix_freq[i][j + 1] * log(score_matrix_freq[i][j + 1],
                                                  2), 3)
    return score_matrix_entropy