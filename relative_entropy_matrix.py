from math import log
from new_zeros_matrix import new_zeros_matrix


def relative_entropy_matrix(motif_width, score_matrix_freq):
    score_matrix_relative_entropy = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_relative_entropy[i][j] = round(
                score_matrix_freq[i][j + 1] * log(
                    score_matrix_freq[i][j + 1] / score_matrix_freq[i][0], 2),
                3)
    return score_matrix_relative_entropy