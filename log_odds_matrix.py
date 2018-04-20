from math import log
from new_zeros_matrix import new_zeros_matrix


def log_odds_matrix(motif_width, score_matrix_odds):
    score_matrix_log_odds = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_log_odds[i][j] = round(
                log(score_matrix_odds[i][j], 2), 3)
    return score_matrix_log_odds