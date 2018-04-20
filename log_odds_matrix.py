from math import log


def log_odds_matrix(motif_width, score_matrix_odds):
    score_matrix_log_odds = []
    for i in range(4):
        score_matrix_log_odds.append([])
        for j in range(motif_width):
            score_matrix_log_odds[i].append(0)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_log_odds[i][j] = round(
                log(score_matrix_odds[i][j], 2), 3)
    return score_matrix_log_odds