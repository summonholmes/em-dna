def odds_matrix(motif_width, score_matrix_freq):
    score_matrix_odds = []
    for i in range(4):
        score_matrix_odds.append([])
        for j in range(motif_width):
            score_matrix_odds[i].append(0)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_odds[i][j] = round(
                score_matrix_freq[i][j + 1] / score_matrix_freq[i][0], 3)
    return score_matrix_odds