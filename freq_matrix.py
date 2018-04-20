from new_zeros_matrix import new_zeros_matrix


def freq_matrix(motif_width, count_background_bases, score_matrix_pseudo):
    score_matrix_freq = new_zeros_matrix(4, motif_width + 1)
    for i in range(4):
        for j in range(motif_width + 1):
            score_matrix_freq[i][j] = score_matrix_pseudo[i][j]
    score_matrix_freq[0][0] = round(
        count_background_bases[0] /
        (count_background_bases[0] + count_background_bases[1] +
         count_background_bases[2] + count_background_bases[3]), 3)
    score_matrix_freq[1][0] = round(
        count_background_bases[1] /
        (count_background_bases[0] + count_background_bases[1] +
         count_background_bases[2] + count_background_bases[3]), 3)
    score_matrix_freq[2][0] = round(
        count_background_bases[2] /
        (count_background_bases[0] + count_background_bases[1] +
         count_background_bases[2] + count_background_bases[3]), 3)
    score_matrix_freq[3][0] = round(
        count_background_bases[3] /
        (count_background_bases[0] + count_background_bases[1] +
         count_background_bases[2] + count_background_bases[3]), 3)
    col_totals = [sum(x) for x in zip(*score_matrix_freq)]
    for i in range(4):
        for j in range(1, motif_width + 1):
            score_matrix_freq[i][j] = score_matrix_freq[i][j] / col_totals[j]
            score_matrix_freq[i][j] = round(score_matrix_freq[i][j], 3)
    return score_matrix_freq