from math import log


def blank_matrix(user_iter):
    blank_matrix = []
    for i in range(user_iter):
        blank_matrix.append([])
    return blank_matrix


def counts_matrix(motif_width, len_list, count_background_bases,
                  normalize_count_all_motif_bases):
    score_matrix = new_zeros_matrix(4, motif_width + 1)
    for i in range(4):
        score_matrix[i][0] = count_background_bases[i]
    for i in range(len_list):
        for j in range(motif_width):
            if normalize_count_all_motif_bases[0][i][j] == j + 1:
                score_matrix[0][j + 1] += 1
            if normalize_count_all_motif_bases[1][i][j] == j + 1:
                score_matrix[1][j + 1] += 1
            if normalize_count_all_motif_bases[2][i][j] == j + 1:
                score_matrix[2][j + 1] += 1
            if normalize_count_all_motif_bases[3][i][j] == j + 1:
                score_matrix[3][j + 1] += 1
    return score_matrix


def entropy_matrix(motif_width, score_matrix_freq):
    score_matrix_entropy = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_entropy[i][j] = round(
                score_matrix_freq[i][j + 1] * log(score_matrix_freq[i][j + 1],
                                                  2), 3)
    return score_matrix_entropy


def freq_matrix(motif_width, count_background_bases, score_matrix_pseudo):
    score_matrix_freq = new_zeros_matrix(4, motif_width + 1)
    for i in range(4):
        for j in range(motif_width + 1):
            score_matrix_freq[i][j] = score_matrix_pseudo[i][j]
    for i in range(4):
        score_matrix_freq[i][0] = round(
            count_background_bases[i] /
            (count_background_bases[0] + count_background_bases[1] +
             count_background_bases[2] + count_background_bases[3]), 3)
    col_totals = [sum(i) for i in zip(*score_matrix_freq)]
    for i in range(4):
        for j in range(1, motif_width + 1):
            score_matrix_freq[i][j] = score_matrix_freq[i][j] / col_totals[j]
            score_matrix_freq[i][j] = round(score_matrix_freq[i][j], 3)
    return score_matrix_freq


def information(score_matrix_entropy_sum):
    score_matrix_information = [
        round(2 - i, 3) for i in score_matrix_entropy_sum
    ]
    return max(score_matrix_information)


def log_odds_matrix(motif_width, score_matrix_odds):
    score_matrix_log_odds = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_log_odds[i][j] = round(
                log(score_matrix_odds[i][j], 2), 3)
    return score_matrix_log_odds


def odds_matrix(motif_width, score_matrix_freq):
    score_matrix_odds = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_odds[i][j] = round(
                score_matrix_freq[i][j + 1] / score_matrix_freq[i][0], 3)
    return score_matrix_odds


def new_zeros_matrix(amount, motif_width):
    zeros_matrix = []
    for i in range(amount):
        zeros_matrix.append([])
        for j in range(motif_width):
            zeros_matrix[i].append(0)
    return zeros_matrix


def new_em_zeros_matrix(len_list, len_seq, motif_width):
    em_seq_motif = []
    for i in range(len_list):
        em_seq_motif.append([])
        for j in range(len_seq[i] - motif_width):
            em_seq_motif[i].append(0)
    return em_seq_motif


def relative_entropy_matrix(motif_width, score_matrix_freq):
    score_matrix_relative_entropy = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_relative_entropy[i][j] = round(
                score_matrix_freq[i][j + 1] * log(
                    score_matrix_freq[i][j + 1] / score_matrix_freq[i][0], 2),
                3)
    return score_matrix_relative_entropy


def sum_entropy(score_matrix_entropy):
    score_matrix_entropy_sum = [sum(i) for i in zip(*score_matrix_entropy)]
    for i in score_matrix_entropy_sum:
        i = round(i, 3) * -1
    return score_matrix_entropy_sum


def sum_relative_entropy(score_matrix_relative_entropy):
    score_matrix_relative_entropy_sum = [
        sum(i) for i in zip(*score_matrix_relative_entropy)
    ]
    for i in score_matrix_relative_entropy_sum:
        i = round(i, 3)
    return score_matrix_relative_entropy_sum