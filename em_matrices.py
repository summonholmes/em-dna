from math import log


def blank_matrix(user_iter):
    blank_matrix = []
    for i in range(user_iter):
        blank_matrix.append([])
    return blank_matrix


def counts_matrix(motif_width, len_list, count_background_bases,
                  normalize_count_all_motif_bases):
    score_matrix = new_zeros_matrix(4, motif_width + 1)
    for i in range(len_list):
        for j in range(motif_width):
            for k in range(4):
                # for l in normalize_count_all_motif_bases.keys():
                score_matrix[k][0] = count_background_bases[k]
                if normalize_count_all_motif_bases[k][i][j] == j + 1:
                    score_matrix[k][j + 1] += 1
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
        score_matrix_freq[i][0] = round(
            count_background_bases[i] /
            (count_background_bases[0] + count_background_bases[1] +
             count_background_bases[2] + count_background_bases[3]), 3)
    score_matrix_freq = get_freq_matrix_score(score_matrix_freq, motif_width)
    return score_matrix_freq


def get_freq_matrix_score(score_matrix_freq, motif_width):
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


def motif_actg_matrix(len_list):
    motif_pos_actg = {
        "motif_pos_a": [],
        "motif_pos_c": [],
        "motif_pos_t": [],
        "motif_pos_g": []
    }
    for i in range(len_list):
        motif_pos_actg["motif_pos_a"].append([])
        motif_pos_actg["motif_pos_c"].append([])
        motif_pos_actg["motif_pos_t"].append([])
        motif_pos_actg["motif_pos_g"].append([])
    return motif_pos_actg


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


def odds_matrix(motif_width, score_matrix_freq):
    score_matrix_odds = new_zeros_matrix(4, motif_width)
    for i in range(4):
        for j in range(motif_width):
            score_matrix_odds[i][j] = round(
                score_matrix_freq[i][j + 1] / score_matrix_freq[i][0], 3)
    return score_matrix_odds


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
    return [
        round(j, 3) * -1 for j in [sum(i) for i in zip(*score_matrix_entropy)]
    ]


def sum_relative_entropy(score_matrix_relative_entropy):
    return [
        round(j, 3)
        for j in [sum(i) for i in zip(*score_matrix_relative_entropy)]
    ]