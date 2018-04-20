from counts import add_pseudocounts
from prep_em_matrices import counts_matrix, entropy_matrix, freq_matrix, information, log_odds_matrix, odds_matrix, sum_entropy


def start_the_em_prep_matrices(motif_width, len_list, count_bases_dict):
    score_matrix = counts_matrix(
        motif_width, len_list, count_bases_dict["count_background_bases"],
        count_bases_dict["normalize_count_motif_bases"])
    score_matrix_pseudo = add_pseudocounts(motif_width, score_matrix)
    score_matrix_freq = freq_matrix(motif_width,
                                    count_bases_dict["count_background_bases"],
                                    score_matrix_pseudo)
    score_matrix_odds = odds_matrix(motif_width, score_matrix_freq)
    score_matrix_log_odds = log_odds_matrix(motif_width, score_matrix_odds)
    score_matrix_entropy = entropy_matrix(motif_width, score_matrix_freq)
    score_matrix_entropy_sum = sum_entropy(score_matrix_entropy)
    score_matrix_information = information(score_matrix_entropy_sum)
    return {
        "score_matrix_information": score_matrix_information,
        "score_matrix_log_odds": score_matrix_log_odds
    }