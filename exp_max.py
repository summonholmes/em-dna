from math import pow
from exp_max_final_max import exp_max_final_max
from exp_max_get_max_pos_score import exp_max_get_max_pos_score
from exp_max_score_log_odds import exp_max_score_log_odds
from exp_max_matrices import exp_max_matrices


def exp_max(user_iter, len_list, len_seq, fasta_file_seq, motif_width,
            score_matrix_log_odds, em_motifs):
    final_max_score_pos_motif = {
        "final_max_score": exp_max_final_max(user_iter),
        "final_max_pos": exp_max_final_max(user_iter),
        "final_max_motif": exp_max_final_max(user_iter)
    }
    final_max_score_pos_motif = exp_max_get_max_pos_score(
        user_iter, len_list, len_seq, motif_width, em_motifs,
        score_matrix_log_odds, final_max_score_pos_motif, fasta_file_seq)
    max_final_score = max(
        final_max_score_pos_motif["final_max_score"][user_iter - 1])
    max_final_sequence = final_max_score_pos_motif["final_max_score"][
        user_iter - 1].index(max_final_score)
    max_final_position = final_max_score_pos_motif["final_max_pos"][
        user_iter - 1][max_final_sequence]
    max_final_motif = final_max_score_pos_motif["final_max_motif"][
        user_iter - 1][max_final_sequence]
    sum_score_max_motif = sum(
        final_max_score_pos_motif["final_max_score"][user_iter - 1])
    return [[
        max_final_motif, max_final_score, max_final_position,
        max_final_sequence, sum_score_max_motif
    ], [
        final_max_score_pos_motif["final_max_pos"][user_iter - 1],
        final_max_score_pos_motif["final_max_score"][user_iter - 1],
        final_max_score_pos_motif["final_max_motif"][user_iter - 1]
    ]]
