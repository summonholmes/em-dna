from math import pow
from prep_em_matrices import blank_matrix
from exp_max import exp_max_get_max_pos_score
from init_motifs import init_scores_pos_motifs


def init_max_final_sco_seq_pos_mot(user_iter, max_score_pos_motif,
                                   scores_pos_motifs):
    max_final_sco_seq_pos_mot = {
        "max_final_score":
        max(max_score_pos_motif["max_scores"][user_iter - 1])
    }
    max_final_sco_seq_pos_mot["max_final_sequence"] = scores_pos_motifs[
        "max_scores"][user_iter - 1].index(
            max_final_sco_seq_pos_mot["max_final_score"])
    max_final_sco_seq_pos_mot["max_final_position"] = scores_pos_motifs[
        "max_pos"][user_iter
                   - 1][max_final_sco_seq_pos_mot["max_final_sequence"]]
    max_final_sco_seq_pos_mot["max_final_motif"] = scores_pos_motifs[
        "max_motifs"][user_iter
                      - 1][max_final_sco_seq_pos_mot["max_final_sequence"]]
    max_final_sco_seq_pos_mot["sum_score_max_motif"] = sum(
        scores_pos_motifs["max_scores"][user_iter - 1])
    return max_final_sco_seq_pos_mot


def start_the_em_run(user_iter, len_list, len_seq, fasta_file_seq, motif_width,
                     score_matrix_log_odds, em_motifs):
    scores_pos_motifs = init_scores_pos_motifs(user_iter)
    max_score_pos_motif = exp_max_get_max_pos_score(
        user_iter, len_list, len_seq, motif_width, em_motifs,
        score_matrix_log_odds, scores_pos_motifs, fasta_file_seq)
    max_final_sco_seq_pos_mot = init_max_final_sco_seq_pos_mot(
        user_iter, max_score_pos_motif, scores_pos_motifs)
    return [[
        max_final_sco_seq_pos_mot["max_final_motif"],
        max_final_sco_seq_pos_mot["max_final_score"],
        max_final_sco_seq_pos_mot["max_final_position"],
        max_final_sco_seq_pos_mot["max_final_sequence"],
        max_final_sco_seq_pos_mot["sum_score_max_motif"],
    ], [
        scores_pos_motifs["max_pos"][user_iter - 1],
        scores_pos_motifs["max_scores"][user_iter - 1],
        scores_pos_motifs["max_motifs"][user_iter - 1]
    ]]