from math import pow
from exp_max_score_log_odds import exp_max_score_log_odds
from exp_max_matrices import exp_max_matrices


def exp_max(user_iter, len_list, len_seq, fasta_file_seq, motif_width, score_matrix_log_odds, em_motifs):
    final_max_score = []
    for i in range(user_iter):
        final_max_score.append([])
    final_max_pos = []
    for i in range(user_iter):
        final_max_pos.append([])
    final_motif = []
    for i in range(user_iter):
        final_motif.append([])

    for i in range(user_iter):
        for j in range(len_list):
            max_score = -1
            max_pos = -1
            last_motif_pos = len_seq[j] - motif_width
            for k in range(last_motif_pos):
                list_em_motif = list(em_motifs[j][k])
                score_em_motif = exp_max_score_log_odds(
                    motif_width, list_em_motif, score_matrix_log_odds)
                sum_score_em_motif = sum(score_em_motif)
                power_score_em_motif = round(pow(2, sum_score_em_motif), 3)
                if power_score_em_motif > max_score:
                    max_pos = k
                    max_score = power_score_em_motif
            final_max_pos[i].append(max_pos)
            final_max_score[i].append(max_score)
            final_motif[i].append(
                fasta_file_seq[j][max_pos:(max_pos + motif_width)])
        score_matrix_log_odds = exp_max_matrices(
            len_list, fasta_file_seq, motif_width, final_max_pos[i])
    max_final_score = max(final_max_score[user_iter-1])
    max_final_sequence = final_max_score[user_iter-1].index(max_final_score)
    max_final_position = final_max_pos[user_iter-1][max_final_sequence]
    max_final_motif = final_motif[user_iter-1][max_final_sequence]
    sum_score_max_motif = sum(final_max_score[user_iter-1])
	
    return [[max_final_motif, max_final_score, max_final_position, max_final_sequence, sum_score_max_motif], [final_max_pos[user_iter-1], final_max_score[user_iter-1], final_motif[user_iter-1]]]
