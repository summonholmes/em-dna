from em_counts import add_pseudocounts, count_all_bases, count_motif_bases, init_background_motif_counts, normalize_counts
from em_motifs import init_motifs
from em_matrices import counts_matrix, freq_matrix, odds_matrix, log_odds_matrix


def exp_max_get_max_pos_score(user_iter, len_list, len_seq, motif_width,
                              em_motifs, score_matrix_log_odds,
                              scores_pos_motifs, fasta_file_seq):
    for i in range(user_iter):
        for j in range(len_list):
            last_motif_pos = len_seq[j] - motif_width
            max_pos_score = exp_max_pos_score_iter(
                last_motif_pos, motif_width, em_motifs, j,
                score_matrix_log_odds, -1, -1)
            scores_pos_motifs["max_pos"][i].append(max_pos_score["max_pos"])
            scores_pos_motifs["max_scores"][i].append(
                max_pos_score["max_score"])
            scores_pos_motifs["max_motifs"][i].append(
                fasta_file_seq[j][max_pos_score["max_pos"]:(
                    max_pos_score["max_pos"] + motif_width)])
        score_matrix_log_odds = exp_max_matrices(
            len_list, fasta_file_seq, motif_width,
            scores_pos_motifs["max_pos"][i])
    return scores_pos_motifs


def exp_max_matrices(len_list, fasta_file_seq, motif_width, final_max_pos):
    count_bases = count_all_bases(fasta_file_seq, len_list)
    motif = init_motifs(motif_width, fasta_file_seq, len_list, final_max_pos)
    count_background_bases = init_background_motif_counts(
        len_list, count_bases, motif)
    count_all_motif_bases = count_motif_bases(motif_width, len_list, motif)
    normalize_count_motif_bases = normalize_counts(motif_width, len_list,
                                                   count_all_motif_bases)
    score_matrix = counts_matrix(motif_width, len_list, count_background_bases,
                                 normalize_count_motif_bases)
    score_matrix_pseudo = add_pseudocounts(motif_width, score_matrix)
    score_matrix_freq = freq_matrix(motif_width, count_background_bases,
                                    score_matrix_pseudo)
    score_matrix_odds = odds_matrix(motif_width, score_matrix_freq)
    score_matrix_log_odds = log_odds_matrix(motif_width, score_matrix_odds)
    return score_matrix_log_odds


def exp_max_pos_score_iter(last_motif_pos, motif_width, em_motifs, j,
                           score_matrix_log_odds, max_score, max_pos):
    for k in range(last_motif_pos):
        score_em_motif = exp_max_score_log_odds(motif_width,
                                                list(em_motifs[j][k]),
                                                score_matrix_log_odds)
        sum_score_em_motif = sum(score_em_motif)
        power_score_em_motif = round(pow(2, sum_score_em_motif), 3)
        if power_score_em_motif > max_score:
            max_pos = k
            max_score = power_score_em_motif
    return {"max_pos": max_pos, "max_score": max_score}


def exp_max_score_log_odds(motif_width, list_em_motif, score_matrix_log_odds):
    score_em_motif = []
    for l in range(motif_width):
        if list_em_motif[l] == 'A':
            score_em_motif.append(score_matrix_log_odds[0][l])
        elif list_em_motif[l] == 'C':
            score_em_motif.append(score_matrix_log_odds[1][l])
        elif list_em_motif[l] == 'T':
            score_em_motif.append(score_matrix_log_odds[2][l])
        elif list_em_motif[l] == 'G':
            score_em_motif.append(score_matrix_log_odds[3][l])
        else:
            print("An error occurred in exp_max")
            exit(0)
    return score_em_motif