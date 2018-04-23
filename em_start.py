from em_counts import add_pseudocounts, count_all_bases, count_motif_bases, init_background_motif_counts, normalize_counts
from em_matrices import blank_matrix, counts_matrix, entropy_matrix, freq_matrix, information, log_odds_matrix, odds_matrix, sum_entropy
from em_motifs import em_motif, init_motifs, init_scores_pos_motifs, start_rand
from math import pow
from em_score import exp_max_get_max_pos_score


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


def start_the_em_prep(max_bit_score_arr, final_record, motif_width, user_align,
                      user_iter, fasta_file_seq):
    for i in range(user_align):
        print("Progress: {:2.1%}".format((i + 1) / user_align), end="\r")
        len_list = len(fasta_file_seq)
        len_seq = [len(seq) for seq in fasta_file_seq]
        count_bases_dict = start_the_em_prep_counts(
            len_list, len_seq, fasta_file_seq, motif_width)
        score_matrix_dict = start_the_em_prep_matrices(motif_width, len_list,
                                                       count_bases_dict)
        max_bit_score_arr.append(score_matrix_dict["score_matrix_information"])
        max_bit_score_arr.sort()
        em_motifs = em_motif(motif_width, fasta_file_seq, len_list, len_seq)
        finish = start_the_em_run(
            user_iter, len_list, len_seq, fasta_file_seq, motif_width,
            score_matrix_dict["score_matrix_log_odds"], em_motifs)
        final_record.append(finish)
    return final_record


def start_the_em_prep_counts(len_list, len_seq, fasta_file_seq, motif_width):
    motif_start_pos = start_rand(len_list, len_seq, motif_width)
    count_bases = count_all_bases(fasta_file_seq, len_list)
    motif = init_motifs(motif_width, fasta_file_seq, len_list, motif_start_pos)
    count_background_bases = init_background_motif_counts(
        len_list, count_bases, motif)
    count_all_motif_bases = count_motif_bases(motif_width, len_list, motif)
    normalize_count_motif_bases = normalize_counts(motif_width, len_list,
                                                   count_all_motif_bases)
    return {
        "count_background_bases": count_background_bases,
        "normalize_count_motif_bases": normalize_count_motif_bases
    }


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


def start_the_em_run(user_iter, len_list, len_seq, fasta_file_seq, motif_width,
                     score_matrix_log_odds, em_motifs):
    scores_pos_motifs = init_scores_pos_motifs(user_iter)
    max_score_pos_motif = exp_max_get_max_pos_score(
        user_iter, len_list, len_seq, motif_width, em_motifs,
        score_matrix_log_odds, scores_pos_motifs, fasta_file_seq)
    max_final_sco_seq_pos_mot = init_max_final_sco_seq_pos_mot(
        user_iter, max_score_pos_motif, scores_pos_motifs)
    for i in scores_pos_motifs.keys():
        scores_pos_motifs[i] = scores_pos_motifs[i][user_iter - 1]
    return {
        "max_final_sco_seq_pos_mot": max_final_sco_seq_pos_mot,
        "scores_pos_motifs": scores_pos_motifs
    }