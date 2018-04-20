from count_all_bases import count_all_bases
from init_motifs import init_motifs
from init_background_motif_counts import init_background_motif_counts
from count_motif_bases import count_motif_bases
from normalize_counts import normalize_counts
from counts_matrix import counts_matrix
from add_pseudocounts import add_pseudocounts
from freq_matrix import freq_matrix
from odds_matrix import odds_matrix
from log_odds_matrix import log_odds_matrix


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