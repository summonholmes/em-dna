from start_rand import start_rand
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
from entropy_matrix import entropy_matrix
from sum_entropy import sum_entropy
from information import information
from em_motif import em_motif
from exp_max import exp_max


def start_the_em(max_bit_score_arr, final_record, motif_width, user_align,
                 user_iter, fasta_file_seq):
    print("\nNow preparing for E-M...")
    print(
        "M = Max Motif\nS = Max Score\nP = Max Position\n# = Max Sequence #"\
        "\nSS = Max Sum Scores\n\nPlease be patient..."
    )
    for i in range(user_align):
        len_list = len(fasta_file_seq)
        len_seq = [len(seq) for seq in fasta_file_seq]
        motif_start_pos = start_rand(len_list, len_seq, motif_width)
        count_bases = count_all_bases(fasta_file_seq, len_list)
        motif = init_motifs(motif_width, fasta_file_seq, len_list,
                            motif_start_pos)
        count_background_bases = init_background_motif_counts(
            len_list, count_bases, motif)
        count_all_motif_bases = count_motif_bases(motif_width, len_list, motif)
        normalize_count_motif_bases = normalize_counts(motif_width, len_list,
                                                       count_all_motif_bases)
        score_matrix = counts_matrix(motif_width, len_list,
                                     count_background_bases,
                                     normalize_count_motif_bases)
        score_matrix_pseudo = add_pseudocounts(motif_width, score_matrix)
        score_matrix_freq = freq_matrix(motif_width, count_background_bases,
                                        score_matrix_pseudo)
        score_matrix_odds = odds_matrix(motif_width, score_matrix_freq)
        score_matrix_log_odds = log_odds_matrix(motif_width, score_matrix_odds)
        score_matrix_entropy = entropy_matrix(motif_width, score_matrix_freq)
        score_matrix_entropy_sum = sum_entropy(score_matrix_entropy)
        score_matrix_information = information(score_matrix_entropy_sum)
        max_bit_score_arr.append(score_matrix_information)
        max_bit_score_arr.sort()
        em_motifs = em_motif(motif_width, fasta_file_seq, len_list, len_seq)
        finish = exp_max(user_iter, len_list, len_seq, fasta_file_seq,
                         motif_width, score_matrix_log_odds, em_motifs)
        final_record.append(finish)
    return final_record