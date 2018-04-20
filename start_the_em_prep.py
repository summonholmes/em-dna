from counts import add_pseudocounts
from init_motifs import em_motif
from start_the_em_prep_counts import start_the_em_prep_counts
from start_the_em_prep_matrices import start_the_em_prep_matrices
from start_the_em_run import start_the_em_run


def start_the_em_prep(max_bit_score_arr, final_record, motif_width, user_align,
                      user_iter, fasta_file_seq):
    print("\nNow preparing for E-M...")
    print(
        "M = Max Motif\nS = Max Score\nP = Max Position\n# = Max Sequence #"\
        "\nSS = Max Sum Scores\n\nPlease be patient..."
    )
    for i in range(user_align):
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