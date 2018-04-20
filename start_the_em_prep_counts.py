from counts import count_all_bases, count_motif_bases, init_background_motif_counts, normalize_counts
from init_motifs import init_motifs, start_rand


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