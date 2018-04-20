from check_if_int import check_if_int
from input_start_align import input_start_align
from input_start_iter import input_start_iter
from input_start_fasta import input_start_fasta
from start_the_em import start_the_em
from pprint import pprint


def main():
    print("Welcome to a Python implementation of Expectation-Maximization")
    max_bit_score_arr = []
    final_record = []
    motif_width = input("Please specify the width of the motif: ")
    motif_width = check_if_int(motif_width)
    user_align = input_start_align()
    user_iter = input_start_iter()
    fasta_file_seq = input_start_fasta()
    em_results = start_the_em(max_bit_score_arr, final_record, motif_width,
                              user_align, user_iter, fasta_file_seq)
    final_motif = max(em_results, key=lambda item: item[0][4])
    print("\nDONE!  First horizontal array: max information.")
    print("Vertical arrays: positions, scores, and motifs. ")
    print("\nMAX: Motif\tScore\t Pos\tSeq\tSumScore")
    pprint(final_motif)


main()