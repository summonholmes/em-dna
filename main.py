#!/usr/bin/env python
from input_start_align import input_start_align
from input_start_iter import input_start_iter
from input_start_fasta import input_start_fasta
from init_param_list import init_param_list
from init_param_seq import init_param_seq
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
from relative_entropy_matrix import relative_entropy_matrix
from sum_relative_entropy import sum_relative_entropy
from information import information
from em_motif import em_motif
from exp_max import exp_max
from pprint import pprint

#          -/oyddmdhs+:.                summonholmes@summonholmes
#      -odNMMMMMMMMNNmhy+-`             OS: Gentoo
#    -yNMMMMMMMMMMMNNNmmdhy+-           Kernel: x86_64 Linux 4.9.6-gentoo-r1
#  `omMMMMMMMMMMMMNmdmmmmddhhy/`        Uptime: 1h 22m
#  omMMMMMMMMMMMNhhyyyohmdddhhhdo`      Packages: 652
# .ydMMMMMMMMMMdhs++so/smdddhhhhdm+`    Shell: bash 4.3.48
#  oyhdmNMMMMMMMNdyooydmddddhhhhyhNd.   WM: OpenBox
#   :oyhhdNNMMMMMMMNNNmmdddhhhhhyymMh   WM Theme: 10x-Orange-G
#     .:+sydNMMMMMNNNmmmdddhhhhhhmMmy   GTK Theme: 10x-Orange-G [GTK2/3]
#        /mMMMMMMNNNmmmdddhhhhhmMNhs:   Icon Theme: 10x-Orange-G
#     `oNMMMMMMMNNNmmmddddhhdmMNhs+`    Font: Sans 10
#   `sNMMMMMMMMNNNmmmdddddmNMmhs/.      CPU: Intel Core i7 CPU L 640 @ 2.134GHz
#  /NMMMMMMMMNNNNmmmdddmNMNdso:`        RAM: 1685MiB / 7781MiB
# +MMMMMMMNNNNNmmmmdmNMNdso/-
# MNNNNNNNmmmmmNNMmhs+/-`
# /hMMNNNNNNNNMNdhs++/-`
# `/ohdmmddhys+++/:.`
#   `-//////:--.
#
# Python 3.4.5 - Biopython (import fasta), pprint (make arrays readable)
# Compiler: GCC 4.9.4 - retrieved using sys.version
# Making executable - Place all files in same directory, ensure the shebang line is at the
# 	top of this file, and run: '$ chmod +x main.py'.  In a terminal, type '$ path/to/file/main.py'


def main():
	print("Welcome to a Python implementation of Expectation-Maximization")

	# This is just my ordering of the process.
	max_bit_score_arr = []
	motif_width = int(input("Please specify the width of the motif: "))
	user_align = input_start_align()
	user_iter = input_start_iter()
	fasta_file_seq = input_start_fasta()
	final_record = []
	# Default 50 rounds
	for i in range(user_align):
		# Gather all necessary counts.
		len_list = init_param_list(fasta_file_seq)
		len_seq = init_param_seq(fasta_file_seq, len_list)
		motif_start_pos = start_rand(len_list, len_seq, motif_width)
		count_bases = count_all_bases(fasta_file_seq, len_list)
		motif = init_motifs(motif_width, fasta_file_seq, len_list, motif_start_pos)
		count_background_bases = init_background_motif_counts(len_list, count_bases, motif)
		count_all_motif_bases = count_motif_bases(motif_width, len_list, motif)
		normalize_count_motif_bases = normalize_counts(motif_width, len_list, count_all_motif_bases)

		# There is redundancy in the creation of additional matrices, because python treats them as "pointers."
		score_matrix = counts_matrix(motif_width, len_list, count_background_bases, normalize_count_motif_bases)
		score_matrix_pseudo = add_pseudocounts(motif_width, score_matrix)
		score_matrix_freq = freq_matrix(motif_width, count_background_bases, score_matrix_pseudo)
		score_matrix_odds = odds_matrix(motif_width, score_matrix_freq)
		score_matrix_log_odds = log_odds_matrix(motif_width, score_matrix_odds)
		score_matrix_entropy = entropy_matrix(motif_width, score_matrix_freq)
		score_matrix_entropy_sum = sum_entropy(motif_width, score_matrix_entropy)
		score_matrix_relative_entropy = relative_entropy_matrix(motif_width, score_matrix_freq)
		score_matrix_relative_entropy_sum = sum_relative_entropy(motif_width, score_matrix_relative_entropy)
		score_matrix_information = information(motif_width, score_matrix_entropy_sum)
		max_bit_score_arr.append(score_matrix_information)
		max_bit_score_arr.sort()
		em_motifs = em_motif(motif_width, fasta_file_seq, len_list, len_seq)
		# Each final iteration from the round is printed below.
		finish = exp_max(user_iter, len_list, len_seq, fasta_file_seq, motif_width, score_matrix_log_odds, em_motifs)
		print("\nROUND", i, "RESULTS: ")
		pprint(finish)
		final_record.append(finish)
	# Print the max values, then the rest of the set
	final_motif = max(final_record, key=lambda item: item[0][4])
	print("\nDONE!  First horizontal array contains the max information related to the max sum of score scores.")
	print("Vertical arrays are all positions, scores, and motifs corresponding to the max sum of scores ")
	print("\nMAX: Motif\tScore\t Pos\tSeq\tSumScore")
	pprint(final_motif)
	return 0
main()
