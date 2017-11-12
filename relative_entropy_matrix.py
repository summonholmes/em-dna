from math import log
# from pprint import pprint


def relative_entropy_matrix(motif_width, score_matrix_freq):
	score_matrix_relative_entropy = []
	for i in range(4):
		score_matrix_relative_entropy.append([])
		for j in range(motif_width):
			score_matrix_relative_entropy[i].append(0)
	for i in range(4):
		for j in range(motif_width):
			score_matrix_relative_entropy[i][j] = \
				round(score_matrix_freq[i][j + 1] * log(score_matrix_freq[i][j + 1]/score_matrix_freq[i][0], 2), 3)
	# print("\nRelative Entropy Matrix: ")
	# pprint(score_matrix_relative_entropy)
	return score_matrix_relative_entropy
