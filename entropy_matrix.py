from math import log
# from pprint import pprint


def entropy_matrix(motif_width, score_matrix_freq):
	score_matrix_entropy = []
	for i in range(4):
		score_matrix_entropy.append([])
		for j in range(motif_width):
			score_matrix_entropy[i].append(0)
	for i in range(4):
		for j in range(motif_width):
			score_matrix_entropy[i][j] = round(score_matrix_freq[i][j+1] * log(score_matrix_freq[i][j+1], 2), 3)
	# print("\nEntropy Matrix: ")
	# pprint(score_matrix_entropy)
	return score_matrix_entropy
