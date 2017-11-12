# from pprint import pprint


def add_pseudocounts(motif_width, score_matrix):
	score_matrix_pseudo = []
	for i in range(4):
		score_matrix_pseudo.append([])
		for j in range(motif_width + 1):
			score_matrix_pseudo[i].append(0)
	for i in range(4):
		for j in range(motif_width + 1):
			score_matrix_pseudo[i][j] = score_matrix[i][j]
	for i in range(4):
		for j in range(1, motif_width+1):
			score_matrix_pseudo[i][j] += 1
	# print("\nScore matrix with pseudocount + 1: ")
	# pprint(score_matrix_pseudo)
	return score_matrix_pseudo
