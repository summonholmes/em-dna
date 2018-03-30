def information(motif_width, score_matrix_entropy_sum):
	score_matrix_information = []
	for i in range(motif_width):
		score_matrix_information.append(round(2 - score_matrix_entropy_sum[i], 3))

	return max(score_matrix_information)
