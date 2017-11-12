# from pprint import pprint


def information(motif_width, score_matrix_entropy_sum):
	score_matrix_information = []
	for i in range(motif_width):
		score_matrix_information.append(round(2 - score_matrix_entropy_sum[i], 3))
	# print("\nInformation: ")
	# pprint(score_matrix_information)
	return max(score_matrix_information)
