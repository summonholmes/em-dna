from pprint import pprint


def init_score_matrix(motif_width):
	new_matrix = []
	for i in range(4):
		new_matrix.append([])
		for j in range(motif_width + 1):
			new_matrix[i].append(0)
	print("\nInitiate score blank matrix")
	for k in range(motif_width + 1):
		print(' ', k, end="")
	print('')
	print('A', 'C', 'G', 'T', pprint(new_matrix))
	return new_matrix
