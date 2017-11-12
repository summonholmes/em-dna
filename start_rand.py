from random import seed
from random import randint


def start_rand(len_list, len_seq, motif_width):
	seed(a=None)
	motif_start_pos = []
	for i in range(len_list):
		motif_start_pos.append(randint(0, (len_seq[i] - motif_width)))
	# print("Motif starting positions for each sequence: ", motif_start_pos)
	return motif_start_pos
