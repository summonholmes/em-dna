def init_motifs(motif_width, fasta_file_seq, len_list, motif_start_pos):
	motif = []
	for i in range(len_list):
		x = fasta_file_seq[i]
		y = motif_start_pos[i]
		z = motif_start_pos[i] + motif_width
		motif.append(x[y:z])
	# print("\nList the motifs from each sequence: ", motif)
	return motif
