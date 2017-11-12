def count_all_bases(fasta_file_seq, len_list):
	count_bases = []
	count_a = 0
	count_c = 0
	count_t = 0
	count_g = 0
	for i in range(len_list):
		count_a += fasta_file_seq[i].count('A')
		count_c += fasta_file_seq[i].count('C')
		count_t += fasta_file_seq[i].count('T')
		count_g += fasta_file_seq[i].count('G')
	# print("\nTotal base counts: ")
	# print('A: ', count_a)
	# print('C: ', count_c)
	# print('T: ', count_t)
	# print('G: ', count_g)
	count_bases.append(count_a)
	count_bases.append(count_c)
	count_bases.append(count_g)
	count_bases.append(count_t)
	return count_bases
