def exp_max_score_log_odds(motif_width, list_em_motif, score_matrix_log_odds):
	score_em_motif = []
	for l in range(motif_width):
		if list_em_motif[l] == 'A':
			score_em_motif.append(score_matrix_log_odds[0][l])
		elif list_em_motif[l] == 'C':
			score_em_motif.append(score_matrix_log_odds[1][l])
		elif list_em_motif[l] == 'G':
			score_em_motif.append(score_matrix_log_odds[2][l])
		elif list_em_motif[l] == 'T':
			score_em_motif.append(score_matrix_log_odds[3][l])
		else:
			print("An error occurred in exp_max")
			exit(0)
	return score_em_motif
