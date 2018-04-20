def counts_matrix(motif_width, len_list, count_background_bases,
                  normalize_count_all_motif_bases):
    score_matrix = []
    for i in range(4):
        score_matrix.append([])
        for j in range(motif_width + 1):
            score_matrix[i].append(0)
    score_matrix[0][0] = count_background_bases[0]
    score_matrix[1][0] = count_background_bases[1]
    score_matrix[2][0] = count_background_bases[2]
    score_matrix[3][0] = count_background_bases[3]
    for i in range(len_list):
        for j in range(motif_width):
            if normalize_count_all_motif_bases[0][i][j] == j + 1:
                score_matrix[0][j + 1] += 1
            if normalize_count_all_motif_bases[1][i][j] == j + 1:
                score_matrix[1][j + 1] += 1
            if normalize_count_all_motif_bases[2][i][j] == j + 1:
                score_matrix[2][j + 1] += 1
            if normalize_count_all_motif_bases[3][i][j] == j + 1:
                score_matrix[3][j + 1] += 1
    return score_matrix