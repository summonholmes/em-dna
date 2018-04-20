def sum_entropy(motif_width, score_matrix_entropy):
    score_matrix_entropy_sum = [sum(x) for x in zip(*score_matrix_entropy)]
    for i in range(motif_width):
        score_matrix_entropy_sum[i] = round(score_matrix_entropy_sum[i],
                                            3) * -1
    return score_matrix_entropy_sum