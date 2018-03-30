def sum_relative_entropy(motif_width, score_matrix_relative_entropy):
    score_matrix_relative_entropy_sum = [
        sum(x) for x in zip(*score_matrix_relative_entropy)]
    for i in range(motif_width):
        score_matrix_relative_entropy_sum[i] = round(
            score_matrix_relative_entropy_sum[i], 3)

    return score_matrix_relative_entropy_sum
