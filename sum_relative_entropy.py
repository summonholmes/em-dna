def sum_relative_entropy(score_matrix_relative_entropy):
    score_matrix_relative_entropy_sum = [
        sum(i) for i in zip(*score_matrix_relative_entropy)
    ]
    for i in score_matrix_relative_entropy_sum:
        i = round(i, 3)
    return score_matrix_relative_entropy_sum