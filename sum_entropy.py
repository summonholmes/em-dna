def sum_entropy(score_matrix_entropy):
    score_matrix_entropy_sum = [sum(i) for i in zip(*score_matrix_entropy)]
    for i in score_matrix_entropy_sum:
        i = round(i, 3) * -1
    return score_matrix_entropy_sum