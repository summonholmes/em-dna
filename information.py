def information(score_matrix_entropy_sum):
    score_matrix_information = [
        round(2 - i, 3) for i in score_matrix_entropy_sum
    ]
    return max(score_matrix_information)