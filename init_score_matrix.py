from pprint import pprint
from new_zeros_matrix import new_zeros_matrix


def init_score_matrix(motif_width):
    new_matrix = new_zeros_matrix(4, motif_width + 1)
    return new_matrix