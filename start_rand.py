from random import seed
from random import randint


def start_rand(len_list, len_seq, motif_width):
    seed(a=None)
    motif_start_pos = [randint(0, (i - motif_width)) for i in len_seq]
    return motif_start_pos