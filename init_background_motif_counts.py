def init_background_motif_counts(len_list, count_bases, motif):
    background_actg = {
        "background_a": count_bases[0],
        "background_c": count_bases[1],
        "background_t": count_bases[2],
        "background_g": count_bases[3]
    }
    for i in range(len_list):
        background_actg["background_a"] -= motif[i].count('A')
        background_actg["background_c"] -= motif[i].count('C')
        background_actg["background_t"] -= motif[i].count('T')
        background_actg["background_g"] -= motif[i].count('G')
    count_background_bases = [
        background_actg["background_a"], background_actg["background_c"],
        background_actg["background_t"], background_actg["background_g"]
    ]
    return count_background_bases