def init_background_motif_counts(len_list, count_bases, motif):
    count_background_bases = []
    background_a = count_bases[0]
    background_c = count_bases[1]
    background_g = count_bases[2]
    background_t = count_bases[3]
    for i in range(len_list):
        background_a -= motif[i].count('A')
        background_c -= motif[i].count('C')
        background_g -= motif[i].count('G')
        background_t -= motif[i].count('T')
    count_background_bases.append(background_a)
    count_background_bases.append(background_c)
    count_background_bases.append(background_g)
    count_background_bases.append(background_t)

    return count_background_bases
