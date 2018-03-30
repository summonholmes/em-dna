def count_motif_bases(motif_width, len_list, motif):
    count_all_motif_bases = []
    motif_pos_a = []
    motif_pos_c = []
    motif_pos_g = []
    motif_pos_t = []
    for i in range(len_list):
        motif_pos_a.append([])
        motif_pos_c.append([])
        motif_pos_g.append([])
        motif_pos_t.append([])
    for i in range(len_list):
        for j in range(motif_width):
            motif_pos_a[i].append(motif[i].find('A', j, j+1))
            motif_pos_c[i].append(motif[i].find('C', j, j+1))
            motif_pos_g[i].append(motif[i].find('G', j, j+1))
            motif_pos_t[i].append(motif[i].find('T', j, j+1))
    count_all_motif_bases.append(motif_pos_a)
    count_all_motif_bases.append(motif_pos_c)
    count_all_motif_bases.append(motif_pos_g)
    count_all_motif_bases.append(motif_pos_t)

    return count_all_motif_bases
