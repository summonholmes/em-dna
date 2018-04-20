def count_all_bases(fasta_file_seq, len_list):
    count_actg = {"count_a": 0, "count_c": 0, "count_t": 0, "count_g": 0}
    for i in range(len_list):
        count_actg["count_a"] += fasta_file_seq[i].count('A')
        count_actg["count_c"] += fasta_file_seq[i].count('C')
        count_actg["count_t"] += fasta_file_seq[i].count('T')
        count_actg["count_g"] += fasta_file_seq[i].count('G')
    count_bases = [
        count_actg["count_a"], count_actg["count_c"], count_actg["count_t"],
        count_actg["count_g"]
    ]
    return count_bases