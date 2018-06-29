from random import seed, randint


class EM_Count:
    def __init__(self, motif_width, fasta_file_seqs):
        self.motif_width = motif_width
        self.fasta_file_seqs = fasta_file_seqs
        self.init_rand_motif_posits()
        self.init_total_base_counts_dict()
        self.count_all_bases()
        self.init_motifs_list()
        self.init_bkgd_motif_counts_dict()
        self.count_all_bkgd_bases()
        self.init_motif_base_posit_freq_dict()
        self.motif_base_posit_freq_dict_populate()

    def init_rand_motif_posits(self):
        seed(a=None)
        self.motif_start_posits = [
            randint(0, (j - self.motif_width))
            for j in [len(i) for i in self.fasta_file_seqs]
        ]

    def init_total_base_counts_dict(self):
        self.total_count_all_bases_dict = {
            "total_count_a": 0,
            "total_count_c": 0,
            "total_count_t": 0,
            "total_count_g": 0
        }

    def count_all_bases(self):
        for each_seq in self.fasta_file_seqs:
            self.total_count_all_bases_dict["total_count_a"] += each_seq.count(
                'A')
            self.total_count_all_bases_dict["total_count_c"] += each_seq.count(
                'C')
            self.total_count_all_bases_dict["total_count_t"] += each_seq.count(
                'T')
            self.total_count_all_bases_dict["total_count_g"] += each_seq.count(
                'G')

    def init_motifs_list(self):
        self.initial_rand_motifs_list = [
            self.fasta_file_seqs[i][self.motif_start_posits[i]:(
                self.motif_start_posits[i] + self.motif_width)]
            for i in range(len(self.fasta_file_seqs))
        ]

    def init_bkgd_motif_counts_dict(self):
        self.count_bkgd_bases_dict = {
            "bkgd_a": self.total_count_all_bases_dict["total_count_a"],
            "bkgd_c": self.total_count_all_bases_dict["total_count_c"],
            "bkgd_t": self.total_count_all_bases_dict["total_count_t"],
            "bkgd_g": self.total_count_all_bases_dict["total_count_g"]
        }

    def count_all_bkgd_bases(self):
        for each_motif in self.initial_rand_motifs_list:
            self.count_bkgd_bases_dict["bkgd_a"] -= each_motif.count('A')
            self.count_bkgd_bases_dict["bkgd_c"] -= each_motif.count('C')
            self.count_bkgd_bases_dict["bkgd_t"] -= each_motif.count('T')
            self.count_bkgd_bases_dict["bkgd_g"] -= each_motif.count('G')

    def init_motif_base_posit_freq_dict(self):
        self.motif_base_posit_freq_dict = {
            "motif_posits_a": [[] for i in range(len(self.fasta_file_seqs))],
            "motif_posits_c": [[] for i in range(len(self.fasta_file_seqs))],
            "motif_posits_t": [[] for i in range(len(self.fasta_file_seqs))],
            "motif_posits_g": [[] for i in range(len(self.fasta_file_seqs))]
        }

    def motif_base_posit_freq_dict_populate(self):  # Flagger fxn
        for i in range(len(self.fasta_file_seqs)):
            for j in range(self.motif_width):  # Add 1 to scale to 0
                self.motif_base_posit_freq_dict["motif_posits_a"][i].append(
                    self.initial_rand_motifs_list[i].find('A', j, j + 1) + 1)
                self.motif_base_posit_freq_dict["motif_posits_c"][i].append(
                    self.initial_rand_motifs_list[i].find('C', j, j + 1) + 1)
                self.motif_base_posit_freq_dict["motif_posits_t"][i].append(
                    self.initial_rand_motifs_list[i].find('T', j, j + 1) + 1)
                self.motif_base_posit_freq_dict["motif_posits_g"][i].append(
                    self.initial_rand_motifs_list[i].find('G', j, j + 1) + 1)