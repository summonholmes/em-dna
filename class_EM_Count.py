from numpy import ones
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
        seed(a=123)
        self.motif_start_posits = [
            randint(0, (j - self.motif_width))
            for j in [len(i) for i in self.fasta_file_seqs]
        ]

    def init_total_base_counts_dict(self):  # First matrix dict
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
        self.initial_rand_motifs_list = list(
            map(lambda x, y: x[y:y + self.motif_width], self.fasta_file_seqs,
                self.motif_start_posits))

    def init_bkgd_motif_counts_dict(self):
        self.count_bkgd_bases_dict = {
            "bkgd_a": self.total_count_all_bases_dict["total_count_a"],
            "bkgd_c": self.total_count_all_bases_dict["total_count_c"],
            "bkgd_t": self.total_count_all_bases_dict["total_count_t"],
            "bkgd_g": self.total_count_all_bases_dict["total_count_g"]
        }

    def count_all_bkgd_bases(self):  # Last integer dict
        for each_motif in self.initial_rand_motifs_list:
            self.count_bkgd_bases_dict["bkgd_a"] -= each_motif.count('A')
            self.count_bkgd_bases_dict["bkgd_c"] -= each_motif.count('C')
            self.count_bkgd_bases_dict["bkgd_t"] -= each_motif.count('T')
            self.count_bkgd_bases_dict["bkgd_g"] -= each_motif.count('G')

    def init_motif_base_posit_freq_dict(self):  # Now move to matrix dicts
        self.motif_base_posit_freq_dict = {  # Ones address pseudocounts
            "motif_posits_a": ones(self.motif_width),
            "motif_posits_c": ones(self.motif_width),
            "motif_posits_t": ones(self.motif_width),
            "motif_posits_g": ones(self.motif_width)
        }

    def motif_base_posit_freq_dict_populate(self):  # Count motif positions
        for each_motif in self.initial_rand_motifs_list:
            for j in range(self.motif_width):
                self.motif_base_posit_freq_dict["motif_posits_a"][
                    j] += each_motif[j].count('A')
                self.motif_base_posit_freq_dict["motif_posits_c"][
                    j] += each_motif[j].count('C')
                self.motif_base_posit_freq_dict["motif_posits_t"][
                    j] += each_motif[j].count('T')
                self.motif_base_posit_freq_dict["motif_posits_g"][
                    j] += each_motif[j].count('G')