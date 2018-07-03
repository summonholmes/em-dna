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
        seed(a=None)
        self.motif_start_posits = [
            randint(0, (j - self.motif_width))
            for j in [len(i) for i in self.fasta_file_seqs]
        ]  # Starting positions are always random

    def init_total_base_counts_dict(self):
        self.total_count_all_bases_dict = {
            "total_count_a": 0,
            "total_count_c": 0,
            "total_count_t": 0,
            "total_count_g": 0
        }

    def count_all_bases(self):
        all_bases = "".join(self.fasta_file_seqs)
        self.total_count_all_bases_dict[
            "total_count_a"] += all_bases.count('A')
        self.total_count_all_bases_dict[
            "total_count_c"] += all_bases.count('C')
        self.total_count_all_bases_dict[
            "total_count_t"] += all_bases.count('T')
        self.total_count_all_bases_dict[
            "total_count_g"] += all_bases.count('G')

    def init_motifs_list(self):
        self.initial_rand_motifs_list = list(  # Isolate motifs
            map(lambda x, y: x[y:y + self.motif_width], self.fasta_file_seqs,
                self.motif_start_posits))

    def init_bkgd_motif_counts_dict(self):
        self.count_bkgd_bases_dict = {
            "bkgd_a": self.total_count_all_bases_dict["total_count_a"],
            "bkgd_c": self.total_count_all_bases_dict["total_count_c"],
            "bkgd_t": self.total_count_all_bases_dict["total_count_t"],
            "bkgd_g": self.total_count_all_bases_dict["total_count_g"]
        }  # All bases minus those contained in any motif

    def count_all_bkgd_bases(self):  # Last integer dict
        self.all_motif_bases = "".join(self.initial_rand_motifs_list)
        self.count_bkgd_bases_dict["bkgd_a"] -= self.all_motif_bases.count('A')
        self.count_bkgd_bases_dict["bkgd_c"] -= self.all_motif_bases.count('C')
        self.count_bkgd_bases_dict["bkgd_t"] -= self.all_motif_bases.count('T')
        self.count_bkgd_bases_dict["bkgd_g"] -= self.all_motif_bases.count('G')

    def init_motif_base_posit_freq_dict(self):  # Now move to matrix dicts
        self.motif_base_posit_freq_dict = {
            "motif_posits_a": ones(self.motif_width),
            "motif_posits_c": ones(self.motif_width),
            "motif_posits_t": ones(self.motif_width),
            "motif_posits_g": ones(self.motif_width)
        }  # 1s address normalization/pseudocounts, exp-fxn does not work with 0s

    def motif_base_posit_freq_dict_populate(self):
        for j in range(self.motif_width):  # Count bases in all motif positions
            self.motif_base_posit_freq_dict["motif_posits_a"][
                j] += self.all_motif_bases[j::self.motif_width].count('A')
            self.motif_base_posit_freq_dict["motif_posits_c"][
                j] += self.all_motif_bases[j::self.motif_width].count('C')
            self.motif_base_posit_freq_dict["motif_posits_t"][
                j] += self.all_motif_bases[j::self.motif_width].count('T')
            self.motif_base_posit_freq_dict["motif_posits_g"][
                j] += self.all_motif_bases[j::self.motif_width].count('G')
