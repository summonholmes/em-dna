from numpy import ones
from random import randint, seed


class EM_Count:
    # Count the data and motifs, ensure randomness
    seed(a=None)

    def __init__(self, motif_width, fasta_file_seqs):
        self.motif_width = motif_width
        self.fasta_file_seqs = fasta_file_seqs
        self.init_rand_motif_posits()
        self.init_motifs_list()
        self.init_total_base_counts_dict()
        self.merge_all_bases()
        self.count_all_bases()
        self.init_bkgd_motif_counts_dict()
        self.count_all_bkgd_bases()
        self.init_motif_base_posit_freq_dict()
        self.motif_base_posit_freq_dict_populate()

    def init_rand_motif_posits(self):
        # Starting positions are always random
        self.motif_start_posits = [
            randint(0, (len_seq - self.motif_width))
            for len_seq in [len(seq) for seq in self.fasta_file_seqs]
        ]

    def init_motifs_list(self):
        # Isolate motifs
        self.initial_rand_motifs_list = [
            seq[pos:pos + self.motif_width]
            for seq, pos in zip(self.fasta_file_seqs, self.motif_start_posits)
        ]

    def init_total_base_counts_dict(self):
        # Initialize at 0
        self.total_count_all_bases_dict = {
            "total_count_a": 0,
            "total_count_c": 0,
            "total_count_t": 0,
            "total_count_g": 0
        }

    def merge_all_bases(self):
        # Merge all the bases & motif bases
        self.all_bases = "".join(self.fasta_file_seqs)
        self.all_motif_bases = "".join(self.initial_rand_motifs_list)

    def count_all_bases(self):
        # All bases
        self.total_count_all_bases_dict[
            "total_count_a"] += self.all_bases.count('A')
        self.total_count_all_bases_dict[
            "total_count_c"] += self.all_bases.count('C')
        self.total_count_all_bases_dict[
            "total_count_t"] += self.all_bases.count('T')
        self.total_count_all_bases_dict[
            "total_count_g"] += self.all_bases.count('G')

    def init_bkgd_motif_counts_dict(self):
        # All bases minus those contained in any motif
        self.count_bkgd_bases_dict = {
            "bkgd_a": self.total_count_all_bases_dict["total_count_a"],
            "bkgd_c": self.total_count_all_bases_dict["total_count_c"],
            "bkgd_t": self.total_count_all_bases_dict["total_count_t"],
            "bkgd_g": self.total_count_all_bases_dict["total_count_g"]
        }

    def count_all_bkgd_bases(self):
        # All bases not in motifs
        self.count_bkgd_bases_dict["bkgd_a"] -= self.all_motif_bases.count('A')
        self.count_bkgd_bases_dict["bkgd_c"] -= self.all_motif_bases.count('C')
        self.count_bkgd_bases_dict["bkgd_t"] -= self.all_motif_bases.count('T')
        self.count_bkgd_bases_dict["bkgd_g"] -= self.all_motif_bases.count('G')

    def init_motif_base_posit_freq_dict(self):
        # Now prep for counts matrix
        self.motif_base_posit_freq_dict = {
            "motif_posits_a": ones(self.motif_width),
            "motif_posits_c": ones(self.motif_width),
            "motif_posits_t": ones(self.motif_width),
            "motif_posits_g": ones(self.motif_width)
        }
        # 1s address pseudocounts, exp-fxn does not work with 0s

    def motif_base_posit_freq_dict_populate(self):
        # Count bases in all motif positions
        for i in range(self.motif_width):
            self.motif_base_posit_freq_dict["motif_posits_a"][
                i] += self.all_motif_bases[i::self.motif_width].count('A')
            self.motif_base_posit_freq_dict["motif_posits_c"][
                i] += self.all_motif_bases[i::self.motif_width].count('C')
            self.motif_base_posit_freq_dict["motif_posits_t"][
                i] += self.all_motif_bases[i::self.motif_width].count('T')
            self.motif_base_posit_freq_dict["motif_posits_g"][
                i] += self.all_motif_bases[i::self.motif_width].count('G')
