from numpy import ones
from random import randint, seed


class EM_Count:
    # Count the data and motifs, ensure randomness
    seed(a=None)

    def __init__(self, motif_width, fasta_file_seqs, total_bases_counts):
        self.motif_width = motif_width
        self.fasta_file_seqs = fasta_file_seqs
        self.total_bases_counts = total_bases_counts
        self.init_rand_posits()
        self.init_em_motifs()
        self.merge_all_motif_bases()
        self.init_bkgd_counts()
        self.count_all_bkgd_bases()
        self.init_motif_posits_freqs()
        self.count_motif_posits_freqs()

    def init_rand_posits(self):
        # Starting positions are always random
        self.motif_start_posits = [
            randint(0, (len_seq - self.motif_width))
            for len_seq in [len(seq) for seq in self.fasta_file_seqs]
        ]

    def init_em_motifs(self):
        # Isolate motifs
        self.em_motifs = [
            seq[pos:pos + self.motif_width]
            for seq, pos in zip(self.fasta_file_seqs, self.motif_start_posits)
        ]

    def merge_all_motif_bases(self):
        # Merge all the bases & motif bases
        self.all_motif_bases = ''.join(self.em_motifs)

    def init_bkgd_counts(self):
        # All bases minus those contained in any motif
        self.total_bkgd_counts = {
            'A': self.total_bases_counts['A'],
            'C': self.total_bases_counts['C'],
            'T': self.total_bases_counts['T'],
            'G': self.total_bases_counts['G']
        }

    def count_all_bkgd_bases(self):
        # All bases not in motifs
        self.total_bkgd_counts['A'] -= self.all_motif_bases.count('A')
        self.total_bkgd_counts['C'] -= self.all_motif_bases.count('C')
        self.total_bkgd_counts['T'] -= self.all_motif_bases.count('T')
        self.total_bkgd_counts['G'] -= self.all_motif_bases.count('G')

    def init_motif_posits_freqs(self):
        # Now prep for counts matrix
        self.motif_posits_freqs = {
            'A': ones(self.motif_width),
            'C': ones(self.motif_width),
            'T': ones(self.motif_width),
            'G': ones(self.motif_width)
        }
        # 1s address pseudocounts, exp-fxn does not work with 0s

    def count_motif_posits_freqs(self):
        # Count bases in all motif positions
        for i in range(self.motif_width):
            self.motif_posits_freqs['A'][i] += self.all_motif_bases[
                i::self.motif_width].count('A')
            self.motif_posits_freqs['C'][i] += self.all_motif_bases[
                i::self.motif_width].count('C')
            self.motif_posits_freqs['T'][i] += self.all_motif_bases[
                i::self.motif_width].count('T')
            self.motif_posits_freqs['G'][i] += self.all_motif_bases[
                i::self.motif_width].count('G')
