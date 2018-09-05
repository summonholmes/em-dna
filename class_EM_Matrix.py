from class_EM_Count import EM_Count
from numpy import empty, fromiter, log2, sum


class EM_Matrix(EM_Count):
    # Matrix operations class
    def __init__(self, motif_width, fasta_file_seqs, count_bkgd_bases_dict,
                 motif_base_posit_freq_dict):
        self.motif_width = motif_width
        self.fasta_file_seqs = fasta_file_seqs
        self.count_bkgd_bases_dict = count_bkgd_bases_dict
        self.motif_base_posit_freq_dict = motif_base_posit_freq_dict
        self.init_em_matrix()
        self.counts_matrix_convert_bgs()
        self.counts_matrix_populate()
        self.freq_matrix_convert_bgs()
        self.freq_col_totals()
        self.freq_matrix_convert_cols()
        self.odds_matrix_populate()
        self.odds_matrix_to_log()

    def init_em_matrix(self):
        # Counts, freq, odds, then log odds
        self.em_log_odds_matrix = empty((4, self.motif_width + 1))
        # Will be overwritten

    def counts_matrix_convert_bgs(self):
        # By Column
        self.em_log_odds_matrix[:, 0] = fromiter(
            self.count_bkgd_bases_dict.values(), dtype=int)

    def counts_matrix_populate(self):
        # Vectorize rather than iterate
        for i, base in enumerate(('a', 'c', 't', 'g')):
            self.em_log_odds_matrix[i, 1:] = self.motif_base_posit_freq_dict[
                "motif_posits_%s" % base]  # By Row

    def freq_matrix_convert_bgs(self):
        # Turn background to freq
        self.em_log_odds_matrix[:, 0] = self.em_log_odds_matrix[:, 0] / sum(
            self.em_log_odds_matrix[:, 0])

    def freq_col_totals(self):
        # Sum columns
        self.col_totals = self.em_log_odds_matrix.sum(axis=0)

    def freq_matrix_convert_cols(self):
        # Turn columns to freq
        for i in range(1, self.motif_width + 1):
            self.em_log_odds_matrix[:, i] = \
                self.em_log_odds_matrix[:, i] / self.col_totals[i]

    def odds_matrix_populate(self):
        # Odds of getting each base in each position
        for i in range(4):
            self.em_log_odds_matrix[i, 1:] = (  # By Row
                self.em_log_odds_matrix[i, 1:] / self.em_log_odds_matrix[i][0])

    def odds_matrix_to_log(self):
        # Log-odds matrix is NOT offset!
        self.em_log_odds_matrix = log2(self.em_log_odds_matrix[:, 1:])
