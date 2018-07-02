from class_EM_Count import EM_Count
from math import log
from numpy import empty, round, vectorize


class EM_Matrix(EM_Count):
    def __init__(self, motif_width, fasta_file_seqs, count_bkgd_bases_dict,
                 motif_base_posit_freq_dict):
        self.motif_width = motif_width
        self.fasta_file_seqs = fasta_file_seqs
        self.count_bkgd_bases_dict = count_bkgd_bases_dict
        self.motif_base_posit_freq_dict = motif_base_posit_freq_dict
        self.init_em_matrix()  # Counts, freq, odds, then log odds
        self.counts_matrix_populate()
        self.freq_matrix_convert_bgs()
        self.freq_matrix_convert_cols()
        self.odds_matrix_populate()
        self.odds_matrix_to_log()

    def init_em_matrix(self):
        self.em_log_odds_matrix = empty(  # Will be overwritten
            (4, self.motif_width + 1))

    def counts_matrix_populate(self):  # Vectorize rather than iterate!
        self.em_log_odds_matrix[:, 0] = list(
            self.count_bkgd_bases_dict.values())  # By Column
        self.em_log_odds_matrix[0, 1:] = self.motif_base_posit_freq_dict[
            "motif_posits_a"]  # By Row
        self.em_log_odds_matrix[1, 1:] = self.motif_base_posit_freq_dict[
            "motif_posits_c"]  # By Row
        self.em_log_odds_matrix[2, 1:] = self.motif_base_posit_freq_dict[
            "motif_posits_t"]  # By Row
        self.em_log_odds_matrix[3, 1:] = self.motif_base_posit_freq_dict[
            "motif_posits_g"]  # By Row

    def freq_matrix_convert_bgs(self):
        self.em_log_odds_matrix[:, 0] = round(  # By Column
            (self.em_log_odds_matrix[:, 0] / sum(
                self.em_log_odds_matrix[:, 0])),
            decimals=3)

    def freq_matrix_convert_cols(self):
        col_totals = self.em_log_odds_matrix.sum(axis=0)
        for i in range(1, self.motif_width + 1):
            self.em_log_odds_matrix[:, i] = round(  # By Column
                (self.em_log_odds_matrix[:, i] / col_totals[i]),
                decimals=3)

    def odds_matrix_populate(self):
        for i in range(4):
            self.em_log_odds_matrix[i, 1:] = round(  # By Row
                (self.em_log_odds_matrix[i, 1:] /
                 self.em_log_odds_matrix[i][0]),
                decimals=3)  # Odds of getting each base in each position

    def odds_matrix_to_log(self):
        logger = vectorize(lambda x: round(log(x, 2), decimals=3))
        self.em_log_odds_matrix = logger(
            self.em_log_odds_matrix[:, 1:])  # log-odds matrix is NOT offset!
