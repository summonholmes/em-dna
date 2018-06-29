from class_EM_Count import EM_Count
from math import log
from numpy import ones, zeros


class EM_Matrix(EM_Count):
    def __init__(self, motif_width, fasta_file_seqs, count_bkgd_bases_dict,
                 motif_base_posit_freq_dict):
        self.motif_width = motif_width
        self.fasta_file_seqs = fasta_file_seqs
        self.count_bkgd_bases_dict = count_bkgd_bases_dict
        self.motif_base_posit_freq_dict = motif_base_posit_freq_dict
        self.init_freq_odds_matrices()
        self.counts_matrix_populate()
        self.freq_matrix_convert_bgs()
        self.freq_matrix_convert_cols()
        self.odds_matrix_populate()
        self.odds_matrix_to_log()

    def init_freq_odds_matrices(self):
        self.em_freq_matrix = ones((4,
                                    self.motif_width + 1))  # For pseudocounts
        self.em_log_odds_matrix = zeros((4, self.motif_width))

    def counts_matrix_populate(self):
        for i in range(len(self.fasta_file_seqs)):
            for j in range(self.motif_width):
                for k in range(4):
                    self.em_freq_matrix[k][0] = list(  # bgs in 1st column
                        self.count_bkgd_bases_dict.values())[k]
                    if list(self.motif_base_posit_freq_dict.values())[k][i][
                            j] == j + 1:  # This is the posit of the element after normalization
                        self.em_freq_matrix[k][
                            j + 1] += 1  # Counter, always offset to avoid bkgd

    def freq_matrix_convert_bgs(self):
        for i in range(4):
            self.em_freq_matrix[i][0] = round(
                list(self.count_bkgd_bases_dict.values())[i] / sum(
                    list(self.count_bkgd_bases_dict.values())), 3)

    def freq_matrix_convert_cols(self):
        col_totals = [sum(i) for i in zip(*self.em_freq_matrix)]
        for i in range(4):
            for j in range(1, self.motif_width + 1):  # Offset again
                self.em_freq_matrix[i][
                    j] = self.em_freq_matrix[i][j] / col_totals[j]
                self.em_freq_matrix[i][j] = round(self.em_freq_matrix[i][j], 3)

    def odds_matrix_populate(self):  # odds matrix is NOT offset!
        for i in range(4):
            for j in range(self.motif_width):
                self.em_log_odds_matrix[i][j] = round(
                    self.em_freq_matrix[i][j + 1] / self.em_freq_matrix[i][0],
                    3)

    def odds_matrix_to_log(self):
        for i in range(4):
            for j in range(self.motif_width):
                self.em_log_odds_matrix[i][j] = round(
                    log(self.em_log_odds_matrix[i][j], 2), 3)
