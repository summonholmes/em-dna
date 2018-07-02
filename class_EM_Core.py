from class_EM_Count import EM_Count
from class_EM_Matrix import EM_Matrix
from class_EM_Run import EM_Run
from math import log
from numpy import round, vectorize


class EM_Core:
    def __init__(self, motif_width, total_rand_aligns, total_em_iters,
                 fasta_file_seqs):
        self.motif_width = motif_width
        self.total_rand_aligns = total_rand_aligns
        self.total_em_iters = total_em_iters
        self.fasta_file_seqs = fasta_file_seqs
        self.total_records = []
        self.iterate_over_total_rand_aligns()

    def iterate_over_total_rand_aligns(self):
        for i in range(
                self.total_rand_aligns):  # Encompasses the rest of the program
            print(
                "Progress: {:2.1%}".format((i + 1) / self.total_rand_aligns),
                end="\r")
            self.em_counter_obj = EM_Count(self.motif_width,
                                           self.fasta_file_seqs)
            self.em_matrix_obj = EM_Matrix(
                self.motif_width, self.fasta_file_seqs,
                self.em_counter_obj.count_bkgd_bases_dict,
                self.em_counter_obj.motif_base_posit_freq_dict)
            self.em_run_obj = EM_Run(self.total_em_iters, self.fasta_file_seqs,
                                     self.motif_width,
                                     self.em_matrix_obj.em_log_odds_matrix)
            self.total_records.append(self.em_run_obj.best_results
                                      )  # Record best results of each round
