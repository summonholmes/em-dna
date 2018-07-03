from class_EM_Matrix import EM_Matrix
from numpy import array, empty, select, power


class EM_Run(EM_Matrix):  # Still working on this class
    def __init__(self, total_em_iters, fasta_file_seqs, motif_width,
                 em_log_odds_matrix):
        self.total_em_iters = total_em_iters
        self.fasta_file_seqs = fasta_file_seqs
        self.motif_width = motif_width
        self.em_log_odds_matrix = em_log_odds_matrix
        self.init_em_motifs()
        self.init_max_likely_dict()
        self.em_start_scoring()
        self.em_finalize()
        self.em_best_result()

    def init_em_motifs(self):  # Every possibility generated
        self.em_motifs = [
            array([  # Numpy character array required
                list(i[j:j + self.motif_width])
                for j in range(len(i) - self.motif_width)
            ]) for i in self.fasta_file_seqs
        ]

    def init_max_likely_dict(self):
        self.max_likely_dict = {  # Preallocate, try and change
            "max_scores_matrix": [[] for i in range(self.total_em_iters)],
            "max_posits_matrix": [[] for i in range(self.total_em_iters)],
            "max_motifs_matrix": [[] for i in range(self.total_em_iters)]
        }  # Deep copies are actually slower

    def em_start_scoring(self):  # Vectorized code is key
        for i in range(self.total_em_iters):
            choice_list = [i for i in self.em_log_odds_matrix]
            for motif_set in self.em_motifs:
                cond_list = [
                    motif_set == 'A', motif_set == 'C', motif_set == 'T',
                    motif_set == 'G'
                ]  # This replaces if/else only iterates per sequence
                score = power(2, select(cond_list, choice_list).sum(axis=1))
                self.update_max_likely_dict(i, motif_set, score)
            self.update_log_odds(self.max_likely_dict["max_posits_matrix"][i])

    def update_max_likely_dict(self, i, motif_set, score):
        self.max_likely_dict["max_scores_matrix"][i].append(max(score))
        self.max_likely_dict["max_posits_matrix"][i].append(score.argmax())
        self.max_likely_dict["max_motifs_matrix"][i].append(''.join(
            motif_set[score.argmax()]))

    def update_log_odds(self, updated_max_posits):  # The ML magic
        self.motif_start_posits = updated_max_posits
        self.init_total_base_counts_dict()
        self.count_all_bases()
        self.init_motifs_list()
        self.init_bkgd_motif_counts_dict()
        self.count_all_bkgd_bases()
        self.init_motif_base_posit_freq_dict()
        self.motif_base_posit_freq_dict_populate()
        self.init_em_matrix()
        self.counts_matrix_populate()
        self.freq_matrix_convert_bgs()
        self.freq_matrix_convert_cols()
        self.odds_matrix_populate()
        self.odds_matrix_to_log()

    def em_finalize(self):
        self.final_scores_dict = {
            "final_score": max(self.max_likely_dict["max_scores_matrix"][-1])
        }  # Update by indexing
        self.final_scores_dict["final_seq"] = self.max_likely_dict[
            "max_scores_matrix"][-1].index(
                self.final_scores_dict["final_score"])
        self.final_scores_dict["final_pos"] = self.max_likely_dict[
            "max_posits_matrix"][-1][self.final_scores_dict["final_seq"]]
        self.final_scores_dict["final_motif"] = self.max_likely_dict[
            "max_motifs_matrix"][-1][self.final_scores_dict["final_seq"]]
        self.final_scores_dict["final_sum_scores"] = sum(
            self.max_likely_dict["max_scores_matrix"][-1])

    def em_best_result(self):  # Still local to em_core for loop
        for i in self.max_likely_dict.keys():
            self.max_likely_dict[i] = self.max_likely_dict[i][-1]
        self.best_results = {
            "final_scores_seqs_posits_motifs": self.final_scores_dict,
            "max_scores_posits_motifs": self.max_likely_dict
        }