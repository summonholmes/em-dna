from class_EM_Matrix import EM_Matrix
from numpy import chararray


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
        self.em_motifs = [[
            i[j:j + self.motif_width]
            for j in range(len(i) - self.motif_width)
        ] for i in self.fasta_file_seqs]

    def init_max_likely_dict(self):
        self.max_likely_dict = {  # Preallocate, try and change
            "max_scores_matrix": [[] for i in range(self.total_em_iters)],
            "max_posits_matrix": [[] for i in range(self.total_em_iters)],
            "max_motifs_matrix": [[] for i in range(self.total_em_iters)]
        }  # Deep copies are actually slower

    def em_start_scoring(self):
        for i in range(self.total_em_iters):
            for j in self.fasta_file_seqs:
                last_motif_posit = len(j) - self.motif_width
                max_likely_iter_dict = self.em_score_posits(  # Temporarily store each iteration
                    self.fasta_file_seqs.index(j), -1, -1,
                    last_motif_posit)  # -1s so scoring works
                self.update_max_likely_dict(max_likely_iter_dict, i, j)
            self.update_log_odds(self.max_likely_dict["max_posits_matrix"][i])

    def em_score_posits(self, j, max_score, max_pos, last_motif_posit):
        for k in range(last_motif_posit):
            score_em_motif = self.em_log_odds(list(self.em_motifs[j][k]), [])
            sum_score_em_motif = sum(score_em_motif)
            power_score_em_motif = round(pow(2, sum_score_em_motif), 3)
            if power_score_em_motif > max_score:  # Local to the current motif
                max_posit = k
                max_score = power_score_em_motif
        return {"max_posit": max_posit, "max_score": max_score}

    def em_log_odds(self, list_em_motif, score_em_motif):
        for l in range(self.motif_width):
            if list_em_motif[l] == 'A':
                score_em_motif.append(self.em_log_odds_matrix[0][l])
            elif list_em_motif[l] == 'C':
                score_em_motif.append(self.em_log_odds_matrix[1][l])
            elif list_em_motif[l] == 'T':
                score_em_motif.append(self.em_log_odds_matrix[2][l])
            elif list_em_motif[l] == 'G':
                score_em_motif.append(self.em_log_odds_matrix[3][l])
            else:
                print("An error occurred in exp_max")
                exit(0)
        return score_em_motif

    def update_max_likely_dict(self, max_likely_iter_dict, i, j):
        self.max_likely_dict["max_scores_matrix"][i].append(
            max_likely_iter_dict["max_score"])
        self.max_likely_dict["max_posits_matrix"][i].append(
            max_likely_iter_dict["max_posit"])
        self.max_likely_dict["max_motifs_matrix"][i].append(
            j[max_likely_iter_dict["max_posit"]:(
                max_likely_iter_dict["max_posit"] + self.motif_width)])

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