from class_EM_Matrix import EM_Matrix
from itertools import tee, islice, chain
from numpy import array, select, power, cumsum


class EM_Run(EM_Matrix):  # The ML class
    def __init__(self, total_em_iters, fasta_file_seqs, motif_width,
                 em_log_odds_matrix):
        self.total_em_iters = total_em_iters
        self.fasta_file_seqs = fasta_file_seqs
        self.motif_width = motif_width
        self.em_log_odds_matrix = em_log_odds_matrix
        self.init_seq_lens()
        self.init_em_motifs()
        self.init_max_likely_dict()
        self.em_start_scoring()
        self.em_finalize()
        self.em_best_result()

    def prev_and_next(self, iterable):
        prevs, items = tee(iterable, 2)
        prevs = chain([None], prevs)
        return zip(prevs, items)  # Thx nosklo!

    def init_seq_lens(self):  # Every contiguous motif
        self.seq_lens = cumsum(  # Use for merging
            [len(i) - self.motif_width for i in self.fasta_file_seqs])

    def init_em_motifs(self):  # Every contiguous motif
        self.em_motifs = array([
            list(i[j:j + self.motif_width]) for i in self.fasta_file_seqs
            for j in range(len(i) - self.motif_width)
        ])

    def init_max_likely_dict(self):
        self.max_likely_dict = {  # Preallocate standard arrays
            "max_scores_matrix": [[] for i in range(self.total_em_iters)],
            "max_posits_matrix": [[] for i in range(self.total_em_iters)],
            "max_motifs_matrix": [[] for i in range(self.total_em_iters)]
        }  # Essentially deep-copy where references are not allowed

    def em_start_scoring(self):  # Vectorized code is key
        for i in range(self.total_em_iters):
            self.score = power(
                2,  # Numpy select replaces for/if/else
                select([
                    self.em_motifs == 'A', self.em_motifs == 'C',
                    self.em_motifs == 'T', self.em_motifs == 'G'
                ], self.em_log_odds_matrix).sum(axis=1))
            self.update_max_likely_dict(i)
            self.update_log_odds(self.max_likely_dict["max_posits_matrix"][i])

    def em_motifs_merge(self):  # Converting to by-sequence format
        self.em_seq_score = [(self.em_motifs[pre:cur], self.score[pre:cur])
                             for pre, cur in self.prev_and_next(self.seq_lens)]

    def update_max_likely_dict(self, i):
        self.em_motifs_merge()  # Tuple results by sequence
        for motif_score in self.em_seq_score:
            self.max_likely_dict["max_scores_matrix"][i].append(
                max(motif_score[1]))
            self.max_likely_dict["max_posits_matrix"][i].append(
                motif_score[1].argmax())
            # Convert motif back to string
            self.max_likely_dict["max_motifs_matrix"][i].append(''.join(
                motif_score[0][motif_score[1].argmax()]))

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
        }  # Update by indexing with last position/last iteration
        self.final_scores_dict["final_seq"] = self.max_likely_dict[
            "max_scores_matrix"][-1].index(
                self.final_scores_dict["final_score"])
        self.final_scores_dict["final_pos"] = self.max_likely_dict[
            "max_posits_matrix"][-1][self.final_scores_dict["final_seq"]]
        self.final_scores_dict["final_motif"] = self.max_likely_dict[
            "max_motifs_matrix"][-1][self.final_scores_dict["final_seq"]]
        self.final_scores_dict["final_sum_scores"] = sum(
            self.max_likely_dict["max_scores_matrix"][-1])

    def em_best_result(self):  # Last fxn within the em_core for loop
        for i in self.max_likely_dict.keys():
            self.max_likely_dict[i] = self.max_likely_dict[i][-1]
        self.best_results = {
            "final_scores_seqs_posits_motifs": self.final_scores_dict,
            "max_scores_posits_motifs": self.max_likely_dict
        }