from class_EM_Matrix import EM_Matrix
from itertools import chain, tee
from numpy import array, cumsum, max, power, select, sum


class EM_Run(EM_Matrix):
    # The ML and scoring class
    def __init__(self, total_em_iters, fasta_file_seqs, motif_width,
                 em_log_odds_matrix):
        self.total_em_iters = total_em_iters
        self.fasta_file_seqs = fasta_file_seqs
        self.motif_width = motif_width
        self.em_log_odds_matrix = em_log_odds_matrix
        self.init_seq_lens()
        self.init_em_motifs()
        self.init_max_likely_dict()
        self.em_iter_over_total_rounds()
        self.em_max_sum()
        self.em_best_result()

    def prev_and_next(self, iterable):
        # Stack overflow iteration technique, Thx nosklo!
        prevs, items = tee(iterable, 2)
        prevs = chain([None], prevs)
        return zip(prevs, items)

    def init_seq_lens(self):
        # Every contiguous motif cumsum for merging, run only once
        self.seq_lens = cumsum(
            [len(i) - self.motif_width for i in self.fasta_file_seqs])

    def init_em_motifs(self):
        # Every contiguous motif, run only once
        self.em_motifs = array([
            list(i[j:j + self.motif_width]) for i in self.fasta_file_seqs
            for j in range(len(i) - self.motif_width)
        ])

    def init_max_likely_dict(self):
        # Preallocate or empty existing
        self.max_likely_dict = {
            "max_posits_matrix": [],
            "max_scores_matrix": [],
            "max_motifs_matrix": []
        }

    def em_iter_over_total_rounds(self):
        # Consolidate looping
        for i in range(self.total_em_iters):
            self.init_max_likely_dict()
            self.em_score()
            self.update_max_likely_dict()
            self.update_log_odds(self.max_likely_dict["max_posits_matrix"])

    def em_score(self):
        # Vectorized code is key, numpy select replaces for/if/else
        self.score = power(
            2,
            select([
                self.em_motifs == 'A', self.em_motifs == 'C',
                self.em_motifs == 'T', self.em_motifs == 'G'
            ], self.em_log_odds_matrix).sum(axis=1))

    def update_max_likely_dict(self):
        # Now simultaneously group scores and sequences, and append to dict
        for pre, cur in self.prev_and_next(self.seq_lens):
            max_pos = self.score[pre:cur].argmax()
            self.max_likely_dict["max_posits_matrix"].append(max_pos)
            self.max_likely_dict["max_scores_matrix"].append(
                self.score[pre:cur][max_pos])
            self.max_likely_dict["max_motifs_matrix"].append(''.join(
                self.em_motifs[pre:cur][max_pos]))

    def update_log_odds(self, updated_max_posits):
        # The ML magic
        self.motif_start_posits = updated_max_posits
        self.init_motifs_list()
        self.init_total_base_counts_dict()
        self.merge_all_bases()
        self.count_all_bases()
        self.init_bkgd_motif_counts_dict()
        self.count_all_bkgd_bases()
        self.init_motif_base_posit_freq_dict()
        self.motif_base_posit_freq_dict_populate()
        self.init_em_matrix()
        self.counts_matrix_convert_bgs()
        self.counts_matrix_populate()
        self.freq_matrix_convert_bgs()
        self.freq_col_totals()
        self.freq_matrix_convert_cols()
        self.odds_matrix_populate()
        self.odds_matrix_to_log()

    def em_max_sum(self):
        # These index the best results below
        self.max_likely_dict["final_score"] = max(
            self.max_likely_dict["max_scores_matrix"])
        self.max_likely_dict["final_sum_scores"] = sum(
            self.max_likely_dict["max_scores_matrix"])

    def em_best_result(self):
        # Last fxn within the em_core for loop
        self.max_likely_dict["final_seq"] = self.max_likely_dict[
            "max_scores_matrix"].index(self.max_likely_dict["final_score"])
        self.max_likely_dict["final_pos"] = self.max_likely_dict[
            "max_posits_matrix"][self.max_likely_dict["final_seq"]]
        self.max_likely_dict["final_motif"] = self.max_likely_dict[
            "max_motifs_matrix"][self.max_likely_dict["final_seq"]]
