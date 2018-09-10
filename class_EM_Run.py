from class_EM_Matrix import EM_Matrix
from itertools import chain, tee
from numpy import max, power, select, sum


class EM_Run(EM_Matrix):
    # The ML and scoring class
    def __init__(self, total_em_iters, fasta_file_seqs, motif_width,
                 em_log_odds, seq_cumsum, contig_motifs, total_bases_counts):
        self.total_em_iters = total_em_iters
        self.fasta_file_seqs = fasta_file_seqs
        self.motif_width = motif_width
        self.em_log_odds = em_log_odds
        self.seq_cumsum = seq_cumsum
        self.contig_motifs = contig_motifs
        self.total_bases_counts = total_bases_counts
        self.em_iter_over_total_rounds()
        self.em_max_sum()
        self.em_best_result()

    def prev_and_next(self, iterable):
        # Stack overflow iteration technique, Thx nosklo!
        prevs, items = tee(iterable, 2)
        prevs = chain([None], prevs)
        return zip(prevs, items)

    def em_iter_over_total_rounds(self):
        # Consolidate looping
        for i in range(self.total_em_iters):
            self.init_max_likely_results()
            self.em_score()
            self.update_max_likely_results()
            self.update_log_odds()

    def init_max_likely_results(self):
        # Preallocate or empty existing
        self.max_likely_results = {
            "max_posits": [],
            "max_scores": [],
            "max_motifs": []
        }

    def em_score(self):
        # Vectorized code is key, numpy select replaces for/if/else
        self.score = power(
            2,
            select([
                self.contig_motifs == 'A', self.contig_motifs == 'C',
                self.contig_motifs == 'T', self.contig_motifs == 'G'
            ], self.em_log_odds).sum(axis=1))

    def update_max_likely_results(self):
        # Now simultaneously group scores and sequences, and append to dict
        for pre, cur in self.prev_and_next(self.seq_cumsum):
            max_pos = self.score[pre:cur].argmax()
            self.max_likely_results["max_posits"].append(max_pos)
            self.max_likely_results["max_scores"].append(
                self.score[pre:cur][max_pos])
            self.max_likely_results["max_motifs"].append(''.join(
                self.contig_motifs[pre:cur][max_pos]))

    def update_log_odds(self):
        # The ML magic
        self.motif_start_posits = self.max_likely_results["max_posits"]
        self.init_em_motifs()
        self.merge_all_motif_bases()
        self.init_bkgd_counts()
        self.count_all_bkgd_bases()
        self.init_motif_posits_freqs()
        self.count_motif_posits_freqs()
        self.init_em_matrix()
        self.counts_matrix_convert_bkgds()
        self.counts_matrix_populate()
        self.freq_matrix_convert_bkgds()
        self.freq_col_totals()
        self.freq_matrix_convert_cols()
        self.odds_matrix_populate()
        self.odds_matrix_to_log()

    def em_max_sum(self):
        # These index the best results below
        self.max_likely_results["final_score"] = max(
            self.max_likely_results["max_scores"])
        self.max_likely_results["final_sum_scores"] = sum(
            self.max_likely_results["max_scores"])

    def em_best_result(self):
        # Last fxn within the em_core for loop
        self.max_likely_results["final_seq"] = self.max_likely_results[
            "max_scores"].index(self.max_likely_results["final_score"])
        self.max_likely_results["final_pos"] = self.max_likely_results[
            "max_posits"][self.max_likely_results["final_seq"]]
        self.max_likely_results["final_motif"] = self.max_likely_results[
            "max_motifs"][self.max_likely_results["final_seq"]]
