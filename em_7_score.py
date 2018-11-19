from em_4_count import em_count
from em_5_matrix import em_matrix
from numpy import power, select


def em_score(self):
    def em_score(self):
        # Vectorized numpy select replaces for/if/else
        self.score = power(
            2,
            select([
                self.contig_motifs == 'A', self.contig_motifs == 'C',
                self.contig_motifs == 'T', self.contig_motifs == 'G'
            ], self.em_log_odds).sum(axis=1))

    def update_max_likely_results(self):
        # Now simultaneously group scores and sequences, and append to dict
        self.max_likely_posits = [
            self.score[pre:cur].argmax()
            for pre, cur in zip(self.seq_cumsum[:-1], self.seq_cumsum[1:])
        ]

    def update_log_odds(self):
        # The ML magic
        self.motif_start_posits = self.max_likely_posits
        em_count(self)
        em_matrix(self)

    def em_iter_over_total_rounds(self):
        # Consolidate looping
        for i in range(self.total_em_iters):
            em_score(self)
            update_max_likely_results(self)
            update_log_odds(self)

    def init_max_likely_results(self):
        # Allocate for final results
        self.max_likely_results = {
            "Final Positions": [],
            "Final Scores": [],
            "Final Motifs": []
        }

    def record_iter_results(self):
        # Only record everything when finished with the ML iters
        for pre, cur, pos in zip(self.seq_cumsum[:-1], self.seq_cumsum[1:],
                                 self.max_likely_posits):
            self.max_likely_results["Final Positions"].append(pos)
            self.max_likely_results["Final Scores"].append(
                self.score[pre:cur][pos])
            self.max_likely_results["Final Motifs"].append(''.join(
                self.contig_motifs[pre:cur][pos]))

    em_iter_over_total_rounds(self)
    init_max_likely_results(self)
    record_iter_results(self)
