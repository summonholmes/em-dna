from em_4_count import em_count
from em_5_matrix import em_matrix
from numpy import power, select


def em_score(self):
    def init_max_likely_results(self):
        # Preallocate or empty existing
        self.max_likely_results = {
            "Final Positions": [],
            "Final Scores": [],
            "Final Motifs": []
        }

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
        for pre, cur in zip(self.seq_cumsum[:-1], self.seq_cumsum[1:]):
            max_pos = self.score[pre:cur].argmax()
            self.max_likely_results["Final Positions"].append(max_pos)
            self.max_likely_results["Final Scores"].append(
                self.score[pre:cur][max_pos])
            self.max_likely_results["Final Motifs"].append(''.join(
                self.contig_motifs[pre:cur][max_pos]))

    def update_log_odds(self):
        # The ML magic
        self.motif_start_posits = self.max_likely_results["Final Positions"]
        em_count(self)
        em_matrix(self)

    def em_iter_over_total_rounds(self):
        # Consolidate looping
        for i in range(self.total_em_iters):
            init_max_likely_results(self)
            em_score(self)
            update_max_likely_results(self)
            update_log_odds(self)

    em_iter_over_total_rounds(self)
