from numpy import empty, fromiter, log2, sum


def em_matrix(self):
    # Matrix operations

    def init_em_matrix(self):
        # Counts, freq, odds, then log odds
        self.em_log_odds = empty((4, self.motif_width + 1))
        # Will be overwritten

    def counts_matrix_convert_bkgds(self):
        # By Column
        self.em_log_odds[:, 0] = fromiter(
            self.total_bkgd_counts.values(), dtype=int)

    def counts_matrix_populate(self):
        # By Row
        for i, base in enumerate(('A', 'C', 'T', 'G')):
            self.em_log_odds[i, 1:] = self.motif_posits_freqs[base]

    def freq_matrix_convert_bkgds(self):
        # Turn background to freq # By Column
        self.em_log_odds[:, 0] /= sum(self.em_log_odds[:, 0])

    def freq_col_totals(self):
        # Sum columns
        self.col_totals = self.em_log_odds.sum(axis=0)

    def freq_matrix_convert_cols(self):
        # Turn columns to freq
        for i in range(1, self.motif_width + 1):  # By Column
            self.em_log_odds[:, i] /= self.col_totals[i]

    def odds_matrix_populate(self):
        # Odds of getting each base in each position
        for i in range(4):  # By Row
            self.em_log_odds[i, 1:] /= self.em_log_odds[i][0]

    def odds_matrix_to_log(self):
        # Log-odds matrix is NOT offset!
        self.em_log_odds = log2(self.em_log_odds[:, 1:])

    init_em_matrix(self)
    counts_matrix_convert_bkgds(self)
    counts_matrix_populate(self)
    freq_matrix_convert_bkgds(self)
    freq_col_totals(self)
    freq_matrix_convert_cols(self)
    odds_matrix_populate(self)
    odds_matrix_to_log(self)
