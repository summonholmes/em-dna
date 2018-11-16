from numpy import array, cumsum


def em_prep(self):
    # Consolidate all consistent and reusable information

    def merge_all_bases(self):
        # Merge all bases to one string
        self.all_bases = ''.join(self.fasta_file_seqs)

    def init_total_bases_dict(self):
        # Initialize at 0
        self.total_bases_counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0}

    def count_all_bases(self):
        # Populate the dict above
        self.total_bases_counts['A'] += self.all_bases.count('A')
        self.total_bases_counts['C'] += self.all_bases.count('C')
        self.total_bases_counts['T'] += self.all_bases.count('T')
        self.total_bases_counts['G'] += self.all_bases.count('G')

    def init_seq_cumsum(self):
        # Every contiguous motif cumsum for merging, run only once
        self.seq_cumsum = [
            len(seq) - self.motif_width for seq in self.fasta_file_seqs
        ]
        self.seq_cumsum.insert(0, 0)
        self.seq_cumsum = cumsum(self.seq_cumsum)

    def init_contig_motifs(self):
        # Every contiguous motif, run only once
        self.contig_motifs = array([
            list(seq[start:start + self.motif_width])
            for seq in self.fasta_file_seqs
            for start in range(len(seq) - self.motif_width)
        ])

    merge_all_bases(self)
    init_total_bases_dict(self)
    count_all_bases(self)
    init_seq_cumsum(self)
    init_contig_motifs(self)
