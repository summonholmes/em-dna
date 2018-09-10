from class_EM_Prep import EM_Prep
from class_EM_Count import EM_Count
from class_EM_Matrix import EM_Matrix
from class_EM_Run import EM_Run
from termcolor import colored


class EM_Core:
    # Runs and contains all EM rounds
    def __init__(self, motif_width, total_rand_aligns, total_em_iters,
                 fasta_file_seqs):
        self.motif_width = motif_width
        self.total_rand_aligns = total_rand_aligns
        self.total_em_iters = total_em_iters
        self.fasta_file_seqs = fasta_file_seqs
        self.em_prep = EM_Prep(fasta_file_seqs, motif_width)
        self.init_em()

    def init_em(self):
        # Get all results in a list comprehension
        self.total_records = [
            self.iter_total_rand_aligns(i)
            for i in range(self.total_rand_aligns)
        ]

    def iter_total_rand_aligns(self, i):
        # Encompasses the rest of the program
        self.print_colored(i)
        self.process_counter_obj()
        self.process_matrix_obj()
        self.process_run_obj()
        return self.em_run_obj.max_likely_results
        # Record best results of each round

    def print_colored(self, i):
        # Color progress meter
        print(
            colored(
                "Progress: {:2.1%}".format((i + 1) / self.total_rand_aligns),
                "blue"),
            end="\r")

    def process_counter_obj(self):
        # Step 1 - inherit Core
        self.em_counter_obj = EM_Count(self.motif_width, self.fasta_file_seqs,
                                       self.em_prep.total_bases_counts)

    def process_matrix_obj(self):
        # Step 2 - inherit Counter
        self.em_matrix_obj = EM_Matrix(self.motif_width, self.fasta_file_seqs,
                                       self.em_counter_obj.total_bkgd_counts,
                                       self.em_counter_obj.motif_posits_freqs)

    def process_run_obj(self):
        # Step 3 - inherit Matrix
        self.em_run_obj = EM_Run(
            self.total_em_iters, self.fasta_file_seqs, self.motif_width,
            self.em_matrix_obj.em_log_odds, self.em_prep.seq_cumsum,
            self.em_prep.contig_motifs, self.em_prep.total_bases_counts)
