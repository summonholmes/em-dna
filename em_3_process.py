from em_4_count import em_count
from em_5_matrix import em_matrix
from em_6_run import em_run
from random import randint
from termcolor import colored


def em_process(self):
    def print_colored(self, i):
        # Color progress meter
        print(
            colored(
                "Progress: {:2.1%}".format((i + 1) / self.total_rand_aligns),
                "blue"),
            end='\r')

    def init_rand_posits(self):
        # Starting positions are always random
        self.motif_start_posits = [
            randint(0, (len_seq - self.motif_width))
            for len_seq in [len(seq) for seq in self.fasta_file_seqs]
        ]

    def iter_total_rand_aligns(self, i):
        # Encompasses the rest of the program
        print_colored(self, i)
        init_rand_posits(self)
        em_count(self)
        em_matrix(self)
        em_run(self)
        return self.max_likely_results
        # Record best results of each round

    def init_em(self):
        # Get all results in a list comprehension
        self.total_records = [
            iter_total_rand_aligns(self, i)
            for i in range(self.total_rand_aligns)
        ]

    def finalize_max_results(self):
        self.final_results = max(
            self.total_records, key=lambda x: x["Max Sum Scores"])

    init_em(self)
    finalize_max_results(self)
