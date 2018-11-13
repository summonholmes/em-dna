from class_EM_Input import EM_Input
from class_EM_Prep import EM_Prep
from class_EM_Count import EM_Count
from class_EM_Matrix import EM_Matrix
from class_EM_Run import EM_Run
from pandas import DataFrame
from termcolor import colored


class EM_Core:
    # Runs and contains all EM rounds
    def __init__(self):
        # Initialize input object
        em_input = EM_Input()
        print("Commencing EM on", colored(em_input.input_fasta_path,
                                          "magenta"))

        # Get params
        self.motif_width = em_input.motif_width
        self.total_rand_aligns = em_input.total_rand_aligns
        self.total_em_iters = em_input.total_em_iters
        self.fasta_file_seqs = em_input.fasta_file_seqs
        self.em_prep = EM_Prep(self.fasta_file_seqs, self.motif_width)

        # Init process
        self.init_em()
        self.finalize_max_results()
        self.gen_dataframe()
        self.clean_final_dict()

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
            end='\r')

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

    def finalize_max_results(self):
        self.final_results = max(
            self.total_records, key=lambda x: x["Max Sum Scores"])

    def color_the_sequences(self):
        # Map the colors
        return map(
            lambda seq, pos, align:
                ' ' * align + seq[:pos] +
                colored(seq[pos:pos + self.motif_width], "red") +
                seq[pos + 1:],
                self.fasta_file_seqs,
                self.final_dataframe["Final Positions"],
                self.final_dataframe["Final Positions"].max(
                ) - self.final_dataframe["Final Positions"])

    def gen_dataframe(self):
        filter_keys = ("Final Scores", "Final Positions", "Final Motifs")
        self.final_dataframe = DataFrame(
            {key: self.final_results[key]
             for key in filter_keys})

    def clean_final_dict(self):
        del self.final_results["Final Positions"]
        del self.final_results["Final Scores"]
        del self.final_results["Final Motifs"]

    def display_results(self):
        # Print final results
        print(colored("\n\nAlignment:", "green"))
        print(*self.color_the_sequences(), sep='\n')
        print(colored("\nResults:", "green"))
        print(self.final_dataframe)
        print(*self.final_results.items(), sep='\n')
