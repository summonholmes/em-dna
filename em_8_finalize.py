from pandas import DataFrame
from termcolor import colored


def em_finalize(self):
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
        print(*color_the_sequences(self), sep='\n')
        print(colored("\nResults:", "green"))
        print(self.final_dataframe)
        print(*self.final_results.items(), sep='\n')

    gen_dataframe(self)
    clean_final_dict(self)
    display_results(self)
