from em_1_input import em_input
from em_2_prep import em_prep
from em_3_process import em_process
from termcolor import colored


class EM_DNA:
    # Runs and contains all EM rounds
    def __init__(self):
        # Initialize input object
        em_input(self)
        em_prep(self)
        print("Commencing EM on", colored(self.input_fasta_path, "magenta"))

        # Init EM process
        em_process(self)
