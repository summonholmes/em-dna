class EM_Input:
    # Record user input.  Commented out interactiveness
    def __init__(self):
        # __init__ always runs everything in the class
        self.input_start_motif_width()
        self.input_start_align()
        self.input_start_iter()
        self.input_start_fasta()
        self.check_if_fasta()

    def input_start_motif_width(self):
        # # Motif width for rest of program
        # self.motif_width = input("Please specify the width of the motif: ")
        self.motif_width = 6

    # def check_if_int(self, int_return):
    #     # Successfully check if integer or close
    #     try:
    #         self.motif_width = int(self.motif_width)
    #         return self.motif_width
    #     except ValueError:
    #         print("You didn't provide an integer ya turkey!")
    #         exit(0)

    def input_start_align(self):
        # # Number of EM rounds from random positions
        # input_rand_starting_aligns = self.check_if_int(
        #     input("Use 50 initial rand starting aligns? "
        #           "('1' for yes or '0' for no): "))
        # if input_rand_starting_aligns == 1:
        #     self.total_rand_aligns = 50
        # elif input_rand_starting_aligns == 0:
        #     self.total_rand_aligns = self.check_if_int(
        #         input("Please specify a value for the initial rand "
        #               "starting aligns: "))
        # else:
        #     print("Unknown option.  Now exiting")
        #     exit(0)
        self.total_rand_aligns = 50

    def input_start_iter(self):
        # # Number of ML iterations per EM round
        # input_total_em_iters = self.check_if_int(
        #     input("Use 50 iters to perform the E-M steps? "
        #           "('1' for yes or '0' for no): "))
        # if input_total_em_iters == 1:
        #     self.total_em_iters = 50
        # elif input_total_em_iters == 0:
        #     self.total_em_iters = self.check_if_int(
        #         input("Please specify a value for the number of iters: "))
        # else:
        #     print("Unknown option.  Now exiting")
        #     exit(0)
        self.total_em_iters = 50

    def input_start_fasta(self):
        # # Fasta file to read
        # self.input_user_fasta_path = input(
        #     "Please specify the path of the fasta file: ")
        # self.check_if_fasta(self.input_user_fasta_path)
        self.input_user_fasta_path = "example.fasta"

    def check_if_fasta(self):
        # Successfully open and process fasta or close
        try:
            with open(self.input_user_fasta_path) as fasta:
                self.process_fasta(fasta)
            fasta.close()
        except IOError:
            print("Unable to open file.  Try again.")
            exit(0)

    def process_fasta(self, fasta):
        # Read, split, and skip every other line
        self.fasta_file_seqs = fasta.read()
        self.fasta_file_seqs = self.fasta_file_seqs.split("\n")
        self.fasta_file_seqs = [seq for seq in self.fasta_file_seqs[1::2]]
