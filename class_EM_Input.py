from Bio import SeqIO


class EM_Input:
    def __init__(self):
        self.fasta_file_seqs = []
        self.input_start_motif_width()
        self.input_start_align()
        self.input_start_iter()
        self.input_start_fasta()

    def input_start_motif_width(self):
        # self.motif_width = self.check_if_int(
        #     input("Please specify the width of the motif: "))
        self.motif_width = 6

    def input_start_align(self):
        # input_rand_starting_aligns = self.check_if_int(
        #     input(
        #         "Use 50 initial rand starting aligns? ('1' for yes or '0' for no): "
        #     ))
        # if input_rand_starting_aligns == 1:
        self.total_rand_aligns = 50
        # elif input_rand_starting_aligns == 0:
        #     self.total_rand_aligns = self.check_if_int(
        #         input(
        #             "Please specify a value for the initial rand starting aligns: "
        #         ))
        # else:
        #     print("Unknown option.  Now exiting")
        #     exit(0)

    def input_start_iter(self):
        # input_total_em_iters = self.check_if_int(
        #     input(
        #         "Use 500 iters to perform the E-M steps? ('1' for yes or '0' for no): "
        #     ))
        # if input_total_em_iters == 1:
        self.total_em_iters = 50
        # elif input_total_em_iters == 0:
        #     self.total_em_iters = self.check_if_int(
        #         input("Please specify a value for the number of iters: "))
        # else:
        #     print("Unknown option.  Now exiting")
        #     exit(0)

    def input_start_fasta(self):
        # input_user_fasta_path = input(
        #     "Please specify the path of the fasta file: ")
        # input_user_fasta_path = "example.fasta"
        self.check_if_fasta("example.fasta")

    def check_if_int(self, int_return):
        try:
            int_return = int(int_return)
            return int_return
        except ValueError:
            print("You didn't provide an integer ya turkey!")
            exit(0)

    def check_if_fasta(self, input_user_fasta_path):
        try:
            for record in SeqIO.parse(input_user_fasta_path, "fasta"):
                self.fasta_file_seqs.append(str(record.seq))
        except IOError:
            print("Unable to open file.  Try again.")
            exit(0)