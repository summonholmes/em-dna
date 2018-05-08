from Bio import SeqIO
from math import pow, log
from pprint import pprint
from random import randint, seed


class EM_Input:
    def __init__(self):
        self.fasta_file_seqs = []
        self.input_start_motif_width()
        self.input_start_align()
        self.input_start_iter()
        self.input_start_fasta()

    def input_start_motif_width(self):
        self.motif_width = self.check_if_int(
            input("Please specify the width of the motif: "))

    def input_start_align(self):
        input_start_align = self.check_if_int(
            input(
                "Use 50 initial random starting alignments? ('1' for yes or '0' for no): "
            ))
        if input_start_align == 1:
            self.total_random_alignments = 50
        elif input_start_align == 0:
            self.total_random_alignments = self.check_if_int(
                input(
                    "Please specify a value for the initial random starting alignments: "
                ))
        else:
            print("Unknown option.  Now exiting")
            exit(0)

    def input_start_iter(self):
        input_start_iter = self.check_if_int(
            input(
                "Use 500 iterations to perform the E-M steps? ('1' for yes or '0' for no): "
            ))
        if input_start_iter == 1:
            self.total_em_iterations = 500
        elif input_start_iter == 0:
            self.total_em_iterations = self.check_if_int(
                input("Please specify a value for the number of iterations: "))
        else:
            print("Unknown option.  Now exiting")
            exit(0)

    def input_start_fasta(self):
        user_fasta_path = input("Please specify the path of the fasta file: ")
        self.check_if_fasta(user_fasta_path)

    def check_if_int(self, int_return):
        try:
            int_return = int(int_return)
            return int_return
        except ValueError:
            print("You didn't provide an integer ya turkey!")
            exit(0)

    def check_if_fasta(self, user_fasta_path):
        try:
            for record in SeqIO.parse(user_fasta_path, "fasta"):
                self.fasta_file_seqs.append(str(record.seq))
        except IOError:
            print("Unable to open file.  Try again.")
            exit(0)


class EM_Core(EM_Input):
    def __init__(self, total_random_alignments, fasta_file_seqs, motif_width,
                 total_em_iterations):
        self.max_bit_score_arr = []
        self.final_record = []
        self.total_random_alignments = total_random_alignments
        self.fasta_file_seqs = fasta_file_seqs
        self.motif_width = motif_width
        self.total_em_iterations = total_em_iterations
        self.iterate_total_random_alignments()

    def iterate_total_random_alignments(self):
        for i in range(self.total_random_alignments):
            print(
                "Progress: {:2.1%}".format(
                    (i + 1) / self.total_random_alignments),
                end="\r")
            self.total_number_fasta_sequences = len(self.fasta_file_seqs)
            self.total_base_pairs_each_seq = [
                len(seq) for seq in self.fasta_file_seqs
            ]
            self.em_counter_obj = EM_Count(self.total_number_fasta_sequences,
                                           self.total_base_pairs_each_seq,
                                           self.motif_width,
                                           self.fasta_file_seqs)
            self.em_matrix_obj = EM_Matrix(
                self.em_counter_obj.motif_width,
                self.em_counter_obj.total_number_fasta_sequences,
                self.em_counter_obj.count_background_and_normalize_dict)
            self.max_bit_score_arr.append(self.em_matrix_obj.score_matrix_dict[
                "score_matrix_information"])
            self.max_bit_score_arr.sort()
            self.em_run_obj = EM_Run(
                self.total_em_iterations,
                self.em_matrix_obj.total_number_fasta_sequences,
                self.em_counter_obj.total_base_pairs_each_seq,
                self.em_matrix_obj.motif_width,
                self.em_matrix_obj.score_matrix_dict["score_matrix_log_odds"],
                self.em_counter_obj.fasta_file_seqs)
            self.final_record.append(self.em_run_obj.finish)


class EM_Count(EM_Core):
    def __init__(self, total_number_fasta_sequences, total_base_pairs_each_seq,
                 motif_width, fasta_file_seqs):
        self.total_number_fasta_sequences = total_number_fasta_sequences
        self.total_base_pairs_each_seq = total_base_pairs_each_seq
        self.motif_width = motif_width
        self.fasta_file_seqs = fasta_file_seqs
        self.init_rand_motif_pos()
        self.count_all_bases()
        self.init_motifs()
        self.init_background_motif_counts()
        self.count_motif_bases()
        self.normalize_counts()
        self.count_background_and_normalize_dict = {
            "count_background_bases": self.count_background_bases,
            "normalize_count_motif_bases": self.normalize_count_motif_bases
        }

    def init_rand_motif_pos(self):
        seed(a=None)
        self.motif_start_pos = [
            randint(0, (i - self.motif_width))
            for i in self.total_base_pairs_each_seq
        ]

    def count_all_bases(self):
        self.total_count_all_bases = {
            "count_a": 0,
            "count_c": 0,
            "count_t": 0,
            "count_g": 0
        }
        for i in range(self.total_number_fasta_sequences):
            self.total_count_all_bases["count_a"] += self.fasta_file_seqs[
                i].count('A')
            self.total_count_all_bases["count_c"] += self.fasta_file_seqs[
                i].count('C')
            self.total_count_all_bases["count_t"] += self.fasta_file_seqs[
                i].count('T')
            self.total_count_all_bases["count_g"] += self.fasta_file_seqs[
                i].count('G')
        return self.total_count_all_bases

    def init_motifs(self):
        self.the_initial_random_motif = []
        for i in range(self.total_number_fasta_sequences):
            x = self.fasta_file_seqs[i]
            y = self.motif_start_pos[i]
            z = self.motif_start_pos[i] + self.motif_width
            self.the_initial_random_motif.append(x[y:z])

    def init_background_motif_counts(self):
        self.count_background_bases = {
            "background_a": self.total_count_all_bases["count_a"],
            "background_c": self.total_count_all_bases["count_c"],
            "background_t": self.total_count_all_bases["count_t"],
            "background_g": self.total_count_all_bases["count_g"]
        }
        for i in range(self.total_number_fasta_sequences):
            self.count_background_bases[
                "background_a"] -= self.the_initial_random_motif[i].count('A')
            self.count_background_bases[
                "background_c"] -= self.the_initial_random_motif[i].count('C')
            self.count_background_bases[
                "background_t"] -= self.the_initial_random_motif[i].count('T')
            self.count_background_bases[
                "background_g"] -= self.the_initial_random_motif[i].count('G')

    def motif_actg_matrix(self):
        self.count_all_motif_bases = {
            "motif_pos_a": [],
            "motif_pos_c": [],
            "motif_pos_t": [],
            "motif_pos_g": []
        }
        for i in range(self.total_number_fasta_sequences):
            self.count_all_motif_bases["motif_pos_a"].append([])
            self.count_all_motif_bases["motif_pos_c"].append([])
            self.count_all_motif_bases["motif_pos_t"].append([])
            self.count_all_motif_bases["motif_pos_g"].append([])

    def count_motif_bases(self):
        self.motif_actg_matrix()
        for i in range(self.total_number_fasta_sequences):
            for j in range(self.motif_width):
                self.count_all_motif_bases["motif_pos_a"][i].append(
                    self.the_initial_random_motif[i].find('A', j, j + 1))
                self.count_all_motif_bases["motif_pos_c"][i].append(
                    self.the_initial_random_motif[i].find('C', j, j + 1))
                self.count_all_motif_bases["motif_pos_t"][i].append(
                    self.the_initial_random_motif[i].find('T', j, j + 1))
                self.count_all_motif_bases["motif_pos_g"][i].append(
                    self.the_initial_random_motif[i].find('G', j, j + 1))

    def normalize_counts(self):
        for i in self.count_all_motif_bases.keys():
            for j in range(self.total_number_fasta_sequences):
                for k in range(self.motif_width):
                    self.count_all_motif_bases[i][j][k] += 1
        self.normalize_count_motif_bases = self.count_all_motif_bases


class EM_Matrix(EM_Count):
    def __init__(self, motif_width, total_number_fasta_sequences,
                 count_background_and_normalize_dict):
        self.motif_width = motif_width
        self.total_number_fasta_sequences = total_number_fasta_sequences
        self.count_background_and_normalize_dict = count_background_and_normalize_dict
        self.counts_matrix()
        self.add_pseudocounts()
        self.freq_matrix()
        self.odds_matrix()
        self.log_odds_matrix()
        self.entropy_matrix()
        self.sum_entropy()
        self.information()
        self.score_matrix_dict = {
            "score_matrix_information": self.score_matrix_information,
            "score_matrix_log_odds": self.score_matrix_log_odds
        }

    def new_zeros_matrix(self, amount, motif_width):
        zeros_matrix = []
        for i in range(amount):
            zeros_matrix.append([])
            for j in range(motif_width):
                zeros_matrix[i].append(0)
        return zeros_matrix

    def counts_matrix(self):
        self.score_matrix = self.new_zeros_matrix(4, self.motif_width + 1)
        for i in range(self.total_number_fasta_sequences):
            for j in range(self.motif_width):
                for k in range(4):
                    self.score_matrix[k][0] = list(
                        self.count_background_and_normalize_dict[
                            "count_background_bases"].values())[k]
                    if list(self.count_background_and_normalize_dict[
                            "normalize_count_motif_bases"].values())[k][i][
                                j] == j + 1:
                        self.score_matrix[k][j + 1] += 1

    def add_pseudocounts(self):
        self.score_matrix_pseudo = self.new_zeros_matrix(
            4, self.motif_width + 1)
        for i in range(4):
            for j in range(self.motif_width + 1):
                self.score_matrix_pseudo[i][j] = self.score_matrix[i][j] + 1

    def freq_matrix(self):
        self.score_matrix_freq = self.new_zeros_matrix(4, self.motif_width + 1)
        for i in range(4):
            for j in range(self.motif_width + 1):
                self.score_matrix_freq[i][j] = self.score_matrix_pseudo[i][j]
            self.score_matrix_freq[i][0] = round(
                list(self.count_background_and_normalize_dict[
                    "count_background_bases"].values())[i] / sum(
                        list(self.count_background_and_normalize_dict[
                            "count_background_bases"].values())), 3)
        self.get_freq_matrix_score()

    def get_freq_matrix_score(self):
        col_totals = [sum(i) for i in zip(*self.score_matrix_freq)]
        for i in range(4):
            for j in range(1, self.motif_width + 1):
                self.score_matrix_freq[i][j] = self.score_matrix_freq[i][
                    j] / col_totals[j]
                self.score_matrix_freq[i][j] = round(
                    self.score_matrix_freq[i][j], 3)

    def odds_matrix(self):
        self.score_matrix_odds = self.new_zeros_matrix(4, self.motif_width)
        for i in range(4):
            for j in range(self.motif_width):
                self.score_matrix_odds[i][j] = round(
                    self.score_matrix_freq[i][j + 1] /
                    self.score_matrix_freq[i][0], 3)

    def log_odds_matrix(self):
        self.score_matrix_log_odds = self.new_zeros_matrix(4, self.motif_width)
        for i in range(4):
            for j in range(self.motif_width):
                self.score_matrix_log_odds[i][j] = round(
                    log(self.score_matrix_odds[i][j], 2), 3)

    def entropy_matrix(self):
        self.score_matrix_entropy = self.new_zeros_matrix(4, self.motif_width)
        for i in range(4):
            for j in range(self.motif_width):
                self.score_matrix_entropy[i][j] = round(
                    self.score_matrix_freq[i][j + 1] * log(
                        self.score_matrix_freq[i][j + 1], 2), 3)
        return self.score_matrix_entropy

    def sum_entropy(self):
        self.score_matrix_entropy_sum = [
            round(j, 3) * -1
            for j in [sum(i) for i in zip(*self.score_matrix_entropy)]
        ]

    def information(self):
        self.score_matrix_information = max(
            [round(2 - i, 3) for i in self.score_matrix_entropy_sum])


class EM_Run(EM_Matrix):
    def __init__(self, total_em_iterations, total_number_fasta_sequences,
                 total_base_pairs_each_seq, motif_width, score_matrix_log_odds,
                 fasta_file_seqs):
        self.total_em_iterations = total_em_iterations
        self.total_number_fasta_sequences = total_number_fasta_sequences
        self.total_base_pairs_each_seq = total_base_pairs_each_seq
        self.motif_width = motif_width
        self.score_matrix_log_odds = score_matrix_log_odds
        self.fasta_file_seqs = fasta_file_seqs
        self.init_em_motifs()
        self.init_max_scores_pos_motifs_matrices()
        self.exp_max_get_max_pos_score()
        self.init_finalize_scores_seqs_pos_motifs()
        self.finish_em()

    def new_em_zeros_matrix(self, total_number_fasta_sequences,
                            total_base_pairs_each_seq, motif_width):
        em_zeros_matrix = []
        for i in range(total_number_fasta_sequences):
            em_zeros_matrix.append([])
            for j in range(total_base_pairs_each_seq[i] - motif_width):
                em_zeros_matrix[i].append(0)
        return em_zeros_matrix

    def blank_matrix(self, total_em_iterations):
        blank_matrix = []
        for i in range(total_em_iterations):
            blank_matrix.append([])
        return blank_matrix

    def init_max_scores_pos_motifs_matrices(self):
        self.max_scores_pos_motifs_matrices = {
            "max_scores_matrix": self.blank_matrix(self.total_em_iterations),
            "max_pos_matrix": self.blank_matrix(self.total_em_iterations),
            "max_motifs_matrix": self.blank_matrix(self.total_em_iterations)
        }

    def init_em_motifs(self):
        self.em_motifs = self.new_em_zeros_matrix(
            self.total_number_fasta_sequences, self.total_base_pairs_each_seq,
            self.motif_width)
        for i in range(self.total_number_fasta_sequences):
            for j in range(
                    self.total_base_pairs_each_seq[i] - self.motif_width):
                x = self.fasta_file_seqs[i]
                z = j + self.motif_width
                self.em_motifs[i][j] = x[j:z]

    def exp_max_get_max_pos_score(self):
        for i in range(self.total_em_iterations):
            for j in range(self.total_number_fasta_sequences):
                self.last_motif_pos = self.total_base_pairs_each_seq[j] - self.motif_width
                self.exp_max_pos_score_iter(j, -1, -1)
                self.max_scores_pos_motifs_matrices["max_pos_matrix"][
                    i].append(self.max_pos_score_dict["max_pos"])
                self.max_scores_pos_motifs_matrices["max_scores_matrix"][
                    i].append(self.max_pos_score_dict["max_score"])
                self.max_scores_pos_motifs_matrices["max_motifs_matrix"][
                    i].append(self.fasta_file_seqs[j][self.max_pos_score_dict[
                        "max_pos"]:(self.max_pos_score_dict["max_pos"] +
                                    self.motif_width)])
            self.exp_max_update_log_odds(
                self.max_scores_pos_motifs_matrices["max_pos_matrix"][i])

    def exp_max_pos_score_iter(self, j, max_score, max_pos):
        for k in range(self.last_motif_pos):
            self.exp_max_score_log_odds(list(self.em_motifs[j][k]))
            self.sum_score_em_motif = sum(self.score_em_motif)
            self.power_score_em_motif = round(
                pow(2, self.sum_score_em_motif), 3)
            if self.power_score_em_motif > max_score:
                max_pos = k
                max_score = self.power_score_em_motif
        self.max_pos_score_dict = {"max_pos": max_pos, "max_score": max_score}

    def exp_max_score_log_odds(self, list_em_motif):
        self.score_em_motif = []
        for l in range(self.motif_width):
            if list_em_motif[l] == 'A':
                self.score_em_motif.append(self.score_matrix_log_odds[0][l])
            elif list_em_motif[l] == 'C':
                self.score_em_motif.append(self.score_matrix_log_odds[1][l])
            elif list_em_motif[l] == 'T':
                self.score_em_motif.append(self.score_matrix_log_odds[2][l])
            elif list_em_motif[l] == 'G':
                self.score_em_motif.append(self.score_matrix_log_odds[3][l])
            else:
                print("An error occurred in exp_max")
                exit(0)

    def exp_max_update_log_odds(self, final_max_pos):
        self.motif_start_pos = final_max_pos
        self.count_all_bases()
        self.init_motifs()
        self.init_background_motif_counts()
        self.count_motif_bases()
        self.normalize_counts()
        self.count_background_and_normalize_dict = {
            "count_background_bases": self.count_background_bases,
            "normalize_count_motif_bases": self.normalize_count_motif_bases
        }
        self.counts_matrix()
        self.add_pseudocounts()
        self.freq_matrix()
        self.odds_matrix()
        self.log_odds_matrix()

    def init_finalize_scores_seqs_pos_motifs(self):
        self.final_scores_seqs_pos_motifs = {
            "max_final_score":
            max(self.max_scores_pos_motifs_matrices["max_scores_matrix"][
                self.total_em_iterations - 1])
        }
        self.final_scores_seqs_pos_motifs[
            "max_final_sequence"] = self.max_scores_pos_motifs_matrices[
                "max_scores_matrix"][self.total_em_iterations - 1].index(
                    self.final_scores_seqs_pos_motifs["max_final_score"])
        self.final_scores_seqs_pos_motifs[
            "max_final_position"] = self.max_scores_pos_motifs_matrices[
                "max_pos_matrix"][self.total_em_iterations - 1][
                    self.final_scores_seqs_pos_motifs["max_final_sequence"]]
        self.final_scores_seqs_pos_motifs[
            "max_final_motif"] = self.max_scores_pos_motifs_matrices[
                "max_motifs_matrix"][self.total_em_iterations - 1][
                    self.final_scores_seqs_pos_motifs["max_final_sequence"]]
        self.final_scores_seqs_pos_motifs["sum_score_max_motif"] = sum(
            self.max_scores_pos_motifs_matrices["max_scores_matrix"][
                self.total_em_iterations - 1])

    def finish_em(self):
        for i in self.max_scores_pos_motifs_matrices.keys():
            self.max_scores_pos_motifs_matrices[
                i] = self.max_scores_pos_motifs_matrices[i][
                    self.total_em_iterations - 1]
        self.finish = {
            "final_scores_seqs_pos_motifs": self.final_scores_seqs_pos_motifs,
            "max_scores_pos_motifs_matrices":
            self.max_scores_pos_motifs_matrices
        }


em_input_obj = EM_Input()
em_core_obj = EM_Core(em_input_obj.total_random_alignments,
                      em_input_obj.fasta_file_seqs, em_input_obj.motif_width,
                      em_input_obj.total_em_iterations)
em_core_obj.final_results = max(
    em_core_obj.final_record,
    key=lambda x: x["final_scores_seqs_pos_motifs"]["sum_score_max_motif"])
pprint(em_core_obj.final_results)
