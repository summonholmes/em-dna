from Bio import SeqIO
from math import log
from pprint import pprint
from random import randint, seed


class EM_Core:
    def __init__(self):
        self.fasta_file_seq = []
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
            self.user_align = 50
        elif input_start_align == 0:
            self.user_align = self.check_if_int(
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
            self.user_iter = 500
        elif input_start_iter == 0:
            self.user_iter = self.check_if_int(
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
                self.fasta_file_seq.append(str(record.seq))
        except IOError:
            print("Unable to open file.  Try again.")
            exit(0)


class EM_Align(EM_Core):
    def __init__(self, user_align, fasta_file_seq, motif_width, user_iter):
        self.max_bit_score_arr = []
        self.final_record = []
        self.user_align = user_align
        self.fasta_file_seq = fasta_file_seq
        self.motif_width = motif_width
        self.user_iter = user_iter
        for i in range(self.user_align):
            print(
                "Progress: {:2.1%}".format((i + 1) / self.user_align),
                end="\r")
            self.len_list = len(self.fasta_file_seq)
            self.len_seq = [len(seq) for seq in self.fasta_file_seq]
            self.em_counter_obj = EM_Count(self.len_seq, self.motif_width,
                                           self.len_list, self.fasta_file_seq)
            self.em_matrix_obj = EM_Matrix(
                self.em_counter_obj.motif_width, self.em_counter_obj.len_list,
                self.em_counter_obj.count_bases_dict)
            self.max_bit_score_arr.append(self.em_matrix_obj.score_matrix_dict[
                "score_matrix_information"])
            self.max_bit_score_arr.sort()
            self.em_run_obj = EM_Run(
                self.em_matrix_obj.len_list, self.em_counter_obj.len_seq,
                self.em_matrix_obj.motif_width,
                self.em_counter_obj.fasta_file_seq, self.user_iter,
                self.em_matrix_obj.score_matrix_log_odds,
                self.em_matrix_obj.count_bases_dict)
            self.final_record.append(self.em_run_obj.finish)
        self.final_results = max(
            self.final_record,
            key=lambda x: x["max_final_sco_seq_pos_mot"]["sum_score_max_motif"]
        )
        pprint(self.final_results)


class EM_Count(EM_Align):
    def __init__(self, len_seq, motif_width, len_list, fasta_file_seq):
        self.len_seq = len_seq
        self.motif_width = motif_width
        self.len_list = len_list
        self.fasta_file_seq = fasta_file_seq
        self.start_rand()
        self.count_all_bases()
        self.init_motifs()
        self.init_background_motif_counts()
        self.count_motif_bases()
        self.normalize_counts()
        self.count_bases_dict = {
            "count_background_bases": self.count_background_bases,
            "normalize_count_motif_bases": self.normalize_count_motif_bases
        }

    def start_rand(self):
        seed(a=None)
        self.motif_start_pos = [
            randint(0, (i - self.motif_width)) for i in self.len_seq
        ]

    def count_all_bases(self):
        self.count_bases = {
            "count_a": 0,
            "count_c": 0,
            "count_t": 0,
            "count_g": 0
        }
        for i in range(self.len_list):
            self.count_bases["count_a"] += self.fasta_file_seq[i].count('A')
            self.count_bases["count_c"] += self.fasta_file_seq[i].count('C')
            self.count_bases["count_t"] += self.fasta_file_seq[i].count('T')
            self.count_bases["count_g"] += self.fasta_file_seq[i].count('G')

    def init_motifs(self):
        self.motif = []
        for i in range(self.len_list):
            x = self.fasta_file_seq[i]
            y = self.motif_start_pos[i]
            z = self.motif_start_pos[i] + self.motif_width
            self.motif.append(x[y:z])

    def init_background_motif_counts(self):
        self.count_background_bases = {
            "background_a": self.count_bases["count_a"],
            "background_c": self.count_bases["count_c"],
            "background_t": self.count_bases["count_t"],
            "background_g": self.count_bases["count_g"]
        }
        for i in range(self.len_list):
            self.count_background_bases["background_a"] -= self.motif[i].count(
                'A')
            self.count_background_bases["background_c"] -= self.motif[i].count(
                'C')
            self.count_background_bases["background_t"] -= self.motif[i].count(
                'T')
            self.count_background_bases["background_g"] -= self.motif[i].count(
                'G')

    def count_motif_bases(self):
        self.motif_actg_matrix()
        for i in range(self.len_list):
            for j in range(self.motif_width):
                self.count_all_motif_bases["motif_pos_a"][i].append(
                    self.motif[i].find('A', j, j + 1))
                self.count_all_motif_bases["motif_pos_c"][i].append(
                    self.motif[i].find('C', j, j + 1))
                self.count_all_motif_bases["motif_pos_t"][i].append(
                    self.motif[i].find('T', j, j + 1))
                self.count_all_motif_bases["motif_pos_g"][i].append(
                    self.motif[i].find('G', j, j + 1))

    def motif_actg_matrix(self):
        self.count_all_motif_bases = {
            "motif_pos_a": [],
            "motif_pos_c": [],
            "motif_pos_t": [],
            "motif_pos_g": []
        }
        for i in range(self.len_list):
            self.count_all_motif_bases["motif_pos_a"].append([])
            self.count_all_motif_bases["motif_pos_c"].append([])
            self.count_all_motif_bases["motif_pos_t"].append([])
            self.count_all_motif_bases["motif_pos_g"].append([])

    def normalize_counts(self):
        for i in self.count_all_motif_bases.keys():
            for j in range(self.len_list):
                for k in range(self.motif_width):
                    self.count_all_motif_bases[i][j][k] += 1
        self.normalize_count_motif_bases = self.count_all_motif_bases


class EM_Matrix(EM_Count):
    def __init__(self, motif_width, len_list, count_bases_dict):
        self.motif_width = motif_width
        self.len_list = len_list
        self.count_bases_dict = count_bases_dict
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

    def counts_matrix(self):
        self.score_matrix = self.new_zeros_matrix(4, self.motif_width + 1)
        for i in range(self.len_list):
            for j in range(self.motif_width):
                for k in range(4):
                    self.score_matrix[k][0] = list(self.count_bases_dict[
                        "count_background_bases"].values())[k]
                    if list(self.count_bases_dict[
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
                list(self.count_bases_dict["count_background_bases"].values())[
                    i] / sum(
                        list(self.count_bases_dict["count_background_bases"]
                             .values())), 3)
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
        return self.score_matrix_log_odds

    def entropy_matrix(self):
        self.score_matrix_entropy = self.new_zeros_matrix(4, self.motif_width)
        for i in range(4):
            for j in range(self.motif_width):
                self.score_matrix_entropy[i][j] = round(
                    self.score_matrix_freq[i][j + 1] * log(
                        self.score_matrix_freq[i][j + 1], 2), 3)

    def sum_entropy(self):
        self.score_matrix_entropy_sum = [
            round(j, 3) * -1
            for j in [sum(i) for i in zip(*self.score_matrix_entropy)]
        ]

    def information(self):
        self.score_matrix_information = max(
            [round(2 - i, 3) for i in self.score_matrix_entropy_sum])

    def new_zeros_matrix(self, amount, motif_width):
        zeros_matrix = []
        for i in range(amount):
            zeros_matrix.append([])
            for j in range(motif_width):
                zeros_matrix[i].append(0)
        return zeros_matrix

    def new_motif_zeros_matrix(self, len_list, len_seq, motif_width):
        motif_zeros_matrix = []
        for i in range(len_list):
            motif_zeros_matrix.append([])
            for j in range(len_seq[i] - motif_width):
                motif_zeros_matrix[i].append(0)
        return motif_zeros_matrix


class EM_Run(EM_Matrix):
    def __init__(self, len_list, len_seq, motif_width, fasta_file_seq,
                 user_iter, score_matrix_log_odds, count_bases_dict):
        self.len_list = len_list
        self.len_seq = len_seq
        self.motif_width = motif_width
        self.fasta_file_seq = fasta_file_seq
        self.user_iter = user_iter
        self.score_matrix_log_odds = score_matrix_log_odds
        self.count_bases_dict = count_bases_dict
        self.init_em_motifs()
        self.init_scores_pos_motifs()
        self.exp_max_get_max_pos_score()
        self.init_max_final_sco_seq_pos_mot()
        for i in self.scores_pos_motifs.keys():
            self.scores_pos_motifs[i] = self.scores_pos_motifs[i][
                self.user_iter - 1]
        self.finish = {
            "max_final_sco_seq_pos_mot": self.max_final_sco_seq_pos_mot,
            "scores_pos_motifs": self.scores_pos_motifs
        }

    def init_em_motifs(self):
        self.em_motifs = self.new_motif_zeros_matrix(
            self.len_list, self.len_seq, self.motif_width)
        for i in range(self.len_list):
            for j in range(self.len_seq[i] - self.motif_width):
                x = self.fasta_file_seq[i]
                z = j + self.motif_width
                self.em_motifs[i][j] = x[j:z]

    def init_scores_pos_motifs(self):
        self.scores_pos_motifs = {
            "max_scores": self.blank_matrix(),
            "max_pos": self.blank_matrix(),
            "max_motifs": self.blank_matrix()
        }

    def exp_max_get_max_pos_score(self):
        self.max_pos_score = {"max_pos": -1, "max_score": -1}
        for i in range(self.user_iter):
            for j in range(self.len_list):
                self.last_motif_pos = self.len_seq[j] - self.motif_width
                self.exp_max_pos_score_iter(j)
                self.scores_pos_motifs["max_pos"][i].append(
                    self.max_pos_score["max_pos"])
                self.scores_pos_motifs["max_scores"][i].append(
                    self.max_pos_score["max_score"])
                self.scores_pos_motifs["max_motifs"][i].append(
                    self.fasta_file_seq[j][self.max_pos_score["max_pos"]:(
                        self.max_pos_score["max_pos"] + self.motif_width)])
            self.exp_max_matrices(i)

    def exp_max_pos_score_iter(self, j):
        for k in range(self.last_motif_pos):
            self.exp_max_score_log_odds(j, k)
            sum_score_em_motif = sum(self.score_em_motif)
            power_score_em_motif = round(pow(2, sum_score_em_motif), 3)
            if power_score_em_motif > self.max_pos_score["max_score"]:
                self.max_pos_score["max_pos"] = k
                self.max_pos_score["max_score"] = power_score_em_motif

    def exp_max_score_log_odds(self, j, k):
        self.score_em_motif = []
        for l in range(self.motif_width):
            if list(self.em_motifs[j][k])[l] == 'A':
                self.score_em_motif.append(self.score_matrix_log_odds[0][l])
            elif list(self.em_motifs[j][k])[l] == 'C':
                self.score_em_motif.append(self.score_matrix_log_odds[1][l])
            elif list(self.em_motifs[j][k])[l] == 'T':
                self.score_em_motif.append(self.score_matrix_log_odds[2][l])
            elif list(self.em_motifs[j][k])[l] == 'G':
                self.score_em_motif.append(self.score_matrix_log_odds[3][l])
            else:
                print("An error occurred in exp_max")
                exit(0)

    def exp_max_matrices(self, i):
        self.count_all_bases()
        self.motif_start_pos = self.scores_pos_motifs["max_pos"][i]
        self.init_motifs()
        self.init_background_motif_counts()
        self.count_motif_bases()
        self.normalize_counts()
        self.counts_matrix()
        self.add_pseudocounts()
        self.freq_matrix()
        self.odds_matrix()
        self.log_odds_matrix()

    def init_max_final_sco_seq_pos_mot(self):
        self.max_final_sco_seq_pos_mot = {
            "max_final_score":
            max(self.scores_pos_motifs["max_scores"][self.user_iter - 1])
        }
        self.max_final_sco_seq_pos_mot[
            "max_final_sequence"] = self.scores_pos_motifs[
                "max_scores"][self.user_iter - 1].index(
                    self.max_final_sco_seq_pos_mot["max_final_score"])
        self.max_final_sco_seq_pos_mot[
            "max_final_position"] = self.scores_pos_motifs["max_pos"][
                self.user_iter
                - 1][self.max_final_sco_seq_pos_mot["max_final_sequence"]]
        self.max_final_sco_seq_pos_mot[
            "max_final_motif"] = self.scores_pos_motifs["max_motifs"][
                self.user_iter
                - 1][self.max_final_sco_seq_pos_mot["max_final_sequence"]]
        self.max_final_sco_seq_pos_mot["sum_score_max_motif"] = sum(
            self.scores_pos_motifs["max_scores"][self.user_iter - 1])

    def blank_matrix(self):
        blank_matrix = []
        for i in range(self.user_iter):
            blank_matrix.append([])
        return blank_matrix


em_core_obj = EM_Core()
em_align_obj = EM_Align(em_core_obj.user_align, em_core_obj.fasta_file_seq,
                        em_core_obj.motif_width, em_core_obj.user_iter)
