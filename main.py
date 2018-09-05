from class_EM_Input import EM_Input
from class_EM_Core import EM_Core
from pprint import pprint
from termcolor import colored


# Map the colors
def color_the_sequences(em_core_obj, final_results):
    return map(lambda x, y: x.replace(x[y:y + em_core_obj.motif_width],
               colored(x[y:y + em_core_obj.motif_width], "red")),
               em_core_obj.fasta_file_seqs, final_results["max_posits_matrix"])


# Initialize input object
em_input_obj = EM_Input()  # Get Input
print("Commencing EM on", colored(em_input_obj.input_user_fasta_path,
                                  "magenta"))

# Initialize core object
em_core_obj = EM_Core(  # Start the nested iterations
    em_input_obj.motif_width, em_input_obj.total_rand_aligns,
    em_input_obj.total_em_iters, em_input_obj.fasta_file_seqs)

# Get the final results
final_results = max(
    em_core_obj.total_records, key=lambda x: x["final_sum_scores"])

# Print final results
print(colored("\n\nInput:", "green"))
print(*color_the_sequences(em_core_obj, final_results), sep="\n")
print(colored("\nResults:", "green"))
pprint(final_results)
