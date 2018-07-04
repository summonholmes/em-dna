from class_EM_Input import EM_Input
from class_EM_Core import EM_Core
from pprint import pprint

em_input_obj = EM_Input()  # Get Input
em_core_obj = EM_Core(  # Start the nested iterations
    em_input_obj.motif_width, em_input_obj.total_rand_aligns,
    em_input_obj.total_em_iters, em_input_obj.fasta_file_seqs)
final_results = max(  # Get the final results
    em_core_obj.total_records,
    key=lambda x: x["final_scores_seqs_posits_motifs"]["final_sum_scores"])
pprint(final_results)  # Make output pretty
