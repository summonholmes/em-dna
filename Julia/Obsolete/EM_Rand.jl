#=
Generate random motifs.

This is required to deal with local
optimization
=#

# Get random positions
motif_start_posits = [rand(1:length(seq) - motif_width) 
    for seq in fasta_file_seqs]

# Get motifs from random positions
current_motif_posits = [seq[pos:pos+motif_width-1] 
    for (seq, pos) in zip(fasta_file_seqs, motif_start_posits)]
