#=
All input operations are handled here.

Once the parameters are set, the file is loaded
with each line as an array element.

Finally, the sequences are filtered out.
=#

# Set params
motif_width = 6
total_rand_aligns = 50
total_em_iters = 50

# Load file
fasta_file_seqs = open("example.fasta") do fasta
    fasta_file_seqs = readlines(fasta)
end

# Process
fasta_file_seqs = fasta_file_seqs[2:2:end]
