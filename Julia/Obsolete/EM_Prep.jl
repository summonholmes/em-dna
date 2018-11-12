#=
Initial counting operations for the E step.

Take counts of bases, motif bases, and 
positional counts.  Run only once.
=#

# Merge all sequences to a single string
all_bases = join(fasta_file_seqs)

# Init total bases dict
bkgd_bases = Dict('A' => 0, 'C' => 0, 'T' => 0, 'G' => 0)

# Create anonymous counting functions
A = A -> A == 'A'
C = C -> C == 'C'
T = T -> T == 'T'
G = G -> G == 'G'

# Count all bases
bkgd_bases['A'] = count(A, all_bases)
bkgd_bases['C'] = count(C, all_bases)
bkgd_bases['T'] = count(T, all_bases)
bkgd_bases['G'] = count(G, all_bases)

# Get cumsum to track scores and motifs
seq_cumsum = cumsum(
    [length(i) - motif_width for i in fasta_file_seqs])

# Get all possible motifs for each sequence
em_motifs = hcat(
    [collect(i[j:j + motif_width - 1]) for i in fasta_file_seqs 
    for j = 1:length(i) - motif_width]...)
