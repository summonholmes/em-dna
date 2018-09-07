#=
All counting operations for the E step.

Take counts of bases, motif bases, and 
positional counts.
=#

# Merge all motifs to a single string
all_motif_bases = join(current_motif_posits)

# Copy all bases to motif_bkgd_bases
motif_bkgd_bases = Dict('A' => bkgd_bases['A'], 
                        'C' => bkgd_bases['C'], 
                        'T' => bkgd_bases['T'], 
                        'G' => bkgd_bases['G'])

# Count all background bases
motif_bkgd_bases['A'] -= count(A, all_motif_bases)
motif_bkgd_bases['C'] -= count(C, all_motif_bases)
motif_bkgd_bases['T'] -= count(T, all_motif_bases)
motif_bkgd_bases['G'] -= count(G, all_motif_bases)

# Init positional pseudocounts
motif_posits_freq = Dict('A' => ones(motif_width), 
                         'C' => ones(motif_width), 
                         'T' => ones(motif_width), 
                         'G' => ones(motif_width))

# Count positional bases
for i = 1:motif_width
    motif_posits_freq['A'][i] = count(
        A, all_motif_bases[i:motif_width:end])
    motif_posits_freq['C'][i] = count(
        C, all_motif_bases[i:motif_width:end])
    motif_posits_freq['T'][i] = count(
        T, all_motif_bases[i:motif_width:end])
    motif_posits_freq['G'][i] = count(
        G, all_motif_bases[i:motif_width:end])
end
