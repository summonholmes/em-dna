#=
Starting at EM-Core for Julia:
The primary iteration of the program.

Loop over total_rand_aligns and record the best
results each time.

This is important to the EM due to the likelihood of
local optimization.
=#

# Include all files
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
seq_cumsum = pushfirst!(cumsum(
    [length(i) - motif_width for i in fasta_file_seqs]), 0)

# Get all possible motifs for each sequence
em_motifs = hcat(
    [collect(i[j:j + motif_width - 1]) for i in fasta_file_seqs 
    for j = 1:length(i) - motif_width]...)

# Use list for final results
final_results = []

# Start the Loops
for i = 1:total_rand_aligns
    #=
    Generate random motifs.

    This is required to deal with local
    optimization
    =#

    # Get random positions
    # motif_start_posits = [rand(1:length(seq) - motif_width - 1) 
    #     for seq in fasta_file_seqs]
    motif_start_posits = [3, 17, 5, 26, 17, 6, 2, 24, 21, 21, 3, 10, 8, 21, 35, 21, 15, 
        10, 0, 27, 5, 24, 4, 0, 20, 28, 6, 28, 2].+1

    # Get motifs from random positions
    current_motif_posits = [seq[pos:pos+motif_width-1] 
        for (seq, pos) in zip(fasta_file_seqs, motif_start_posits)]

    for j = 1:total_em_iters
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
            motif_posits_freq['A'][i] += count(
                A, all_motif_bases[i:motif_width:end])
            motif_posits_freq['C'][i] += count(
                C, all_motif_bases[i:motif_width:end])
            motif_posits_freq['T'][i] += count(
                T, all_motif_bases[i:motif_width:end])
            motif_posits_freq['G'][i] += count(
                G, all_motif_bases[i:motif_width:end])
        end

        #=
        All matrix operations for the E step.
        
        Builds on the counting script and
        constructs the log-odds matrix
        for scoring
        =#
        
        # Initialize matrix
        em_log_odds = zeros(4, motif_width + 1)
        
        # Read in count dicts to complete counts matrix
        for (i, base) in enumerate(('A', 'C', 'T', 'G'))
            em_log_odds[i, 1] = motif_bkgd_bases[base]
            em_log_odds[i, 2:end] = motif_posits_freq[base]
        end

        # Get col totals
        col_totals = sum(em_log_odds[:, 2:end], dims=1)

        # Convert background to freq, complete freq matrix
        em_log_odds[:, 1] /= sum(em_log_odds[:, 1])

        # Divide by col_totals to complete freq matrix
        for i = 2:motif_width + 1
            em_log_odds[:, i:i] /= col_totals[i - 1]
        end

        # Complete odds matrix
        for i = 1:4
            em_log_odds[i:i, 2:end] /= em_log_odds[i:i, 1][1]
        end

        # Finish log-odds matrix by taking log2 of odds
        em_log_odds = log2.(em_log_odds[:, 2:end])


        #=
        The scoring operations for the M step.
        
        Score against the log odds, and take the max score
        and assign its new positions to perform ML
        =#
        
        # Init score array
        score_motifs = Array{Any,2}(em_motifs)
        
        # Score each motif position
        for i = 1:motif_width
            score_motifs[i, :] = replace(
                score_motifs[i, :], 'A' => em_log_odds[1, i])
            score_motifs[i, :] = replace(
                score_motifs[i, :], 'C' => em_log_odds[2, i])
            score_motifs[i, :] = replace(
                score_motifs[i, :], 'T' => em_log_odds[3, i])
            score_motifs[i, :] = replace(
                score_motifs[i, :], 'G' => em_log_odds[4, i])
        end

        # Sum the columns and compute 2^x
        score_motifs = exp2.(sum(score_motifs, dims=(1)))

        # Initialize max_likely_dict for each iter
        global max_likely_dict = Dict("max_posits" => Int[], 
                                      "max_scores" => Float64[], 
                                      "max_motifs" => String[])

        # Populate the max_likely_dict
        for (pre, cur) in zip(seq_cumsum[1:end-1].+1, seq_cumsum[2:end])
            max_pos = findmax(score_motifs[pre:cur])[2]
            push!(max_likely_dict["max_posits"], max_pos)
            push!(max_likely_dict["max_scores"], 
                score_motifs[pre:cur][max_pos])
            push!(max_likely_dict["max_motifs"], 
                join(em_motifs[:, pre:cur][:, max_pos]))
        end

        # Finally, update positions
        motif_start_posits = max_likely_dict["max_posits"]
        current_motif_posits = [seq[pos:pos+motif_width-1] 
            for (seq, pos) in zip(
                fasta_file_seqs, motif_start_posits)]
    end

    # Record alignment results after ML iterations
    push!(final_results, max_likely_dict["max_posits"])
    push!(final_results, max_likely_dict["max_scores"])
    push!(final_results, max_likely_dict["max_motifs"])
end

# Display final results
max_sum_results = [sum(i) for i in final_results[2:3:end]]
max_ind = findmax(max_sum_results)
print(final_results[(max_ind[2] * 3) - 2])
print(final_results[(max_ind[2] * 3) - 1])
print(final_results[max_ind[2] * 3])
