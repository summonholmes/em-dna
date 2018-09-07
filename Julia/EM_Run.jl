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
max_likely_dict = Dict("max_posits" => Int[], 
                       "max_scores" => Float64[], 
                       "max_motifs" => String[])

# Populate the max_likely_dict
for (pre, cur) in zip(seq_cumsum[1:end-1], seq_cumsum[2:end])
    max_pos = findmax(score_motifs[pre:cur])[2]
    push!(max_likely_dict["max_posits"], max_pos)
    push!(max_likely_dict["max_scores"], 
        score_motifs[pre:cur][max_pos])
    push!(max_likely_dict["max_motifs"], 
        join(em_motifs[:, pre:cur][:, max_pos]))
end
