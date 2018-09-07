#=
The primary iteration of the program.

Loop over total_rand_aligns and record the best
results each time.

This is important to the EM due to the likelihood of
local optimization.
=#

# Include all files
include("EM_Input.jl")
include("EM_Prep.jl")

# Use list for final results
final_results = []

# Start the Loops
for i = 1:total_rand_aligns
    include("EM_Rand.jl")
    for j = 1:total_em_iters
        include("EM_Count.jl")
        include("EM_Matrix.jl")
        include("EM_Run.jl")
        current_motif_posits = max_likely_dict["max_posits"]
    end
    push!(final_results, max_likely_dict["max_posits"])
    push!(final_results, max_likely_dict["max_scores"])
    push!(final_results, max_likely_dict["max_motifs"])
end
