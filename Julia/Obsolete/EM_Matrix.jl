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
