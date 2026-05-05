"""
    compute_frequencies(vals; skip_rigid=true)

Compute frequencies from eigenvalues.

- Negative eigenvalues → `missing` (imaginary modes)
- By default skips first 6 rigid-body modes (optional)
"""
function compute_frequencies(vals; skip_rigid=false)
    freqs = Vector{Union{Missing, Float64}}(undef, length(vals))

    for (i, λ) in enumerate(vals)
        if λ < 0
            freqs[i] = missing
        else
            freqs[i] = sqrt(λ)
        end
    end

    return freqs
end