"""
    frequencies_thz(vals, masses; unit=:thz)

Convert NMA eigenvalues into frequencies.

Assumptions:
- `vals` come from a mass-weighted Hessian diagonalization
- Hessian units are consistent with:
    energy  = kJ/mol
    length  = Å
    mass    = amu

Supported units:
- `:thz`  → terahertz
- `:cm1`  → wavenumbers (cm^-1)
- `:hz`   → hertz
- `:radps` → angular frequency in rad/s
- `:arb`  → raw sqrt(λ), arbitrary units

Negative eigenvalues are returned as imaginary frequencies.

Note:
`masses` is included for API consistency and validation workflows.
It is not used directly here if `vals` are already from the mass-weighted Hessian.
"""
function frequencies_thz(vals, masses; unit=:thz)
    isempty(vals) && return ComplexF64[]

    # Physical constants
    NA   = 6.02214076e23          # mol^-1
    amu  = 1.66053906660e-27      # kg
    ang  = 1.0e-10                # m
    c_cm = 2.99792458e10          # cm/s

    # Conversion:
    # sqrt(kJ/mol / (amu * Å^2)) -> rad/s
    internal_to_radps = sqrt((1000.0 / NA) / (amu * ang^2))

    factor = if unit == :radps
        internal_to_radps
    elseif unit == :hz
        internal_to_radps / (2π)
    elseif unit == :thz
        internal_to_radps / (2π) / 1.0e12
    elseif unit == :cm1
        internal_to_radps / (2π * c_cm)
    elseif unit == :arb
        1.0
    else
        throw(ArgumentError("Unsupported unit: $unit. Use :thz, :cm1, :hz, :radps, or :arb"))
    end

    freqs = Vector{ComplexF64}(undef, length(vals))

    for (i, λ) in enumerate(vals)
        if λ < 0
            freqs[i] = complex(0.0, sqrt(abs(λ)) * factor)
        else
            freqs[i] = complex(sqrt(λ) * factor, 0.0)
        end
    end

    return freqs
end