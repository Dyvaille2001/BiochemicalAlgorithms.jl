############################################################
# nma_report.jl  -  small table helpers (console + Markdown + TSV)
############################################################
using Printf

fmt(x::AbstractString) = String(x)
fmt(x::Integer)        = string(x)
fmt(x::Bool)           = x ? "yes" : "no"
fmt(::Nothing)         = ""
fmt(::Missing)         = ""
function fmt(x::AbstractFloat)
    isnan(x) && return "nan"
    isinf(x) && return x > 0 ? "inf" : "-inf"
    # integer-valued floats (counts that were promoted to Float64 inside a numeric row)
    (x == round(x) && abs(x) < 1e9) && return string(Int(x))
    (x != 0 && (abs(x) < 1e-3 || abs(x) >= 1e5)) && return @sprintf("%.3e", x)
    return @sprintf("%.4f", x)
end
fmt(x) = string(x)

"Print a Markdown-style aligned table."
function print_table(io::IO, title::AbstractString, header::Vector{String}, rows)
    println(io, "\n", title)
    cells = [[fmt(c) for c in r] for r in rows]
    w = [maximum(vcat(length(header[j]), [length(c[j]) for c in cells]))
         for j in eachindex(header)]
    line(r) = "| " * join([rpad(r[j], w[j]) for j in eachindex(r)], " | ") * " |"
    println(io, line(header))
    println(io, "|" * join(["-"^(w[j] + 2) for j in eachindex(header)], "|") * "|")
    for c in cells
        println(io, line(c))
    end
    return nothing
end

"Write a tab-separated table."
function save_tsv(path::AbstractString, header::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(header, '\t'))
        for r in rows
            println(io, join([fmt(c) for c in r], '\t'))
        end
    end
    return path
end

"""
    report_table(title, header, rows; md=nothing, tsv=nothing)

Print to the console and optionally append to a Markdown file / write a TSV.
"""
function report_table(title::AbstractString, header::Vector{String}, rows;
                      md=nothing, tsv=nothing)
    print_table(stdout, title, header, rows)
    md  === nothing || open(io -> print_table(io, title, header, rows), md, "a")
    tsv === nothing || save_tsv(tsv, header, rows)
    return nothing
end