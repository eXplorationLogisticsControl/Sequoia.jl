# Sparse SDPA (`.dat-s`) input. The problem is the SDPA dual form
#
#     minimize cᵀx  subject to  Σᵢ Fᵢ xᵢ − F₀ ≽ 0
#
# A positive block is symmetric; a negative block size `-k` is a diagonal block of `k`
# linear inequalities. `*INTEGERS*` and `*<index>` mark 1-based integer variables.

using CodecZlib: GzipDecompressorStream
using LinearAlgebra: Symmetric, dot

"""
    PSDBlock(coef)

One semidefinite block. `coef[1]` is `F₀` and `coef[i+1]` is `Fᵢ`, each an `n × n`
symmetric matrix.
"""
struct PSDBlock
    coef::Vector{Matrix{Float64}}
end

"""
    DiagonalBlock(coef)

One diagonal block, a list of linear inequalities `Σᵢ Fᵢ xᵢ − F₀ ≥ 0` componentwise.
`coef[1]` is `F₀` and `coef[i+1]` is `Fᵢ`, each a vector.
"""
struct DiagonalBlock
    coef::Vector{Vector{Float64}}
end

"""
    SDPAProblem(c, blocks, integers)

A sparse-SDPA problem: `minimize cᵀx` subject to `Σᵢ Fᵢ xᵢ − F₀ ≽ 0`. `integers` holds
the 1-based indices of the integer variables, in increasing order.

`blocks` mixes [`PSDBlock`](@ref) and [`DiagonalBlock`](@ref). A problem may have no
variables: the constraints are then constant.
"""
struct SDPAProblem
    c::Vector{Float64}
    blocks::Vector{Union{PSDBlock,DiagonalBlock}}
    integers::Vector{Int}
end

"""
    read_dat_s(path::AbstractString) -> SDPAProblem
    read_dat_s(io::IO) -> SDPAProblem

Read a sparse SDPA file (`.dat-s`). A path ending in `.gz` is decompressed.

Lines starting with `*` or `"` are comments. `*INTEGERS*`, followed by `*<index>` lines,
names the integer variables. Header lines may carry a trailing note (`28 = m`), and only
their leading numbers are read. Repeated entries at the same position are summed. An
off-diagonal stored on one side only is mirrored. Both sides may be stored; if they
disagree, an `ArgumentError` is thrown.
"""
function read_dat_s(path::AbstractString)
    if endswith(lowercase(path), ".gz")
        return open(GzipDecompressorStream, path) do io
            read_dat_s(io)
        end
    end
    return open(read_dat_s, path)
end

function read_dat_s(io::IO)
    lines = readlines(io)
    pos = Ref(1)

    function nextline()
        while pos[] <= length(lines)
            line = strip(lines[pos[]])
            pos[] += 1
            isempty(line) || return line
        end
        return nothing
    end

    # Next non-comment line, or `:integers` when `*INTEGERS*` is reached.
    function next_data()
        while true
            line = nextline()
            line === nothing && return nothing
            line == "*INTEGERS*" && return :integers
            (startswith(line, '*') || startswith(line, '"')) && continue
            return line
        end
    end

    function require_line(what)
        line = next_data()
        (line === nothing || line === :integers) &&
            throw(ArgumentError("truncated SDPA file: missing $what"))
        return line
    end

    # Leading numeric tokens of successive lines, until `n` values are collected.
    function read_numbers(n::Int, parseone, what)
        vals = Float64[]
        while length(vals) < n
            line = require_line(what)
            took = 0
            for tok in split(line)
                v = parseone(tok)
                v === nothing && break
                push!(vals, v)
                took += 1
                length(vals) == n && break
            end
            took == 0 && throw(ArgumentError("expected $what, got: $line"))
        end
        return vals
    end

    m = Int(read_numbers(1, tryparse_int, "the variable count")[1])
    nblocks = Int(read_numbers(1, tryparse_int, "the block count")[1])
    m >= 0 && nblocks >= 0 ||
        throw(ArgumentError("variable and block counts must be nonnegative, got $m and $nblocks"))
    sizes = Int.(read_numbers(nblocks, tryparse_int, "block sizes"))
    c = read_numbers(m, tryparse_float, "the objective")

    acc, psd, order = allocate_blocks(sizes, m)
    saw_integers = false
    while true
        line = next_data()
        if line === nothing
            break
        elseif line === :integers
            saw_integers = true
            break
        end
        add_entry!(acc, psd, order, m, line)
    end

    integers = saw_integers ? read_integers!(nextline, m) : Int[]
    blocks = build_blocks(acc, psd)
    return SDPAProblem(c, blocks, integers)
end

"Parse a base-10 integer token, or `nothing` when `tok` is not one."
tryparse_int(tok) = tryparse(Int, tok)

"Parse a floating-point token, or `nothing` when `tok` is not one."
tryparse_float(tok) = tryparse(Float64, tok)

function allocate_blocks(sizes, m)
    nblocks = length(sizes)
    psd = Vector{Bool}(undef, nblocks)
    order = Vector{Int}(undef, nblocks)
    acc = Vector{Any}(undef, nblocks)
    for (b, s) in enumerate(sizes)
        s == 0 && throw(ArgumentError("block $b has size 0"))
        psd[b] = s > 0
        order[b] = abs(s)
        n = order[b]
        acc[b] = if psd[b]
            [(zeros(n, n), falses(n, n)) for _ in 0:m]
        else
            [zeros(n) for _ in 0:m]
        end
    end
    return acc, psd, order
end

function add_entry!(acc, psd, order, m, line)
    toks = split(line)
    length(toks) >= 5 || throw(ArgumentError("SDPA entry needs 5 numbers, got: $line"))
    mat = tryparse_int(toks[1])
    blk = tryparse_int(toks[2])
    i = tryparse_int(toks[3])
    j = tryparse_int(toks[4])
    val = tryparse_float(toks[5])
    (mat === nothing || blk === nothing || i === nothing || j === nothing || val === nothing) &&
        throw(ArgumentError("SDPA entry needs 5 numbers, got: $line"))
    (0 <= mat <= m) || throw(ArgumentError("matrix index $mat is outside 0:$m"))
    (1 <= blk <= length(acc)) || throw(ArgumentError("block index $blk is outside 1:$(length(acc))"))
    n = order[blk]
    (1 <= i <= n && 1 <= j <= n) ||
        throw(ArgumentError("position ($i,$j) is outside block $blk of order $n"))
    if psd[blk]
        vals, written = acc[blk][mat + 1]
        vals[i, j] += val
        written[i, j] = true
    else
        i == j || throw(ArgumentError("diagonal block $blk has an off-diagonal entry ($i,$j)"))
        acc[blk][mat + 1][i] += val
    end
    return
end

function read_integers!(nextline, m)
    found = Int[]
    while true
        line = nextline()
        line === nothing && break
        if startswith(line, '"')
            continue
        end
        matched = match(r"^\*\s*(\d+)\s*$", line)
        if matched !== nothing
            push!(found, parse(Int, matched.captures[1]))
            continue
        end
        startswith(line, '*') && continue
        throw(ArgumentError("unexpected line after *INTEGERS*: $line"))
    end
    ints = sort!(unique(found))
    for k in ints
        1 <= k <= m || throw(ArgumentError("integer index $k is outside 1:$m"))
    end
    return ints
end

function build_blocks(acc, psd)
    blocks = Union{PSDBlock,DiagonalBlock}[]
    for b in eachindex(acc)
        if psd[b]
            coef = Matrix{Float64}[]
            for (vals, written) in acc[b]
                symmetrize!(vals, written, b)
                push!(coef, vals)
            end
            push!(blocks, PSDBlock(coef))
        else
            push!(blocks, DiagonalBlock(acc[b]))
        end
    end
    return blocks
end

"Mirror a one-sided off-diagonal. Both sides may be present; they must agree."
function symmetrize!(vals::Matrix{Float64}, written::BitMatrix, block::Int)
    n = size(vals, 1)
    for i in 1:n, j in (i + 1):n
        wij, wji = written[i, j], written[j, i]
        if wij && wji
            isapprox(vals[i, j], vals[j, i]; atol=1e-8, rtol=1e-8) ||
                throw(ArgumentError(
                    "block $block has inconsistent entries at ($i,$j) and ($j,$i)"))
        elseif wij
            vals[j, i] = vals[i, j]
        elseif wji
            vals[i, j] = vals[j, i]
        end
    end
    return vals
end

"""
    jump_model(problem::SDPAProblem) -> JuMP.Model

JuMP model of `problem`, minimized. Variables are created in file order. Those listed in
`problem.integers` are integer; the file gives them no bounds, so any bounds come from
the diagonal blocks. [`solve`](@ref) on this model returns the integer variables in
increasing index order.
"""
function jump_model(problem::SDPAProblem)
    m = length(problem.c)
    model = JuMP.Model()
    x = JuMP.@variable(model, [1:m])
    for k in problem.integers
        JuMP.set_integer(x[k])
    end
    # An empty dot product is not an expression JuMP will accept as an objective.
    m == 0 ? JuMP.@objective(model, Min, 0.0) : JuMP.@objective(model, Min, dot(problem.c, x))
    for block in problem.blocks
        add_block!(model, x, block)
    end
    return model
end

"Constant affine expression. A bare number is not a PSD-cone argument."
function affine_constant(value::Real)
    expr = zero(JuMP.AffExpr)
    JuMP.add_to_expression!(expr, value)
    return expr
end

function add_block!(model, x, block::PSDBlock)
    F = block.coef
    m = length(x)
    if m == 0
        # The empty sum is a numeric matrix, and `@constraint` rejects that in a PSD cone.
        n = size(F[1], 1)
        expr = [affine_constant(-F[1][i, j]) for i in 1:n, j in 1:n]
        JuMP.@constraint(model, Symmetric(expr) in JuMP.PSDCone())
    else
        JuMP.@constraint(model, Symmetric(sum(x[i] * F[i + 1] for i in 1:m) - F[1]) in JuMP.PSDCone())
    end
    return model
end

function add_block!(model, x, block::DiagonalBlock)
    F = block.coef
    m = length(x)
    for j in eachindex(F[1])
        if m == 0
            row = affine_constant(-F[1][j])
            JuMP.@constraint(model, row >= 0)
        else
            JuMP.@constraint(model, sum(x[i] * F[i + 1][j] for i in 1:m) - F[1][j] >= 0)
        end
    end
    return model
end
