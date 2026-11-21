"""
    Options(; kwargs...)
    Options(o::Options; kwargs...)

Settings for [`solve`](@ref). The second form copies `o` with some fields replaced, e.g.
`Options(o; time_limit=60)`.

# Fields

- `time_limit::Float64 = Inf`: wall-clock limit in seconds, checked between nodes.
- `node_limit::Int = typemax(Int)`: stop once this many relaxations have been solved. A
  branching step solves two children, so the count can exceed the limit by one.
- `rel_gap::Float64 = 1e-6`, `abs_gap::Float64 = 1e-8`: a node with lower bound `b` is
  pruned once `b ≥ inc - max(abs_gap, rel_gap * |inc|)`, where `inc` is the incumbent
  value. With the default values, `:optimal` means optimal within these tolerances.
- `int_tol::Float64 = 1e-6`: a relaxed value within `int_tol` of an integer counts as
  integral.
- `integral_objective::Bool = false`: declare that the objective takes integer values at
  every feasible point. Lower bounds are then rounded up before pruning, to
  `ceil(b - max(int_tol, abs_gap, rel_gap·|b|))`.
- `node_selection::Symbol = :best`: `:best` (best-bound first) or `:depth` (depth first,
  last-in first-out).
- `log_every::Float64 = 5.0`: seconds between progress lines (`Inf` for silent).
"""
Base.@kwdef struct Options
    time_limit::Float64 = Inf
    node_limit::Int = typemax(Int)
    rel_gap::Float64 = 1e-6
    abs_gap::Float64 = 1e-8
    int_tol::Float64 = 1e-6
    integral_objective::Bool = false
    node_selection::Symbol = :best
    log_every::Float64 = 5.0
end

function Options(o::Options; kwargs...)
    fields = (; (f => getfield(o, f) for f in fieldnames(Options))...)
    return Options(; merge(fields, values(kwargs))...)
end

const NODE_SELECTIONS = (:best, :depth)

"Throw an `ArgumentError` if `o` holds a value [`solve`](@ref) cannot work with."
function validate(o::Options)
    o.node_selection in NODE_SELECTIONS ||
        throw(ArgumentError("node_selection must be one of $NODE_SELECTIONS, got :$(o.node_selection)"))
    0 <= o.int_tol < 0.5 || throw(ArgumentError("int_tol must be in [0, 0.5), got $(o.int_tol)"))
    o.rel_gap >= 0 && o.abs_gap >= 0 ||
        throw(ArgumentError("rel_gap and abs_gap must be nonnegative, got $(o.rel_gap) and $(o.abs_gap)"))
    o.time_limit >= 0 || throw(ArgumentError("time_limit must be nonnegative, got $(o.time_limit)"))
    o.node_limit >= 0 || throw(ArgumentError("node_limit must be nonnegative, got $(o.node_limit)"))
    o.log_every > 0 || throw(ArgumentError("log_every must be positive (Inf for silent), got $(o.log_every)"))
    return o
end
