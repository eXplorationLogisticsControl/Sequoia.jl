# Primal heuristics that plug into `solve(...; heuristic)`.

"""
    rounding_heuristic(evaluate; rule=:random, rounds=1, rng=nothing, every=1, cache=true)
        -> heuristic(x, lb, ub)

Rounding with a fixed-integer re-solve ("round, fix and solve"), the pattern of SCIP-SDP's
`sdprand` and `sdpfracround` heuristics. The node's relaxed point `x` is rounded to an
integer point `xint` inside the node box `[lb, ub]`. Then `evaluate(xint)` solves the problem
with every integer variable fixed, so the continuous variables are re-optimized; rounding
the relaxed continuous values as well would rarely be feasible.

- `evaluate(xint) -> value`, `Inf` (infeasible) or `nothing` (undecided): the exact objective
  with the integers fixed at `xint`, the same contract as `polish` in [`solve`](@ref).
  [`fixed_oracle`](@ref) builds one from the node oracle.
- `rule`:
  - `:nearest`: round each fractional value to the nearest integer. SCIP-SDP's fractional
    rounding fixes in order of increasing fractionality and propagates after every fixing;
    without propagation the order does not matter and this is plain rounding.
  - `:random`: randomized rounding (Raghavan & Thompson, 1987): round `xᵢ` up with
    probability `xᵢ - ⌊xᵢ⌋`, so `E[xint] = x`.
- `rounds`: random draws per call (`:random` only); the best feasible point is returned.
- `rng`: random number generator for `:random` (e.g. `Random.Xoshiro(1)` for reproducible
  runs); `nothing` uses the global one.
- `every`: run on every `every`-th call only. Each call costs at least one fixed-integer solve,
  which can cost as much as the node relaxation itself.
- `cache`: remember evaluated points, so a point reached again costs nothing.

The returned function has the `heuristic` signature of [`solve`](@ref): it returns
`(value, xint)` for the best feasible point found, or `nothing`.
"""
function rounding_heuristic(evaluate; rule::Symbol=:random, rounds::Integer=1, rng=nothing,
                            every::Integer=1, cache::Bool=true)
    rule in (:nearest, :random) || throw(ArgumentError("rule must be :nearest or :random, got :$rule"))
    rounds >= 1 && every >= 1 || throw(ArgumentError("rounds and every must be positive"))
    seen = Dict{Vector{Float64},Union{Nothing,Float64}}()
    calls = Ref(0)
    draw() = rng === nothing ? rand() : rand(rng)
    score(xint) = cache ? get!(() -> _value(evaluate(xint)), seen, xint) : _value(evaluate(xint))
    return function (x, lb, ub)
        calls[] += 1
        (calls[] - 1) % every == 0 || return nothing
        best = nothing
        for _ in 1:(rule === :nearest ? 1 : rounds)
            xint = similar(x)
            for i in eachindex(x)
                f = floor(x[i])
                up = rule === :nearest ? x[i] - f >= 0.5 : draw() < x[i] - f
                xint[i] = clamp(up ? f + 1 : f, lb[i], ub[i])
            end
            v = score(xint)
            v !== nothing && isfinite(v) && (best === nothing || v < best[1]) && (best = (v, xint))
        end
        return best
    end
end

_value(v::Nothing) = nothing
_value(v::Real) = isnan(v) ? nothing : Float64(v)

"""
    fixed_oracle(relax) -> evaluate(xint)

Exact objective at an integer point from the node oracle: `relax(xint, xint)` is the problem
with every integer variable fixed. Returns the primal value, `Inf` if infeasible, or `nothing`
if the solve failed. Use with [`rounding_heuristic`](@ref) or as `polish`.
"""
fixed_oracle(relax) = function (xint)
    r = relax(copy(xint), copy(xint))
    return r.status === :optimal ? r.value : r.status === :infeasible ? Inf : nothing
end
