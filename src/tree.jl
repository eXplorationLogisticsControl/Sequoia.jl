"""
    Tree

Mutable state of one branch-and-bound run (minimization). Users never see it: it is
created and consumed by [`solve`](@ref).

Every point of the original box lies in exactly one of these sets, and each comes with a
lower bound on its points. Their minimum is the global bound:

- `open`: nodes still to be processed, each with its own `bound`.
- `pruned`: discarded subtrees. Their smallest (rounded) bound is kept, so that a gap
  tolerance never turns into an overstated final bound.
- `lost`: fully fixed nodes that could not be resolved. Their smallest bound is kept.
- closed subtrees whose best point is known (infeasible, or integral relaxation). These are
  bounded by the incumbent.
"""
mutable struct Tree{R,H,P,B}
    const relax::R
    const heuristic::H
    const polish::P
    const branch::B
    const options::Options
    const nvars::Int
    const t0::Float64
    incumbent::Float64
    x_incumbent::Union{Nothing,Vector{Float64}}
    nodes::Int
    unknown::Int
    const open::Vector{Node}
    lost::Float64
    pruned::Float64
    const trace::Vector{Snapshot}
    last_log::Float64
end

function Tree(relax, heuristic, polish, branch, options::Options, nvars::Integer,
              incumbent::Real, x_incumbent)
    return Tree(relax, heuristic, polish, branch, options, Int(nvars), time(),
                Float64(incumbent), x_incumbent, 0, 0, Node[], Inf, Inf, Snapshot[], -Inf)
end

"""
Bound used for pruning. When the objective is known to be integral it is rounded up, after
subtracting a slack of `max(int_tol, abs_gap, rel_gap·|b|)`. The slack keeps a dual bound
that overshoots by the solver's (relative) accuracy from jumping a whole unit.
"""
function rounded_bound(t::Tree, b::Real)
    o = t.options
    (o.integral_objective && isfinite(b)) || return Float64(b)
    return ceil(b - max(o.int_tol, o.abs_gap, o.rel_gap * abs(b)))
end

"Pruning threshold: a node with rounded bound at least this cannot improve the incumbent."
function prune_threshold(t::Tree)
    inc, o = t.incumbent, t.options
    return isfinite(inc) ? inc - max(o.abs_gap, o.rel_gap * abs(inc)) : inc
end

"Whether a node with lower bound `b` cannot improve the incumbent (beyond the gap tolerance)."
is_closed(t::Tree, b::Real) = rounded_bound(t, b) >= prune_threshold(t)

"""
    discard!(t, b) -> Bool

If a node with bound `b` is closed, record its bound for the final global bound and return
`true`. Otherwise return `false`.
"""
function discard!(t::Tree, b::Real)
    is_closed(t, b) || return false
    t.pruned = min(t.pruned, rounded_bound(t, b))
    return true
end

"""
    improve!(t, value, x) -> Bool

Make `(value, round.(x))` the incumbent if `value` is finite and better. Then drop the open
nodes that the new incumbent closes. Returns whether the incumbent changed.
"""
function improve!(t::Tree, value::Real, x)
    (isfinite(value) && value < t.incumbent) || return false
    length(x) == t.nvars ||
        throw(DimensionMismatch("integer point has length $(length(x)), expected $(t.nvars)"))
    t.incumbent = value
    t.x_incumbent = round.(Float64.(x))
    filter!(n -> !discard!(t, n.bound), t.open)
    return true
end

"Current global lower bound."
function global_bound(t::Tree)
    b = min(t.incumbent, t.lost, t.pruned)
    for n in t.open
        b = min(b, n.bound)
    end
    return b
end

"""
    evaluate!(t, lb, ub, parent_bound) -> Float64

Solve the relaxation of the node `[lb, ub]` and act on the result:

- infeasible: close the node.
- failed (`:unknown`, or a relaxed `x` with `NaN`/`Inf` entries): keep the node with its
  parent's bound, to be branched on.
- integral: score the integer point (`polish`, or the relaxed value) and update the
  incumbent. The subtree is closed only if the node's bound then closes it. Otherwise the
  point was rejected or left undecided, or its score lies above the node's bound, and the
  rest of the subtree may hold a better point, so it stays open. A fully fixed node holds
  only this point, so it is closed once scored.
- fractional: keep the node unless its bound closes it.

The oracle and the callbacks get copies, so they cannot corrupt the tree by mutating or
reusing vectors. Returns the relaxation's own bound: `Inf` if infeasible, `-Inf` if failed
or unavailable.
"""
function evaluate!(t::Tree, lb::Vector{Float64}, ub::Vector{Float64}, parent_bound::Float64)
    r = t.relax(copy(lb), copy(ub))
    r isa NodeResult || throw(ArgumentError("the node oracle must return a NodeResult, got $(typeof(r))"))
    t.nodes += 1
    r.status === :infeasible && return Inf
    r.status in (:optimal, :unknown) ||
        throw(ArgumentError("NodeResult status must be :optimal, :infeasible or :unknown, got :$(r.status)"))
    if r.status === :optimal
        r.x === nothing && throw(ArgumentError("a NodeResult with status :optimal must carry x"))
        length(r.x) == t.nvars ||
            throw(DimensionMismatch("NodeResult x has length $(length(r.x)), expected $(t.nvars)"))
    end
    if r.status === :unknown || !all(isfinite, r.x)
        t.unknown += 1
        discard!(t, parent_bound) || push!(t.open, Node(lb, ub, parent_bound, nothing))
        return -Inf
    end
    x = copy(r.x)   # the oracle may reuse its buffer
    own = isnan(r.bound) ? -Inf : r.bound
    b = max(own, parent_bound)

    if t.branch(copy(x), copy(lb), copy(ub), t.options.int_tol) === nothing
        # integral relaxation: the integer point closest to x is the candidate
        xi = round.(x)
        v = t.polish === nothing ? r.value : t.polish(copy(xi))
        decided = v !== nothing && !isnan(v)
        decided && isfinite(v) && improve!(t, v, xi)
        if decided && isfinite(v) && discard!(t, b)
            # bound meets incumbent: nothing better in the subtree
        elseif has_free(lb, ub)
            # rejected, undecided, or scored above the bound: the rest stays open
            discard!(t, b) || push!(t.open, Node(lb, ub, b, nothing))
        elseif !decided
            t.lost = min(t.lost, b)   # fully fixed and undecided
        end                           # fully fixed and scored or rejected: closed
    else
        discard!(t, b) || push!(t.open, Node(lb, ub, b, x))
    end
    return own
end

"Record a progress snapshot, and print it if `log_every` seconds have passed (or `force`)."
function log!(t::Tree, force::Bool=false)
    now = time() - t.t0
    snap = (t=now, nodes=t.nodes, bound=global_bound(t), incumbent=t.incumbent, open=length(t.open))
    push!(t.trace, snap)
    if isfinite(t.options.log_every) && (force || now - t.last_log >= t.options.log_every)
        t.last_log = now
        @printf("  nodes=%6d  t=%7.1fs  inc=%.6g  LB=%.6g  open=%d\n",
                t.nodes, now, t.incumbent, snap.bound, length(t.open))
    end
    return snap
end
