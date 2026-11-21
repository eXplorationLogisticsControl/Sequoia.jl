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
    const log_sign::Float64
end

function Tree(relax, heuristic, polish, branch, options::Options, nvars::Integer,
              incumbent::Real, x_incumbent, log_sign::Real=1.0)
    return Tree(relax, heuristic, polish, branch, options, Int(nvars), time(),
                Float64(incumbent), x_incumbent, 0, 0, Node[], Inf, Inf, Snapshot[], -Inf,
                Float64(log_sign))
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

const PROGRESS_HEADER =
    "    Nodes    |          Objective Bounds          |  Work   \n" *
    " Expl Unexpl |     Incumbent        BestBd    Gap |     Time\n"

"Right-aligned objective, or `-` when the value is not finite."
format_log_objective(x::Real, width::Int) = isfinite(x) ? @sprintf("%*.*e", width, 5, x) : lpad("-", width)

"Relative gap as a percent, or `-` when the incumbent or the bound is not finite."
function format_log_gap(incumbent::Real, bound::Real)
    (isfinite(incumbent) && isfinite(bound)) || return lpad("-", 6)
    percent = 100 * abs(incumbent - bound) / max(abs(incumbent), 1e-12)
    return @sprintf("%5.2f%%", percent)
end

"One progress row. `log_sign` maps the internal minimization values onto the displayed sense."
function progress_row(nodes::Integer, open::Integer, incumbent::Real, bound::Real, t::Real;
                      log_sign::Real=1.0)
    inc, bd = log_sign * incumbent, log_sign * bound
    return @sprintf("%6d %5d |%s %s %s | %7.1fs", nodes, open,
                    format_log_objective(inc, 14), format_log_objective(bd, 13),
                    format_log_gap(inc, bd), t)
end

"Objective in the closing summary: twelve digits, or `-` when it is not finite."
format_summary_objective(x::Real) = isfinite(x) ? @sprintf("%.12e", x) : "-"

"Relative gap in the closing summary, with four decimals, or `-`."
function format_summary_gap(incumbent::Real, bound::Real)
    (isfinite(incumbent) && isfinite(bound)) || return "-"
    percent = 100 * abs(incumbent - bound) / max(abs(incumbent), 1e-12)
    return @sprintf("%.4f%%", percent)
end

"Status sentence of the closing summary."
function summary_status(status::Symbol, rel_gap::Real)
    status === :optimal &&
        return "Optimal solution found (tolerance $(@sprintf("%.2e", rel_gap)))"
    status === :infeasible && return "Model is infeasible"
    status === :time_limit && return "Time limit reached"
    status === :node_limit && return "Node limit reached"
    status === :incomplete && return "Optimization incomplete"
    return "Stopped with status :$status"
end

"""
    termination_summary(status, nodes, seconds, objective, bound, rel_gap; log_sign=1, unknown=0)

Gurobi-style closing summary. `objective` and `bound` are minimization values; `log_sign`
maps them onto the displayed sense. `unknown` counts relaxations that failed.
"""
function termination_summary(status::Symbol, nodes::Integer, seconds::Real,
                             objective::Real, bound::Real, rel_gap::Real;
                             log_sign::Real=1.0, unknown::Integer=0)
    noun = nodes == 1 ? "node" : "nodes"
    failed = unknown > 0 ? " ($unknown failed)" : ""
    lines = ["Explored $nodes $noun$failed in $(@sprintf("%.2f", seconds)) seconds",
             "",
             summary_status(status, rel_gap)]
    if status !== :infeasible
        inc, bd = log_sign * objective, log_sign * bound
        push!(lines, "Best objective $(format_summary_objective(inc)), " *
                     "best bound $(format_summary_objective(bd)), " *
                     "gap $(format_summary_gap(inc, bd))")
    end
    return join(lines, "\n")
end

"Print the closing summary when progress logging is on."
function log_summary!(t::Tree, status::Symbol, bound::Real, elapsed::Real)
    isfinite(t.options.log_every) || return nothing
    println()
    println(termination_summary(status, t.nodes, elapsed, t.incumbent, bound, t.options.rel_gap;
                                log_sign=t.log_sign, unknown=t.unknown))
    return nothing
end

"Record a progress snapshot, and print it if `log_every` seconds have passed (or `force`)."
function log!(t::Tree, force::Bool=false)
    now = time() - t.t0
    snap = (t=now, nodes=t.nodes, bound=global_bound(t), incumbent=t.incumbent, open=length(t.open))
    push!(t.trace, snap)
    if isfinite(t.options.log_every) && (force || now - t.last_log >= t.options.log_every)
        isfinite(t.last_log) || print(PROGRESS_HEADER)
        t.last_log = now
        println(progress_row(t.nodes, snap.open, t.incumbent, snap.bound, now; log_sign=t.log_sign))
    end
    return snap
end
