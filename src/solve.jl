"""
    solve(relax, lb, ub; options=Options(), incumbent=Inf, x_incumbent=nothing,
          heuristic=nothing, polish=nothing, branch=most_fractional) -> Result

Branch and bound over integer variables `x ∈ [lb, ub] ∩ ℤⁿ` (minimization). Bounds are
rounded inward to integers, and may be infinite.

# Arguments

- `relax(lb, ub) -> NodeResult`: the node oracle. It solves the continuous relaxation with
  the integer variables restricted to the box `[lb, ub]`. See [`NodeResult`](@ref).
- `incumbent`, `x_incumbent`: optional warm start (objective value and integer point).
  `incumbent = Inf` means none, and `x_incumbent` is then ignored. `-Inf` is an error.
- `heuristic(x, lb, ub) -> (value, xint)` or `nothing`: called on every node taken from the
  queue that has a relaxed solution `x`. It looks for a feasible integer point; `value`
  must be its exact objective value.
- `polish(xint) -> value`, `Inf`, or `nothing`: exact objective at an integer point, e.g.
  a re-solve with the integers fixed. If given, it scores relaxations that are integral
  only within `int_tol`, instead of trusting the relaxed value. `Inf` means the point is
  infeasible and `nothing` undecided. Either way, a node that is not fully fixed is
  branched further. A fully fixed node is closed on `Inf` and left unresolved on
  `nothing` (status `:incomplete`).
- `branch(x, lb, ub, int_tol) -> index or nothing`: branching rule (see
  [`most_fractional`](@ref)). If it returns `nothing` for a node that still has free
  variables, [`first_free`](@ref) is used.

All callbacks receive copies of Sequoia's vectors, so they may modify or keep them.

See [`Options`](@ref) for limits and tolerances and [`Result`](@ref) for the output.
"""
function solve(relax, lb::AbstractVector, ub::AbstractVector;
               options::Options=Options(), incumbent::Real=Inf, x_incumbent=nothing,
               heuristic=nothing, polish=nothing, branch=most_fractional)
    o = validate(options)
    length(lb) == length(ub) ||
        throw(DimensionMismatch("lb and ub have lengths $(length(lb)) and $(length(ub))"))
    (any(isnan, lb) || any(isnan, ub)) && throw(ArgumentError("lb and ub must not contain NaN"))
    (isnan(incumbent) || incumbent == -Inf) &&
        throw(ArgumentError("incumbent must be finite, or Inf for none; got $incumbent"))
    n = length(lb)
    isfinite(incumbent) || (x_incumbent = nothing)
    if x_incumbent !== nothing
        length(x_incumbent) == n ||
            throw(DimensionMismatch("x_incumbent has length $(length(x_incumbent)), expected $n"))
        x_incumbent = Vector{Float64}(x_incumbent)
    end
    lb0 = Float64[ceil(l - o.int_tol) for l in lb]
    ub0 = Float64[floor(u + o.int_tol) for u in ub]

    t = Tree(relax, heuristic, polish, branch, o, n, incumbent, x_incumbent)
    if any(lb0 .> ub0)   # no integer point: infeasible whatever the warm start claims
        return Result(:infeasible, Inf, Inf, nothing, 0, 0, Inf, time() - t.t0, [log!(t, true)])
    end
    root_bound = evaluate!(t, lb0, ub0, -Inf)
    status = :optimal
    while !isempty(t.open)
        if time() - t.t0 > o.time_limit
            status = :time_limit
            break
        elseif t.nodes >= o.node_limit
            status = :node_limit
            break
        end
        node = select_node!(t.open, o.node_selection)
        discard!(t, node.bound) && continue

        if node.x !== nothing && heuristic !== nothing
            h = heuristic(copy(node.x), copy(node.lb), copy(node.ub))
            h === nothing || improve!(t, h[1], h[2])
            discard!(t, node.bound) && continue
        end
        k = branch(node.x === nothing ? nothing : copy(node.x), copy(node.lb), copy(node.ub), o.int_tol)
        if k === nothing   # the rule sees nothing to branch on, but variables may be free
            k = first_free(node.x, node.lb, node.ub, o.int_tol)
        end
        if k === nothing   # fully fixed, and the relaxation failed or was rejected
            t.lost = min(t.lost, node.bound)
            continue
        end
        node.lb[k] < node.ub[k] ||
            throw(ArgumentError("the branching rule returned variable $k, which is fixed in this node"))
        v = split_point(node.x === nothing ? nothing : node.x[k], node.lb[k], node.ub[k])
        for (lbc, ubc) in child_boxes(node, k, v)
            evaluate!(t, lbc, ubc, node.bound)
        end
        log!(t)
    end
    log!(t, true)
    if status === :optimal
        status = discard!(t, t.lost) ? (isfinite(t.incumbent) ? :optimal : :infeasible) : :incomplete
    end
    bound = status in (:optimal, :infeasible) ? min(t.incumbent, t.pruned) : global_bound(t)
    return Result(status, t.incumbent, bound, t.x_incumbent, t.nodes, t.unknown, root_bound,
                  time() - t.t0, t.trace)
end
