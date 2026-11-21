const SOLVED_TERMINATION = (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED)
const SOLVED_PRIMAL = (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)

"""
    jump_result(model, xvars) -> NodeResult

Read the result of an optimized JuMP model into a [`NodeResult`](@ref), in the
minimization sense (values are negated for a `Max` model). Useful when writing your own
node oracle.

- solved: `OPTIMAL`, `LOCALLY_SOLVED`, `ALMOST_OPTIMAL` or `ALMOST_LOCALLY_SOLVED` with a
  (nearly) feasible primal point, or `SLOW_PROGRESS` (e.g. MOSEK stalling at its last digit)
  with *feasible* primal and dual points.
- infeasible: `INFEASIBLE`, `LOCALLY_INFEASIBLE`, or `SLOW_PROGRESS` with a dual
  infeasibility certificate.
- unknown: anything else, including unbounded relaxations.

The bound is the dual objective value if the solver reports a feasible dual point. Otherwise
it is `-Inf`, and the node inherits its parent's bound. A solver configured not to return
duals therefore gets no pruning by bound. `x` holds the values of `xvars`, and plain numbers
in `xvars` are passed through.

!!! note
    A dual objective value is a valid bound only up to the solver's feasibility tolerances.
    `LOCALLY_SOLVED` and `LOCALLY_INFEASIBLE` are trusted like their global counterparts,
    which is correct only for convex relaxations.
"""
function jump_result(model::JuMP.Model, xvars)
    ts = JuMP.termination_status(model)
    ps = JuMP.primal_status(model)
    ds = JuMP.dual_status(model)
    solved = (ts in SOLVED_TERMINATION && ps in SOLVED_PRIMAL) ||
             (ts == MOI.SLOW_PROGRESS && ps == MOI.FEASIBLE_POINT && ds == MOI.FEASIBLE_POINT)
    if solved
        sgn = JuMP.objective_sense(model) == MOI.MAX_SENSE ? -1.0 : 1.0
        bound = ds == MOI.FEASIBLE_POINT ? sgn * JuMP.dual_objective_value(model) : -Inf
        x = Float64[v isa Number ? v : JuMP.value(v) for v in xvars]
        return NodeResult(:optimal, bound, sgn * JuMP.objective_value(model), x)
    elseif ts in (MOI.INFEASIBLE, MOI.LOCALLY_INFEASIBLE) ||
           (ts == MOI.SLOW_PROGRESS && ds == MOI.INFEASIBILITY_CERTIFICATE)
        return NodeResult(:infeasible)
    end
    return NodeResult(:unknown)
end

"Constraint sets that `relax_integrality` plus variable-bound branching cannot handle."
function unsupported_set(S)
    return S <: Union{MOI.SOS1,MOI.SOS2,MOI.Semicontinuous,MOI.Semiinteger,MOI.Indicator}
end

"Integer bounds of variable `v` as given in the model (binaries are clipped to [0, 1])."
function integer_bounds(v::JuMP.VariableRef)
    JuMP.is_fixed(v) && return (JuMP.fix_value(v), JuMP.fix_value(v))
    l = JuMP.has_lower_bound(v) ? JuMP.lower_bound(v) : -Inf
    u = JuMP.has_upper_bound(v) ? JuMP.upper_bound(v) : Inf
    return JuMP.is_binary(v) ? (max(l, 0.0), min(u, 1.0)) : (l, u)
end

"""
    move_interval_bounds!(relaxed, ints, lb, ub)

Turn `VariableRef`-in-`Interval` constraints on the integer variables `ints` of the relaxed
copy into entries of the integer box `[lb, ub]`, and delete the constraints. Branching then
controls these bounds like any other, and they cannot clash with the bounds that
`relax_integrality` and branching set.
"""
function move_interval_bounds!(relaxed::JuMP.Model, ints, lb, ub)
    position = Dict(v => i for (i, v) in enumerate(ints))
    for c in JuMP.all_constraints(relaxed, JuMP.VariableRef, MOI.Interval{Float64})
        obj = JuMP.constraint_object(c)
        i = get(position, obj.func, nothing)
        i === nothing && continue
        lb[i] = max(lb[i], obj.set.lower)
        ub[i] = min(ub[i], obj.set.upper)
        JuMP.delete(relaxed, c)
    end
    return lb, ub
end

"Restrict the relaxed variable `v` to `[l, u]`, replacing any previous bounds or fixing."
function set_box!(v::JuMP.VariableRef, l::Float64, u::Float64)
    JuMP.is_fixed(v) && JuMP.unfix(v)
    if l == u
        JuMP.fix(v, l; force=true)
        return v
    end
    if isfinite(l)
        JuMP.set_lower_bound(v, l)
    elseif JuMP.has_lower_bound(v)
        JuMP.delete_lower_bound(v)
    end
    if isfinite(u)
        JuMP.set_upper_bound(v, u)
    elseif JuMP.has_upper_bound(v)
        JuMP.delete_upper_bound(v)
    end
    return v
end

"""
    solve(model::JuMP.Model; optimizer, kwargs...) -> Result

Branch and bound on a JuMP model with `Bin`/`Int` variables. The model itself is not
modified: Sequoia relaxes a copy and solves each node with `optimizer` (e.g.
`Clarabel.Optimizer`, or `optimizer_with_attributes(...)`). Each node's relaxation is read
with [`jump_result`](@ref).

Other keyword arguments are as in the core [`solve`](@ref solve(::Any, ::AbstractVector, ::AbstractVector)).
`x` and the vectors passed to `heuristic`/`polish` are ordered like the integer variables in
`JuMP.all_variables(model)`. `incumbent`, the values returned by `heuristic` and `polish`,
and the `objective`, `bound`, `root_bound` and `trace` of the result are all in the model's
own objective sense. An infinite `incumbent` (of either sign) means no warm start.

Bounds of the integer variables may be given as variable bounds, by fixing, or as
`x in MOI.Interval(l, u)` constraints.

SOS, semicontinuous, semi-integer and indicator constraints are not supported.
"""
function solve(model::JuMP.Model; optimizer, kwargs...)
    for (F, S) in JuMP.list_of_constraint_types(model)
        unsupported_set(S) && throw(ArgumentError("constraints of type $F-in-$S are not supported"))
    end
    ints0 = [v for v in JuMP.all_variables(model) if JuMP.is_binary(v) || JuMP.is_integer(v)]
    boxes = integer_bounds.(ints0)
    lb0, ub0 = first.(boxes), last.(boxes)

    relaxed, refmap = JuMP.copy_model(model)
    ints = [refmap[v] for v in ints0]
    move_interval_bounds!(relaxed, ints, lb0, ub0)
    JuMP.relax_integrality(relaxed)
    JuMP.set_optimizer(relaxed, optimizer)

    function relax(lb, ub)
        for (v, l, u) in zip(ints, lb, ub)
            set_box!(v, l, u)
        end
        JuMP.optimize!(relaxed)
        return jump_result(relaxed, ints)
    end

    sgn = JuMP.objective_sense(model) == MOI.MAX_SENSE ? -1.0 : 1.0
    sgn == 1 && return solve(relax, lb0, ub0; kwargs...)

    # maximization: the core minimizes -f, and user callbacks speak in the model's sense.
    # An infinite incumbent of either sign means "none".
    kw = Dict{Symbol,Any}(kwargs)
    haskey(kw, :incumbent) && (kw[:incumbent] = isinf(kw[:incumbent]) ? Inf : -kw[:incumbent])
    if get(kw, :heuristic, nothing) !== nothing
        h0 = kw[:heuristic]
        kw[:heuristic] = (x, l, u) -> (h = h0(x, l, u); h === nothing ? h : (-h[1], h[2]))
    end
    if get(kw, :polish, nothing) !== nothing
        p0 = kw[:polish]
        kw[:polish] = x -> (v = p0(x); v === nothing ? v : -v)
    end
    r = solve(relax, lb0, ub0; kw...)
    trace = Snapshot[(; s..., bound=-s.bound, incumbent=-s.incumbent) for s in r.trace]
    return Result(r.status, -r.objective, -r.bound, r.x, r.nodes, r.unknown, -r.root_bound,
                  r.time, trace)
end
