# Linear presolve for the JuMP interface. Singleton rows on integer variables become
# bounds. What remains is propagated, at the root and at every node, by activity: a
# row `lo ≤ aᵀx ≤ hi` cannot hold unless each variable stays inside the bounds implied
# by the others. Continuous bounds are tracked only here; they are never written back
# to the model. Integer bounds are returned and become the node's box.

"One linear row, `lo ≤ Σ coef[t] · x[idx[t]] ≤ hi`, with the constant folded into the sides."
struct LinearRow
    idx::Vector{Int}
    coef::Vector{Float64}
    lo::Float64
    hi::Float64
end

"""
    LinearRows

Linear part of a relaxed model, indexed like `JuMP.all_variables`.

- `lb`, `ub`: static bounds. Integer entries are replaced by the node box during
  propagation; continuous entries are tightened in place by singleton rows and, during a
  call, by activity.
- `integers`: position of each integer variable, in the box order.
- `integer_index`: the reverse map, `0` for a continuous variable.
- `rows`: affine rows still to propagate. Integer singleton rows are removed once they
  have been turned into bounds.
- `cons`: the JuMP constraint of each row, parallel to `rows` until singletons are
  extracted.
"""
struct LinearRows
    lb::Vector{Float64}
    ub::Vector{Float64}
    integers::Vector{Int}
    integer_index::Vector{Int}
    rows::Vector{LinearRow}
    cons::Vector{Any}
end

const LINEAR_SET = Union{MOI.GreaterThan{Float64},MOI.LessThan{Float64},
                         MOI.EqualTo{Float64},MOI.Interval{Float64}}

"Relative change below which a continuous bound is left as it is, so passes terminate."
const CONTINUOUS_TOL = 1e-3

"Activity may overshoot a right-hand side by this fraction before the row is infeasible."
const ACTIVITY_TOL = 1e-6

"""
    linear_rows(model, ints) -> LinearRows

Affine rows of `model` whose set is `GreaterThan`, `LessThan`, `EqualTo` or `Interval`.
`ints` are the integer variables, in box order. Variable bounds and fixings are stored as
static bounds, not as rows.
"""
function linear_rows(model::JuMP.Model, ints)
    vars = JuMP.all_variables(model)
    vindex = Dict(v => i for (i, v) in enumerate(vars))
    viindex = Dict(JuMP.index(v) => i for (i, v) in enumerate(vars))
    integer_index = zeros(Int, length(vars))
    integers = Vector{Int}(undef, length(ints))
    for (k, v) in enumerate(ints)
        j = vindex[v]
        integer_index[j] = k
        integers[k] = j
    end
    lb, ub = static_bounds(model, vars, vindex)
    rows = LinearRow[]
    cons = Any[]
    for (F, S) in JuMP.list_of_constraint_types(model)
        S <: LINEAR_SET || continue
        for c in JuMP.all_constraints(model, F, S)
            parsed = linear_terms(JuMP.constraint_object(c).func, vindex, viindex)
            parsed === nothing && continue
            acc, constant = parsed
            lo, hi = shifted_bounds(JuMP.constraint_object(c).set, constant)
            push!(rows, linear_row(acc, lo, hi))
            push!(cons, c)
        end
    end
    return LinearRows(lb, ub, integers, integer_index, rows, cons)
end

"Variable bounds, fixings, and `VariableRef`-in-`Interval` constraints."
function static_bounds(model, vars, vindex)
    n = length(vars)
    lb = fill(-Inf, n)
    ub = fill(Inf, n)
    for (i, v) in enumerate(vars)
        if JuMP.is_fixed(v)
            lb[i] = ub[i] = JuMP.fix_value(v)
        else
            JuMP.has_lower_bound(v) && (lb[i] = JuMP.lower_bound(v))
            JuMP.has_upper_bound(v) && (ub[i] = JuMP.upper_bound(v))
        end
    end
    for c in JuMP.all_constraints(model, JuMP.VariableRef, MOI.Interval{Float64})
        obj = JuMP.constraint_object(c)
        i = vindex[obj.func]
        lb[i] = max(lb[i], obj.set.lower)
        ub[i] = min(ub[i], obj.set.upper)
    end
    return lb, ub
end

function add_coef!(acc::Dict{Int,Float64}, j::Int, a::Real)
    iszero(a) && return acc
    acc[j] = get(acc, j, 0.0) + a
    iszero(acc[j]) && delete!(acc, j)
    return acc
end

"Affine terms as a coefficient map, or `nothing` when `func` is not an affine row."
linear_terms(::Any, ::Any, ::Any) = nothing

function linear_terms(func::JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, vindex, ::Any)
    acc = Dict{Int,Float64}()
    for (v, a) in func.terms
        add_coef!(acc, vindex[v], a)
    end
    return acc, Float64(func.constant)
end

function linear_terms(func::MOI.ScalarAffineFunction{Float64}, ::Any, viindex)
    acc = Dict{Int,Float64}()
    for term in func.terms
        add_coef!(acc, viindex[term.variable], term.coefficient)
    end
    return acc, Float64(func.constant)
end

function linear_row(acc::Dict{Int,Float64}, lo::Float64, hi::Float64)
    idx = Vector{Int}(undef, length(acc))
    coef = Vector{Float64}(undef, length(acc))
    for (t, (j, a)) in enumerate(acc)
        idx[t] = j
        coef[t] = a
    end
    return LinearRow(idx, coef, lo, hi)
end

shifted_bounds(set::MOI.GreaterThan{Float64}, constant) = (set.lower - constant, Inf)
shifted_bounds(set::MOI.LessThan{Float64}, constant) = (-Inf, set.upper - constant)
function shifted_bounds(set::MOI.EqualTo{Float64}, constant)
    value = set.value - constant
    return (value, value)
end
function shifted_bounds(set::MOI.Interval{Float64}, constant)
    return (set.lower - constant, set.upper - constant)
end

"Bounds on `x` implied by `lo ≤ a·x ≤ hi`."
function implied_bounds(a::Float64, lo::Float64, hi::Float64)
    if a > 0
        l = isfinite(lo) ? lo / a : -Inf
        u = isfinite(hi) ? hi / a : Inf
    else
        l = isfinite(hi) ? hi / a : -Inf
        u = isfinite(lo) ? lo / a : Inf
    end
    return l, u
end

"""
    extract_singletons!(relaxed, rows, ints, lb, ub)

Turn a row on a single integer variable into a tightening of `[lb, ub]` and delete it from
`relaxed`. A row on a single continuous variable tightens that variable's static bound and
stays in the model. `ints` are the integer variables in box order.
"""
function extract_singletons!(relaxed::JuMP.Model, rows::LinearRows, ints, lb, ub)
    (length(ints) == length(lb) == length(ub) == length(rows.integers)) ||
        throw(DimensionMismatch("integer box and integer variables differ in length"))
    keep_rows = LinearRow[]
    keep_cons = Any[]
    for (row, con) in zip(rows.rows, rows.cons)
        if length(row.idx) == 1
            j = row.idx[1]
            l, u = implied_bounds(row.coef[1], row.lo, row.hi)
            k = rows.integer_index[j]
            if k != 0
                lb[k] = max(lb[k], l)
                ub[k] = min(ub[k], u)
                JuMP.delete(relaxed, con)
                continue
            end
            rows.lb[j] = max(rows.lb[j], l)
            rows.ub[j] = min(rows.ub[j], u)
        end
        push!(keep_rows, row)
        push!(keep_cons, con)
    end
    empty!(rows.rows)
    append!(rows.rows, keep_rows)
    empty!(rows.cons)
    append!(rows.cons, keep_cons)
    return rows
end

min_infinite(a, l, u) = a > 0 ? !isfinite(l) : !isfinite(u)
max_infinite(a, l, u) = a > 0 ? !isfinite(u) : !isfinite(l)
min_contrib(a, l, u) = a > 0 ? a * l : a * u
max_contrib(a, l, u) = a > 0 ? a * u : a * l

"Finite part of the row activity, and how many infinite terms were left out of it."
function row_activity(row::LinearRow, lb, ub)
    minact = 0.0
    maxact = 0.0
    nmininf = 0
    nmaxinf = 0
    for (j, a) in zip(row.idx, row.coef)
        l, u = lb[j], ub[j]
        if min_infinite(a, l, u)
            nmininf += 1
        else
            minact += min_contrib(a, l, u)
        end
        if max_infinite(a, l, u)
            nmaxinf += 1
        else
            maxact += max_contrib(a, l, u)
        end
    end
    return minact, maxact, nmininf, nmaxinf
end

"Whether `act` (from `ninf` infinite terms) falls on the infeasible side of `rhs`."
function activity_infeasible(act, ninf, rhs, ::Val{:above})
    (ninf == 0 && isfinite(rhs)) || return false
    return act > rhs + ACTIVITY_TOL * max(1.0, abs(rhs))
end
function activity_infeasible(act, ninf, rhs, ::Val{:below})
    (ninf == 0 && isfinite(rhs)) || return false
    return act < rhs - ACTIVITY_TOL * max(1.0, abs(rhs))
end

"""
Bounds on `x` implied by one term `a·x` of a row, given the activity of the whole row.
`-Inf` or `Inf` means that side gives no bound.
"""
function inferred_bounds(a, l, u, row::LinearRow, minact, maxact, nmininf, nmaxinf)
    new_lb = -Inf
    new_ub = Inf
    if isfinite(row.hi) && (nmininf == 0 || (nmininf == 1 && min_infinite(a, l, u)))
        resid = nmininf == 0 ? minact - min_contrib(a, l, u) : minact
        rhs = (row.hi - resid) / a
        a > 0 ? (new_ub = rhs) : (new_lb = rhs)
    end
    if isfinite(row.lo) && (nmaxinf == 0 || (nmaxinf == 1 && max_infinite(a, l, u)))
        resid = nmaxinf == 0 ? maxact - max_contrib(a, l, u) : maxact
        rhs = (row.lo - resid) / a
        a > 0 ? (new_lb = rhs) : (new_ub = rhs)
    end
    return new_lb, new_ub
end

function accept_bound(old, new, continuous, tighter)
    (isfinite(new) && tighter(new, old)) || return false
    continuous || return true
    isfinite(old) || return true
    return abs(new - old) > CONTINUOUS_TOL * max(1.0, abs(old))
end

"Tighten every variable of `row` from its activity. `nothing` means the row is infeasible."
function tighten_row!(row::LinearRow, rows::LinearRows, full_lb, full_ub, out_lb, out_ub, int_tol)
    minact, maxact, nmininf, nmaxinf = row_activity(row, full_lb, full_ub)
    activity_infeasible(minact, nmininf, row.hi, Val(:above)) && return nothing
    activity_infeasible(maxact, nmaxinf, row.lo, Val(:below)) && return nothing
    changed = false
    for (j, a) in zip(row.idx, row.coef)
        l, u = full_lb[j], full_ub[j]
        cand_lb, cand_ub = inferred_bounds(a, l, u, row, minact, maxact, nmininf, nmaxinf)
        k = rows.integer_index[j]
        if k != 0
            if isfinite(cand_lb)
                rounded = ceil(cand_lb - int_tol)
                if rounded > out_lb[k]
                    out_lb[k] = rounded
                    changed = true
                end
            end
            if isfinite(cand_ub)
                rounded = floor(cand_ub + int_tol)
                if rounded < out_ub[k]
                    out_ub[k] = rounded
                    changed = true
                end
            end
            out_lb[k] > out_ub[k] && return nothing
            full_lb[j] = out_lb[k]
            full_ub[j] = out_ub[k]
        else
            if accept_bound(l, cand_lb, true, >)
                full_lb[j] = cand_lb
                changed = true
            end
            if accept_bound(u, cand_ub, true, <)
                full_ub[j] = cand_ub
                changed = true
            end
            full_lb[j] > full_ub[j] && return nothing
        end
    end
    return changed
end

"""
    propagate_linear(rows, lb, ub; int_tol, max_passes=10) -> (lb, ub) or nothing

Propagate `rows` over the integer box `[lb, ub]`. Integer bounds are rounded inward with
`int_tol`; continuous bounds move only inside `rows`, and only when the change exceeds a
relative tolerance, so that the passes stop. Returns `nothing` when the box is empty.
The input vectors are not modified.
"""
function propagate_linear(rows::LinearRows, lb, ub; int_tol::Real, max_passes::Int=10)
    nint = length(rows.integers)
    (length(lb) == nint && length(ub) == nint) ||
        throw(DimensionMismatch("integer box has lengths $(length(lb)) and $(length(ub)), expected $nint"))
    full_lb = copy(rows.lb)
    full_ub = copy(rows.ub)
    out_lb = Vector{Float64}(undef, nint)
    out_ub = Vector{Float64}(undef, nint)
    for (k, j) in enumerate(rows.integers)
        out_lb[k] = Float64(lb[k])
        out_ub[k] = Float64(ub[k])
        out_lb[k] > out_ub[k] && return nothing
        full_lb[j] = out_lb[k]
        full_ub[j] = out_ub[k]
    end
    for j in eachindex(full_lb)
        rows.integer_index[j] == 0 && full_lb[j] > full_ub[j] && return nothing
    end
    for _ in 1:max_passes
        changed = false
        for row in rows.rows
            step = tighten_row!(row, rows, full_lb, full_ub, out_lb, out_ub, Float64(int_tol))
            step === nothing && return nothing
            changed |= step
        end
        changed || break
    end
    return (out_lb, out_ub)
end

"""
    integral_objective(model, ints) -> Bool

Whether the objective of `model` is affine with an integer constant and integer
coefficients, all of them on the integer variables `ints`. Such an objective takes integer
values at every integer-feasible point.
"""
function integral_objective(model::JuMP.Model, ints)
    return _integral_objective(JuMP.objective_function(model), Set(ints))
end

_integral_objective(value::Real, ::Set) = isinteger(value)
_integral_objective(v::JuMP.VariableRef, allowed::Set) = v in allowed

function _integral_objective(func::JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, allowed::Set)
    isinteger(func.constant) || return false
    for (v, a) in func.terms
        iszero(a) && continue
        (isinteger(a) && v in allowed) || return false
    end
    return true
end

_integral_objective(::Any, ::Set) = false
