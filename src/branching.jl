# A branching rule is any callable `rule(x, lb, ub, int_tol) -> Union{Int,Nothing}`. It
# returns the index of a variable with `lb[i] < ub[i]` to branch on, or `nothing` when `x`
# is integral (within `int_tol`) on all free variables. `x === nothing` means the node's
# relaxation failed or was rejected, and the rule must still pick a free variable if there
# is one.

"""
    most_fractional(x, lb, ub, int_tol) -> Union{Int,Nothing}

Branching rule: pick the free variable whose relaxed value is farthest from an integer
(ties go to the lowest index). Returns `nothing` if all free variables are integral within
`int_tol`. If `x === nothing`, falls back to [`first_free`](@ref).
"""
function most_fractional(x, lb, ub, int_tol)
    x === nothing && return first_free(x, lb, ub, int_tol)
    best, idx = float(int_tol), nothing
    for i in eachindex(x, lb, ub)
        f = abs(x[i] - round(x[i]))
        if f > best && lb[i] < ub[i]
            best, idx = f, i
        end
    end
    return idx
end

"""
    first_free(x, lb, ub, int_tol) -> Union{Int,Nothing}

Branching rule: pick the first variable that is not fixed (`lb[i] < ub[i]`), whatever its
relaxed value. This is the fallback when a node has no relaxed solution.
"""
first_free(x, lb, ub, int_tol) = findfirst(i -> lb[i] < ub[i], eachindex(lb, ub))

"""
    split_point(xk, l, u) -> Float64

Non-integer split value `v` for a variable with integer bounds `l < u` and relaxed value
`xk` (`nothing` if there is none). The children get `x ≤ floor(v)` and `x ≥ ceil(v)`, so
both are non-empty and both exclude `xk` when it is fractional. Throws an `ArgumentError`
if `|v| ≥ 2^53`, where consecutive integers are no longer representable as `Float64`.
"""
function split_point(xk, l, u)
    v = if xk !== nothing && isfinite(xk)
        clamp(xk, l, u)
    elseif isfinite(l) && isfinite(u)
        (l + u) / 2
    elseif isfinite(l)
        l
    elseif isfinite(u)
        u
    else
        0.0
    end
    abs(v) < maxintfloat(Float64) ||
        throw(ArgumentError("cannot branch on a value of magnitude ≥ 2^53 (got $v); bound the integer variables"))
    return isinteger(v) ? (v < u ? v + 0.5 : v - 0.5) : float(v)
end
