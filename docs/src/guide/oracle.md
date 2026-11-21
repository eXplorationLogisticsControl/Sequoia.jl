# Custom node oracles

The core entry point is

```julia
solve(relax, lb, ub; options, incumbent, x_incumbent, heuristic, polish, branch)
```

Sequoia only knows the integer variables `x ∈ [lb, ub] ∩ ℤⁿ`. Everything else about the
problem lives in your **node oracle** `relax(lb, ub)`. It must:

1. solve the continuous relaxation with the integer variables restricted to `[lb, ub]`;
2. return a [`NodeResult`](@ref) (in the **minimization** sense):
   - `NodeResult(:optimal, bound, value, x)`: `bound` is a valid lower bound (a dual
     objective value is the safe choice), `value` the primal objective, and `x` the relaxed
     values of the integer variables;
   - `NodeResult(:infeasible)`;
   - `NodeResult(:unknown)` if the solver failed. Sequoia never prunes such a node: it
     keeps its parent's bound and is branched on. An `x` containing `NaN` or `Inf` is
     treated the same way.

The oracle receives copies of the bounds, and Sequoia copies `x`, so you may reuse a
solver's work vectors.

Because the oracle gets the bounds, it may build a **different model per node**. For
example, it can drop the constraint blocks of indicators fixed to zero, or reject a box
that violates a cardinality budget without calling a solver at all.

## A first example

Minimize a separable convex quadratic over integers. The relaxation over a box has a
closed form (projection), so the bound is exact:

```@example oracle
using Sequoia

t = [0.3, 0.8, 0.45, 2.6, -1.3]
function relax(lb, ub)
    x = clamp.(t, lb, ub)
    f = sum(abs2, x .- t)
    return NodeResult(:optimal, f, f, x)
end

result = solve(relax, [0, 0, 0, -5, -5], [1, 1, 1, 5, 5]; options=Options(log_every=Inf))
result.x
```

## Oracles built on JuMP

[`jump_result`](@ref) converts an optimized JuMP model into a `NodeResult`. It maps solver
statuses conservatively (e.g. a stalled solve counts only if both primal and dual points
are feasible) and takes the dual objective as the bound:

```@example oracle
using JuMP, Clarabel

A = [3.0 1 2; 1 4 1]
b = [4.0, 3]
c = [-5.0, -4, -3]

function lp_relax(lb, ub)
    m = Model(optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false))
    @variable(m, lb[i] <= x[i=1:3] <= ub[i])
    @constraint(m, A * x .<= b)
    @objective(m, Min, c' * x)
    optimize!(m)
    return jump_result(m, x)
end

result = solve(lp_relax, zeros(3), fill(3.0, 3); options=Options(log_every=Inf))
(result.objective, result.x)
```

For speed, build the model once outside the oracle and only update the variable bounds
inside it, as [`solve(model::JuMP.Model)`](@ref solve(::JuMP.Model)) does.

## Bound propagation

`propagate(lb, ub)` is called before every relaxation, including the root. It returns a
tighter box `(lb, ub)`, or `nothing` when the box holds no feasible point. `nothing`, or a
box with no integer point, closes the node without a solve and without counting it.
Sequoia intersects the returned box with the node's box, so a bound can only get tighter,
and then rounds it inward the same way as the root bounds. The oracle, the heuristic and
the branching rule all see the tightened box.

```julia
# a cardinality cut the relaxation does not know about
propagate(lb, ub) = sum(ub) < 2 ? nothing : (lb, min.(ub, 1.0))
solve(relax, lb, ub; propagate)
```

## Heuristics and polishing

- `heuristic(x, lb, ub) -> (value, xint)` or `nothing` is called on each node taken from the
  queue, with its relaxed solution. Rounding plus a feasibility repair is typical. `value`
  must be the *exact* objective value of the feasible integer point `xint`; Sequoia trusts
  it.
- `polish(xint) -> value | Inf | nothing` gives the exact objective of an integer point,
  usually by re-solving with the integers fixed. With `polish`, a relaxation that is
  integral only within `int_tol` no longer yields its own (approximate) value as an
  incumbent. `Inf` rejects the point (infeasible), and `nothing` means "could not decide".
  In both cases the node is branched further if anything is still free. A fully fixed,
  undecided node makes the final status `:incomplete` instead of `:optimal`.

```julia
polish(xint) = begin
    fix_and_solve!(model, xint)           # your code
    r = jump_result(model, [])
    r.status === :optimal ? r.value : r.status === :infeasible ? Inf : nothing
end
```

## Custom branching rules

A branching rule is a function `branch(x, lb, ub, int_tol) -> Union{Int,Nothing}`. It
returns the index of a free variable (`lb[i] < ub[i]`) to branch on, or `nothing` when `x`
is integral. `x === nothing` means the node has no usable relaxed solution; the rule must
then still pick a free variable, and Sequoia falls back to [`first_free`](@ref) if it does
not. The built-in rules are [`most_fractional`](@ref) (the default) and `first_free`.

```julia
# branch on the fractional variable with the largest objective weight
weighted(x, lb, ub, tol) = x === nothing ? first_free(x, lb, ub, tol) :
    let cand = [i for i in eachindex(x) if lb[i] < ub[i] && abs(x[i] - round(x[i])) > tol]
        isempty(cand) ? nothing : cand[argmax(abs.(c[cand]))]
    end
solve(relax, lb, ub; branch=weighted)
```
