"""
    NodeResult(status, bound, value, x)
    NodeResult(status)

Result of one node relaxation, returned by a node oracle (see [`solve`](@ref)).

- `status`: `:optimal`, `:infeasible`, or `:unknown`. `:unknown` means the solver
  failed. A failed node is never pruned: it keeps its parent's bound and is branched on.
- `bound`: a valid lower bound on the node relaxation, e.g. the dual objective value. Use
  `-Inf` (or `NaN`) if none is available; the node then inherits its parent's bound.
- `value`: the primal objective of the relaxed solution. It becomes an incumbent value
  when `x` is integral and no `polish` function is given.
- `x`: relaxed values of the integer variables, in the order of `lb` and `ub`.

The one-argument form is for `:infeasible` and `:unknown` results.
"""
struct NodeResult
    status::Symbol
    bound::Float64
    value::Float64
    x::Union{Nothing,Vector{Float64}}
end
NodeResult(status::Symbol) = NodeResult(status, NaN, NaN, nothing)

"One progress snapshot: time, nodes solved, global lower bound, incumbent, open nodes."
const Snapshot = NamedTuple{(:t, :nodes, :bound, :incumbent, :open),Tuple{Float64,Int,Float64,Float64,Int}}

"""
    Result

Outcome of [`solve`](@ref).

- `status`:
  - `:optimal`: the incumbent is optimal within `rel_gap`/`abs_gap`.
  - `:infeasible`: the tree was exhausted without a feasible point.
  - `:time_limit`, `:node_limit`: a limit was reached first.
  - `:incomplete`: the tree was exhausted, but a fully fixed node could not be resolved
    (its relaxation failed, or `polish` returned `nothing`) and might hold a better
    point than the incumbent.
- `objective`, `x`: the incumbent value and integer point (`Inf` and `nothing` if there is
  none). `x` can be `nothing` while `objective` is finite if the warm-start `incumbent`
  was given without `x_incumbent` and never improved.
- `bound`: a global lower bound on the optimal value. Once the search finishes it can lie
  below `objective` by up to the gap tolerance.
- `nodes`: relaxations solved. `unknown`: relaxations that returned `:unknown`.
- `root_bound`: the root relaxation's bound (`-Inf` if it failed, `Inf` if infeasible).
- `time`: wall-clock seconds.
- `trace`: `(t, nodes, bound, incumbent, open)` snapshots, for gap-vs-time plots.

The JuMP interface reports `objective`, `bound`, `root_bound` and `trace` in the model's own
objective sense.
"""
struct Result
    status::Symbol
    objective::Float64
    bound::Float64
    x::Union{Nothing,Vector{Float64}}
    nodes::Int
    unknown::Int
    root_bound::Float64
    time::Float64
    trace::Vector{Snapshot}
end

"""
    gap(r::Result) -> Float64

Relative optimality gap `|objective - bound| / max(|objective|, 1e-12)`. Works for either
objective sense. `Inf` if there is no incumbent.

An `:optimal` result satisfies `|objective - bound| ≤ max(abs_gap, rel_gap·|objective|)`.
When the objective is (near) zero, `abs_gap` decides, and `gap` can be large.
"""
gap(r::Result) = isfinite(r.objective) ? abs(r.objective - r.bound) / max(abs(r.objective), 1e-12) : Inf

function Base.show(io::IO, ::MIME"text/plain", r::Result)
    println(io, "Sequoia.Result")
    println(io, "  status:     ", r.status)
    println(io, "  objective:  ", r.objective)
    println(io, "  bound:      ", r.bound, isfinite(r.objective) ? "  (gap $(gap(r)))" : "")
    println(io, "  nodes:      ", r.nodes, r.unknown > 0 ? " ($(r.unknown) failed)" : "")
    println(io, "  root bound: ", r.root_bound)
    print(io, "  time:       ", round(r.time; sigdigits=4), " s")
end
