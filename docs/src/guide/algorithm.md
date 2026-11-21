# Algorithm

Sequoia runs a textbook LP-style branch and bound over the integer variables. Every node is
a box `[lb, ub]` of integer bounds. Internally, everything is in the minimization sense.

## Main loop

1. **Root.** Round the bounds inward (`ceil(lb - int_tol)`, `floor(ub + int_tol)`) and
   solve the root relaxation.
2. **Select** an open node: best bound first (`node_selection = :best`, ties go to the
   newest node) or depth first (`:depth`, last in, first out).
3. **Prune** it if its bound cannot improve the incumbent, i.e. if
   `b ≥ inc − max(abs_gap, rel_gap·|inc|)`. With `integral_objective = true`, `b` is first
   rounded up to `ceil(b − max(int_tol, abs_gap, rel_gap·|b|))`. The slack keeps a dual
   bound that overshoots by the solver's accuracy from jumping a whole unit.
4. **Heuristic.** If the node has a relaxed solution, call `heuristic` and update the
   incumbent. A better incumbent immediately removes every open node it closes.
5. **Branch** on the variable `k` chosen by the branching rule. With split value `v`, the
   children get `x_k ≤ floor(v)` and `x_k ≥ ceil(v)`. Each child's relaxation is solved
   right away (eager evaluation), and the child is stored only if it still needs work.

The loop stops when no open node is left, or when `time_limit` or `node_limit` is reached.

## What happens to a solved child

| Relaxation result | Action |
|:--|:--|
| `:infeasible` | closed |
| `:unknown` (solver failure), or `x` containing `NaN`/`Inf` | kept with the **parent's** bound and no relaxed point; it is branched with [`first_free`](@ref) |
| `:optimal`, fractional `x` | kept with bound `max(own bound, parent bound)` |
| `:optimal`, integral `x` | the integer point `round.(x)` is scored by `polish`, or by the relaxed `value` if there is no `polish`, and may become the incumbent (see below) |

After an integral point has been scored, its subtree is closed only if the node's lower
bound `max(own bound, parent bound)` now meets the pruning test. When the relaxation is
solved accurately, the bound and the score agree, and the subtree closes at once. The node
is instead kept and branched further, while any variable is free, when:

- `polish` rejects the point (`Inf`) or cannot decide (`nothing`), or the relaxed `value`
  is `NaN`;
- the score lies above the node's bound by more than the gap tolerance. This happens with
  a stalled solve whose primal value is far above its dual bound, a missing dual bound
  (`-Inf`), or a `polish` value that disagrees with the relaxation. Other points of the
  subtree may then be better.

A fully fixed node holds only its one point, so it is closed once scored. If it cannot be
decided, it is recorded as *unresolved*.

The oracle, `heuristic`, `polish` and the branching rule all receive copies of Sequoia's
vectors, so callbacks that modify or reuse vectors cannot corrupt the tree.

## Statuses and guarantees

| Status | Meaning |
|:--|:--|
| `:optimal` | Tree exhausted. `objective − bound ≤ max(abs_gap, rel_gap·\|objective\|)` (within the oracle's accuracy). When the optimum is (near) zero, `abs_gap` decides, and the *relative* [`gap`](@ref) can be large. |
| `:infeasible` | Tree exhausted with no feasible point (and no warm start). |
| `:incomplete` | Tree exhausted, but an unresolved fully fixed node has a bound below the incumbent. The incumbent is feasible but not certified optimal. `bound` accounts for the unresolved node. |
| `:time_limit`, `:node_limit` | A limit was reached. `bound` is the smallest bound of the open, unresolved and pruned nodes. |

The reported `bound` is always a valid global lower bound (upper bound for maximization),
provided the oracle's bounds are valid. When nodes are pruned within the gap tolerance,
their bounds are kept, so `bound` may lie slightly below `objective` even at `:optimal`.
[`gap`](@ref) reports the remaining relative gap.

## Accuracy caveats

Branch and bound is only as exact as its node solver:

- Conic interior-point solvers return dual objective values that are valid bounds only up
  to their feasibility tolerances. Keep `rel_gap` above the solver tolerance.
- Without `polish`, an integral relaxation is scored by its own primal value. That value
  may be slightly infeasible (`ALMOST_OPTIMAL`, `NEARLY_FEASIBLE_POINT`). Use `polish` when
  certification matters.
- Unbounded relaxations are reported as `:unknown` by [`jump_result`](@ref). With infinite
  integer bounds, the tree may then keep growing until a limit is hit. Give integer
  variables finite bounds where possible.
