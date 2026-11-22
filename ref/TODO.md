# Speed

The node relaxation dominates both CPU time and peak memory. `solve(::JuMP.Model)` keeps one relaxed copy and, at every node, rewrites the integer bounds in `set_box!` (`src/jump.jl`) and calls `JuMP.optimize!`. Clarabel's MathOptInterface wrapper implements `copy_to` and does not implement `modify`, so a bound change makes JuMP build a new `Clarabel.Solver`.

Where a node's time goes (TrussTopology, Clarabel 0.11.1, chordal decomposition off):

| Instance | Node | Clarabel `solve!` | Rebuild (`copy_to` + `setup!`) |
| --- | ---: | ---: | ---: |
| 2x4_2scen_3bars | 12.5 ms | 10.5 ms, 20.5 iterations | 1.7 ms |
| 2x5_1scen_12bars | 37.6 ms | 25.9 ms, 25 iterations | 11 ms |
| 5x5_1bar | 1.5 s | 1.48 s, 22 iterations | 31 ms |

The interior-point iterations dominate at every size. Sequoia's own work per node (bounds, propagation, bookkeeping) is under 0.5 ms.

## One persistent solver (done)

`ext/SequoiaClarabelExt.jl`, on by default (`reuse_solver = true`). It does not keep the iterate: Clarabel's `solve!` always starts from `solver_default_start!`, and the numeric factorization runs on every iteration anyway. Only the copy, the equilibration and the symbolic factorization are saved.

| Instance | Rebuild every node | Reused solver |
| --- | ---: | ---: |
| 2x4_2scen_3bars (400 nodes) | 12.4 ms | 10.0 ms |
| 3x3_3scen (200 nodes) | 21.0 ms | 18.2 ms |
| 2x5_1scen_12bars (60 nodes) | 35.0 ms | 26.9 ms |
| 5x5_1bar (10 nodes) | 1.47 s | 1.44 s |
| bridge_3x9_2bars (4 nodes) | 2.73 s | 2.64 s |

Iteration counts are the same in both modes.

- The extension builds its own cache and bridge layers around a wrapped Clarabel optimizer. The wrapper accepts a constant change on a `Nonnegatives` row through `Clarabel.update_b!`, which scales the new values with the stored equilibration. Any other change raises a `NotAllowedError`, and MOI's caching layer rebuilds the solver as before. That covers adding or deleting a row, and a solver whose presolve or chordal decomposition reshaped the problem.
- MOI's variable `VectorizeBridge` is removed from that bridge layer. It adds a bounded variable as `x = y + l` with `y ≥ 0`, so a new `l` changes every row holding `x`, and MOI refuses the change.
- `set_box!` keeps two bound constraints on a fixed integer instead of calling `JuMP.fix`, so the model's structure does not change.
- Only the bounds that changed since the last solve are written. Consecutive solves differ in a median of 4 bounds (at most 25), not one: best-first search jumps across the tree, and propagation tightens more than the branched variable.
- An infinite integer bound has no stable row. A placeholder above `1e20` turns on Clarabel's presolve, which blocks updates. When such a bound becomes finite, the solver is rebuilt.
- Clarabel 0.9 has no `is_data_update_allowed`, so there the solver is always rebuilt.

## Match the barrier tolerance to the gap

Clarabel's default `tol_gap_rel` is `1e-8`. Sequoia's default `rel_gap` is `1e-6`. A dual bound is valid only up to the solver's feasibility tolerance. Loosening `tol_gap_rel` and `tol_gap_abs` to `1e-6` saves about one iteration in 20 on 2x4_2scen_3bars, with the root bound unchanged to ten digits. Leave `tol_feas` tight enough that the dual objective is still safe to prune on.

## Fewer nodes

- 19 of the 39 TrussTopology runs end the 60 s with no incumbent, so bound pruning never starts. A diving or rounding `heuristic` on the relaxed solution produces a feasible point and lets `discard!` drop open nodes. The callback is already part of `solve`.
- Linear presolve does not see semidefinite rows. Propagate the PSD blocks as well, and close a box before calling `optimize!`.

## Parallel nodes

The relaxed model is updated in place, so two nodes cannot share it. Keep one relaxed copy per thread, each with its own persistent solver. Memory grows with the number of copies: each holds the SDP data and a KKT factorization. Parallel evaluation in the core `solve` must be opt-in, because a custom node oracle may share state between calls.

## Once a node is cheap

These sit under the node solve at the current rates.

- Best-bound selection scans the whole open list (`select_node!` in `src/nodes.jl`).
- The branching rule runs twice per node: once in `evaluate!`, once in the main loop.
- Each node copies the integer box and the relaxed solution several times (`evaluate!`, the branching rule, `child_boxes`, `propagate_linear`).

# Keep stable

Downstream code with its own node oracle depends on these:

- `solve(relax, lb, ub; options, heuristic, polish, incumbent, x_incumbent)`, and `heuristic` running on every node taken from the queue.
- `NodeResult(:optimal, bound, value, x)`, `NodeResult(:infeasible)`, `NodeResult(:unknown)`.
- `jump_result(model, xvars)` with an empty `xvars` or plain numbers in it, and its status rules (MOSEK's `SLOW_PROGRESS` included).
- `Options(o; kwargs...)`, and the `Result` fields `status`, `objective`, `bound`, `x`, `nodes`, `unknown`, `root_bound`, `time` and `trace` (`t`, `nodes`, `bound`, `incumbent`, `open`).
