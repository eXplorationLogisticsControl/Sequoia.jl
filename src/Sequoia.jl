"""
    Sequoia

Small, generic branch and bound for mixed-integer convex problems (MISDP, MISOCP,
MILP, ...). Each node solves a continuous relaxation with a convex solver of the user's
choice.

There are two entry points:

- [`solve(relax, lb, ub; ...)`](@ref solve(::Any, ::AbstractVector, ::AbstractVector)):
  the core. You supply a *node oracle* `relax(lb, ub) -> NodeResult` that solves the
  relaxation with the integer variables restricted to the box `[lb, ub]`. You may
  rebuild or simplify the model at every node (e.g. eliminate blocks whose indicator is
  fixed to 0).
- [`solve(model::JuMP.Model; optimizer, ...)`](@ref solve(::JuMP.Model)): convenience
  wrapper for a JuMP model with `Bin`/`Int` variables. It relaxes a copy of the model and
  changes variable bounds at every node.
- [`read_dat_s`](@ref) reads a sparse SDPA file into an [`SDPAProblem`](@ref), and
  [`jump_model`](@ref) turns that into a JuMP model [`solve`](@ref) can branch on.

The core works in the minimization sense, and the oracle's `bound` must be a valid lower
bound on the node relaxation (a dual objective value is the safe choice).
"""
module Sequoia

import JuMP
using JuMP: MOI
using Printf: @printf, @sprintf

export NodeResult, Options, Result, solve, gap, most_fractional, first_free, jump_result,
       read_dat_s, SDPAProblem, jump_model

include("options.jl")     # Options: tolerances, limits, node selection
include("results.jl")     # NodeResult (oracle output) and Result (solver output)
include("branching.jl")   # branching rules and split points
include("nodes.jl")       # Node, node selection, child creation
include("tree.jl")        # search-tree state: incumbent, pruning, node evaluation, logging
include("solve.jl")       # core user-level API: solve(relax, lb, ub)
include("presolve.jl")    # linear presolve: singleton bounds, propagation
include("jump.jl")        # JuMP interface: jump_result, solve(model)
include("io.jl")          # sparse SDPA reader

end # module
