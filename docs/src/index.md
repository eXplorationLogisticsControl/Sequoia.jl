```@meta
CurrentModule = Sequoia
```

# Sequoia.jl

Sequoia is a small, generic **branch-and-bound** solver for mixed-integer *convex*
problems: mixed-integer semidefinite (MISDP), second-order cone (MISOCP), and linear
(MILP) programs. Each node of the search tree solves a continuous relaxation with a convex
solver of your choice, such as [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl),
MOSEK, or SCS.

It is built for problems where generic MI-conic solvers struggle, and where you know
something about the structure:

- **Bring your own node oracle.** At each node, Sequoia only asks for "the relaxation over
  this box of integer bounds". You may rebuild a smaller model per node, for instance
  dropping the constraint blocks of indicators fixed to zero.
- **Robust to solver failures.** Conic solvers sometimes stall. A failed relaxation is
  never used to prune: the node keeps its parent's bound and is branched on, and the final
  status says when a result cannot be certified.
- **Exact scoring of integer points.** An optional `polish` callback re-evaluates
  candidate integer points exactly, so a relaxation that is integral only within a
  tolerance cannot yield a wrong incumbent.

## Installation

Sequoia is not registered. Install it from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/eXplorationLogisticsControl/Sequoia.jl")
```

## Quick start

A 0/1 knapsack problem, written in JuMP and solved with Clarabel at the nodes:

```@example quickstart
using JuMP, Clarabel, Sequoia

w = [12, 7, 11, 8, 9]   # weights
v = [24, 13, 23, 15, 16]  # values

model = Model()
@variable(model, z[1:5], Bin)
@constraint(model, w' * z <= 26)
@objective(model, Max, v' * z)

optimizer = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
result = Sequoia.solve(model; optimizer, options=Options(log_every=Inf))
```

```@example quickstart
result.x   # values of the integer variables, in the order of all_variables(model)
```

## Where to go next

- [Solving JuMP models](@ref): the convenience interface used above.
- [Custom node oracles](@ref): the core interface, for structure-exploiting relaxations.
- [Algorithm](@ref): how nodes are selected, pruned and branched, and what each status
  guarantees.
- [API reference](@ref).
