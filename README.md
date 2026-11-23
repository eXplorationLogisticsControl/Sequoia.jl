# Sequoia.jl

[![Dev docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://eXplorationLogisticsControl.github.io/Sequoia.jl/dev/)
[![CI](https://github.com/eXplorationLogisticsControl/Sequoia.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eXplorationLogisticsControl/Sequoia.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A small, generic branch-and-bound solver for mixed-integer convex problems (MISDP, MISOCP,
MILP, ...). Each node solves a continuous relaxation with the convex solver of your
choice.

- **Bring your own node oracle**: `relax(lb, ub) -> NodeResult` may rebuild a smaller
  model per node to exploit problem structure.
- **Robust to solver failures**: a failed relaxation is never used to prune, and the final
  status tells you when a result cannot be certified.
- **JuMP convenience interface** for models with `Bin`/`Int` variables.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/eXplorationLogisticsControl/Sequoia.jl")
```

## Quick start

```julia
using JuMP, Clarabel, Sequoia

model = Model()
@variable(model, z[1:5], Bin)
@constraint(model, [12, 7, 11, 8, 9]' * z <= 26)
@objective(model, Max, [24, 13, 23, 15, 16]' * z)

result = Sequoia.solve(model; optimizer=Clarabel.Optimizer)
result.status, result.objective, result.x
```

With your own node oracle (minimization):

```julia
t = [0.3, 0.8, 2.6]
relax(lb, ub) = (x = clamp.(t, lb, ub); f = sum(abs2, x .- t); NodeResult(:optimal, f, f, x))
result = solve(relax, [0, 0, -5], [1, 1, 5])
```

See the [documentation](https://eXplorationLogisticsControl.github.io/Sequoia.jl/dev/) for
the algorithm, statuses and guarantees, and the full API.

## Development

```sh
julia --project -e 'using Pkg; Pkg.test()'                                  # tests
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()' # docs env, once
julia --project=docs docs/make.jl                                           # build docs
```

The Documentation workflow publishes the docs to GitHub Pages on every push to `main`. It
needs Settings → Pages → Source set to "GitHub Actions".
