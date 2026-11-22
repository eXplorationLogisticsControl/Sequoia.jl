# Solving JuMP models

[`solve(model::JuMP.Model; optimizer, ...)`](@ref solve(::JuMP.Model)) takes a JuMP model
with `Bin` and `Int` variables and any constraints your continuous solver supports. Sequoia
copies the model, relaxes integrality on the copy, and at every node only changes the
bounds of the integer variables before re-solving. Your model is never modified.

```@example jump
using JuMP, Clarabel, Sequoia, LinearAlgebra
const OPT = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
nothing # hide
```

## A mixed-integer SDP

Choose a cheapest subset of PSD "resources" `Aᵢ` whose sum dominates a requirement `B` in
the semidefinite order:

```math
\min_{z \in \{0,1\}^n} \; c^\top z
\quad \text{s.t.} \quad \sum_i z_i A_i \succeq B .
```

```@example jump
using Random
Random.seed!(2)
n, d = 7, 3
A = [(G = randn(d, 2); G * G') for _ in 1:n]
B = 1.8 * sum(A) / n
c = rand(n) .+ 0.5

model = Model()
@variable(model, z[1:n], Bin)
@constraint(model, Symmetric(sum(z[i] * A[i] for i in 1:n) - B) in PSDCone())
@objective(model, Min, c' * z)

result = Sequoia.solve(model; optimizer=OPT, options=Options(log_every=Inf))
```

## General integers and continuous variables

Integer variables may have any bounds, including infinite ones. Fractional bounds are
rounded inward. Continuous variables are left to the relaxation:

```@example jump
model = Model()
@variable(model, -5 <= x[1:2] <= 5, Int)
@variable(model, t)
@constraint(model, [t; x[1] - 2.6; x[2] + 1.3] in SecondOrderCone())
@objective(model, Min, t)

result = Sequoia.solve(model; optimizer=OPT, options=Options(log_every=Inf))
result.x
```

Only the integer variables are returned in `result.x`. To recover the continuous ones, fix
the integers at `result.x` in your model and solve once more.

## Options and warm starts

All keyword arguments of the core solver are available; see [`Options`](@ref):

```julia
opts = Options(time_limit=600, rel_gap=1e-4, node_selection=:depth, log_every=10)
Sequoia.solve(model; optimizer=OPT, options=opts,
              incumbent=12.0, x_incumbent=z0)     # warm start (model's own sense)
```

If every feasible point has an integer objective value (e.g. a cardinality objective),
set `integral_objective=true`. Lower bounds are then rounded up before pruning, which can
shrink the tree a lot. With `presolve` left on, Sequoia sets this itself when the objective
is affine and integer-valued on the integer variables.

## Presolve

`presolve = true` (the default) reads the linear rows of the relaxed copy before the
search:

- a row on a single integer variable becomes part of the integer box, and the row is
  dropped from the copy;
- the remaining rows are propagated at every node, including the root. Fixing one variable
  of an at-most-one row, for example, forces the others, and a row that cannot be satisfied
  closes the node without a solve.

Pass `presolve = false` to branch on the model exactly as written. Your own `propagate`
callback still runs; the linear propagation runs first.

## Reusing the solver

Clarabel's MathOptInterface wrapper builds a new solver whenever the model changes. With
`reuse_solver = true` (the default) and Clarabel loaded, Sequoia keeps one solver for the
whole search and writes only the integer bounds that changed since the last node. Each
bound is a row of its own, so a new value changes one entry of the right-hand side, and the
solver skips the copy, equilibration and symbolic factorization. The interior-point method
still starts from scratch. On small SDPs, where building the solver is a noticeable part of
a node, this saved 13–23% of the time per node on the TrussTopology instances.

The solver is rebuilt when a bound becomes infinite or finite (a row appears or disappears),
or when Clarabel's presolve or chordal decomposition has reshaped the problem. Other
optimizers are attached as usual. Pass `reuse_solver = false` to rebuild at every node.

## Maximization

For `Max` models, Sequoia minimizes the negated objective internally. Everything you pass
in (`incumbent`, the values returned by `heuristic` and `polish`) and everything you get
back (`objective`, `bound`, `root_bound`, `trace`) is in the model's own sense.

## Limitations

- SOS1/SOS2, semicontinuous, semi-integer and indicator constraints raise an
  `ArgumentError`. Reformulate them with binaries (e.g. big-M or perspective
  constraints).
- The relaxations must be convex for the bounds to be valid.
- The node relaxation is re-solved from scratch at every node; no warm starting of the
  continuous solver is attempted. Only Clarabel's setup is reused (see above).
