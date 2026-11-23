# API reference

```@meta
CurrentModule = Sequoia
```

```@docs
Sequoia
```

## Solving

```@docs
solve(::Any, ::AbstractVector, ::AbstractVector)
solve(::JuMP.Model)
Options
Result
gap
```

## Input

```@docs
read_dat_s
SDPAProblem
jump_model
Sequoia.PSDBlock
Sequoia.DiagonalBlock
```

## Node oracles

```@docs
NodeResult
jump_result
```

## Primal heuristics

```@docs
rounding_heuristic
fixed_oracle
```

## Branching rules

```@docs
most_fractional
first_free
```

## Internals

These are not part of the public API and may change without notice.

```@docs
Sequoia.Node
Sequoia.Tree
Sequoia.split_point
Sequoia.select_node!
Sequoia.child_boxes
Sequoia.evaluate!
Sequoia.discard!
Sequoia.improve!
Sequoia.validate
```
