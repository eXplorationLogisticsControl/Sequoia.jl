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

## Node oracles

```@docs
NodeResult
jump_result
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
