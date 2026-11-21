# Benchmarks

Solve the sparse SDPA groups in `ref/instances` and write one performance table per group.

```sh
julia --project=benchmark benchmark/run.jl [TrussTopology ...]
```

With no group names, every directory under `ref/instances` is solved. Each instance stops at 60 seconds. Call `benchmark(groups; time_limit=Inf)` from this project to run without a cap. Clarabel's chordal decomposition is off: it crashes on these sparse PSD blocks, and the blocks are small enough for the dense factorization.

- [TrussTopology](TrussTopology/README.md)
