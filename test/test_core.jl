@testset "analytic oracle, both node selections" begin
    t = [0.3, 0.8, 0.45, 2.6, -1.3]
    for sel in (:best, :depth)
        r = solve(projection_oracle(t), [0, 0, 0, -5, -5], [1, 1, 1, 5, 5];
                  options=Options(QUIET; node_selection=sel))
        @test r.status === :optimal
        @test r.x == [0, 1, 0, 3, -1]
        @test r.objective ≈ sum(abs2, r.x .- t)
        check_invariants(r)
    end
end

@testset "infinite and fractional bounds" begin
    seen = Vector{Float64}[]
    oracle = projection_oracle([2.6, -1.3])
    r = solve((l, u) -> (push!(seen, copy(l), copy(u)); oracle(l, u)), [-Inf, -Inf], [Inf, Inf];
              options=QUIET)
    @test r.status === :optimal && r.x == [3, -1]

    # bounds are rounded inward before the root is solved
    empty!(seen)
    r = solve((l, u) -> (push!(seen, copy(l), copy(u)); oracle(l, u)), [0.2, -2.5], [2.7, 4.0];
              options=QUIET)
    @test seen[1] == [1, -2] && seen[2] == [2, 4]
    @test r.status === :optimal && r.x == [2, -1]
    # tolerance: a bound within int_tol of an integer is that integer
    empty!(seen)
    solve((l, u) -> (push!(seen, copy(l), copy(u)); oracle(l, u)), [1 - 1e-9, 0], [3 + 1e-9, 1];
          options=QUIET)
    @test seen[1] == [1, 0] && seen[2] == [3, 1]
end

@testset "empty box and zero variables" begin
    calls = Ref(0)
    r = solve(projection_oracle([0.5]; calls), [2.0], [1.5]; options=QUIET)
    @test r.status === :infeasible && r.nodes == 0 && calls[] == 0 && r.root_bound == Inf
    r = solve((l, u) -> NodeResult(:optimal, 3, 3, Float64[]), Float64[], Float64[]; options=QUIET)
    @test r.status === :optimal && r.objective == 3 && r.x == Float64[]
end

@testset "a failing oracle never prunes on failure" begin
    t = [0.3, 0.8, 0.45]
    calls = Ref(0)
    f = projection_oracle(t)
    relax(lb, ub) = (calls[] += 1; calls[] == 1 ? NodeResult(:unknown) : f(lb, ub))
    r = solve(relax, zeros(3), ones(3); options=QUIET)
    @test r.status === :optimal && r.x == [0, 1, 0] && r.unknown == 1
    @test r.root_bound == -Inf

    # every relaxation fails: the tree is enumerated and nothing can be certified
    r = solve((l, u) -> NodeResult(:unknown), zeros(2), ones(2); options=QUIET)
    @test r.status === :incomplete && r.unknown == r.nodes == 7 && r.bound == -Inf
end

@testset "NaN bound or value from the oracle" begin
    f = projection_oracle([0.2, 0.9])
    nanbound(l, u) = (r = f(l, u); NodeResult(r.status, NaN, r.value, r.x))
    r = solve(nanbound, zeros(2), ones(2); options=QUIET)
    @test r.status === :optimal && r.x == [0, 1]

    # integral root with an undefined value: the subtree must not be lost
    t = [0.0, 1.0]
    g = projection_oracle(t)
    first = Ref(true)
    nanroot(l, u) = first[] ? (first[] = false; NodeResult(:optimal, 0.0, NaN, copy(t))) : g(l, u)
    r = solve(nanroot, zeros(2), ones(2); options=QUIET)
    @test r.status === :optimal && r.x == [0, 1] && r.objective == 0
end

@testset "polish semantics" begin
    # integral node: the heuristic is not used to score it
    r = solve(lp_oracle, [0.0, 0], [1.0, 1]; options=QUIET, heuristic=(x, l, u) -> nothing)
    @test r.status === :optimal && r.objective ≈ 1 && r.x == [1, 0]
    # polish undecided at the only integral point that matters => :incomplete, never wrong
    r = solve(lp_oracle, [0.0, 0], [1.0, 1]; options=QUIET, polish=x -> nothing)
    @test r.status === :incomplete && r.bound <= 1
    # polish rejects (1,0): branching must still find (0,1) and certify it
    r = solve(lp_oracle, [0.0, 0], [1.0, 1]; options=QUIET,
              polish=x -> x == [1, 0] ? Inf : x[1] + 2x[2])
    @test r.status === :optimal && r.objective ≈ 2 && r.x == [0, 1]
    # polish cannot decide (1,0): (0,1) is found, but the result is not certified
    r = solve(lp_oracle, [0.0, 0], [1.0, 1]; options=QUIET,
              polish=x -> x == [1, 0] ? nothing : x[1] + 2x[2])
    @test r.status === :incomplete && r.objective ≈ 2 && r.bound <= 1 + 1e-9
    r = solve(lp_oracle, [0.0, 0], [1.0, 1]; options=QUIET, polish=x -> x[1] + 2x[2])
    @test r.status === :optimal && r.objective ≈ 1
    # polish scores (1,0) above the node bound: the rest of the subtree is still searched
    r = solve(lp_oracle, [0.0, 0], [1.0, 1]; options=QUIET,
              polish=x -> x == [1, 0] ? 3.0 : x[1] + 2x[2])
    @test r.status === :optimal && r.objective ≈ 2 && r.x == [0, 1]
end

@testset "heuristic and warm start" begin
    t = [0.3, 0.8, 0.45]
    f = projection_oracle(t)
    obj(x) = sum(abs2, x .- t)
    r = solve(f, zeros(3), ones(3); options=QUIET, heuristic=(x, l, u) -> (z = round.(x); (obj(z), z)))
    @test r.status === :optimal && r.x == [0, 1, 0]
    @test_throws DimensionMismatch solve(f, zeros(3), ones(3); options=QUIET,
                                         heuristic=(x, l, u) -> (0.0, [1.0]))
    # warm start at the optimum: kept, and certified
    r = solve(f, zeros(3), ones(3); options=QUIET, incumbent=obj([0, 1, 0]), x_incumbent=[0, 1, 0])
    @test r.status === :optimal && r.x == [0, 1, 0]
    # warm start value without a point
    r = solve(f, zeros(3), ones(3); options=QUIET, incumbent=obj([0, 1, 0]))
    @test r.status === :optimal && r.objective ≈ obj([0, 1, 0])
    @test_throws DimensionMismatch solve(f, zeros(3), ones(3); incumbent=1.0, x_incumbent=[0, 1])
    @test_throws ArgumentError solve(f, zeros(3), ones(3); incumbent=NaN)
end

@testset "limits" begin
    f = projection_oracle(fill(0.5, 8) .+ 0.01 .* (1:8))
    r = solve(f, zeros(8), ones(8); options=Options(QUIET; node_limit=5))
    @test r.status === :node_limit && 5 <= r.nodes <= 6
    @test r.bound <= r.objective
    r = solve(f, zeros(8), ones(8); options=Options(QUIET; time_limit=0.0))
    @test r.status === :time_limit && r.nodes == 1
    @test isfinite(r.bound) && r.bound <= minimum(sum(abs2, z .- (fill(0.5, 8) .+ 0.01 .* (1:8)))
                                                  for z in binvecs(8))
end

@testset "final bound honours the gap tolerance" begin
    # root bound 1.0 (fractional); heuristic finds 1.5; rel_gap 0.5 closes the root
    root = (l, u) -> NodeResult(:optimal, 1.0, 1.0, [0.5])
    r = solve(root, [0.0], [1.0]; options=Options(QUIET; rel_gap=0.5),
              heuristic=(x, l, u) -> (1.5, [1.0]))
    @test r.status === :optimal && r.objective == 1.5 && r.bound == 1.0
    @test gap(r) ≈ 1 / 3
end

@testset "invalid oracle output and branching rules" begin
    @test_throws ArgumentError solve((l, u) -> 1.0, [0], [1])
    @test_throws ArgumentError solve((l, u) -> NodeResult(:feasible), [0], [1])
    @test_throws ArgumentError solve((l, u) -> NodeResult(:optimal, 0, 0, nothing), [0], [1])
    @test_throws DimensionMismatch solve((l, u) -> NodeResult(:optimal, 0, 0, [0.5, 0.5]), [0], [1])
    @test_throws DimensionMismatch solve(lp_oracle, [0, 0], [1])
    @test_throws ArgumentError solve(lp_oracle, [0, NaN], [1, 1])

    f = projection_oracle([0.3, 0.6])
    # a rule that picks a fixed variable is an error, not an infinite loop
    @test_throws ArgumentError solve(f, [0, 0], [1, 1]; options=QUIET,
                                     branch=(x, l, u, tol) -> x === nothing ? nothing : 1)
    # a rule that cannot handle x === nothing falls back to first_free
    first = Ref(true)
    flaky(l, u) = first[] ? (first[] = false; NodeResult(:unknown)) : f(l, u)
    r = solve(flaky, [0, 0], [1, 1]; options=QUIET,
              branch=(x, l, u, tol) -> x === nothing ? nothing : most_fractional(x, l, u, tol))
    @test r.status === :optimal && r.x == [0, 1]
end
