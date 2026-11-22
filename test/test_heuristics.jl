# rounding_heuristic / fixed_oracle

@testset "nearest rounding re-solves with integers fixed" begin
    calls = Ref(0)
    oracle = projection_oracle([0.3, 0.8, 0.45, 2.6]; calls)
    h = rounding_heuristic(fixed_oracle(oracle); rule=:nearest)
    v, xint = h([0.3, 0.8, 0.45, 2.6], zeros(4), fill(5.0, 4))
    @test xint == [0, 1, 0, 3]
    @test v ≈ 0.3^2 + 0.2^2 + 0.45^2 + 0.4^2
    # node box is respected, and a point seen before is not re-solved
    n = calls[]
    v2, x2 = h([0.3, 0.8, 0.45, 2.6], zeros(4), fill(5.0, 4))
    @test x2 == xint && calls[] == n
    _, x3 = h([0.3, 0.8, 0.45, 2.6], [1, 0, 0, 0], [1, 1, 1, 2])
    @test x3 == [1, 1, 0, 2]
end

@testset "randomized rounding: unbiased, reproducible, best of rounds" begin
    seen = Vector{Float64}[]
    ev = x -> (push!(seen, copy(x)); sum(x))           # every point feasible, value = Σx
    h = rounding_heuristic(ev; rule=:random, rounds=2000, rng=Random.Xoshiro(3), cache=false)
    v, xint = h([0.25, 0.5, 0.75], zeros(3), ones(3))
    @test length(seen) == 2000
    @test isapprox(sum(seen) ./ length(seen), [0.25, 0.5, 0.75]; atol=0.05)
    @test v == 0 && xint == zeros(3)                    # the best draw is returned
    h1 = rounding_heuristic(ev; rule=:random, rng=Random.Xoshiro(7))
    h2 = rounding_heuristic(ev; rule=:random, rng=Random.Xoshiro(7))
    @test all(h1([0.5, 0.5, 0.5, 0.5], zeros(4), ones(4)) == h2([0.5, 0.5, 0.5, 0.5], zeros(4), ones(4)) for _ in 1:5)
    @test h1([0.0, 1.0], zeros(2), ones(2))[2] == [0, 1]   # integral values are kept
end

@testset "infeasible / undecided points, frequency, argument checks" begin
    @test rounding_heuristic(x -> Inf; rule=:nearest)([0.4], [0.0], [1.0]) === nothing
    @test rounding_heuristic(x -> nothing; rule=:nearest)([0.4], [0.0], [1.0]) === nothing
    @test rounding_heuristic(x -> NaN; rule=:nearest)([0.4], [0.0], [1.0]) === nothing
    @test fixed_oracle(lp_oracle)([0.0, 0.0]) == Inf
    @test fixed_oracle((l, u) -> NodeResult(:unknown))([0.0]) === nothing
    n = Ref(0)
    h = rounding_heuristic(x -> (n[] += 1; 1.0); rule=:nearest, every=3, cache=false)
    foreach(_ -> h([0.4], [0.0], [1.0]), 1:7)
    @test n[] == 3                                      # calls 1, 4, 7
    @test_throws ArgumentError rounding_heuristic(identity; rule=:up)
    @test_throws ArgumentError rounding_heuristic(identity; rounds=0)
end

@testset "inside solve: finds the optimum earlier on a knapsack-like MILP" begin
    # min -Σ cᵢ xᵢ s.t. Σ wᵢ xᵢ ≤ W, x ∈ {0,1}: the LP relaxation (greedy by ratio) has
    # one fractional item; nearest rounding of the relaxation gives a feasible point
    c, w, W = [10.0, 9, 8, 7, 6, 5, 4, 3], [5.0, 5, 5, 5, 5, 5, 5, 5], 22.0
    function knap(lb, ub)
        x, cap = copy(lb), W - sum(w .* lb)
        cap < -1e-12 && return NodeResult(:infeasible)
        for i in sortperm(c ./ w; rev=true)
            lb[i] == ub[i] && continue
            x[i] = clamp(cap / w[i], 0, 1); cap -= w[i] * x[i]
        end
        f = -sum(c .* x)
        return NodeResult(:optimal, f, f, x)
    end
    ev(x) = sum(w .* x) <= W + 1e-9 ? -sum(c .* x) : Inf
    r0 = solve(knap, zeros(8), ones(8); options=QUIET)
    r1 = solve(knap, zeros(8), ones(8); options=QUIET, heuristic=rounding_heuristic(ev; rule=:nearest))
    @test r0.status === r1.status === :optimal
    @test r1.objective ≈ r0.objective ≈ -(10 + 9 + 8 + 7)
    @test r1.nodes <= r0.nodes
    check_invariants(r1)
end
