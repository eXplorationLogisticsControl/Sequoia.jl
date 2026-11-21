@testset "knapsack (MILP), max sense, warm start, depth first" begin
    Random.seed!(1)
    n = 10
    w, v, cap = rand(1:20, n), rand(1:30, n), 40
    best = maximum(dot(v, z) for z in binvecs(n) if dot(w, z) <= cap)

    model = Model()
    @variable(model, z[1:n], Bin)
    @constraint(model, dot(w, z) <= cap)
    @objective(model, Max, dot(v, z))
    for opts in (QUIET, Options(QUIET; node_selection=:depth))
        r = solve(model; optimizer=OPT, options=opts)
        @test r.status === :optimal
        @test r.objective ≈ best atol = 1e-6
        @test dot(v, r.x) ≈ best atol = 1e-6
        @test dot(w, r.x) <= cap
        @test r.bound ≈ best atol = 1e-6
        @test r.root_bound >= best - 1e-6     # max sense: root bound is an upper bound
        check_invariants(r; minimize=false)
    end
    @test is_binary(z[1]) && !has_lower_bound(z[1])   # user model untouched
    r = solve(model; optimizer=OPT, options=QUIET, incumbent=best)
    @test r.status === :optimal && r.objective ≈ best
end

@testset "general integers (MISOCP)" begin
    model = Model()
    @variable(model, -5 <= x[1:2] <= 5, Int)
    @variable(model, t)
    @constraint(model, [t; x[1] - 2.6; x[2] + 1.3] in SecondOrderCone())
    @objective(model, Min, t)
    r = solve(model; optimizer=OPT, options=QUIET)
    @test r.status === :optimal
    @test r.x ≈ [3, -1]
    @test r.objective ≈ hypot(0.4, 0.3) atol = 1e-6
    @test lower_bound(x[1]) == -5 && is_integer(x[1])
end

@testset "small MISDP: min-cost cover Σ zᵢAᵢ ⪰ B" begin
    Random.seed!(2)
    n, d = 7, 3
    A = [(G = randn(d, 2); G * G') for _ in 1:n]
    B = 0.6 * sum(A) / n * 3
    c = rand(n) .+ 0.5
    feas(z) = eigmin(Symmetric(sum(z[i] * A[i] for i in 1:n) - B)) >= -1e-9
    best = minimum((dot(c, z) for z in binvecs(n) if feas(z)); init=Inf)
    @test isfinite(best)

    model = Model()
    @variable(model, z[1:n], Bin)
    @constraint(model, Symmetric(sum(z[i] * A[i] for i in 1:n) - B) in PSDCone())
    @objective(model, Min, dot(c, z))
    r = solve(model; optimizer=OPT, options=QUIET)
    @test r.status === :optimal
    @test r.objective ≈ best rtol = 1e-6
    @test feas(r.x)
    check_invariants(r)
end

@testset "infeasible and integral objective" begin
    model = Model()
    @variable(model, z[1:3], Bin)
    @constraint(model, sum(z) >= 2.5)
    @constraint(model, z[1] + z[2] <= 1.2)
    @objective(model, Min, sum(z))
    r = solve(model; optimizer=OPT, options=Options(QUIET; integral_objective=true))
    @test r.status === :infeasible && r.x === nothing && r.objective == Inf

    model = Model()
    @variable(model, z[1:4], Bin)
    @constraint(model, 2z[1] + 2z[2] + 2z[3] + 2z[4] >= 3)   # LP bound 1.5, optimum 2
    @objective(model, Min, sum(z))
    r = solve(model; optimizer=OPT, options=Options(QUIET; integral_objective=true))
    @test r.status === :optimal && r.objective ≈ 2
    r2 = solve(model; optimizer=OPT, options=QUIET)
    @test r2.status === :optimal && r2.objective ≈ 2 && r.nodes <= r2.nodes
end

@testset "fixed variables, fractional bounds, max-sense callbacks" begin
    model = Model()
    @variable(model, z[1:3], Bin)
    fix(z[1], 1)
    @constraint(model, z[2] + z[3] >= 1)
    @objective(model, Min, sum(z))
    r = solve(model; optimizer=OPT, options=QUIET)
    @test r.status === :optimal && r.objective ≈ 2 && r.x[1] == 1

    # an integer variable with fractional bounds is rounded inward
    model = Model()
    @variable(model, 0.5 <= y <= 2.5, Int)
    @objective(model, Min, (y - 0.1) * 1.0)
    r = solve(model; optimizer=OPT, options=QUIET)
    @test r.status === :optimal && r.x == [1]

    # max sense: heuristic values are given in the model's own sense
    model = Model()
    @variable(model, z[1:3], Bin)
    @constraint(model, sum(z) <= 2.5)
    @objective(model, Max, sum(z))
    r = solve(model; optimizer=OPT, options=QUIET, heuristic=(x, l, u) -> (1.0, [1.0, 0, 0]))
    @test r.status === :optimal && r.objective ≈ 2 && gap(r) < 1e-6
    r = solve(model; optimizer=OPT, options=QUIET, polish=x -> sum(x))
    @test r.status === :optimal && r.objective ≈ 2

    # root_bound is the root relaxation's value, even if the root is integral
    model = Model()
    @variable(model, z[1:2], Bin)
    @constraint(model, z[1] + z[2] >= 1)
    @objective(model, Min, z[1] + 2z[2])
    r = solve(model; optimizer=OPT, options=QUIET, incumbent=5.0)
    @test r.root_bound ≈ 1 atol = 1e-6
end

@testset "mixed-integer with continuous variables" begin
    model = Model()
    @variable(model, z[1:3], Bin)
    @variable(model, 0 <= y[1:3] <= 10)
    @constraint(model, [i in 1:3], y[i] <= 10z[i])
    @constraint(model, sum(y) >= 12)
    @objective(model, Min, sum(5z) + sum(y))
    r = solve(model; optimizer=OPT, options=QUIET)
    @test r.status === :optimal && sum(r.x) == 2
    @test r.objective ≈ 22 atol = 1e-5
end

@testset "unsupported constraints" begin
    model = Model()
    @variable(model, x[1:2] >= 0)
    @constraint(model, x in SOS1())
    @test_throws ArgumentError solve(model; optimizer=OPT)
    model = Model()
    @variable(model, z, Bin)
    @variable(model, x)
    @constraint(model, z --> {x >= 1})
    @test_throws ArgumentError solve(model; optimizer=OPT)
end

@testset "jump_result" begin
    model = Model(OPT)
    @variable(model, x >= 1)
    @objective(model, Min, 2x)
    optimize!(model)
    r = jump_result(model, [x, 3])
    @test r.status === :optimal
    @test r.bound ≈ 2 atol = 1e-6
    @test r.value ≈ 2 atol = 1e-6
    @test r.x ≈ [1, 3] atol = 1e-6
    @test jump_result(model, []).x == Float64[]
    @objective(model, Max, -2x)
    optimize!(model)
    r = jump_result(model, [x])
    @test r.value ≈ 2 atol = 1e-6     # minimization sense
    @test r.bound ≈ 2 atol = 1e-6

    @constraint(model, x <= 0)
    optimize!(model)
    @test jump_result(model, [x]).status === :infeasible

    # statuses a convex solver rarely produces on demand, via MOI's mock optimizer
    function mocked(ts, ps, ds; obj=2.0, dobj=1.9)
        mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}();
                                               eval_objective_value=false, eval_dual_objective_value=false)
        m = direct_model(mock)
        @variable(m, y)
        @objective(m, Min, y)
        optimize!(m)
        MOI.set(mock, MOI.TerminationStatus(), ts)
        MOI.set(mock, MOI.ResultCount(), 1)
        MOI.set(mock, MOI.PrimalStatus(), ps)
        MOI.set(mock, MOI.DualStatus(), ds)
        MOI.set(mock, MOI.ObjectiveValue(), obj)
        MOI.set(mock, MOI.DualObjectiveValue(), dobj)
        MOI.set(mock, MOI.VariablePrimal(), index(y), obj)
        return jump_result(m, [y])
    end
    r = mocked(MOI.SLOW_PROGRESS, MOI.FEASIBLE_POINT, MOI.FEASIBLE_POINT)
    @test r.status === :optimal && r.bound == 1.9 && r.value == 2.0
    @test mocked(MOI.SLOW_PROGRESS, MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT).status === :unknown
    @test mocked(MOI.SLOW_PROGRESS, MOI.NO_SOLUTION, MOI.INFEASIBILITY_CERTIFICATE).status === :infeasible
    r = mocked(MOI.ALMOST_OPTIMAL, MOI.NEARLY_FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    @test r.status === :optimal && r.bound == -Inf
    @test mocked(MOI.DUAL_INFEASIBLE, MOI.NO_SOLUTION, MOI.NO_SOLUTION).status === :unknown
    @test mocked(MOI.NUMERICAL_ERROR, MOI.NO_SOLUTION, MOI.NO_SOLUTION).status === :unknown
    @test mocked(MOI.LOCALLY_INFEASIBLE, MOI.NO_SOLUTION, MOI.NO_SOLUTION).status === :infeasible
end
