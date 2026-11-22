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

@testset "linear presolve" begin
    # an affine bound is not a variable bound, so only presolve puts it in the box
    model = Model()
    @variable(model, x, Int)
    @constraint(model, 2x >= -1)
    @constraint(model, 2x <= 3)
    @objective(model, Min, x)
    relaxed, refmap = copy_model(model)
    rows = Sequoia.linear_rows(relaxed, [refmap[x]])
    lb, ub = [-Inf], [Inf]
    Sequoia.extract_singletons!(relaxed, rows, [refmap[x]], lb, ub)
    @test lb[1] ≈ -0.5 && ub[1] ≈ 1.5
    @test isempty(rows.rows)
    @test num_constraints(relaxed, AffExpr, MOI.GreaterThan{Float64}) == 0
    @test num_constraints(model, AffExpr, MOI.GreaterThan{Float64}) == 1   # model untouched

    seen = Ref{Any}(nothing)
    r = solve(model; optimizer=OPT, options=QUIET,
              propagate=(l, u) -> (seen[] = (copy(l), copy(u)); (l, u)))
    @test r.status === :optimal && r.x == [0]
    @test seen[][1] == [0.0] && seen[][2] == [1.0]   # rounded inward

    # a continuous singleton tightens the static bound and stays in the model
    model = Model()
    @variable(model, y)
    @constraint(model, 2y <= 4)
    rows = Sequoia.linear_rows(model, JuMP.VariableRef[])
    Sequoia.extract_singletons!(model, rows, JuMP.VariableRef[], Float64[], Float64[])
    @test rows.ub[1] == 2
    @test length(rows.rows) == 1
    @test num_constraints(model, AffExpr, MOI.LessThan{Float64}) == 1

    # at most one: fixing the first variable forces the other two off
    model = Model()
    @variable(model, z[1:3], Bin)
    @constraint(model, sum(z) <= 1)
    rows = Sequoia.linear_rows(model, z)
    fixed = Sequoia.propagate_linear(rows, [1.0, 0, 0], [1.0, 1, 1]; int_tol=1e-6)
    @test fixed isa Tuple && fixed[1] == [1.0, 0, 0] && fixed[2] == [1.0, 0, 0]

    # the same row, inside the search: same optimum, and the infeasible child is not solved
    model = Model()
    @variable(model, z[1:3], Bin)
    @constraint(model, sum(z) <= 1)
    @constraint(model, 2 * sum(z) >= 1.5)
    @objective(model, Min, 0.5 * sum(z))
    r = solve(model; optimizer=OPT, options=QUIET)
    r2 = solve(model; optimizer=OPT, options=QUIET, presolve=false)
    @test r.status === :optimal && r2.status === :optimal
    @test r.objective ≈ 0.5 atol = 1e-6
    @test r2.objective ≈ 0.5 atol = 1e-6
    @test sum(r.x) == 1
    @test r.nodes < r2.nodes

    # propagation proves the root empty, so no relaxation is solved
    model = Model()
    @variable(model, z[1:2], Bin)
    @constraint(model, sum(z) >= 3)
    @objective(model, Min, sum(z))
    r = solve(model; optimizer=OPT, options=QUIET)
    @test r.status === :infeasible && r.nodes == 0 && r.x === nothing

    # an integer objective is detected; a fractional coefficient or a continuous term is not
    model = Model()
    @variable(model, z[1:2], Bin)
    @variable(model, y)
    @objective(model, Min, z[1] + 2z[2])
    @test Sequoia.integral_objective(model, z)
    @objective(model, Min, 0.5 * z[1] + 0.5 * z[2])
    @test !Sequoia.integral_objective(model, z)
    @objective(model, Min, y + z[1])
    @test !Sequoia.integral_objective(model, z)
    @objective(model, Min, 4)
    @test Sequoia.integral_objective(model, z)

    model = Model()
    @variable(model, z[1:4], Bin)
    @constraint(model, 2 * sum(z) >= 3)
    @objective(model, Min, sum(z))
    r = solve(model; optimizer=OPT, options=QUIET)
    r2 = solve(model; optimizer=OPT, options=QUIET, presolve=false)
    @test r.status === :optimal && r.objective ≈ 2
    @test r2.objective ≈ 2 && r.nodes <= r2.nodes

    model = Model()
    @variable(model, z[1:2], Bin)
    @constraint(model, sum(z) >= 1)
    @objective(model, Min, 0.5 * sum(z))
    r = solve(model; optimizer=OPT, options=QUIET)
    @test r.status === :optimal && sum(r.x) == 1
    @test r.objective ≈ 0.5 atol = 1e-6
end

@testset "reused Clarabel solver" begin
    ext = Base.get_extension(Sequoia, :SequoiaClarabelExt)
    @test ext !== nothing

    function conic_model(reuse; huge_row=false)
        cm = Model()
        Sequoia.set_node_optimizer!(cm, OPT, reuse)
        @variable(cm, cx[1:3])
        @constraint(cm, Symmetric([cx[1] 0.3; 0.3 cx[2]]) in PSDCone())
        @constraint(cm, [cx[3]; cx[1] - 0.2; cx[2] - 0.7] in SecondOrderCone())
        huge_row && @constraint(cm, cx[1] + cx[2] <= 1e30)   # Clarabel's presolve drops it
        @objective(cm, Min, cx[1] + 2cx[2] + cx[3])
        return cm, cx
    end
    boxes = [([0.0, 0, 0], [1.0, 1, 1]), ([0.5, 0, 0], [1.0, 0.4, 1]),
             ([0.0, 0.6, 0], [0.2, 0.6, 1]), ([1.0, 1, 0], [1.0, 1, 0.1])]
    clarabel_solver(m) = unsafe_backend(m).inner.solver

    for huge_row in (false, true)
        m, x = conic_model(true; huge_row)
        solvers = []
        for (l, u) in boxes
            foreach(Sequoia.set_box!, x, l, u)
            optimize!(m)
            push!(solvers, clarabel_solver(m))
            fresh, fx = conic_model(false; huge_row)
            foreach(Sequoia.set_box!, fx, l, u)
            optimize!(fresh)
            r, rf = jump_result(m, x), jump_result(fresh, fx)
            @test r.status === rf.status
            if r.status === :optimal
                @test r.value ≈ rf.value atol = 1e-6
                @test r.x ≈ rf.x atol = 1e-5
            end
        end
        # new bound values update the solver in place, unless presolve has dropped rows
        @test allunique(solvers) == huge_row
        @test allequal(solvers) == !huge_row

        # dropping a bound deletes its row, so the solver is rebuilt
        Sequoia.set_box!(x[3], 0.0, Inf)
        optimize!(m)
        @test clarabel_solver(m) !== last(solvers)
        @test jump_result(m, x).status === :optimal
    end

    # branching from infinite bounds adds rows; results match a solver rebuilt at every node
    model = Model()
    @variable(model, x[1:2], Int)
    @variable(model, t)
    @constraint(model, [t; x[1] - 2.6; x[2] + 1.3] in SecondOrderCone())
    @objective(model, Min, t)
    r = solve(model; optimizer=OPT, options=QUIET)
    r2 = solve(model; optimizer=OPT, options=QUIET, reuse_solver=false)
    @test r.status === r2.status === :optimal
    @test r.x == r2.x == [3, -1]
    @test r.objective ≈ r2.objective atol = 1e-7

    Random.seed!(2)
    n, d = 7, 3
    A = [(G = randn(d, 2); G * G') for _ in 1:n]
    B = 0.6 * sum(A) / n * 3
    c = rand(n) .+ 0.5
    model = Model()
    @variable(model, z[1:n], Bin)
    @constraint(model, Symmetric(sum(z[i] * A[i] for i in 1:n) - B) in PSDCone())
    @objective(model, Min, dot(c, z))
    r = solve(model; optimizer=OPT, options=QUIET)
    r2 = solve(model; optimizer=OPT, options=QUIET, reuse_solver=false)
    @test r.status === r2.status === :optimal
    @test r.x == r2.x
    @test r.objective ≈ r2.objective rtol = 1e-7
    check_invariants(r)
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
