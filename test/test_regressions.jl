# Regression tests for defects found in the adversarial review. Each testset names the
# failure it guards against.

@testset "integral node with a weak bound is not closed on its primal value" begin
    # one binary, f(0) = 3, f(1) = 5. The root reports a valid bound 0 but an unconverged
    # integral primal point x = 1 with value 5; the optimum 3 must still be found.
    f(x) = x == 0 ? 3.0 : 5.0
    relax(lb, ub) = lb[1] == ub[1] ? NodeResult(:optimal, f(lb[1]), f(lb[1]), [lb[1]]) :
                    NodeResult(:optimal, 0.0, 5.0, [1.0])
    r = solve(relax, [0.0], [1.0]; options=QUIET)
    @test r.status === :optimal && r.objective == 3 && r.x == [0] && r.bound <= 3
    # same with no dual bound at all
    relax2(lb, ub) = lb[1] == ub[1] ? relax(lb, ub) : NodeResult(:optimal, -Inf, 5.0, [1.0])
    r = solve(relax2, [0.0], [1.0]; options=QUIET)
    @test r.status === :optimal && r.objective == 3
end

@testset "integral_objective tolerates a tiny relative bound overshoot" begin
    c = [100000.0, 100001.0]
    function lp(lb, ub)
        (lb[1] + lb[2] > 1 || ub[1] + ub[2] < 1) && return NodeResult(:infeasible)
        x1 = clamp(1 - lb[2], lb[1], ub[1])
        x = lb == [0, 0] && ub == [1, 1] ? [0.5, 0.5] : [x1, 1 - x1]
        v = x == [0.5, 0.5] ? c[1] : c[1] * x[1] + c[2] * x[2]
        return NodeResult(:optimal, v * (1 + 1e-10), v, x)
    end
    for io in (false, true)
        r = solve(lp, [0.0, 0], [1.0, 1]; options=Options(QUIET; integral_objective=io),
                  heuristic=(x, l, u) -> (c[2], [0.0, 1.0]))
        @test r.status === :optimal && r.objective == c[1] && r.x == [1, 0]
        @test r.bound <= c[1] + 1e-9 * c[1]
    end
end

@testset "non-finite relaxed x is a failure, never an incumbent" begin
    r = solve((l, u) -> NodeResult(:optimal, 1, 1, [NaN]), [0.0], [1.0]; options=QUIET)
    @test r.status === :incomplete && r.x === nothing && r.unknown == r.nodes
    r = solve((l, u) -> NodeResult(:optimal, 1, 1, [Inf]), [0.0], [Inf]; options=Options(QUIET; node_limit=20))
    @test r.x === nothing && r.status !== :optimal
end

@testset "incumbent = ±Inf" begin
    f = projection_oracle([0.3])
    @test_throws ArgumentError solve(f, [0], [1]; incumbent=-Inf)
    # Inf means "none" and x_incumbent is then ignored
    r = solve((l, u) -> NodeResult(:infeasible), [0.0], [1.0]; options=QUIET, x_incumbent=[1.0])
    @test r.status === :infeasible && r.x === nothing
    # max-sense JuMP model: an infinite incumbent of either sign means none
    model = Model()
    @variable(model, z[1:2], Bin)
    @objective(model, Max, z[1] + 2z[2])
    for inc in (Inf, -Inf)
        r = solve(model; optimizer=OPT, options=QUIET, incumbent=inc)
        @test r.status === :optimal
        @test r.objective ≈ 3 atol = 1e-6
    end
    r = solve(model; optimizer=OPT, options=QUIET, incumbent=3.0, x_incumbent=[1, 1])
    @test r.status === :optimal && r.x == [1, 1]
end

@testset "Interval bounds on integer variables (JuMP)" begin
    for kind in (:int, :bin)
        model = Model()
        x = kind === :int ? @variable(model, integer = true) : @variable(model, binary = true)
        @constraint(model, x in MOI.Interval(kind === :int ? -1.5 : 0.0, kind === :int ? 3.5 : 1.0))
        @objective(model, Min, 7.4 - x)   # wants x large: the interval binds
        r = solve(model; optimizer=OPT, options=QUIET)
        @test r.status === :optimal && r.x == [kind === :int ? 3 : 1]
        @test num_constraints(model, VariableRef, MOI.Interval{Float64}) == 1   # model untouched
    end
end

@testset "LOCALLY_SOLVED counts as solved" begin
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}();
                                       eval_objective_value=false, eval_dual_objective_value=false)
    m = direct_model(mock)
    @variable(m, y)
    @objective(m, Min, y)
    optimize!(m)
    for (ts, ps) in ((MOI.LOCALLY_SOLVED, MOI.FEASIBLE_POINT),
                     (MOI.ALMOST_LOCALLY_SOLVED, MOI.NEARLY_FEASIBLE_POINT))
        MOI.set(mock, MOI.TerminationStatus(), ts)
        MOI.set(mock, MOI.ResultCount(), 1)
        MOI.set(mock, MOI.PrimalStatus(), ps)
        MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
        MOI.set(mock, MOI.ObjectiveValue(), 2.0)
        MOI.set(mock, MOI.DualObjectiveValue(), 2.0)
        MOI.set(mock, MOI.VariablePrimal(), index(y), 2.0)
        @test jump_result(m, [y]).status === :optimal
    end
end

@testset "oracle reusing its x buffer cannot corrupt the tree" begin
    t = [0.3, 0.8, 0.45, 0.6]
    buf = zeros(4)
    function bufo(lb, ub)
        buf .= clamp.(t, lb, ub)
        f = sum(abs2, buf .- t)
        return NodeResult(:optimal, f, f, buf)
    end
    for sel in (:best, :depth)
        r = solve(bufo, zeros(4), ones(4); options=Options(QUIET; node_selection=sel))
        @test r.status === :optimal && r.x == round.(t) && r.objective ≈ sum(abs2, round.(t) .- t)
    end
end

@testset "callbacks mutating their arguments cannot corrupt the tree" begin
    t = [0.45, 0.4]
    proj = projection_oracle(t)
    for h in ((x, lb, ub) -> (lb[2] = 1.0; nothing),       # in-place "diving"
              (x, lb, ub) -> (x .= round.(x); nothing),    # in-place rounding
              (x, lb, ub) -> (ub .= lb; nothing))
        r = solve(proj, [0.0, 0.0], [1.0, 1.0]; options=QUIET, heuristic=h)
        @test r.status === :optimal && r.x == [0, 0] && r.objective ≈ sum(abs2, t)
    end
    # a polish that mutates its input must not change the stored incumbent
    r = solve(lp_oracle, [0.0, 0], [1.0, 1]; options=QUIET,
              polish=x -> (v = x[1] + 2x[2]; x .= 7; v))
    @test r.status === :optimal && r.x == [1, 0]
    # a mutating oracle
    r = solve((l, u) -> (res = proj(l, u); l .= 0; u .= 1; res), [0.0, 0.0], [1.0, 1.0]; options=QUIET)
    @test r.status === :optimal && r.x == [0, 0]
end

@testset "branching rule returning nothing with free variables" begin
    # a stateful rule that gives up after its first call: Sequoia falls back to first_free
    calls = Ref(0)
    rule(x, l, u, tol) = (calls[] += 1; calls[] <= 1 ? most_fractional(x, l, u, tol) : nothing)
    f = projection_oracle([0.3, 0.6, 0.45])
    r = solve(f, zeros(3), ones(3); options=QUIET, branch=rule)
    @test r.status === :optimal && r.x == [0, 1, 0]
end

@testset "empty box ignores the warm start" begin
    r = solve(projection_oracle([0.2, 0.2]), [0.0, 2.0], [1.0, 1.0]; options=QUIET,
              incumbent=0.0, x_incumbent=[0.0, 0.0])
    @test r.status === :infeasible && r.x === nothing && r.objective == Inf && r.nodes == 0
end

@testset "huge integer values are refused, not looped on" begin
    @test_throws ArgumentError Sequoia.split_point(nothing, 2.0^53, Inf)
    @test_throws ArgumentError solve((l, u) -> NodeResult(:unknown), [2.0^53], [Inf];
                                     options=Options(QUIET; node_limit=50))
end
