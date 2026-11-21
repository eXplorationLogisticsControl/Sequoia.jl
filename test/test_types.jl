@testset "Options defaults, copy and validation" begin
    o = Options()
    @test o.node_selection === :best && o.rel_gap == 1e-6 && isinf(o.time_limit)
    o2 = Options(o; time_limit=60, node_selection=:depth)
    @test o2.time_limit == 60 && o2.node_selection === :depth && o2.rel_gap == o.rel_gap
    @test o.time_limit == Inf   # original untouched

    f = (lb, ub) -> NodeResult(:infeasible)
    for bad in (Options(node_selection=:breadth), Options(int_tol=0.5), Options(rel_gap=-1),
                Options(time_limit=-1), Options(node_limit=-1), Options(log_every=0))
        @test_throws ArgumentError solve(f, [0], [1]; options=bad)
    end
end

@testset "NodeResult" begin
    r = NodeResult(:unknown)
    @test r.status === :unknown && isnan(r.bound) && r.x === nothing
    r = NodeResult(:optimal, 1, 2, [1, 0])
    @test r.bound === 1.0 && r.x == [1.0, 0.0]
end

@testset "gap and show" begin
    r = Result(:optimal, 10.0, 9.0, [1.0], 3, 0, 8.0, 0.1, Sequoia.Snapshot[])
    @test gap(r) ≈ 0.1
    @test gap(Result(:optimal, -10.0, -11.0, nothing, 0, 0, 0.0, 0.0, Sequoia.Snapshot[])) ≈ 0.1
    @test gap(Result(:infeasible, Inf, Inf, nothing, 0, 0, Inf, 0.0, Sequoia.Snapshot[])) == Inf
    s = sprint(show, MIME"text/plain"(), r)
    @test occursin("optimal", s) && occursin("gap", s)
end
