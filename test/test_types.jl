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

@testset "progress log" begin
    @test Sequoia.PROGRESS_HEADER ==
          "    Nodes    |          Objective Bounds          |  Work   \n" *
          " Expl Unexpl |     Incumbent        BestBd    Gap |     Time\n"
    # Internal values are minimization; log_sign = -1 displays a Max model.
    missing = Sequoia.progress_row(3, 2, Inf, -52.5714, 0.1; log_sign=-1)
    done = Sequoia.progress_row(15, 0, -51.0, -51.0, 0.1; log_sign=-1)
    @test missing == "     3     2 |             -   5.25714e+01      - |     0.1s"
    @test done == "    15     0 |   5.10000e+01   5.10000e+01  0.00% |     0.1s"

    t = Sequoia.Tree(identity, nothing, nothing, nothing, Options(log_every=1e-6), 0, Inf, nothing, -1.0)
    t.nodes = 3
    push!(t.open, Sequoia.Node(Float64[], Float64[], -52.5714, nothing),
          Sequoia.Node(Float64[], Float64[], 0.0, nothing))
    pipe = Pipe()
    Base.link_pipe!(pipe; reader_supports_async=true, writer_supports_async=true)
    redirect_stdout(pipe) do
        Sequoia.log!(t, true)
        Sequoia.log!(t, true)
    end
    close(pipe.in)
    text = read(pipe, String)
    @test count("    Nodes    |", text) == 1
    @test count("     3     2 |             -   5.25714e+01      - |", text) == 2

    optimal = Sequoia.termination_summary(:optimal, 15, 0.12, -51.0, -51.0, 1e-6; log_sign=-1)
    @test optimal == """
    Explored 15 nodes in 0.12 seconds

    Optimal solution found (tolerance 1.00e-06)
    Best objective 5.100000000000e+01, best bound 5.100000000000e+01, gap 0.0000%"""
    infeasible = Sequoia.termination_summary(:infeasible, 1, 0.0, Inf, Inf, 1e-6)
    @test infeasible == "Explored 1 node in 0.00 seconds\n\nModel is infeasible"
    limited = Sequoia.termination_summary(:time_limit, 3, 5.0, Inf, -52.5714, 1e-4;
                                           log_sign=-1, unknown=2)
    @test limited == """
    Explored 3 nodes (2 failed) in 5.00 seconds

    Time limit reached
    Best objective -, best bound 5.257140000000e+01, gap -"""
    @test occursin("Node limit reached",
                   Sequoia.termination_summary(:node_limit, 4, 1.0, -2.0, -3.0, 1e-6; log_sign=-1))
    @test occursin("Optimization incomplete",
                   Sequoia.termination_summary(:incomplete, 4, 1.0, -2.0, -3.0, 1e-6))

    pipe = Pipe()
    Base.link_pipe!(pipe; reader_supports_async=true, writer_supports_async=true)
    redirect_stdout(pipe) do
        solve((lb, ub) -> NodeResult(:infeasible), [0.0], [-1.0]; options=Options(log_every=1e-9))
    end
    close(pipe.in)
    logged = read(pipe, String)
    @test count("Model is infeasible", logged) == 1
    @test occursin("Explored 0 nodes", logged)

    pipe = Pipe()
    Base.link_pipe!(pipe; reader_supports_async=true, writer_supports_async=true)
    redirect_stdout(pipe) do
        solve((lb, ub) -> NodeResult(:infeasible), [0.0], [-1.0]; options=Options(log_every=Inf))
    end
    close(pipe.in)
    @test isempty(read(pipe, String))
end

@testset "gap and show" begin
    r = Result(:optimal, 10.0, 9.0, [1.0], 3, 0, 8.0, 0.1, Sequoia.Snapshot[])
    @test gap(r) ≈ 0.1
    @test gap(Result(:optimal, -10.0, -11.0, nothing, 0, 0, 0.0, 0.0, Sequoia.Snapshot[])) ≈ 0.1
    @test gap(Result(:infeasible, Inf, Inf, nothing, 0, 0, Inf, 0.0, Sequoia.Snapshot[])) == Inf
    s = sprint(show, MIME"text/plain"(), r)
    @test occursin("optimal", s) && occursin("gap", s)
end
