# min x s.t. x I ≽ I, x integer. Optimum is x = 1.
const PSD_DAT = """
1
1
2
1.0
1 1 1 1 1.0
1 1 2 2 1.0
0 1 1 1 1.0
0 1 2 2 1.0
*INTEGERS*
*1
"""

# min -x s.t. x ≤ 2.5, x integer. The relaxation is 2.5; the optimum is x = 2.
const DIAG_DAT = """
1
1
-1
-1.0
1 1 1 1 -1.0
0 1 1 1 -2.5
*INTEGERS*
*1
"""

function dat_problem(text::AbstractString)
    return read_dat_s(IOBuffer(text))
end

@testset "PSD block solves" begin
    p = dat_problem(PSD_DAT)
    @test p.c == [1.0]
    @test p.integers == [1]
    F1 = p.blocks[1].coef[2]
    @test F1 == [1.0 0.0; 0.0 1.0]
    r = solve(jump_model(p); optimizer=OPT, options=QUIET)
    @test r.status === :optimal
    @test r.objective ≈ 1.0 atol = 1e-6
    @test r.x ≈ [1.0] atol = 1e-6
end

@testset "diagonal block branches" begin
    p = dat_problem(DIAG_DAT)
    @test p.blocks[1] isa Sequoia.DiagonalBlock
    @test p.blocks[1].coef[1] == [-2.5]
    @test p.blocks[1].coef[2] == [-1.0]
    r = solve(jump_model(p); optimizer=OPT, options=QUIET)
    @test r.status === :optimal
    @test r.objective ≈ -2.0 atol = 1e-6
    @test r.x ≈ [2.0] atol = 1e-6
end

# min x s.t. 0 ≤ x ≤ 1, with both bounds stored in a diagonal block. Optimum is x = 0.
const BOUNDED_DAT = """
1
1
-2
1.0
1 1 1 1 1.0
0 1 1 1 0.0
1 1 2 2 -1.0
0 1 2 2 -1.0
*INTEGERS*
*1
"""

@testset "diagonal block bounds are finite" begin
    p = dat_problem(BOUNDED_DAT)
    seen = Ref{Any}(nothing)
    r = solve(jump_model(p); optimizer=OPT, options=QUIET,
              propagate=(l, u) -> (seen[] = (copy(l), copy(u)); (l, u)))
    @test r.status === :optimal
    @test r.objective ≈ 0 atol = 1e-6
    @test r.x ≈ [0] atol = 1e-6
    @test seen[][1] == [0.0] && seen[][2] == [1.0]
end

@testset "off-diagonal symmetry" begin
    one = dat_problem("""
    1
    1
    2
    0.0
    1 1 1 2 0.5
    """)
    F = one.blocks[1].coef[2]
    @test F[1, 2] == 0.5 && F[2, 1] == 0.5

    both = dat_problem("""
    1
    1
    2
    0.0
    1 1 1 2 0.5
    1 1 2 1 0.5
    """)
    @test both.blocks[1].coef[2] == [0.0 0.5; 0.5 0.0]

    @test_throws ArgumentError dat_problem("""
    1
    1
    2
    0.0
    1 1 1 2 0.5
    1 1 2 1 0.4
    """)
end

@testset "no variables" begin
    # Constant constraints. An empty sum must not collapse the PSD cone into a bare matrix.
    psd = dat_problem("0\n1\n1\n0 1 1 1 -1.0\n")
    diag = dat_problem("0\n1\n-1\n0 1 1 1 -2.0\n")
    for p in (psd, diag)
        r = solve(jump_model(p); optimizer=OPT, options=QUIET)
        @test r.status === :optimal
        @test r.objective ≈ 0 atol = 1e-6
        @test r.x === nothing || isempty(r.x)
    end
end

@testset "truss topology file" begin
    path = joinpath(@__DIR__, "..", "ref", "instances", "TrussTopology",
                    "demonstsmall_1bar_4scen.dat-s.gz")
    p = read_dat_s(path)
    @test length(p.c) == 28
    @test length(p.blocks) == 2
    @test p.blocks[1] isa Sequoia.PSDBlock
    @test size(p.blocks[1].coef[1]) == (22, 22)
    @test all(issymmetric, p.blocks[1].coef)
    @test p.blocks[2] isa Sequoia.DiagonalBlock
    @test length(p.blocks[2].coef[1]) == 82
    @test p.integers == collect(2:28)
end
