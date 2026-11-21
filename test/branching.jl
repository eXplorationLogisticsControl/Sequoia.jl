@testset "most_fractional" begin
    lb, ub = zeros(4), ones(4)
    @test most_fractional([0.0, 0.3, 0.5, 1.0], lb, ub, 1e-6) == 3
    @test most_fractional([0.0, 1e-7, 1.0, 1 - 1e-7], lb, ub, 1e-6) === nothing
    @test most_fractional([0.5, 0.5, 0.0, 0.0], lb, ub, 1e-6) == 1          # ties: lowest index
    @test most_fractional([0.5, 0.4, 0, 0], [0.0, 0, 0, 0], [0.0, 1, 1, 1], 1e-6) == 2  # skips fixed
    @test most_fractional([2.5, -1.25], [-5.0, -5], [5.0, 5], 1e-6) == 1   # general integers
    @test most_fractional(nothing, [1.0, 0, 0], [1.0, 1, 1], 1e-6) == 2    # fallback
end

@testset "first_free" begin
    @test first_free(nothing, [1.0, 0, 0], [1.0, 0, 1], 1e-6) == 3
    @test first_free([0.5], [1.0], [1.0], 1e-6) === nothing
    @test first_free(nothing, Float64[], Float64[], 1e-6) === nothing
end

@testset "split_point" begin
    sp = Sequoia.split_point
    @test sp(0.3, 0.0, 1.0) == 0.3
    @test sp(1.7, 0.0, 1.0) == 0.5          # clamped to u, then made non-integer
    @test sp(nothing, 0.0, 4.0) == 2.5      # midpoint 2 is integer, shifted up
    @test sp(nothing, 0.0, 3.0) == 1.5
    @test sp(nothing, -Inf, 3.0) == 2.5
    @test sp(nothing, 2.0, Inf) == 2.5
    @test sp(nothing, -Inf, Inf) == 0.5
    @test sp(NaN, 0.0, 1.0) == 0.5
    # both children non-empty for every case
    for (xk, l, u) in ((0.3, 0.0, 1.0), (nothing, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 0.0, 1.0),
                       (nothing, -Inf, Inf), (5.0, -Inf, 5.0), (-3.0, -3.0, Inf))
        v = sp(xk, l, u)
        @test !isinteger(v) && l <= floor(v) && ceil(v) <= u
    end
end
