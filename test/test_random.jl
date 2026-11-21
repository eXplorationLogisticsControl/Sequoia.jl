"All integer points of the box `[lb, ub]` (finite, small)."
boxpoints(lb, ub) = (collect(Float64, p) for p in Iterators.product((l:u for (l, u) in zip(lb, ub))...))

@testset "random separable QPs with flaky oracle" begin
    rng = Xoshiro(42)
    for trial in 1:60
        n = rand(rng, 1:4)
        lb = rand(rng, -3:0, n)
        ub = lb .+ rand(rng, 0:3, n)
        tgt = 4 .* randn(rng, n)
        w = rand(rng, n) .+ 0.1
        f(x) = sum(w .* (x .- tgt) .^ 2)
        best = minimum(f, boxpoints(lb, ub))
        pfail = trial <= 30 ? 0.0 : 0.3
        relax = function (l, u)
            rand(rng) < pfail && return NodeResult(:unknown)
            x = clamp.(tgt, l, u)
            return NodeResult(:optimal, f(x), f(x), x)
        end
        for sel in (:best, :depth)
            r = solve(relax, lb, ub; options=Options(QUIET; node_selection=sel))
            check_invariants(r)
            if pfail == 0
                @test r.status === :optimal
            end
            @test r.status in (:optimal, :incomplete)
            r.status === :optimal && @test r.objective ≈ best rtol = 1e-9
            @test r.bound <= best + 1e-9
            r.x === nothing || @test f(r.x) ≈ r.objective
        end
    end
end

@testset "random pure-integer LPs (JuMP + Clarabel) vs enumeration" begin
    rng = Xoshiro(7)
    solved = 0
    for trial in 1:12
        n, m = 3, 3
        A = rand(rng, -4:4, m, n)
        b = rand(rng, 0:6, m)
        c = rand(rng, -5:5, n)
        sense = isodd(trial) ? MIN_SENSE : MAX_SENSE
        pts = [p for p in boxpoints(fill(-2, n), fill(2, n)) if all(A * p .<= b)]
        model = Model()
        @variable(model, -2 <= x[1:n] <= 2, Int)
        @constraint(model, A * x .<= b)
        @objective(model, sense, dot(c, x))
        r = solve(model; optimizer=OPT, options=Options(QUIET; integral_objective=true))
        if isempty(pts)
            @test r.status === :infeasible
        else
            best = (sense == MIN_SENSE ? minimum : maximum)(dot(c, p) for p in pts)
            @test r.status === :optimal
            @test r.objective ≈ best atol = 1e-6
            @test all(A * r.x .<= b)
            @test dot(c, r.x) ≈ best atol = 1e-6
            check_invariants(r; minimize=sense == MIN_SENSE)
            solved += 1
        end
    end
    @test solved >= 6
end
