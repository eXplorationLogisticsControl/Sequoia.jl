# Shared fixtures for the test suite.

const OPT = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
const QUIET = Sequoia.Options(log_every=Inf)

"All 0/1 vectors of length `n`."
binvecs(n) = (digits(i; base=2, pad=n) for i in 0:(2^n - 1))

"""
Analytic oracle for `min Σ (xᵢ - tᵢ)²` over a box: the relaxation optimum is the
projection of `t` onto the box, and its value is exact (bound == value).
"""
function projection_oracle(t; calls=Ref(0))
    return function (lb, ub)
        calls[] += 1
        x = clamp.(t, lb, ub)
        f = sum(abs2, x .- t)
        return Sequoia.NodeResult(:optimal, f, f, x)
    end
end

"""
Analytic oracle for the LP `min x₁ + 2x₂  s.t. x₁ + x₂ ≥ 1` over the box `[lb, ub]`:
greedy fill of the cheaper variable first.
"""
function lp_oracle(lb, ub)
    x1 = clamp(1 - lb[2], lb[1], ub[1])
    x2 = clamp(1 - x1, lb[2], ub[2])
    x1 + x2 < 1 - 1e-12 && return Sequoia.NodeResult(:infeasible)
    f = x1 + 2x2
    return Sequoia.NodeResult(:optimal, f, f, [x1, x2])
end

"Check the invariants every finished `Result` must satisfy."
function check_invariants(r; minimize=true)
    s = minimize ? 1 : -1
    @test r.nodes >= 0 && r.unknown <= r.nodes
    @test s * r.bound <= s * r.objective + 1e-9
    bounds = [s * snap.bound for snap in r.trace]
    incs = [s * snap.incumbent for snap in r.trace]
    @test all(incs[i+1] <= incs[i] for i in 1:length(incs)-1)   # incumbent never gets worse
    @test all(bounds .<= incs .+ 1e-9)            # bound below incumbent at all times
    if r.status === :optimal   # default tolerances
        @test abs(r.objective - r.bound) <= max(1e-8, 1e-6 * abs(r.objective)) + 1e-12
    elseif r.status === :infeasible
        @test r.objective == Inf && r.x === nothing
    end
end
