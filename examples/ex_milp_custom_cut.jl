# julia --project=examples examples/ex_milp_custom_cut.jl

"""Example solving a 0-1 knapsack MILP with a custom node oracle that adds cover cuts"""

import Pkg
Pkg.activate(@__DIR__)

using JuMP, HiGHS
include(joinpath(@__DIR__, "..", "src", "Sequoia.jl"))

# maximize p'z subject to w'z <= W, z binary
w = [23, 31, 29, 44, 53, 38, 63, 85, 89, 82, 47, 52, 61, 35, 28, 71, 57, 43, 66, 39]
p = w .+ 10
W = 518
n = length(w)

"""
    separate_cover(x) -> Union{Vector{Int},Nothing}

Look for a cover `C` (a set of items with `sum(w[C]) > W`) whose inequality
`sum(z[C]) <= |C| - 1` is violated by the relaxed point `x`. Items are added greedily by
increasing `(1 - x[i]) / w[i]`, and the cover is then made minimal, which only increases
the violation. Returns `nothing` if no violated cover is found.
"""
function separate_cover(x; tol=1e-6)
    C, weight = Int[], 0
    for i in sortperm((1 .- x) ./ w)
        push!(C, i)
        weight += w[i]
        weight > W && break
    end
    weight > W || return nothing
    for i in sort(C; by=i -> x[i])
        weight - w[i] > W || continue
        filter!(!=(i), C)
        weight -= w[i]
    end
    return sum(x[C]) > length(C) - 1 + tol ? C : nothing
end

"""
    make_oracle(; cuts=true, max_rounds=10) -> (relax, pool)

Node oracle `relax(lb, ub)` for `Sequoia.solve`. The closure keeps one LP relaxation and
the pool of cover cuts across calls. Cover cuts are valid for every integer point, so
they stay in the model for all later nodes.
"""
function make_oracle(; cuts=true, max_rounds=10)
    m = Model(optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
    @variable(m, 0 <= z[1:n] <= 1)
    @constraint(m, w' * z <= W)
    @objective(m, Max, p' * z)
    pool = Vector{Int}[]

    function relax(lb, ub)
        set_lower_bound.(z, lb)
        set_upper_bound.(z, ub)
        optimize!(m)
        r = Sequoia.jump_result(m, z)
        for _ in 1:max_rounds
            (cuts && r.status === :optimal) || break
            C = separate_cover(r.x)
            C === nothing && break
            push!(pool, C)
            @constraint(m, sum(z[C]) <= length(C) - 1)
            optimize!(m)
            r = Sequoia.jump_result(m, z)
        end
        return r
    end
    return relax, pool
end

opts = Sequoia.Options(log_every=Inf, integral_objective=true)
for cuts in (false, true)
    relax, pool = make_oracle(; cuts)
    # jump_result reports in the minimization sense, so the core minimizes -p'z
    result = Sequoia.solve(relax, zeros(n), ones(n); options=opts, log_sign=-1)
    println("cuts=$cuts: status=$(result.status), objective=$(-result.objective), ",
            "root bound=$(-result.root_bound), nodes=$(result.nodes), cuts added=$(length(pool))")
end

# reference solution from HiGHS
ref = Model(optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
@variable(ref, z[1:n], Bin)
@constraint(ref, w' * z <= W)
@objective(ref, Max, p' * z)
optimize!(ref)
println("HiGHS MIP objective: $(objective_value(ref))")
