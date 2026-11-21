# julia --project=examples examples/ex_milp.jl

"""Example solving a MILP"""

import Pkg
Pkg.activate(@__DIR__)

using JuMP, HiGHS
include(joinpath(@__DIR__, "..", "src", "Sequoia.jl"))

model = Model()
@variable(model, z[1:5], Bin)
@constraint(model, [12, 7, 11, 8, 9]' * z <= 26)
@objective(model, Max, [24, 13, 23, 15, 16]' * z)
set_silent(model)

opt = optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=> false)
result = Sequoia.solve(model; optimizer=opt)
result.status, result.objective, result.x
