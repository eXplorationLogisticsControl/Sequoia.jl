# julia --project=examples examples/ex_misdp.jl

"""Example solving a mixed integer SDP: E-optimal experiment design"""

import Pkg
Pkg.activate(@__DIR__)

using JuMP, Clarabel, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "Sequoia.jl"))

# choose measurement directions within a cost budget to maximize the smallest eigenvalue
# of the information matrix
V = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 1.0 0.0; 0.0 1.0 1.0; 1.0 0.0 1.0]
cost = [3, 3, 3, 2, 2, 2]
n = size(V, 1)

model = Model()
@variable(model, z[1:n], Bin)
@variable(model, t)
@constraint(model, cost' * z <= 6)
@constraint(model, Symmetric(sum(z[i] * V[i, :] * V[i, :]' for i in 1:n) - t * I(3)) in PSDCone())
@objective(model, Max, t)
set_silent(model)

opt = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
result = Sequoia.solve(model; optimizer=opt)
result.status, result.objective, result.x
