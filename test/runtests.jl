using Test
using Sequoia
using JuMP
using Clarabel
using LinearAlgebra
using Random

include("helpers.jl")

@testset "Sequoia.jl" begin
    @testset "code quality (Aqua)" include("aqua.jl")
    @testset "options and results" include("types.jl")
    @testset "branching rules" include("branching.jl")
    @testset "node handling" include("nodes.jl")
    @testset "core solve" include("core.jl")
    @testset "JuMP interface" include("jump.jl")
    @testset "randomized vs brute force" include("random.jl")
    @testset "review regressions" include("regressions.jl")
end
