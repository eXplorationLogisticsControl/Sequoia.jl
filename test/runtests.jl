using Test
using Sequoia
using JuMP
using Clarabel
using LinearAlgebra
using Random

include("helpers.jl")

@testset "Sequoia.jl" begin
    @testset "code quality (Aqua)" include("test_aqua.jl")
    @testset "options and results" include("test_types.jl")
    @testset "branching rules" include("test_branching.jl")
    @testset "node handling" include("test_nodes.jl")
    @testset "core solve" include("test_core.jl")
    @testset "primal heuristics" include("test_heuristics.jl")
    @testset "JuMP interface" include("test_jump.jl")
    @testset "sparse SDPA input" include("test_io.jl")
    @testset "randomized vs brute force" include("test_random.jl")
    @testset "review regressions" include("test_regressions.jl")
end
