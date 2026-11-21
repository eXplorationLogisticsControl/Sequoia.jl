const Node = Sequoia.Node

@testset "select_node!" begin
    mk(b) = Node([0.0], [1.0], b, nothing)
    open = [mk(3.0), mk(1.0), mk(2.0), mk(1.0)]
    last_best = open[4]
    n = Sequoia.select_node!(open, :best)
    @test n === last_best
    @test length(open) == 3 && count(m -> m.bound == 1.0, open) == 1
    # ties go to the most recent node
    open = [mk(1.0), mk(5.0), mk(1.0)]
    third = open[3]
    @test Sequoia.select_node!(open, :best) === third
    # depth first is LIFO
    open = [mk(1.0), mk(5.0)]
    @test Sequoia.select_node!(open, :depth).bound == 5.0
end

@testset "child_boxes" begin
    node = Node([0.0, 0.0], [3.0, 1.0], 0.0, [1.4, 0.0])
    (dl, du), (ul, uu) = Sequoia.child_boxes(node, 1, 1.4)
    @test dl == [0, 0] && du == [1, 1] && ul == [2, 0] && uu == [3, 1]
    @test node.lb == [0, 0] && node.ub == [3, 1]     # parent untouched
    @test !(dl === node.lb) && !(uu === node.ub)     # children own their bounds
end

@testset "has_free" begin
    @test Sequoia.has_free([0.0, 1.0], [0.0, 2.0])
    @test !Sequoia.has_free([0.0, 1.0], [0.0, 1.0])
    @test !Sequoia.has_free(Float64[], Float64[])
end
