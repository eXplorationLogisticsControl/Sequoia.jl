"""
    Node

An open node of the search tree. The integer variables are restricted to `[lb, ub]`.

- `bound`: a lower bound on every point in the node. It is the maximum of the node's own
  relaxation bound and its parent's bound.
- `x`: the relaxed solution, or `nothing` if the relaxation failed or its integral solution
  was rejected. Without `x`, the node is branched with the [`first_free`](@ref) fallback.
"""
struct Node
    lb::Vector{Float64}
    ub::Vector{Float64}
    bound::Float64
    x::Union{Nothing,Vector{Float64}}
end

"Remove and return the next node to process from the nonempty list `open`."
function select_node!(open::Vector{Node}, rule::Symbol)
    rule === :depth && return pop!(open)
    # best bound first; ties go to the most recently created node, which dives
    best, i = Inf, lastindex(open)
    for j in lastindex(open):-1:firstindex(open)
        if open[j].bound < best
            best, i = open[j].bound, j
        end
    end
    return popat!(open, i)
end

"Bounds of the two children of `node` when branching on variable `k` at split value `v`."
function child_boxes(node::Node, k::Integer, v::Real)
    down_ub = copy(node.ub)
    down_ub[k] = floor(v)
    up_lb = copy(node.lb)
    up_lb[k] = ceil(v)
    return ((copy(node.lb), down_ub), (up_lb, copy(node.ub)))
end

"Whether some variable in the box `[lb, ub]` is not fixed."
has_free(lb, ub) = any(i -> lb[i] < ub[i], eachindex(lb, ub))
