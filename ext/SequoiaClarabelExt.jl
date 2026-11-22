module SequoiaClarabelExt

import Clarabel
import Sequoia
using Sequoia: MOI

"""
    ReusableOptimizer(inner::Clarabel.Optimizer)

A Clarabel optimizer that keeps its solver when only the constant of a `Nonnegatives` row
changes. Variable bounds reach Clarabel as such rows, so a branch-and-bound node that moves
bounds re-solves without `copy_to` and `setup!`. Every other change raises a
`NotAllowedError`, and MOI's caching layer then rebuilds the solver.
"""
struct ReusableOptimizer{T} <: MOI.AbstractOptimizer
    inner::Clarabel.Optimizer{T}
end

function Sequoia.reusable_optimizer(inner::Clarabel.Optimizer{T}) where {T}
    cache = MOI.Utilities.CachingOptimizer(MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}()),
                                           ReusableOptimizer(inner))
    bridged = MOI.Bridges.full_bridge_optimizer(cache, T)
    # This bridge adds a bounded variable as `x = y + l` with `y ≥ 0`, so a new `l` would
    # change every row that holds `x`. Without it, a bound is a row of its own.
    MOI.Bridges.remove_bridge(bridged, MOI.Bridges.Variable.VectorizeBridge{T})
    return bridged
end

# Clarabel 0.9 has no `is_data_update_allowed`; there the solver is always rebuilt.
function update_allowed(solver)
    solver === nothing && return false
    isdefined(Clarabel, :is_data_update_allowed) || return false
    return Clarabel.is_data_update_allowed(solver)
end

function MOI.modify(o::ReusableOptimizer{T},
                    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{T},MOI.Nonnegatives},
                    change::MOI.VectorConstantChange{T}) where {T}
    solver = o.inner.solver
    update_allowed(solver) || throw(MOI.ModifyConstraintNotAllowed(ci, change))
    # Clarabel's rows are `Ax + s = b` with `b` equal to the function's constants.
    rows = collect(Clarabel.DefaultInt, o.inner.rowranges[ci.value])
    Clarabel.update_b!(solver, zip(rows, change.new_constant))
    return
end

MOI.empty!(o::ReusableOptimizer) = MOI.empty!(o.inner)
MOI.is_empty(o::ReusableOptimizer) = MOI.is_empty(o.inner)
MOI.optimize!(o::ReusableOptimizer) = MOI.optimize!(o.inner)
MOI.copy_to(o::ReusableOptimizer, src::MOI.ModelLike) = MOI.copy_to(o.inner, src)

function MOI.supports_constraint(o::ReusableOptimizer, F::Type{<:MOI.AbstractFunction},
                                 S::Type{<:MOI.AbstractSet})
    return MOI.supports_constraint(o.inner, F, S)
end

const SolverAttribute = Union{MOI.AbstractOptimizerAttribute,MOI.AbstractModelAttribute}

MOI.supports(o::ReusableOptimizer, attr::SolverAttribute) = MOI.supports(o.inner, attr)
MOI.get(o::ReusableOptimizer, attr::SolverAttribute) = MOI.get(o.inner, attr)
MOI.set(o::ReusableOptimizer, attr::SolverAttribute, value) = MOI.set(o.inner, attr, value)

function MOI.supports(o::ReusableOptimizer, attr::MOI.AbstractVariableAttribute,
                      ::Type{MOI.VariableIndex})
    return MOI.supports(o.inner, attr, MOI.VariableIndex)
end
function MOI.get(o::ReusableOptimizer, attr::MOI.AbstractVariableAttribute, vi::MOI.VariableIndex)
    return MOI.get(o.inner, attr, vi)
end
function MOI.get(o::ReusableOptimizer, attr::MOI.AbstractConstraintAttribute,
                 ci::MOI.ConstraintIndex)
    return MOI.get(o.inner, attr, ci)
end

end # module
