# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

# Constraint

using FEMBase

"""
    SparseMatrixFEM

A variant of IJV/COO matrix format, where matrix indices are stored, instead of
integer vectors I and J, to tuples (nid,dof) where nid is node id and dof is
name of degree of freedom. This type can be converted to a SparseMatrixCOO with
the help of DOFMap, detemining the actual indices where values should be added.
"""
struct SparseMatrixFEM{Tv,Ti<:Tuple{Int,Symbol}} <: AbstractSparseMatrix{Tv,Ti}
    I :: Vector{Ti}
    J :: Vector{Ti}
    V :: Vector{Tv}
end

function SparseMatrixFEM()
    I = Tuple{Int64,Symbol}[]
    J = Tuple{Int64,Symbol}[]
    V = Float64[]
    return SparseMatrixFEM(I,J,V)
end

function FEMBase.add!(A::SparseMatrixFEM{Tv,Ti}, dof1::Ti, dof2::Ti, value::Tv) where {Tv,Ti}
    push!(A.I, dof1)
    push!(A.J, dof2)
    push!(A.V, value)
    return nothing
end

####

struct SparseVectorFEM{Tv,Ti<:Tuple{Int,Symbol}} <: AbstractSparseVector{Tv,Ti}
    I :: Vector{Ti}
    V :: Vector{Tv}
end

function SparseVectorFEM()
    I = Tuple{Int64,Symbol}[]
    V = Float64[]
    return SparseVectorFEM(I,V)
end

function FEMBase.add!(A::SparseVectorFEM{Tv,Ti}, dof::Ti, value::Tv) where {Tv,Ti}
    push!(A.I, dof1)
    push!(A.V, value)
    return nothing
end

function Base.display(A::SparseMatrixFEM)
    # TODO
    return nothing
end

function Base.display(A::SparseVectorFEM)
    # TODO
    return nothing
end

A = SparseMatrixFEM()
b = SparseVectorFEM()
add!(A, (1,:u1), (1,:u1), 1.0)

"""
    Constraint - General type nodal constraint

Constraints are in a vector of tuples, where each tuple is
    ((node_id, dof_name, multiplier), coefficient)

Constraints are assembled in preprocess before actual assemble using DOFMap to
get actual matrix rows and columns. Thus type has two presentations, one for
user input and another one ready for assembly procedure.

Internal presentation: constraints are in a vector of tuples, where each tuple is
    (row, gdofs, multipliers, coefficient)

For example, equation `1.5*N1[:u1] + 1.0*N2[:u2] == 0.5` is given

```julia
input_data = ((1, :u1, 1.5), (2, :u2, 1.0), 0.5)
bc = Constraint()
push!(bc.constraints, input_data)
```
"""
struct Constraint <: BoundaryProblem
    C1 :: SparseMatrixFEM
    C2 :: SparseMatrixFEM
    g :: SparseVectorFEM
end

function Constraint()
    C1 = SparseMatrixFEM()
    C2 = SparseMatrixFEM()
    g = SparseVectorFEM()
    return Constraint(C1, C2, g)
end

function add_constraint!(problem, node_id::Int64, dof_name::Union{Colon, Symbol}, value)
    info("add spc for node $node_id, dofs: $dof_name, value: $value")
    p = problem.properties
    # The main difference to standard add!(C1, gdof, gdof, 1.0) comes here.
    # Instead of giving gdof as explicit integer corresponding to some matrix
    # location, give a tuple containing node number and dof name, which is
    # later on determined by dofmap (after initialization)
    gdof = (node_id, dof_name)
    push!(p.C1, gdof, gdof, 1.0)
    push!(p.C2, gdof, gdof, 1.0)
    push!(p.g, gdof, value)
    return nothing
end

function add_constraint!(problem, node_ids, dof_name, value)
    for node_id in node_ids
        add_constraint!(problem, node_id, dof_name, value)
    end
    return nothing
end

function add_constraint!(problem, node_ids, dof_names::Tuple, value)
    for dof_name in dof_names
        add_constraint!(problem, node_ids, dof_name, value)
    end
    return nothing
end

function add_constraint!(problem, multipliers::NTuple{N,T}, coefficient) where {N,T<:Tuple{Int64,Symbol,Float64}}
    info("adding multipoint constraint with $N multipliers")
    info(T)
    C1 = []
    C2 = []
    for (node_id, dof_name, multiplier) in multipliers
        push!(C1, (:C1, node_id, dof_name, multiplier))
        push!(C2, (:C2, node_id, dof_name, multiplier))
    end
    push!()
end

problem4 = Problem(Constraint, "nodal constraints", 0, "N/A")
add_constraint!(problem4, 2, :T, 100.0)
add_constraint!(problem4, (4, 5), :, 0.0)
add_constraint!(problem4, 6, (:u1, :u2, :u6, :T), 0.0)
add_constraint!(problem4, ((1, :u1, 1.0), (3, :u1, 1.0)), 0.0)
add_constraint!(problem4, ((1, :u2, 1.0), (3, :u2, 1.0)), 0.0)
