# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

using Base.Test

#using FEMStatics

using FEMBase

abstract type AbstractDOFMap end

mutable struct SimpleDOFMap <: AbstractDOFMap
    local_dof_indices :: Dict{Symbol, Int64}
    permutation :: Dict{Int64, Int64}
end

function SimpleDOFMap()
    return SimpleDOFMap(Dict(), Dict())
end

function FEMBase.get_gdofs(dofmap::SimpleDOFMap, nodes, dof_names)
    ldi = dofmap.local_dof_indices
    perm = dofmap.permutation
    max_dim = length(ldi)
    return (max_dim*(perm[j]-1)+ldi[n] for j in nodes for n in dof_names)
end

mutable struct Static <: AbstractAnalysis
    t0 :: Float64
    t1 :: Float64
    initial_increment_time :: Float64
    min_increment_time :: Float64
    max_increment_time :: Float64
    max_iterations :: Int64
    dofs :: AbstractDOFMap
end

function Static()
    return Static(0.0, 1.0, Inf, 1.0e-5, Inf, 20, SimpleDOFMap())
end

function initialize_dofs!(dofmap, problems)
    info("Initializing DOFMap.")
    j = 1
    ldi = dofmap.local_dof_indices
    for problem in problems
        dof_abbv = get_dof_abbv(problem)
        info("Problem $(problem.name) dofs = $dof_abbv")
        for dof in dof_abbv
            haskey(ldi, dof) && continue
            while j in values(ldi) j += 1 end
            info("Adding dof $dof to local dof indices as index $j")
            ldi[dof] = j
        end
    end
    info("Dofmap initialized. ldi = $ldi")
    node_ids = Set{Int64}()
    for problem in problems
        for element in get_elements(problem)
            union!(node_ids, get_connectivity(element))
        end
    end
    nnodes = length(node_ids)
    info("Number of nodes in analysis: $nnodes. Maximum node id: $(maximum(node_ids)).")
    # TODO: Add NodeNumbering.jl here.
    sorted_node_ids = sort(collect(node_ids))
    dofmap.permutation = perm = Dict(j=>i for (i, j) in enumerate(sorted_node_ids))
    info("Node permutation: $perm.")
    for problem in problems
        dof_names = get_dof_abbv(problem)
        for element in get_elements(problem)
            connectivity = get_connectivity(element)
            gdofs = get_gdofs(dofmap, connectivity, dof_names)
            set_gdofs!(problem, element, collect(gdofs))
        end
    end
    info("Analysis dofmap initialized.")
end

function initialize_dofs(analysis)
    dofmap = SimpleDOFMap()
    initialize_dofs!(dofmap, get_problems(analysis))
    return dofmap
end

function run!(analysis::Analysis{Static})

    initialize_dofs(analysis)

    p = analysis.properties
    dt = p.initial_increment_time
    time = p.t0
    n_increment = 0
    max_iterations = 2
    done = false
    while !done
        # determine dt
        dt = clamp(dt, p.min_increment_time, p.max_increment_time)
        dt = min(dt, p.t1-time)
        n_increment += 1
        time = min(time+dt, p.t1)
        converged = false
        info("Starting increment $n_increment at time $time.")
        local i
        for i=1:p.max_iterations
            info("Increment # $n_increment. Iteration # $i. Time increment $dt.")
            converged = true
            if converged
                info("Increment # $n_increment converged in $i iterations.")
                break
            end
        end
        if i == p.max_iterations && !converged
            error("Did not converge in $(p.max_iterations) iterations.")
        end
        if isapprox(time, p.t1)
            info("Step completed!")
            done = true
        end
    end

    return nothing
end

struct Dummy1 <: FieldProblem end
struct Dummy2 <: FieldProblem end
struct Dummy3 <: FieldProblem end

function get_dof_abbv(::Problem{Dummy1})
    return (:u1, :u2)
end

function get_dof_abbv(::Problem{Dummy2})
    return (:u1, :u2, :ur3)
end

function get_dof_abbv(::Problem{Dummy3})
    return (:u1, :u2, :u3, :T)
end

problem1 = Problem(Dummy1, "test problem 1", 2)
problem2 = Problem(Dummy2, "test problem 2", 3)
problem3 = Problem(Dummy3, "test problem 3", 4)

element1 = Element(Seg2, [1, 3])
element2 = Element(Seg2, [2, 1])
element3 = Element(Seg2, [1, 5])
add_elements!(problem1, [element1])
add_elements!(problem2, [element1, element2])
add_elements!(problem3, [element1, element3])

analysis = Analysis(Static, "static analysis")
add_problems!(analysis, [problem1, problem2, problem3])
run!(analysis)

info("element 1 gdofs in problem 1: ", get_gdofs(problem1, element1))
info("element 1 gdofs in problem 2: ", get_gdofs(problem2, element1))
info("element 1 gdofs in problem 3: ", get_gdofs(problem3, element1))
info("element 2 gdofs in problem 1: ", get_gdofs(problem1, element2))
info("element 2 gdofs in problem 2: ", get_gdofs(problem2, element2))
info("element 2 gdofs in problem 3: ", get_gdofs(problem3, element2))
info("element 3 gdofs in problem 1: ", get_gdofs(problem1, element3))
info("element 3 gdofs in problem 2: ", get_gdofs(problem2, element3))
info("element 3 gdofs in problem 3: ", get_gdofs(problem3, element3))
