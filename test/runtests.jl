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

"""
    Static

Static analysis runs analysis from `t0` to `t1` in several increments. In each
increment, force equilibrium is achieved by Newton-Rhapson iterations.
"""
mutable struct Static <: AbstractAnalysis
    t0 :: Float64
    t1 :: Float64
    initial_increment_time :: Float64
    min_increment_time :: Float64
    max_increment_time :: Float64
    max_iterations :: Int64
    max_increments :: Int64
    dofmap :: AbstractDOFMap
end

function Static()
    return Static(0.0, 1.0, Inf, 1.0e-5, Inf, 20, 10000, SimpleDOFMap())
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
    dofmap.permutation = Dict(j=>i for (i, j) in enumerate(sorted_node_ids))
    info("Analysis dofmap initialized.")
end

struct LUSolver <: AbstractLinearSystemSolver end

"""
    solve!(ls::LinearSystem, solver::LUSolver)

Solve linear system using LU factorization. If linear system has zero rows,
1.0 is added to diagonal to make matrix non-singular.
"""
function FEMBase.solve!(ls::LinearSystem, solver::LUSolver)

    info("Solving linear system using LUSolver")

    A = [ls.K ls.C1; ls.C2 ls.D]
    b = [ls.f; ls.g]

    # add 1.0 to diagonal for any zero rows in system
    p = ones(2*ls.dim)
    p[unique(rowvals(A))] = 0.0
    A += spdiagm(p)

    # solve A*x = b using LU factorization and update solution vectors
    x = lufact(A) \ full(b)
    ls.u = x[1:ls.dim]
    ls.la = x[ls.dim+1:end]

    return nothing
end

function run!(analysis::Analysis{Static})

    dofmap = SimpleDOFMap()
    initialize_dofs!(dofmap, get_problems(analysis))

    # set global dofs for each problem according to the dofmap
    perm = dofmap.permutation
    info("Node permutation: $perm.")
    for problem in get_problems(analysis)
        dof_names = get_dof_abbv(problem)
        for element in get_elements(problem)
            connectivity = get_connectivity(element)
            gdofs = get_gdofs(dofmap, connectivity, dof_names)
            set_gdofs!(problem, element, collect(gdofs))
        end
    end

    N = length(dofmap.permutation)*length(dofmap.local_dof_indices)
    info("Maximum number of dofs in assembly: $N.")

    p = analysis.properties
    dt = p.initial_increment_time
    time = p.t0
    ls = LinearSystem(N)

    for n=1:p.max_increments
        # determine dt
        dt = clamp(dt, p.min_increment_time, p.max_increment_time)
        dt = min(dt, p.t1-time)
        time = min(time+dt, p.t1)

        info("Starting increment $n at time $time.")
        for i=1:p.max_iterations
            info("Increment # $n. Iteration # $i. Time increment $dt.")

            # Assemble matrices for all problems
            for problem in get_problems(analysis)
                assemble!(problem, time)
            end

            # Collect everything to LinearSystem
            for problem in get_problems(analysis)
                ls.K += sparse(problem.assembly.K, N, N)
                ls.Kg += sparse(problem.assembly.Kg, N, N)
                ls.f += sparse(problem.assembly.f, N, 1)
                ls.C1 += sparse(problem.assembly.C1, N, N)
                ls.C2 += sparse(problem.assembly.C2, N, N)
                ls.D += sparse(problem.assembly.D, N, N)
                ls.g += sparse(problem.assembly.g, N, 1)
            end

            # solve linear system
            solvers = [LUSolver()]
            solve!(ls, solvers)

            # update solution back to elements
            for problem in get_problems(analysis)
                field_name = get_unknown_field_name(problem)
                field_dim = get_unknown_field_dimension(problem)

            end

            converged = true
            if converged
                info("Increment # $n converged in $i iterations.")
                break
            end
            if i == p.max_iterations
                error("Solution of nonlinear system did not converge in $i ",
                      "iterations.")
            end
        end
        if isapprox(time, p.t1)
            info("Step completed in $n increments.")
            break
        end
        if n == p.max_increments
            error("Step did not finish in $n increments.")
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
