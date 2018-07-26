# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

using Base.Test

#using FEMStatics

using FEMBase

abstract type AbstractDOFMap end

mutable struct DOFMap <: AbstractDOFMap
    local_dof_indices :: Dict{Symbol, Int64}
    permutation :: Dict{Int64, Int64}
    fields :: Dict{String, Vector{Symbol}}
end

function DOFMap()
    return DOFMap(Dict(), Dict(), Dict())
end

function FEMBase.get_gdofs(dofmap::DOFMap, nodes, dof_names)
    ldi = dofmap.local_dof_indices
    perm = dofmap.permutation
    max_dim = length(ldi)
    return (max_dim*(perm[j]-1)+ldi[n] for j in nodes for n in dof_names)
end

function extract_nodal_field_scalar(dofmap, field_name, u)
    sol = Dict{Int64, Float64}()
    ldi = dofmap.local_dof_indices
    max_dim = length(ldi)
    k = ldi[first(dofmap.fields[field_name])]
    for (i, j) in dofmap.permutation
        dof = max_dim*(j-1)+k
        sol[i] = u[dof]
    end
    return sol
end

function extract_nodal_field_vector(dofmap, field_name, u)
    sol = Dict{Int64, Vector{Float64}}()
    ldi = dofmap.local_dof_indices
    dof_names = dofmap.fields[field_name]
    max_dim = length(ldi)
    for (i, j) in dofmap.permutation
        sol[i] = [u[max_dim*(j-1)+ldi[n]] for n in dof_names]
    end
    return sol
end

"""
    extract_nodal_field(dofmap, field_name, u)

Given a solution vector `u`, create dictionary where key is node id and value
is a solution for that node. The type of value depends from the field type,
for scalar field return type is `Dict{Int64, Float64}` and for vector field the
type is `Dict{Int64, Vector{Float64}}`.
"""
function extract_nodal_field(dofmap, field_name, u)
    if length(dofmap.fields[field_name]) == 1
        return extract_nodal_field_scalar(dofmap, field_name, u)
    else
        return extract_nodal_field_vector(dofmap, field_name, u)
    end
end

dofmap = DOFMap()
dofmap.local_dof_indices = Dict(:u1 => 1, :u2 => 2, :u6 => 3, :T => 4)
dofmap.permutation = Dict(1 => 2, 2 => 1)
dofmap.fields = Dict("displacement" => [:u1, :u2, :u6], "temperature" => [:T])
u = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
temperature = extract_nodal_field(dofmap, "temperature", u)
using Base.Test
@test temperature[1] == 8.0
@test temperature[2] == 4.0
displacement = extract_nodal_field(dofmap, "displacement", u)
@test displacement[1] == [5.0, 6.0, 7.0]
@test displacement[2] == [1.0, 2.0, 3.0]

"""
    Static

Run (quasi)static analysis for a set of problems.

Time span (`t0`, `t1`) is divided into number of increments. Time stepping is
controlled using `initial_increment_time`, `min_increment_time` and
`max_increment_time`. In each increment, force equilibrium is achieved by
Newton-Rhapson iterations controlled using `max_iterations` and
`convergence_tolerance`. If system is known to be linear, boolean `linear` can
be used to limit number of iterations to 1.
"""
mutable struct Static <: AbstractAnalysis
    t0 :: Float64
    t1 :: Float64
    initial_increment_time :: Float64
    min_increment_time :: Float64
    max_increment_time :: Float64
    max_iterations :: Int64
    max_increments :: Int64
    convergence_tolerance :: Float64
    linear :: Bool
    dofmap :: AbstractDOFMap
    linear_system_solvers :: Vector{AbstractLinearSystemSolver}
    u :: Vector{Float64}
    la :: Vector{Float64}
end

function Static()
    t0 = 0.0
    t1 = 1.0
    initial_increment_time = Inf
    min_increment_time = 1.0e-5
    max_increment_time = Inf
    max_iterations = 20
    max_increments = 10000#Int((t1-t0)/min_increment_time)
    convergence_tolerance = 1.0e-5
    linear = false
    dofmap = DOFMap()
    linear_system_solvers = [LUSolver()]
    u = zeros(0)
    la = zeros(0)
    return Static(t0, t1, initial_increment_time, min_increment_time,
                  max_increment_time, max_iterations, max_increments,
                  convergence_tolerance, linear, dofmap, linear_system_solvers,
                  u, la)
end

function initialize_dofs!(dofmap, problems; verbose=false)

    info("Setting up degrees of freedom map.")

    # Determine node permutation
    node_ids = Set{Int64}()
    for problem in problems
        for element in get_elements(problem)
            union!(node_ids, get_connectivity(element))
        end
    end
    nnodes = length(node_ids)
    max_node_id = maximum(node_ids)
    # TODO: Here we should implement algorithm to calculate permutation mimizing
    # the bandwidth of assembled matrices e.g. using reverse Cuthill-McKee algorithm
    # For now, just sort nodes according to their index number.
    info("Node permutation algorithm: sort by node id.")
    permutation = Dict(j=>i for (i, j) in enumerate(sort(collect(node_ids))))
    dofmap.permutation = permutation
    if nnodes < 20 || verbose
        info("Node permutation for node i is perm[i], where perm[i] is")
        for j=1:nnodes
            info("perm[$j] = $(permutation[j])")
        end
    end

    # Define "local dof index", determining how local dofs are ordered in node
    j = 1
    ldi = dofmap.local_dof_indices
    info("Global degrees of freedom for node i are (perm[i]-1)+ldi[j], where ldi[j] is")
    for problem in problems
        is_field_problem(problem) || continue
        if !haskey(dofmap.fields, field_name(problem))
            dofmap.fields[field_name(problem)] = []
        end
        for dof in map(first, field_dofs(problem))
            haskey(ldi, dof) && continue
            while j in values(ldi) j += 1 end
            ldi[dof] = j
            info("ldi[$j] => $dof ($(field_name(problem)))")
            push!(dofmap.fields[field_name(problem)], dof)
        end
    end

    nfields = length(dofmap.fields)
    info("DOFMap defines $nfields fields:")
    for fn in keys(dofmap.fields)
        dofmap.fields[fn] = sort(dofmap.fields[fn])
        field_dim = length(dofmap.fields[fn])
        if field_dim == 1
            s = first(dofmap.fields[fn])
            info("$fn (scalar field) => $s")
        else
            s = join(dofmap.fields[fn], ", ")
            info("$fn (vector field) => ($s)")
        end
    end

    info("DOFMap summary:")
    dofs_per_node = length(ldi)
    max_dofs = nnodes*dofs_per_node
    info("Number of nodes in analysis    $nnodes")
    info("Number of dofs per node        $dofs_per_node")
    info("Maximum node id                $max_node_id")
    info("Maximum number of dofs         $max_dofs")
    return dofmap
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

    p = analysis.properties

    dofmap = initialize_dofs!(p.dofmap, get_problems(analysis))

    # set global dofs for each element according to the dofmap
    perm = dofmap.permutation
    for problem in get_problems(analysis)
        is_field_problem(problem) || continue
        dof_names = map(first, field_dofs(problem))
        for element in get_elements(problem)
            connectivity = get_connectivity(element)
            gdofs = get_gdofs(dofmap, connectivity, dof_names)
            set_gdofs!(problem, element, collect(gdofs))
        end
    end

    N = length(dofmap.permutation)*length(dofmap.local_dof_indices)
    info("Maximum number of dofs in assembly: $N.")

    dt = p.initial_increment_time
    time = p.t0
    ls = LinearSystem(N)
    p.u = zeros(N)
    p.la = zeros(N)

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

            # solve linear system and update solution:
            # u[i+1] = u[i] + Δu
            # la[i+1] = la[i]
            solve!(ls, p.linear_system_solvers)
            p.u = p.u + ls.u
            p.la = ls.la

            # update solution (u,la) back to elements
            # for problem in get_problems(analysis)
            #     for element in get_elements(problem)
            #         update_element!(problem, element, u, la, time)
            #     end
            # end

            converged = p.linear || norm(ls.u) < p.convergence_tolerance
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

"""
    update_element!(problem, element, u, la, time)

Given problem and element, update solution (u,la) to element fields.
"""
function update_element!(problem, element, u, la, time)
    gdofs = get_gdofs(problem, element)
    connectivity = get_connectivity(element)
    nnodes = length(element)
    ue = reshape(u[gdofs], field_dim(problem), nnodes)
    lae = reshape(la[gdofs], field_dim(problem), nnodes)
    ued = Dict(j => ue[:,i] for (i,j) in enumerate(connectivity))
    lad = Dict(j => lae[:,i] for (i,j) in enumerate(connectivity))
    update!(element, field_name(problem), time => ued)
    update!(element, "lambda", time => lad)
    return nothing
end

struct Elasticity <: FieldProblem end
struct Beam <: FieldProblem end
struct Heat <: FieldProblem end
struct Dirichlet <: BoundaryProblem end
struct MPC <: BoundaryProblem end

field_name(::Problem{Elasticity}) = "displacement"
field_dofs(::Problem{Elasticity}) = (:u1 => 1, :u2 => 2)

field_name(::Problem{Beam}) = "displacement"
field_dofs(::Problem{Beam}) = (:u1 => 1, :u2 => 2, :u6 => 3)

field_name(::Problem{Heat}) = "temperature"
field_dofs(::Problem{Heat}) = (:T => 1,)

field_dim(problem) = length(field_dofs(problem))

problem1 = Problem(Elasticity, "elasticity problem", 2)
problem2 = Problem(Beam, "beam problem", 3)
problem3 = Problem(Heat, "heat problem", 1)
problem4 = Problem(Dirichlet, "fixed displacement", 2, "displacement")
problem5 = Problem(Dirichlet, "fixed temperature", 1, "temperature")
problem6 = Problem(MPC, "multi-body constraint problem", 3, "displacement")

element1 = Element(Quad4, [5, 1, 2, 6])
element2 = Element(Seg2, [3, 4])
element3 = Element(Poi1, [2])

element4 = Element(Seg2, [6, 5])
element5 = Element(Poi1, [4])
element6 = Element(Poi1, [1])
element7 = Element(Poi1, [3])

add_elements!(problem1, [element1])
add_elements!(problem2, [element2])
add_elements!(problem3, [element1, element2, element3])
add_elements!(problem4, [element4, element5])
add_elements!(problem5, [element4, element5])
add_elements!(problem6, [element6, element7])

analysis = Analysis(Static, "static analysis")
add_problems!(analysis, [problem1, problem2, problem3, problem4, problem5, problem6])
run!(analysis)
