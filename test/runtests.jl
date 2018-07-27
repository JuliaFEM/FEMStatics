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

    info("Setting global dofs for elements")
    perm = dofmap.permutation
    for problem in problems
        is_field_problem(problem) || continue
        dof_names = map(first, field_dofs(problem))
        for element in get_elements(problem)
            connectivity = get_connectivity(element)
            gdofs = get_gdofs(dofmap, connectivity, dof_names)
            set_gdofs!(problem, element, collect(gdofs))
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
            field_names = keys(p.dofmap.fields)
            fields = Dict(fn => extract_nodal_field(p.dofmap, fn, p.u) for fn in field_names)
            fields_la = Dict(fn => extract_nodal_field(p.dofmap, fn, p.la) for fn in field_names)
            for problem in get_problems(analysis)
                for (fn, fv) in fields
                    update!(problem, fn, time => fv)
                    update!(problem, "$fn (lambda)", time => fv)
                end
            end

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

X = Dict(
    1 => [2.0, 0.0],
    2 => [2.0, 2.0],
    3 => [2.0, 0.0],
    4 => [4.0, 0.0],
    5 => [0.0, 0.0],
    6 => [0.0, 2.0])

elements = [element1, element2, element3, element4, element5, element6, element7]
update!(elements, "geometry", X)

function FEMBase.assemble_elements!(problem::Problem{Elasticity}, assembly::Assembly,
                            elements::Vector{Element{E}}, time) where E
    # local stiffness matrix calculated for E=288, nu=1/3 element
    ke = [
    144.0   54.0  -90.0    0.0  -72.0  -54.0   18.0    0.0
     54.0  144.0    0.0   18.0  -54.0  -72.0    0.0  -90.0
    -90.0    0.0  144.0  -54.0   18.0    0.0  -72.0   54.0
      0.0   18.0  -54.0  144.0    0.0  -90.0   54.0  -72.0
    -72.0  -54.0   18.0    0.0  144.0   54.0  -90.0    0.0
    -54.0  -72.0    0.0  -90.0   54.0  144.0    0.0   18.0
     18.0    0.0  -72.0   54.0  -90.0    0.0  144.0  -54.0
      0.0  -90.0   54.0  -72.0    0.0   18.0  -54.0  144.0
    ]
    element = first(elements)
    gdofs = get_gdofs(problem, element)
    info("elasticity gdofs = $(collect(gdofs))")
    add!(assembly.K, gdofs, gdofs, ke)
end

function FEMBase.assemble_elements!(problem::Problem{Heat}, assembly::Assembly,
                            elements::Vector{Element{E}}, time) where E <: Seg2
    # local conductance matrix calculated for k = 6.0, A = 2.0, L = 2.0
    ke = [
    6.0  -6.0
   -6.0   6.0
    ]
    element = first(elements)
    gdofs = get_gdofs(problem, element)
    info("heat seg2 gdofs = $(collect(gdofs))")
    add!(assembly.K, gdofs, gdofs, ke)
end

function FEMBase.assemble_elements!(problem::Problem{Heat}, assembly::Assembly,
                            elements::Vector{Element{E}}, time) where E <: Quad4
    # local conductance matrix calculated for k = 6.0
    ke = [
    4.0  -1.0  -2.0  -1.0
   -1.0   4.0  -1.0  -2.0
   -2.0  -1.0   4.0  -1.0
   -1.0  -2.0  -1.0   4.0
    ]
    element = first(elements)
    gdofs = get_gdofs(problem, element)
    info("heat quad4 gdofs = $(collect(gdofs))")
    add!(assembly.K, gdofs, gdofs, ke)
end

function FEMBase.assemble_elements!(problem::Problem{Beam}, assembly::Assembly,
                            elements::Vector{Element{E}}, time) where E
    # local stiffness matrix calculated for L = 2.0, E = 10.0, I = 1.0, A = 2.0
    ke = [
    10.0    0.0    0.0  -10.0    0.0    0.0
     0.0   15.0   15.0    0.0  -15.0   15.0
     0.0   15.0   20.0    0.0  -15.0   10.0
   -10.0    0.0    0.0   10.0    0.0    0.0
     0.0  -15.0  -15.0    0.0   15.0  -15.0
     0.0   15.0   10.0    0.0  -15.0   20.0
    ]
    # local force vector calculated for fy = 3.0, L = 2.0
    fe = [0.0, 6.0, 4.0, 0.0, 6.0, -4.0]
    element = first(elements)
    gdofs = get_gdofs(problem, element)
    info("beam gdofs = $(collect(gdofs))")
    add!(assembly.K, gdofs, gdofs, ke)
    add!(assembly.f, gdofs, fe)
end

add_elements!(problem1, [element1])
add_elements!(problem2, [element2])
add_elements!(problem3, [element1, element2, element3])
add_elements!(problem4, [element4, element5])
add_elements!(problem5, [element4, element5])
add_elements!(problem6, [element6, element7])

time = 0.0
dofmap = initialize_dofs!(DOFMap(), [problem1, problem2, problem3])
assemble!(problem1, time)
assemble!(problem2, time)
assemble!(problem3, time)
N = 24
K1 = full(problem1.assembly.K, N, N)
K2 = full(problem2.assembly.K, N, N)
K3 = full(problem3.assembly.K, N, N)
f1 = full(problem1.assembly.f, N, 1)
f2 = full(problem2.assembly.f, N, 1)
f3 = full(problem3.assembly.f, N, 1)
display(round.(Int, K1))
display(round.(Int, K2))
display(round.(Int, K3))
K = K1+K2+K3
f = f1+f2+f3

function diagc(K)
    p = ones(size(K,2))
    p[unique(rowvals(K))] = 0.0
    return dropzeros(spdiagm(p))
end
diagc(K::Matrix) = diagc(sparse(K))

display(diagc(K))

display(round.(Int, K))
display(round.(Int, f))
free_dofs = 1:4*3
u = zeros(f)
u[free_dofs] = (K+diagc(K))[free_dofs, free_dofs] \ f[free_dofs]
display(u)
info("displacements:")
for (k,v) in extract_nodal_field(dofmap, "displacement", u)
    println("$k => $v")
end
info("temperature:")
for (k,v) in extract_nodal_field(dofmap, "temperature", u)
    println("$k => $v")
end
# analysis = Analysis(Static, "static analysis")
# add_problems!(analysis, [problem1, problem2, problem3, problem4, problem5, problem6])
# run!(analysis)
