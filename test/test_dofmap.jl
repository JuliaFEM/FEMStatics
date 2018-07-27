# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

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
