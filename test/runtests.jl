# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

using Base.Test
using Revise

using FEMStatics

revise()

mutable struct Dummy <: FieldProblem end

problem = Problem(Dummy, "test", 3)
analysis = Analysis(Static, "static analysis")
add_problems!(analysis, [problem])
run!(analysis)
