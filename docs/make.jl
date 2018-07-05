# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

using Documenter, FEMStatics

makedocs(modules=[FEMStatics],
         format = :html,
         checkdocs = :all,
         sitename = "FEMStatics.jl",
         pages = ["index.md"]
        )
