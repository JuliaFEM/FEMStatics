# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

using JuliaFEM

problem = Problem(Elasticity, "elasticity problem", 2)
element = Element(Quad4, [1, 2, 3, 4])
X = Dict(
    1 => [0.0, 0.0],
    2 => [2.0, 0.0],
    3 => [2.0, 2.0],
    4 => [0.0, 2.0])
update!(element, "geometry", X)

problem.properties.formulation = :plane_stress
add_elements!(problem, [element])
update!(element, "youngs modulus", 288.0)
update!(element, "poissons ratio", 1/3)
time = 0.0
assemble!(problem, time)
k = full(problem.assembly.K)
k[abs.(k).<1.0e-9] = 0
display(k)

L = 2.0
E = 10.0
I = 1.0
k1 = E*I/L^3 * [
      0.0    0.0    0.0      0.0    0.0    0.0
      0.0    12.0   6.0*L    0.0  -12.0    6.0*L
      0.0    6.0*L  4.0*L^2  0.0   -6.0*L  2.0*L^2
      0.0    0.0    0.0      0.0    0.0    0.0
      0.0  -12.0   -6.0*L    0.0   12.0   -6.0*L
      0.0    6.0*L  2.0*L^2  0.0   -6.0*L  4.0*L^2
    ]
k2 = zeros(6,6)
b = [1, 4]
A = 2.0
k2[b,b] = E*A/L * [1.0 -1.0; -1.0 1.0]
display(k1+k2)

problem = Problem(PlaneHeat, "heat problem", 1)
element = Element(Quad4, [1, 2, 3, 4])
X = Dict(
    1 => [0.0, 0.0],
    2 => [2.0, 0.0],
    3 => [2.0, 2.0],
    4 => [0.0, 2.0])
update!(element, "geometry", X)
update!(element, "thermal conductivity", 6.0)
add_elements!(problem, [element])
time = 0.0
assemble!(problem, time)
k = full(problem.assembly.K)
k[abs.(k).<1.0e-9] = 0
display(k)

k = 6*A/L * [1.0 -1.0; -1.0 1.0]
display(k)

fy = 3.0
L = 2.0
fe = [0.0, fy*L, fy*L^2/3, 0.0, fy*L, -fy*L^2/3]
display(fe)
