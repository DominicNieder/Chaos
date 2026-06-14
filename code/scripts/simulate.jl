using JLD2
using LinearAlgebra
include("../models/henon_heiles.jl")



# --- initial position & energy
E       =  0.005:0.005:0.1
x0      = -0.4:0.05:0.4
y0      = -0.4:0.05:0.4

# calculate max momentum from energy restraint
v       = m*w^2/2*(x^2 + y^2) + a*(x^2*y - y^3/3)
py_max    = norm(2*e*m - v)

px0     =  0.0:py_max/10:py_max
py0     = sqrt(2E - y0^2 - x0^2 - 2x0^2*y0 + 2y0^3/3)
u0      = [x0, y0, px0, 0.0]

println("Simulating Hénon-Heiles at E = $E ...")
sol = solve_trajectory(u0; 
    tspan=(0.0, num_int[1]), dt=num_int[2], p=param
    )

outfile = joinpath(DATA_DIR, "trajectory_E$(E).jld2")
jldsave(outfile; t=sol.t, u=sol.u)
println("Saved to $outfile")
