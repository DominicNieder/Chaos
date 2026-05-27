using JLD2

include("../models/henon_heiles.jl")

E    = 1/8
x0   = 0.0
y0   = 0.0
px0  = sqrt(2E - y0^2 - x0^2 - 2x0^2*y0 + 2y0^3/3)
u0   = [x0, y0, px0, 0.0]

println("Simulating Hénon-Heiles at E = $E ...")
sol = solve_trajectory(u0, tspan=(0.0, 1000.0))

outfile = joinpath(DATA_DIR, "trajectory_E$(E).jld2")
jldsave(outfile; t=sol.t, u=sol.u)
println("Saved to $outfile")
