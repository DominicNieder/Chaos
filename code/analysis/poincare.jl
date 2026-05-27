using OrdinaryDiffEq

include("../models/henon_heiles.jl")

# Collect Poincaré section at y=0, py>0
function poincare_section(sol)
    xs, pxs = Float64[], Float64[]
    for i in 2:length(sol.t)
        y_prev = sol[2, i-1]
        y_curr = sol[2, i]
        if y_prev < 0 && y_curr >= 0
            # linear interpolation to y=0 crossing
            α = -y_prev / (y_curr - y_prev)
            x  = sol[1, i-1] + α * (sol[1, i]  - sol[1, i-1])
            px = sol[3, i-1] + α * (sol[3, i]  - sol[3, i-1])
            push!(xs, x)
            push!(pxs, px)
        end
    end
    xs, pxs
end
