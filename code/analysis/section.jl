using OrdinaryDiffEq

include("../models/henon_heiles.jl")

# Collect section at y=input, py>0
function surface_section_y0(sol; y_cut=0)
    xs, pxs = Float64[], Float64[]
    for i in 2:length(sol.t)
        y_prev = sol[2, i-1]
        y_curr = sol[2, i]
        if y_prev < y_cut && y_curr >= y_cut
            # linear interpolation to y=y_cut crossing
            α = -y_prev / (y_curr - y_prev)
            x  = sol[1, i-1] + α * (sol[1, i]  - sol[1, i-1])
            px = sol[3, i-1] + α * (sol[3, i]  - sol[3, i-1])
            push!(xs, x)
            push!(pxs, px)
        end
    end
    xs, pxs
end

# Collect section at x=0, px>0
function surface_section_x0(sol; x_cut=0)
    ys, pys = Float64[], Float64[]
    for i in 2:length(sol.t)
        x_prev = sol[1, i-1]
        x_curr = sol[1, i]
        if x_prev < x_cut && x_curr >= x_cut
            # linear interpolation to x=0 crossing
            α = -x_prev / (x_curr - x_prev)
            y  = sol[2, i-1] + α * (sol[2, i]  - sol[2, i-1])
            py = sol[4, i-1] + α * (sol[4, i]  - sol[4, i-1])
            push!(ys, y)
            push!(pys, py)
        end
    end
    ys, pys
end