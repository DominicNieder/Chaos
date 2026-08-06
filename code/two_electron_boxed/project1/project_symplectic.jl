#=
Event-driven symplectic integrator for two charges in two disjoint boxes.

Free flight is advanced with a fixed step `dt` using velocity Verlet (or a
4th-order Yoshida composition of it).  When a step would carry a particle out
of its box, the collision time inside that step is located by bisection, the
particle is placed exactly on the wall, its velocity is flipped, and a fresh
constant-dt grid starts from there.

Velocities, not momenta, are the state here -- that is what Verlet wants.
=#
import Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
using GLMakie
quarto_theme=joinpath(@__DIR__,"../../styles/makie_theme.jl")
include(quarto_theme)
set_theme!(QUARTO_THEME)
# ----------------------------------------------------------------------
# physics
# ----------------------------------------------------------------------

function accel(x, p)
    d = x[1] - x[2]
    f = p.C * sign(d) / d^2          # force on particle 1
    (f / p.m1, -f / p.m2)
end

energy(x, v, p) = 0.5p.m1 * v[1]^2 + 0.5p.m2 * v[2]^2 + p.C / abs(x[1] - x[2])

# ----------------------------------------------------------------------
# one step
# ----------------------------------------------------------------------

"velocity Verlet, 2nd order, out of place"
function verlet(x, v, h, p)
    a  = accel(x, p)
    xn = (x[1] + h * v[1] + 0.5h^2 * a[1],
          x[2] + h * v[2] + 0.5h^2 * a[2])
    an = accel(xn, p)
    vn = (v[1] + 0.5h * (a[1] + an[1]),
          v[2] + 0.5h * (a[2] + an[2]))
    xn, vn
end

# Yoshida's 4th-order composition of three Verlet steps
const Y1 = 1 / (2 - 2^(1 / 3))
const Y0 = -2^(1 / 3) / (2 - 2^(1 / 3))

function yoshida4(x, v, h, p)
    x, v = verlet(x, v, Y1 * h, p)
    x, v = verlet(x, v, Y0 * h, p)
    x, v = verlet(x, v, Y1 * h, p)
    x, v
end

# ----------------------------------------------------------------------
# walls
# ----------------------------------------------------------------------

"""
    which_wall(x, p) -> (i, s)

Returns the particle and side of the first violated wall, or (0, 0) if the
configuration is legal. s = 1 is the left wall, s = 2 the right one.
"""
function which_wall(x, p)
    for i in 1:2, s in 1:2
        w = p.boxes[i][s]
        if (s == 1 && x[i] < w) || (s == 2 && x[i] > w)
            return (i, s)
        end
    end
    (0, 0)
end

"""
    collision_time(x, v, h, p, step) -> τ

Smallest τ ∈ (0, h] for which a step of length τ from (x, v) ends outside a
box. Plain bisection on the step length: robust, and the cost is only paid
once per collision. `step` is the integrator used for free flight, so the
landing point is consistent with the rest of the trajectory.
"""
function collision_time(x, v, h, p, step; iters = 60)
    lo, hi = 0.0, h                  # lo: known legal, hi: known illegal
    for _ in 1:iters
        mid = 0.5 * (lo + hi)
        xm, _ = step(x, v, mid, p)
        which_wall(xm, p) == (0, 0) ? (lo = mid) : (hi = mid)
    end
    hi
end

reflect(v, i) = i == 1 ? (-v[1], v[2]) : (v[1], -v[2])
place(x, i, w) = i == 1 ? (w, x[2]) : (x[1], w)

# ----------------------------------------------------------------------
# the driver
# ----------------------------------------------------------------------

"""
    evolve(p, x0, v0, T; dt, step, section, sample)

Integrate to time `T`. Records the other particle's (x, v) each time the wall
named by `section = (particle, side)` is hit. `sample` is the interval at
which the trajectory is stored for plotting.

Returns (ts, xs, Es, pts).
"""
function evolve(p, x0, v0, T; dt = 1e-3, step = yoshida4,
                section = (1, 2), sample = 0.05)

    x, v, t = Tuple(x0), Tuple(v0), 0.0
    pts  = Tuple{Float64,Float64}[]
    ts   = Float64[]
    xs   = Tuple{Float64,Float64}[]
    Es   = Float64[]
    tnext = 0.0

    while t < T
        h = min(dt, T - t)
        xn, vn = step(x, v, h, p)
        i, s = which_wall(xn, p)

        if i == 0                                   # ordinary free-flight step
            x, v, t = xn, vn, t + h
        else                                        # collision inside this step
            τ = collision_time(x, v, h, p, step)
            x, v = step(x, v, τ, p)
            t += τ
            x = place(x, i, p.boxes[i][s])          # snap exactly onto the wall
            v = reflect(v, i)
            if (i, s) == section
                j = 3 - i
                push!(pts, (x[j], v[j]))
            end
        end

        if t ≥ tnext                                # thin out the trajectory
            push!(ts, t); push!(xs, x); push!(Es, energy(x, v, p))
            tnext += sample
        end
    end
    ts, xs, Es, pts
end

# ----------------------------------------------------------------------
# example
# ----------------------------------------------------------------------

p = (C = 1.0, m1 = 1.0, m2 = 1.0,
     boxes = ((-2.0, -1.0), (1.0, 2.0)))

x0 = (-1.5, 1.3)
v0 = (0.4, -0.2)                     # velocities, since m = 1 these equal p
T  = 200_000.0

ts, xs, Es, pts = evolve(p, x0, v0, T; dt = 1e-3, step = yoshida4)

E0 = energy(x0, v0, p)
println("collisions recorded: ", length(pts),
        "   max |ΔE/E₀|: ", maximum(abs.(Es .- E0)) / abs(E0))

# ----------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------

fig = Figure(size = (1000, 700))

ax1 = Axis(fig[1, 1:2]; xlabel = "t", ylabel = "x", title = "trajectories")
n = searchsortedfirst(ts, 100.0)
lines!(ax1, ts[1:n], [x[1] for x in xs[1:n]]; label = "particle 1")
lines!(ax1, ts[1:n], [x[2] for x in xs[1:n]]; label = "particle 2")
for (left, right) in p.boxes
    hlines!(ax1, [left, right]; linestyle = :dash)
end
axislegend(ax1; position = :rc)

ax2 = Axis(fig[2, 1]; xlabel = "x₂", ylabel = "v₂",
           title = "section: particle 1 at its right wall")
scatter!(ax2, first.(pts), last.(pts); markersize = 3)

ax3 = Axis(fig[2, 2]; xlabel = "t", ylabel = "|ΔE / E₀|", yscale = log10,
           title = "energy drift")
lines!(ax3, ts, abs.(Es .- E0) ./ abs(E0) .+ eps())

name    = "time_predict_alg01.png"
savename = joinpath("../../../figures/1D-boxed-chr-particles/explore", name)
outfig  = joinpath(@__DIR__, savename)
mkpath(dirname(outfig))
save(outfig, fig, px_per_unit = 1)
close(save_screen)
display(GLMakie.Screen(), fig)