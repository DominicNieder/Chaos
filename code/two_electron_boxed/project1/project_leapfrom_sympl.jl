import Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
using GLMakie
quarto_theme = joinpath(@__DIR__, "../../styles/makie_theme.jl")
include(quarto_theme)
set_theme!(QUARTO_THEME)

#=
Leapfrog is not adaptive, so wall events can't be handed to a callback the
way Vern9's can: a step that would cross a wall is bisected in TIME (not in
space) until the crossing is found to machine precision, the particle is
placed exactly on the wall, and its velocity is flipped. A fresh leapfrog
sequence then starts from there with the original step size dt.

State here is (x, v), not (x, p) -- leapfrog is naturally written in
velocities. With m1 = m2 = 1 below, v and p are numerically identical, so
the section plot is comparable to the Vern9 file's (x₂, p₂) panel.
=#

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
# integrators (both symplectic; leapfrog is the one you asked for)
# ----------------------------------------------------------------------

"leapfrog / velocity Verlet, 2nd order, kick-drift-kick"
function leapfrog(x, v, h, p)
    a  = accel(x, p)
    xn = (x[1] + h * v[1] + 0.5h^2 * a[1],
          x[2] + h * v[2] + 0.5h^2 * a[2])
    an = accel(xn, p)
    vn = (v[1] + 0.5h * (a[1] + an[1]),
          v[2] + 0.5h * (a[2] + an[2]))
    xn, vn
end

# Yoshida's 4th-order composition of three leapfrog steps -- drop-in upgrade
# if 2nd order isn't tight enough; pass step = yoshida4 to evolve().
const Y1 = 1 / (2 - 2^(1 / 3))
const Y0 = -2^(1 / 3) / (2 - 2^(1 / 3))

function yoshida4(x, v, h, p)
    x, v = leapfrog(x, v, Y1 * h, p)
    x, v = leapfrog(x, v, Y0 * h, p)
    x, v = leapfrog(x, v, Y1 * h, p)
    x, v
end

# ----------------------------------------------------------------------
# walls
# ----------------------------------------------------------------------

function which_wall(x, p)
    for i in 1:2, s in 1:2
        w = p.boxes[i][s]
        if (s == 1 && x[i] < w) || (s == 2 && x[i] > w)
            return (i, s)
        end
    end
    (0, 0)
end

"bisect the step length to find the moment a wall is crossed"
function collision_time(x, v, h, p, step; iters = 60)
    lo, hi = 0.0, h
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
# driver
# ----------------------------------------------------------------------

function evolve(p, x0, v0, T; dt = 1e-3, step = leapfrog,
                section = (1, 2), sample = 0.02)
    x, v, t = Tuple(x0), Tuple(v0), 0.0
    pts = Tuple{Float64,Float64}[]
    ts  = Float64[]
    xs  = Tuple{Float64,Float64}[]
    Es  = Float64[]
    tnext = 0.0

    while t < T
        h = min(dt, T - t)
        xn, vn = step(x, v, h, p)
        i, s = which_wall(xn, p)

        if i == 0
            x, v, t = xn, vn, t + h
        else
            τ = collision_time(x, v, h, p, step)
            x, v = step(x, v, τ, p)
            t += τ
            x = place(x, i, p.boxes[i][s])
            v = reflect(v, i)
            if (i, s) == section
                j = 3 - i
                push!(pts, (x[j], v[j]))
            end
        end

        if t ≥ tnext
            push!(ts, t); push!(xs, x); push!(Es, energy(x, v, p))
            tnext += sample
        end
    end
    ts, xs, Es, pts
end

# ----------------------------------------------------------------------
# run  (same p, initial state, T as the Vern9 file, for a fair comparison)
# ----------------------------------------------------------------------

p = (C     = 1.0,
     m1    = 1.0,
     m2    = 1.0,
     boxes = ((-2.0, -1.0), (1.0, 2.0)))

x0 = (-1.5, 1.3)
v0 = (0.4, -0.2)                     # = (p1, p2) numerically, since m = 1
T  = 200_000.0

ts, xs, Es, pts = evolve(p, x0, v0, T; dt = 1e-3, step = leapfrog, section = (1, 2))

E0 = energy(x0, v0, p)
dE = maximum(abs.(Es .- E0)) / abs(E0)
println("[leapfrog] crossings: ", length(pts), "   relative energy drift: ", dE)

# ----------------------------------------------------------------------
# plots  (same layout as the Vern9 file; no sol(t), so slice the sampled
# trajectory directly instead of calling a continuous interpolant)
# ----------------------------------------------------------------------

fig = Figure(size = (1000, 700))

ax1 = Axis(fig[1, 1:2]; xlabel = "t", ylabel = "x", title = "trajectories")
n = searchsortedlast(ts, min(T, 100.0))
lines!(ax1, ts[1:n], [x[1] for x in xs[1:n]]; label = "particle 1")
lines!(ax1, ts[1:n], [x[2] for x in xs[1:n]]; label = "particle 2")
for (left, right) in p.boxes
    hlines!(ax1, [left, right]; color = C_CREAM, linestyle = :dash)
end
axislegend(ax1; position = :rc)

ax2 = Axis(fig[2, 1]; xlabel = "x₂", ylabel = "p₂",
           title = "section: particle 1 at its right wall")
scatter!(ax2, first.(pts), last.(pts); markersize = 5, color = C_RED)

ax3 = Axis(fig[2, 2]; xlabel = "t", ylabel = "|ΔE / E₀|", yscale = log10,
           title = "energy drift")
drift = abs.(Es .- E0) ./ abs(E0) .+ eps()
lines!(ax3, ts, drift)

save_screen = GLMakie.Screen(; visible = false)
display(save_screen, fig)

name    = "leapfrog01.png"
savename = joinpath("../../../figures/1D-boxed-chr-particles/explore", name)
outfig  = joinpath(@__DIR__, savename)
mkpath(dirname(outfig))
save(outfig, fig, px_per_unit = 1)
close(save_screen)
display(GLMakie.Screen(), fig)