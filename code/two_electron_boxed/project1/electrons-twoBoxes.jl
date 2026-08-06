import Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
using OrdinaryDiffEq
using GLMakie
using DiffEqCallbacks, NonlinearSolve, ADTypes
quarto_theme=joinpath(@__DIR__,"../../styles/makie_theme.jl")
include(quarto_theme)
set_theme!(QUARTO_THEME)

# ----------------------------------------------------------------------
# equations of motion
# ----------------------------------------------------------------------
 
function eom!(du, u, p, t)
    x1, x2, p1, p2 = u
    d = x1 - x2
    f = p.C * sign(d) / d^2        # force on particle 1
    du[1] = p1 / p.m1
    du[2] = p2 / p.m2
    du[3] =  f
    du[4] = -f
end

V_int(u,p)= p.C / abs(u[1] - u[2])  # interaction potential

Kin(u,p)= u[3]^2 / (2p.m1) + u[4]^2 / (2p.m2)  # total kinetic energy

energy(u, p) = Kin(u,p) + V_int(u,p)  # total energy

# ----------------------------------------------------------------------
# all possible init conditions on manningfold of energy()=E
# ----------------------------------------------------------------------

"""
    Determening the momentum p2 (second particle)

u0 = [x1,x2,p1, _]
"""
function init_u0(u0, E, p)
    K = E - V_int(u0, p) - u0[3]^2/(2p.m1)
    K ≥ 0 || error("no real p₂: energy E is below the potential + p₁ contribution")
    [u0[1], u0[2], u0[3], sqrt(2p.m2 * K)]
end


function poincare_boundary(E,p)
    
end

# ----------------------------------------------------------------------
# hard walls
# ----------------------------------------------------------------------
 
"""
    wall_cb(p, i, s, section, pts)
 
One elastic wall: particle `i`, side `s` (1 = left, 2 = right). The wall
position and the outward direction are baked in by the closure, so no event
index has to be decoded.
"""
function wall_cb(p, i, s, section, pts)
    wall    = p.boxes[i][s]
    outward = s == 1 ? -1.0 : 1.0
    mom     = 2 + i                  # momentum of particle i lives at u[2+i]
 
    condition(u, t, integrator) = u[i] - wall
 
    function affect!(integrator)
        # reflect only if the particle is actually leaving; after a bounce it
        # sits exactly on the wall, so the event can fire a second time
        if sign(integrator.u[mom]) == outward
            integrator.u[mom] = -integrator.u[mom]
            if (i, s) == section
                j = 3 - i            # the other particle
                push!(pts, (integrator.u[j], integrator.u[2+j]))
            end
        end
    end
 
    ContinuousCallback(condition, affect!; interp_points = 20)
end
 
"""
    wall_callback(p; section = (1, 2))
 
Elastic reflection at all four walls. Every time the wall given by `section`
is hit, the *other* particle's (x, p) is recorded. Returns `(cb, pts)`.
"""
function wall_callback(p; section = (1, 2))
    pts = Tuple{Float64,Float64}[]
    cbs = [wall_cb(p, i, s, section, pts) for i in 1:2 for s in 1:2]
    return CallbackSet(cbs...), pts
end
 

"""
    energy_projection(p, u0; kw...)

Projects the state back onto the surface of constant energy E₀ = energy(u0, p).
"""
function energy_projection(p, u0; kw...)
    E0 = energy(u0, p)
    g!(resid, u, pp, t) = (resid[1] = energy(u, pp) - E0)
    ManifoldProjection(g!; autodiff = AutoForwardDiff(),
                       resid_prototype = zeros(1), kw...)
end

# ----------------------------------------------------------------------
# example
# ----------------------------------------------------------------------
 
p = (C     = 2.0,
     m1    = 1.0,
     m2    = 1.0,
     boxes = ((-2.0, -1.0), (1.0, 2.0)))

u1    = [-1.5, 1.3, 0.4, -0.2]              # x1, x2, p1, p2
E0 = energy(u1, p)

u0 = u1

tspan = (0.0, 2_000_000.0)

poincare_sec_at = (1, 2)  # record when particle 1 hits its right wall

cbwalls, section_pts = wall_callback(p; section = poincare_sec_at)   
proj = energy_projection(p, u0; nlopts = Dict(:ftol => 1e-14))

prob    = ODEProblem(eom!, u0, tspan, p)
sol     = solve(prob, Vern9(lazy=Val{false}()); callback = CallbackSet(cbwalls, proj), abstol = 1e-12, reltol = 1e-12)  # for not saving all positions: , saveat = range(tspan..., length = 5000)
 

dE = maximum(abs(energy(u, p) - E0) for u in sol.u) / abs(E0)
println("crossings: ", length(section_pts), "   relative energy drift: ", dE)
 
# ----------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------

fig = Figure(size = (1000, 700))
    
# --- trajectories over a short window, with the walls drawn in ---------
ax1 = Axis(fig[1, 1:2]; xlabel = "t", ylabel = "x", title = "trajectories")
tt = range(tspan[1], min(tspan[2], 1000.0); length = 2000)
lines!(ax1, tt, [sol(t)[1] for t in tt]; label = "particle 1")
lines!(ax1, tt, [sol(t)[2] for t in tt]; label = "particle 2")
for (left, right) in p.boxes
    hlines!(ax1, [left, right]; color = C_CREAM, linestyle = :dash)
end
axislegend(ax1; position = :rc)
    
# --- Poincaré section --------------------------------------------------
ax2 = Axis(fig[2, 1]; xlabel = "x₂", ylabel = "p₂",
        title = "section: particle 1 at its right wall")
scatter!(ax2, first.(section_pts), last.(section_pts); markersize = 5, color = C_RED)
    
# --- energy conservation ----------------------------------------------
ax3 = Axis(fig[2, 2]; xlabel = "t", ylabel = "|ΔE / E₀|", yscale = log10,
        title = "energy drift")
drift = [abs(energy(u, p) - E0) / abs(E0) + eps() for u in sol.u]
lines!(ax3, sol.t, drift)
save_screen = GLMakie.Screen(; visible=false)
display(save_screen, fig)

name = "interaction$(p.C)06.png"
savename=joinpath("../../../figures/1D-boxed-chr-particles/explore", name)
outfig=joinpath(@__DIR__, savename)
save(outfig, fig, px_per_unit=1)
close(save_screen)
display(GLMakie.Screen(), fig)