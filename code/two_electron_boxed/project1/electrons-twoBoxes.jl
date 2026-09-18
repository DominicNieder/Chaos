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
