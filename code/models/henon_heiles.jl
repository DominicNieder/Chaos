using OrdinaryDiffEq
using JSON3  # or JSON
using Polynomials
# ---------------------------------------------------
#               Simulation Settings
# ---------------------------------------------------
const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const FIGURES_DIR = joinpath(@__DIR__, "../../figures")
const DATA_DIR    = joinpath(@__DIR__, "../../data")

cfg = JSON3.read(read(CONFIG_DIR, String))

# system parameters
param       =   (cfg.a.value, cfg.m.value, cfg.w.value ) 
# running variables
init_var    =   (cfg.E.value, cfg.x0.value, cfg.y0.value, cfg.py0.value)
# integration variables
num_int     =   (cfg.T.value, cfg.timestep.value)

# ---------------------------------------------------
#                   System equations
# ---------------------------------------------------

struct HenonHeiles
    E::Float64
end

"""
Equation of motion for Hénon-Heiles potential. 
"""
function equations!(du, u, p, t)
    a, m, w = p 
    x, y, px, py = u
    du[1] =  px
    du[2] =  py
    du[3] = -m*w^2*x - 2a*x*y
    du[4] = -m*w^2*y - a*(x^2 - y^2)
end


# Energies, param= (a, m, w)
potential(x::Real, y::Real; p=param) = p[2]*p[3]^2/2 *(x^2+y^2) + p[1]*(x^2*y - y^3/3)

kinetic(px, py; m=param[2]) = (px^2 + py^2)/(2*m)

hamiltonian(x, y, px, py; p=param) = kinetic(px, py, m=p[2]) + potential(x,y, p=p)


"""
makes a function that filters the surface of section values whilst solving the ODE
"""
function section_callback(section_q, section_p)
    condition(u, t, integrator) = u[1]
    affect!(integrator) = begin
        if integrator.u[3] >= 0
            push!(section_q, integrator.u[2])
            push!(section_p, integrator.u[4])
        end
    end
    ContinuousCallback(condition, affect!)
end

"""
u0 -> initial conditions
tspan   =   (tmin, tmax)
dt      =   timestep size
p       =   parameters of the henon heiles modle

Uses Vern9 algorithem: 
"Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy 9th order interpolant)"
"""
function solve_trajectory(
    u0; 
    tspan=(0.0, num_int[1]), 
    dt=num_int[2], 
    p=param, 
    callback=nothing
    )
    prob = ODEProblem(
        equations!, u0, tspan, p
        )
    solve(
        prob, 
        Vern9(), 
        abstol=1e-10, 
        reltol=1e-10, 
        saveat=dt, 
        callback=callback
        )
end


"""
calculates the limit for y0, for surface section x=0.
    Beware that the expression is set at
    alpha   =   1
    m       =   1
    w       =   1
"""
function limit_of_initial_y0(E::Float64)
    p=Polynomial([2E, 0,-1, 2/3])
    roots(p)
end

function limit_of_initial_py0(y::Real, E::Float64)
   sqrt(max(0.0, 2*(E - potential(0.0, y))))
end

function limit_of_initial_py0(y::AbstractVector{<:Real}, E::Float64)
    map(yi -> sqrt( max(0.0, 2*(E - potential(0.0, yi) )) ), y)
end