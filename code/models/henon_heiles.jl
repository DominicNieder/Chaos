using OrdinaryDiffEq
using JSON3  # or JSON

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
    m, a, w = p 
    x, y, px, py = u
    du[1] =  px
    du[2] =  py
    du[3] = -m*w^2*x - 2a*x*y
    du[4] = -m*w^2*y - a*(x^2 - y^2)
end


# Energies, param= (a, m, w)
potential(x, y; p=param) = p[2]*p[3]^2/2 *(x^2+y^2) + p[1]*(x^2*y - y^3/3)

kinetic(px, py; m=param[1]) = (px^2 + py^2)/(2*m)

hamiltonian(x, y, px, py; p=param) = kinetic(px, py, m=p[2]) + potential(x,y, p=p)


# numerical integration, num_int = (T, dt)
function solve_trajectory(u0; tspan=(0.0, num_int[1]), dt=num_int[2], p=param)
    prob = ODEProblem(equations!, u0, tspan, p)
    solve(prob, RK4(), abstol=1e-10, reltol=1e-10, saveat=dt)
end


