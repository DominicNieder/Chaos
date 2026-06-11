using OrdinaryDiffEq
using JSON3  # or JSON

# ---------------------------------------------------
#               Simulation Settings
# ---------------------------------------------------

cfg = JSON3.read(read("../sim_config/henon_heiles.json", String))
# system parameters
a       =   cfg.a.value   
m       =   cfg.m.value   
w       =   cfg.w.value   
param   =   (a,m,w) 
# running variables
E0      =   cfg.E.value  
x0      =   cfg.x0.value
y0      =   cfg.y0.value
py0     =   cfg.py0.value
var     =   (E0, x0, y0, py0)

const FIGURES_DIR = joinpath(@__DIR__, "../../figures")
const DATA_DIR    = joinpath(@__DIR__, "../../data")

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


# Energies
potential(x, y; m=m, w=w, a=a) = m*w^2/2 *(x^2+y^2) + a*(x^2*y - y^3/3)

kinetic(px, py; m=m, w=w) = (px^2 + py^2)/(2*m)

hamiltonian(x, y, px, py; m=m, w=w, a=a) = kinetic(px, py, m=m, w=w) + potential(x,y, m=m, w=w, a=a)


# numerical integration
function solve_trajectory(u0; tspan=(0.0, 500.0), dt=0.1, p=param)
    prob = ODEProblem(equations!, u0, tspan, p)
    solve(prob, Vern9(), abstol=1e-10, reltol=1e-10, saveat=dt)
end


