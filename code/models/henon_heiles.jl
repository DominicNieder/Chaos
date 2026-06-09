using OrdinaryDiffEq

const FIGURES_DIR = joinpath(@__DIR__, "../../figures")
const DATA_DIR    = joinpath(@__DIR__, "../../data")

struct HenonHeiles
    E::Float64
end


function equations!(du, u, p, t)
    m = 1
    a = 1
    w = 1
    x, y, px, py = u
    du[1] =  px
    du[2] =  py
    du[3] = -m*w^2*x - 2a*x*y
    du[4] = -m*w^2*y - a*(x^2 + y^2)
end

potential(x, y; m=1, w=1, a=1) = m*w^2/2 (x^2+y^2) + a*(x^2*y - y^3/3)

kinetic(px, py; m=1, w=1, a=1) = (px^2 + py^2)/2

hamiltonian(x, y, px, py; m=1, w=1, a=1) = kinetic(px, py, m=1, w=1, a=1) + potential(x,y, m=1, w=1, a=1)

function solve_trajectory(u0; tspan=(0.0, 500.0), dt=0.01)
    prob = ODEProblem(equations!, u0, tspan)
    solve(prob, Vern9(), abstol=1e-10, reltol=1e-10, saveat=dt)
end


