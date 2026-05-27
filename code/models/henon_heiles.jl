using OrdinaryDiffEq

const FIGURES_DIR = joinpath(@__DIR__, "../../figures")
const DATA_DIR    = joinpath(@__DIR__, "../../data")

struct HenonHeiles
    E::Float64
end

function equations!(du, u, p, t)
    x, y, px, py = u
    du[1] =  px
    du[2] =  py
    du[3] = -x - 2x*y
    du[4] = -y - x^2 + y^2
end

energy(x, y, px, py) = 0.5(px^2 + py^2) + 0.5(x^2 + y^2) + x^2*y - y^3/3

function solve_trajectory(u0; tspan=(0.0, 500.0), dt=0.01)
    prob = ODEProblem(equations!, u0, tspan)
    solve(prob, Vern9(), abstol=1e-10, reltol=1e-10, saveat=dt)
end
