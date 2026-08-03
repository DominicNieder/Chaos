import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3
include("../styles/makie_theme.jl")
include("../models/henon_heiles.jl")
using .HenonHeiles

CONFIG_DIR      = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
DATA_DIR        = joinpath(@__DIR__, "../../section_data/henon-heiles/simulation")
FIG_DIR         = joinpath(@__DIR__, "../../figures/henon-heiles/simulation/08-07-sectionANDtrajectory/")

FIG_ORIENTATION  = joinpath(FIG_DIR, "orientation.json")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")

configurations      = JSON3.read(read(CONFIG_DIR, String))
cfg                 = configurations.explore
param               = [Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value)]

T                   = Float64(cfg.T.value)
dt                  = Float64(cfg.dt.value)
x0                  = Float64(cfg.x0.value)

# --- initial position & energy ---






plane = (1, 0.0)                      # x = 0

px_from_E(y, py, E, p) = sqrt(2p[2]*(E - potential(0.0, y, p)) - py^2)

function Tn(v, n, E, p)               # v = [y, py] on the section
    px = px_from_E(v[1], v[2], E, p)
    reinit!(pmap, [0.0, v[1], px, v[2]])
    step!(pmap, n)                    # exactly n intersections, then stops
    u = current_state(pmap)
    return u
end

F(v, n, E, p) = Tn(v, n, E, p) .- v            # solve F(v) = 0, e.g. with NLsolve.jl

"Is the section point v = [y, py] inside the energy boundary?"
in_section(v, E, p) = !isnan(px_from_energy(v[1], v[2], E, p))


E0 = 0.1127
prez = 10^(-10)
sec_boundry = HenonHeiles.section_boundary_ranges(E0, param, 120)
y_max, py_max = sec_boundry[1], sec_boundry[2]
boundary = HenonHeiles.section_boundary(y_max, py_max)



y0  = 0
py0 = 0
u0= [x0, y0, px0, py0]
v = [u0[2], u0[4]]
px0 = px_from_E(y0, py0, E0, param)


damping = 1

ds = CoupledODEs(equations!, u0, param;
                 diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9))
pmap  = PoincareMap(ds, plane; direction = +1)   # dx/dt > 0 ⟺ px > 0

# --- preparing minimization arrays ---
mem_1 =  F(v, py0, E0, param)
mem_2ab = zeros(Float64, 2, 2)
diff_ab = zeros(Float64, 2, 2)
JT      = zeros(Float64, 2, 2)

N_max  = 50
pmap_values = zeros(Float64, 2, N_max)
fig = Figure(size=(1920 , 1200))
ax  = Axis(fig[1,1], xlabel="y", ylabel="p_y")

# --- check for in section ---
in_section(v, E, p) || (v = v, converged = false, iterations = 0, resnorm = Inf)

for i in 1:N_max
    pmap_values[:,i] = v
    v_newy = [v[1]+d, v[2]]
    v_newpy = [v[1], v[2]+d]
    mem_2ab[:,1] = F(v_newy, py0, E0, param)
    mem_2ab[:,2] = F(v_newpy, py0, E0, param)
    diff_ab[:,1] = mem_2ab[:,1] - mem_1
    diff_ab[:,2] = mem_2ab[:,2] - mem_1 
    JT[:,1]      = diff_ab[:,1] ./d
    JT[:,2]      = diff_ab[:,2] ./d

    v = v .- damping .* (JT \ F)

end


scatter!(ax, pmap_values)