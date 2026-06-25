import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ProgressMeter
using JLD2
using PolygonOps
using LinearAlgebra
include("../models/henon_heiles.jl")
using .HenonHeiles

cfg         = load_config(joinpath(@__DIR__, "../sim_config/henon_heiles.json"))
DATA_DIR    = joinpath(@__DIR__, "../../data/henon-heiles/simulation")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")

param = (cfg.a, cfg.m, cfg.w)


function section_grid(e, param, resolution; n_grid=60)
    y_range, py_range =section_boundary_ranges(e, param, resolution)
    boundary = section_boundary(y_range, py_range)
    y_min, y_max = extrema(y0_range)
    p_max        = maximum(py0_range)

    ys  = range(y_min, y_max, length=n)
    pys = range(-p_max, p_max, length=n)

    grid_y, grid_py = Float64[], Float64[]
    for y in ys, py in pys
        if inpolygon((y, py), boundary) != 0   # 0=outside, 1=inside, -1=on edge
            push!(grid_y, y)
            push!(grid_py, py)
        end
    end
    grid_y, grid_py
end

# --- initial position & energy
E0      =  cfg.E    
x0      =  cfg.x0
resolution = 10

yroots= limit_of_initial_y0(E0, param)
y_min, y_max = yroots[1], yroots[2]

y0, py0 = collect(range(y_min, y_max, resolution)), 0.0
# --- determening px ---
v       = map(yi -> potential(x0, yi, param), y0)  # 
psquare = map(v_i -> 2*param[2]*(E0 - v_i), v)  # p_y^2 if p_x=0 
pdiff   = map(psquare_i -> (psquare_i - py0^2), psquare)  
px0     = map(pdiff_i -> sqrt(max(0.0, pdiff_i)), pdiff)


u_init     = [[x0, y0[i], px0[i], py0] for i in 1:resolution]
simconfigs = [SimConfig(
    cfg.a, cfg.m, cfg.w, 
    E0, x0, y0[i], py0, 
    cfg.T, cfg.dt
    ) for i in 1:resolution
    ]

println("Simulating $resolution trajectories at E = $E0 ...")
io_lock = ReentrantLock()
progress = Progress(resolution; desc="Simulating: ")
Threads.@threads for i in 1:resolution
    config = simconfigs[i]
    u0     = u_init[i]

    y_sec, py_sec = Float64[], Float64[]
    cb = section_callback(y_sec, py_sec)
    s  = solve_trajectory(u0, param, (0.0, config.T), config.dt; callback=cb)

    fname   = "E$(round(E0,digits=4))-x0$(x0)-y0$(round(y0[i],digits=4))-i$(i).jld2"
    outfile = joinpath(DATA_DIR, fname)
    jldsave(outfile; t=s.t, u=s.u, sec_y=y_sec, sec_py=py_sec)

    lock(io_lock) do
        save_config(DATA_ORIENTATION, fname, config)
    end
    next!(progress)
end