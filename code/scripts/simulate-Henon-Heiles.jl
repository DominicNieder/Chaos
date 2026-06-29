import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using ProgressMeter
using JLD2
using PolygonOps
using LinearAlgebra
model = joinpath(@__DIR__, "../models/henon_heiles.jl")
include(model)
using .HenonHeiles
using GLMakie
jl_style = joinpath(@__DIR__, "../styles/makie_theme.jl")
include(jl_style)


cfg         = HenonHeiles.load_config(joinpath(@__DIR__, "../sim_config/henon_heiles.json"))
DATA_DIR    = joinpath(@__DIR__, "../../data/henon-heiles/simulation")
FIG_DIR     = joinpath(@__DIR__, "../../figures/henon-heiles/simulation")

FIG_ORIENTATION  = joinpath(FIG_DIR, "orientation.json")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")
param = (cfg.a, cfg.m, cfg.w)


function section_grid(e, param, resolution; n_grid=60)
    y_range, py_range = HenonHeiles.section_boundary_ranges(e, param, resolution)
    boundary = HenonHeiles.section_boundary(y_range, py_range)
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

"""

parameters: e, param, n
piocare section boundary
-> (boundary, ylim, pylim)
"""
function section_set_plot_lim(e,param, n)
    sec_boundry = HenonHeiles.section_boundary_ranges(e, param, n)
    y_max, py_max = sec_boundry[1], sec_boundry[2]
    boundary = HenonHeiles.section_boundary(y_max, py_max)
    ylim = extrema(y_max)
    ydiff = (ylim[2]-ylim[1])/50
    ylim = ylim[1]-ydiff, ylim[2]+ydiff 
    p_max_value = maximum(py_max)
    pydiff = p_max_value/50
    pylim = (-p_max_value-pydiff, p_max_value+pydiff)
    (boundary, ylim, pylim)
end

# --- initial position & energy
E0      =  cfg.E    
x0      =  cfg.x0
resolution = 256

yroots= HenonHeiles.limit_of_initial_y0(E0, param)
y_min, y_max = yroots[1], yroots[2]

y0, py0 = collect(range(y_min, y_max, resolution)), 0.0
println("yo\n$y0")
# --- determening px ---
v       = map(yi -> HenonHeiles.potential(x0, yi, param), y0)  # 
psquare = map(v_i -> 2*param[2]*(E0 - v_i), v)  # p_y^2 if p_x=0 
pdiff   = map(psquare_i -> (psquare_i - py0^2), psquare)  
px0     = map(pdiff_i -> sqrt(max(0.0, pdiff_i)), pdiff)


u_init     = [[x0, y0[i], px0[i], py0] for i in 1:resolution]
simconfigs = [HenonHeiles.SimConfig(
    cfg.a, cfg.m, cfg.w, 
    E0, x0, y0[i], py0, 
    cfg.T, cfg.dt
    ) for i in 1:resolution
    ]
println(length(u_init),length(simconfigs))
println("Simulating $resolution trajectories at E = $E0 ...")
# io_lock = ReentrantLock()
progress = Progress(resolution; desc="Simulating: ")

results = Vector{NamedTuple}(undef, resolution)

Threads.@threads for i in 1:resolution
    config = simconfigs[i]
    u0     = u_init[i]

    y_sec, py_sec = Float64[], Float64[]

    cb = HenonHeiles.section_callback(y_sec, py_sec)
    s  = HenonHeiles.solve_trajectory(u0, param, (0.0, config.T), config.dt; callback=cb)

    results[i] = (sec_y=y_sec, sec_py=py_sec)

    next!(progress)
end

b, ylim, pylim  = section_set_plot_lim(E0, param,120)
fname   = "E$(round(E0,digits=4))-T$(cfg.T)-py$py0-n$resolution.png"
outfig = joinpath(FIG_DIR, fname)

with_theme(QUARTO_THEME) do
    fig1 = Figure(size=(1200,1200))
    ax  = Axis(fig1[1,1], xlabel="y", ylabel="p_y")
    scatter!(ax, b)
    colors = cgrad(:viridis, resolution)[range(0, 1, length=resolution)]
    for i in eachindex(results)
        println("plotting initial condition \ny0: $(simconfigs[i].y0)")
        scatter!(ax, results[i].sec_y, results[i].sec_py, color=colors[i], markersize=3)
    end
    ylims!(ax, pylim...)
    xlims!(ax, ylim...)
    # --- save figure ---
    save_screen = GLMakie.Screen(; visible=false)
    display(save_screen, fig1)
    save(outfig, fig1, px_per_unit=2)
    close(save_screen)

    display(GLMakie.Screen(), fig1)
end

fname   = "E$(round(E0,digits=4))-T$(cfg.T)-py$py0-n$resolution.jld2"
outfile = joinpath(DATA_DIR, fname)
jldsave(outfile; results=results)