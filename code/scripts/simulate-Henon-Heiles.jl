import Pkg 
Pkg.activate(joinpath(@__DIR__, ".."))
using ProgressMeter, JLD2, JSON3, PolygonOps, LinearAlgebra, GLMakie, BenchmarkTools

model = joinpath(@__DIR__, "../models/henon_heiles.jl")
jl_style = joinpath(@__DIR__, "../styles/makie_theme.jl")
include(jl_style)
include(model)
using .HenonHeiles


CONFIG_DIR      = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
DATA_DIR        = joinpath(@__DIR__, "../../data/henon-heiles/simulation")
FIG_DIR         = joinpath(@__DIR__, "../../figures/henon-heiles/simulation")

FIG_ORIENTATION  = joinpath(FIG_DIR, "orientation.json")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")

configurations      = JSON3.read(read(CONFIG_DIR, String))
cfg                 = configurations.sim
param               = (Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value))

sample_resolution   = Int64(cfg.y0.resolution)
T                   = Float64(cfg.T.value)
dt                  = Float64[] # Float64(cfg.dt.value)
x0                  = Float64(cfg.x0.value)

# --- initial position & energy ---
e_min = Float64(cfg.E.range.min)                             # 0.01
e_max =  Float64(cfg.E.range.max)                            # 0.167
energy_resolution = Int64(cfg.E.range.resolution)          # 10
    # 128
# e0 = range(e_min,e_max, energy_resolution)          # Energies
e0 = vcat(
    logrange(e_min-0.0005, 0.1, energy_resolution),          # coarse at low E
    range(0.1, e_max, energy_resolution)               # linear, denser at high E
) |> unique |> sort
E0= e0[1:2:end]


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



function initializing_y0_sampling(E0::Float64, x0::Float64, sample_resolution::Int64, param)
    yroots= HenonHeiles.limit_of_initial_y0(E0, param)
    y_min, y_max = yroots[1], yroots[2]
    y0, py0 = collect(range(y_min, y_max, sample_resolution)), 0.0

    # --- determening px ---
    v       = map(yi -> HenonHeiles.potential(x0, yi, param), y0)  # 
    psquare = map(v_i -> 2*param[2]*(E0 - v_i), v)  # p_y^2 if p_x=0 
    pdiff   = map(psquare_i -> (psquare_i - py0^2), psquare)  
     #println(typeof(pdiff),typeof(yroots))
    px0     = map(pdiff_i -> sqrt(max(0.0, pdiff_i)), pdiff)
    # --- output ---
    u_init     = [[x0, y0[i], px0[i], py0] for i in 1:sample_resolution]
    data_e0 = Vector{typeof((sec_y=Float64[], sec_py=Float64[]))}(undef, sample_resolution)    
    u_init, data_e0
end


function simulate_y0_lattice(E0::Float64, T::Float64, dt::Vector{Float64}, x0::Float64, sampling::Int64, param::NTuple{3,Float64}, DATA_DIR::String)

    u_init, data = initializing_y0_sampling(E0, x0, sampling, param)

    println("Simulating $sampling trajectories at E = $E0 ")
    progress = Progress(sampling; desc="Simulating: ")
    Threads.@threads for i in 1:sampling
        u0     = u_init[i]

        y_sec, py_sec = Float64[], Float64[]

        cb = HenonHeiles.section_callback(y_sec, py_sec)
        HenonHeiles.solve_trajectory(u0, param, (0.0, T), dt; callback=cb)

        data[i] = (sec_y=Float32.(y_sec), sec_py=Float32.(py_sec))

        next!(progress)
    end
    fname   = "E$(round(E0,digits=4))-T$(T)-py0.0-n$sampling.jld2"
    outfile = joinpath(DATA_DIR, fname)
    jldsave(outfile; results=data)
end



# for (i, e) in enumerate(E0[end-4:end]) 
#     simulate_y0_lattice(e, T, dt, x0, sample_resolution, param, DATA_DIR)
# end





b, ylim, pylim  = section_set_plot_lim(E0, param, 120)
fname   = "E$(round(E0,digits=4))-T$(cfg.T.value)-py$py0-n$sample_resolution.png"
outfig = joinpath(FIG_DIR, fname)

with_theme(QUARTO_THEME) do
    fig1 = Figure(size=(1200,1200))
    ax  = Axis(fig1[1,1], xlabel="y", ylabel="p_y")
    scatter!(ax, b)
    colors = cgrad(:viridis, sample_resolution)[range(0, 1, length=sample_resolution)]    
    for i in eachindex(data)
        println("plotting initial condition \ny0: $(simconfigs[i].y0)")
        scatter!(ax, data[i].sec_y, data[i].sec_py, color=colors[i], markersize=3)
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

