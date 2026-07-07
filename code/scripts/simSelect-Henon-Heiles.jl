import Pkg 
Pkg.activate(joinpath(@__DIR__, ".."))
using ProgressMeter, JLD2, JSON3, PolygonOps, LinearAlgebra, GLMakie, BenchmarkTools

model = joinpath(@__DIR__, "../models/henon_heiles.jl")
jl_style = joinpath(@__DIR__, "../styles/makie_theme.jl")
include(jl_style)
include(model)
using .HenonHeiles


CONFIG_DIR      = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
DATA_DIR        = joinpath(@__DIR__, "../../section_data/henon-heiles/simulation")
FIG_DIR         = joinpath(@__DIR__, "../../figures/henon-heiles/simulation")

FIG_ORIENTATION  = joinpath(FIG_DIR, "orientation.json")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")

configurations      = JSON3.read(read(CONFIG_DIR, String))
cfg                 = configurations.SimSelect
param               = (Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value))

T                   = Float64(cfg.T.value)
dt                  = Float64(cfg.dt.value)
x0                  = Float64(cfg.x0.value)

# --- initial position & energy ---
E0 = cfg.E.value

y0 = Float64.(cfg.y0.values)
py0 = Float64.(cfg.py0.values)



println("length of y0: $(length(y0)),\npy0: $(length(py0))")
l_sampling = length(y0)

# --- initialising u0 ---
v       = map(yi -> HenonHeiles.potential(x0, yi, param), y0)  # 
psquare = map(v_i -> 2*param[2]*(E0 - v_i), v)  # p_y^2 if p_x=0 
pdiff   = psquare - py0.^2  
 #println(typeof(pdiff),typeof(yroots))
px0     = map(pdiff_i -> sqrt(max(0.0, pdiff_i)), pdiff)
# --- output ---
u_init     = [[x0, y0[i], px0[i], py0[i]] for i in 1:l_sampling]




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


section_data = Vector{typeof((sec_y=Float64[], sec_py=Float64[]))}(undef, l_sampling)    
traj_data = Vector{typeof((traj_x=Float64[], traj_y=Float64[]))}(undef, l_sampling)

Threads.@threads for i in 1:l_sampling
    u0     = u_init[i]
    y_sec, py_sec = Float64[], Float64[]
    cb = HenonHeiles.section_callback(y_sec, py_sec)
    s = HenonHeiles.solve_trajectory(u0, param, (0.0, T), dt; callback=cb)
    section_data[i] = (sec_y=Float32.(y_sec), sec_py=Float32.(py_sec))
    traj_data[i]    = (traj_x=Float32.(s[1, :]), traj_y=Float32.(s[2,:]))
end





b, ylim, pylim = section_set_plot_lim(E0, param, 120)
fname_sections = "sections-E$(round(E0,digits=4))-T$(cfg.T.value)-py$py0-n$l_sampling.png"

# --- potential contour ---
r      = range(-1.0, 1.0, length=120)
levels = logrange(5.0*0.009, 6.9*0.089, 7)
epot   = [HenonHeiles.potential(x, y, param) for x in r, y in r]

with_theme(QUARTO_THEME) do
    save_screen = GLMakie.Screen(; visible=false)

    # --- fig1: combined section ---
    fig1 = Figure(size=(1920, 1200))
    ax   = Axis(fig1[1,1], xlabel="y", ylabel="p_y")
    scatter!(ax, b)
    colors = cgrad(:viridis, l_sampling)[range(0, 1, length=l_sampling)]

    for i in eachindex(y0)
        scatter!(ax, section_data[i].sec_y, section_data[i].sec_py,
                 color=colors[i], markersize=4)

        # --- fig2: trajectory for IC i ---
        fig2    = Figure(size=(1920, 1200))
        ax_traj = Axis(fig2[1,1], xlabel="x", ylabel="y")
        contour!(ax_traj, r, r, epot, labels=true, levels=levels, colormap=:hsv)
        lines!(ax_traj, traj_data[i].traj_x, traj_data[i].traj_y, linewidth=1.0)

        # --- fig3: single-IC section ---
        fig3 = Figure(size=(1920, 1200))
        ax3  = Axis(fig3[1,1], xlabel="y", ylabel="p_y")
        scatter!(ax3, b)
        scatter!(ax3, section_data[i].sec_y, section_data[i].sec_py,
                 color=colors[i], markersize=4)
        xlims!(ax3, ylim...); ylims!(ax3, pylim...)

        # --- save per-IC figures ---
        fname_traj    = "trajectory-E$(round(E0,digits=4))-T$(cfg.T.value)-y$(y0[i])-py$(py0[i]).png"
        fname_section = "section-E$(round(E0,digits=4))-T$(cfg.T.value)-y$(y0[i])-py$(py0[i]).png"
        display(save_screen, fig2)
        save(joinpath(FIG_DIR, fname_traj), fig2, px_per_unit=1)
        display(save_screen, fig3)
        save(joinpath(FIG_DIR, fname_section), fig3, px_per_unit=1)
    end

    xlims!(ax, ylim...); ylims!(ax, pylim...)

    display(save_screen, fig1)
    save(joinpath(FIG_DIR, fname_sections), fig1, px_per_unit=1)
    close(save_screen)

    display(GLMakie.Screen(), fig1)
end