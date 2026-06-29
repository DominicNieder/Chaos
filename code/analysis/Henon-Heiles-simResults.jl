import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using GLMakie, JLD2
include("../models/henon_heiles.jl")
using .HenonHeiles

SIM_DATA_DIR=joinpath(@__DIR__,"../data/henon-heiles/simulation")
EXPL_DATA_DOR = joinpath(@__DIR__, "data/henon-heiles/explore")
FIGURE_DIR = joinpath(@__DIR__, "../figures/henon-heiles/simulation")

println(isdir(SIM_DATA_DIR), isdir(EXPL_DATA_DOR), isdir(FIGURE_DIR))

fsize=(1200,1200)

# --- contour for E=0.1243 ---
e = 0.1243
param = (1,1,1)
n=200
function section_set_plot_lim(e,param, n)
    sec_boundry = section_boundary_ranges(e, param, n)
    y_max, py_max = sec_boundry[1], sec_boundry[2]
    boundary = section_boundary(y_max, py_max)
    ylim = extrema(y_max)
    ydiff = (ylim[2]-ylim[1])/50
    ylim = ylim[1]-ydiff, ylim[2]+ydiff 
    p_max_value = maximum(py_max)
    pydiff = p_max_value/50
    pylim = (-p_max_value-pydiff, p_max_value+pydiff)
    (boundary, ylim, pylim)
end


names1 =[
    "surf_of_sec-T100000-E0.1243-x00.0-y0-0.4-py00.0.jld2",
]
files1 = joinpath(EXPL_DATA_DOR, "surf_of_sec-T100000-E0.1243-x00.0-y0-0.4-py00.0.jld2")
name1 = "phase_map_E$e.png"
fname1 = joinpath(FIGURE_DIR, name1)
println(ispath(file1))
sec1  = load(file1)
py1, y1 = sec1[:py], sec1[:y]
b1, ylim1, pylim1 = section_set_plot_lim(e,param,n)

with_theme(QUARTO_THEME) do
    fig1 = Figure(size=fsize)
    ax  = Axis(fig1[1,1], xlabel="y", ylabel="p_y")
    scatter!(b1)
    scatter!(r[:"y"], r[:"py"])
    ylims!(ax, pylim1...)
    xlims!(ax, ylim1...)
    save(fname1, fig1, px_per_unit=2)
end
